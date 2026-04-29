//! Phase 6c-D — write Phase-6c-A-substituted boundary curves into
//! BRep edges so they propagate to the cap-vs-fillet pcurve fits.
//!
//! Without this wiring, the substituted curves in `BlendBoundary`
//! never reach the BRep edges that the cap-outer-boundary metric
//! reads. `apply_three_quad_decomposition_to_blend` calls this
//! helper after the poke (which creates seam edges with `None` 3-D
//! curves) and before `reproject_loop_pcurves_onto_gregory` (which
//! samples `edge.curve` to fit pcurves on the new Gregory surface).
//!
//! Activation: same env-gate as the substitution itself
//! (`POPKERN_SPHERE_FIT_BOUNDARY=1`). Falls back silently when no
//! match is found between BRep half-edges and boundary segments.
//!
//! # Match heuristic
//!
//! For each cap sub-face's outer-loop half-edge:
//! * Look up the half-edge's origin / destination vertex 3-D
//!   positions.
//! * Pick the boundary segment whose endpoints (in either
//!   direction) minimise the sum of distances to the half-edge's
//!   vertices.
//! * Accept the match if the total distance is below a slack
//!   tolerance (5 cm — wide enough for the rolling-ball /
//!   sphere-fit endpoint coincidence on uniform-radius fixtures,
//!   tight enough to avoid false matches across the apex).

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{CurveGeometry, Point3, fitting};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::blend::boundary::{BlendBoundary, BoundarySegment};

/// For every cap sub-face, replace the 3-D curve on each outer-
/// loop edge that matches a `BlendBoundary` segment with that
/// segment's curve. Returns the number of edges rewritten — useful
/// for diagnostic logging and tests.
pub(crate) fn wire_substituted_boundary_to_edges(
    brep: &mut BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> usize {
    // Match tolerance: rolling-ball footprint endpoints coincide
    // with sphere-fit termination endpoints to ~1e-12 on the
    // uniform-radius cube. Variable-radius cases either drift
    // beyond 5 cm (substitution wouldn't have fired anyway) or
    // sit within the rolling-ball overshoot region.
    //
    // Phase 6c-D Stage 1: optional 30 cm tolerance + vertex
    // relocation under `POPKERN_GREGORY_CAP_REBUILD=1`. Closes the
    // BRep ↔ substituted-boundary alignment gap for fixtures where
    // Path B reroute can't run (pyramid topology — neighbour
    // cylinder axes don't pass through per-fillet sphere centres).
    // Default off so cube_mixed/extreme's Path B path is undisturbed.
    let do_relocate = std::env::var("POPKERN_GREGORY_CAP_REBUILD").is_ok_and(|v| v == "1");
    let tol = if do_relocate { 0.30 } else { 0.05 };

    let segment_curves: Vec<&NurbsCurve> = boundary
        .segments
        .iter()
        .map(|seg| match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        })
        .collect();
    let segment_endpoints: Vec<[popkern_core::geometry::Point3; 2]> = segment_curves
        .iter()
        .map(|c| {
            let (t0, t1) = c.parameter_range();
            [c.point_at(t0), c.point_at(t1)]
        })
        .collect();

    // Phase A — assignment scan: for each cap sub-face's outer-loop
    // half-edge, find the best matching segment (proximity, with
    // the wider tolerance). Build a worklist of
    // `(he, seg_idx, reversed, target_p_origin, target_p_dest)`.
    // We do NOT mutate during the scan so the matching sees
    // pre-relocation positions consistently.
    let mut consumed = vec![false; segment_curves.len()];
    type Job = (HalfEdgeHandle, usize, bool, Point3, Point3);
    let mut worklist: Vec<Job> = Vec::new();

    for &sub_face in sub_faces {
        let Some(face_data) = brep.faces.get(sub_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;

        let hes = collect_loop_half_edges(brep, start_he);

        for he_handle in hes {
            let Some(he) = brep.half_edges.get(he_handle) else {
                continue;
            };
            let v_origin = he.origin;
            let next_he = he.next;
            let Some(next_he_data) = brep.half_edges.get(next_he) else {
                continue;
            };
            let v_dest = next_he_data.origin;

            let p_origin = brep.vertices[v_origin].position;
            let p_dest = brep.vertices[v_dest].position;

            let mut best: Option<(usize, f64, bool)> = None;
            for (i, eps) in segment_endpoints.iter().enumerate() {
                if consumed[i] {
                    continue;
                }
                let d_fwd = p_origin.distance_to(&eps[0]) + p_dest.distance_to(&eps[1]);
                let d_rev = p_origin.distance_to(&eps[1]) + p_dest.distance_to(&eps[0]);
                let (d, reversed) = if d_fwd <= d_rev {
                    (d_fwd, false)
                } else {
                    (d_rev, true)
                };
                if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                    best = Some((i, d, reversed));
                }
            }
            let Some((seg_idx, d, reversed)) = best else {
                continue;
            };
            if d > tol {
                continue;
            }
            consumed[seg_idx] = true;

            let eps = segment_endpoints[seg_idx];
            let (target_p_origin, target_p_dest) = if reversed {
                (eps[1], eps[0])
            } else {
                (eps[0], eps[1])
            };
            worklist.push((he_handle, seg_idx, reversed, target_p_origin, target_p_dest));
        }
    }

    // Phase B — execute. For each job, snap both endpoint vertices
    // to the segment endpoints (idempotent if already coincident),
    // write the curve to the edge (oriented to match the edge's
    // primary half-edge direction), and refit pcurves on both sides.
    let mut wired = 0usize;
    for (he_handle, seg_idx, _reversed, target_p_origin, target_p_dest) in worklist {
        let Some(he) = brep.half_edges.get(he_handle) else {
            continue;
        };
        let v_origin = he.origin;
        let edge_h = he.edge;
        let next_he = he.next;
        let Some(next_he_data) = brep.half_edges.get(next_he) else {
            continue;
        };
        let v_dest = next_he_data.origin;

        // Snap vertices first when relocation is enabled (pyramid
        // path — no Path B). Keeps loop closure invariants when
        // the curve below is written. snap_vertex_to syncs *all*
        // incident half-edges' edge curves CP[0]/CP[last] to the
        // new vertex position, so sibling edges' curves still join
        // cleanly at the moved vertex.
        if do_relocate {
            snap_vertex_to(brep, v_origin, target_p_origin);
            if v_dest != v_origin {
                snap_vertex_to(brep, v_dest, target_p_dest);
            }
        }

        // Orient the curve to match the edge's primary half-edge
        // direction: edge.curve's CP[0] should equal primary
        // half-edge's origin position. If our `he_handle` IS the
        // primary, curve.CP[0] := target_p_origin; if it's the twin,
        // curve.CP[0] := target_p_dest.
        let edge_primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let our_is_primary = edge_primary == Some(he_handle);
        let primary_origin_target = if our_is_primary {
            target_p_origin
        } else {
            target_p_dest
        };

        let segment = segment_curves[seg_idx];
        let (s_t0, s_t1) = segment.parameter_range();
        let seg_cp0 = segment.point_at(s_t0);
        let seg_cp_last = segment.point_at(s_t1);
        let curve = if seg_cp0.distance_to(&primary_origin_target)
            <= seg_cp_last.distance_to(&primary_origin_target)
        {
            segment.clone()
        } else {
            reverse_nurbs(segment)
        };
        let (t0, t1) = curve.parameter_range();

        // Capture face surfaces for pcurve refit (this side + twin)
        // before we mutate the edge curve.
        let twin_he = brep.half_edges[he_handle].twin;
        let cap_face = brep
            .loops
            .get(brep.half_edges[he_handle].loop_)
            .map(|l| l.face);
        let other_face = brep
            .loops
            .get(brep.half_edges[twin_he].loop_)
            .map(|l| l.face);
        let cap_surf = cap_face.and_then(|f| {
            brep.faces
                .get(f)
                .and_then(|fd| fd.surface.as_ref())
                .map(|s| s.geometry.clone())
        });
        let other_surf = other_face.and_then(|f| {
            brep.faces
                .get(f)
                .and_then(|fd| fd.surface.as_ref())
                .map(|s| s.geometry.clone())
        });

        if let Some(edge_data) = brep.edges.get_mut(edge_h) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve.clone()),
                t_start: t0,
                t_end: t1,
            });
            wired += 1;
        }

        // Re-fit both half-edge pcurves from the new 3-D curve.
        if let Some(cap_surf) = cap_surf.as_ref()
            && let Some(pc) = refit_pcurve_on_surface(&curve, cap_surf, t0, t1)
        {
            brep.half_edges[he_handle].pcurve = Some(pc.into());
        }
        if let Some(other_surf) = other_surf.as_ref()
            && let Some(pc) = refit_pcurve_on_surface(&curve, other_surf, t0, t1)
        {
            brep.half_edges[twin_he].pcurve = Some(pc.into());
        }
    }
    wired
}

/// Phase 6c-D: reverse a NurbsCurve by flipping CPs, weights, and
/// reflecting the knot vector. Equivalent to
/// `popkern_core::geometry::edge_triple::analytic_helpers::reverse_nurbs`
/// (kept private there); duplicated locally for cohesion.
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cps = curve.control_points.clone();
    reversed_cps.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cps,
        weights: reversed_weights,
    }
}

/// Phase 6c-D: snap a vertex to `new_pos`, syncing every incident
/// half-edge's edge curve CP[0] (primary side) or CP[last] (twin
/// side) to keep loop closure intact. No-op when already within
/// 1e-9 of `new_pos`.
///
/// This is a drop-in equivalent of
/// `var_radius::update_incident_edges_for_vertex_move` — kept
/// local to the wire pass for cohesion. Both functions produce
/// identical results on identical inputs.
fn snap_vertex_to(brep: &mut BRep, v: VertexHandle, new_pos: Point3) {
    let Some(vd) = brep.vertices.get(v) else {
        return;
    };
    if vd.position.distance_to(&new_pos) < 1e-9 {
        return;
    }

    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        let primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let is_primary = primary == Some(he_h);
        let Some(edge) = brep.edges.get_mut(edge_h) else {
            continue;
        };
        let Some(coe) = edge.curve.as_mut() else {
            continue;
        };
        if let CurveGeometry::Nurbs(nc) = &mut coe.geometry
            && !nc.control_points.is_empty()
        {
            let n = nc.control_points.len();
            if is_primary {
                nc.control_points[0] = new_pos;
            } else {
                nc.control_points[n - 1] = new_pos;
            }
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
}

/// Sample `curve_3d` at uniform parameter values, project each
/// sample to `surface`'s UV, and fit a degree-3 NURBS pcurve
/// through those (u, v, 0) points. Returns `None` if any
/// projection fails or the fit raises an error.
fn refit_pcurve_on_surface(
    curve_3d: &NurbsCurve,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
) -> Option<NurbsCurve> {
    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let t = t0 + s * (t1 - t0);
        let p3 = curve_3d.point_at(t);
        let (u, v) = project_point_to_uv(surface, &p3)?;
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut hes = Vec::new();
    let mut cursor = start;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    hes
}
