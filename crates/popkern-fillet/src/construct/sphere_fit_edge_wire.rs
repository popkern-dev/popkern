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
use popkern_core::geometry::pcurve::Pcurve;
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
    //
    // Phase 6c-D Stage 4: empirical finding — variable-radius
    // pyramid fixtures need the wider tolerance to match ANY
    // boundary segment to a cap sub-edge (`wire_substituted… wrote
    // 0 of 8` at 5cm tolerance on pyramid4_mixed). The wider
    // tolerance + relocation `do_relocate` couples two effects:
    // wider matching AND vertex snap. The vertex snap breaks
    // cube_mixed by collapsing the apex topology. Decouple so we
    // can match-without-snap on injection-active fixtures by
    // default — analytic boundary curve gets written to BRep
    // edges, vertices stay at their existing positions. Edge
    // endpoints may not exactly equal curve(t_start)/curve(t_end),
    // but downstream pcurve refit / measurement tolerates this.
    let do_relocate = std::env::var("POPKERN_GREGORY_CAP_REBUILD").is_ok_and(|v| v == "1");
    let stage4_wide_tol = boundary.sphere_fit_substituted
        && std::env::var("POPKERN_STAGE4_WIDE_WIRE_TOL").is_ok_and(|v| v != "0");
    let tol = if do_relocate || stage4_wide_tol {
        0.30
    } else {
        0.05
    };

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
    //
    // Phase 6c-D Stage 4: when injection is active, allow multiple
    // sub-edges to share a single segment (`stage4_wide_tol`). The
    // cap face's outer loop after sliver-absorb topology surgery
    // can have 5+ sub-edges along a single logical boundary
    // (FilletEnd[k] / FaceEdge[k]); the original "consume each
    // segment once" policy left short sub-edges unwired, holding
    // worst-case `cap_outer_g0` at ~10 cm. Letting multiple
    // sub-edges share a segment as their 3D curve writes the
    // analytic boundary to ALL of them — endpoint mismatch
    // tolerated since downstream consumers (validation, pcurve
    // refit) operate on the curve geometry, not the vertex-to-
    // curve-endpoint identity.
    let allow_multi_sub_edge = stage4_wide_tol;
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
                if consumed[i] && !allow_multi_sub_edge {
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
            if !allow_multi_sub_edge {
                consumed[seg_idx] = true;
            }

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
    for (he_handle, seg_idx, reversed, target_p_origin, target_p_dest) in worklist {
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
        //
        // Phase 6c-D Stage 4: under injection, the cap-side pcurve
        // refit (`refit_pcurve_on_surface` → `project_point_to_uv` →
        // `fit_curve`) produces UVs OUTSIDE the Gregory patch's
        // [0, 1]² parameter domain (e.g. (2.25, 2.22)). The Gregory
        // surface evaluated at these invalid UVs extrapolates the
        // bidegree-9 polynomial, returning points far from the
        // analytic boundary — `cap_outer_g0` stays high.
        //
        // Bypass: write an analytic `Linear` pcurve on the cap side,
        // determined by the segment index (parity tells us whether
        // we're on the row 0 or col 9 boundary):
        //
        // * Even seg_idx (FilletEnd[k]) → cap row 0 (u=0): pcurve
        //   `Linear (0,0) → (0,1)` with parameter range matching
        //   the edge curve's `[t0, t1]`. Reverse if the edge orient-
        //   ation flipped.
        // * Odd seg_idx (FaceEdge[k]) → cap col 9 (v=1): pcurve
        //   `Linear (0,1) → (1,1)`.
        //
        // The cylinder/plane side keeps the NURBS refit — it's an
        // analytic surface where `project_point_to_uv` is reliable,
        // so the existing pcurve refit produces valid UVs there.
        let analytic_cap_pcurve = if stage4_wide_tol {
            cap_pcurve_for_segment(seg_idx, reversed, t0, t1)
        } else {
            None
        };
        if let Some(pc) = analytic_cap_pcurve {
            brep.half_edges[he_handle].pcurve = Some(pc);
        } else if let Some(cap_surf) = cap_surf.as_ref()
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
/// Phase 6c-D Stage 4 — analytic `Linear` cap-side pcurve for an
/// edge wired with a sphere-fit-substituted boundary segment.
///
/// Segment index parity tells us which cap boundary the edge is
/// on (per the layout `sphere_fit_boundary::apply_sphere_fit_substitution`
/// emits — `segments[2k] = FilletEnd[k]`,
/// `segments[2k+1] = FaceEdge[k]`):
///
/// * `seg_idx % 2 == 0` (FilletEnd) → cap **row 0** (`u = 0`):
///   pcurve traces `(0, v)` from `v = 0` to `v = 1`.
/// * `seg_idx % 2 == 1` (FaceEdge) → cap **col 9** (`v = 1`):
///   pcurve traces `(u, 1)` from `u = 0` to `u = 1`.
///
/// `reversed` flips the start/end UVs so the pcurve advances in
/// the same direction as the underlying 3D curve. Parameter
/// range `[t_start, t_end]` matches the edge curve so downstream
/// pcurve evaluation at any `t` in that range maps to a valid
/// Gregory `[0, 1]² UV.
fn cap_pcurve_for_segment(
    seg_idx: usize,
    reversed: bool,
    t_start: f64,
    t_end: f64,
) -> Option<Pcurve> {
    let (start, end) = if seg_idx.is_multiple_of(2) {
        // Row 0: u = 0, v varies 0 → 1.
        ((0.0, 0.0), (0.0, 1.0))
    } else {
        // Col 9: v = 1, u varies 0 → 1.
        ((0.0, 1.0), (1.0, 1.0))
    };
    let (start, end) = if reversed { (end, start) } else { (start, end) };
    Some(Pcurve::linear(start, end, t_start, t_end))
}

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
