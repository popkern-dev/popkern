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
use popkern_core::geometry::{Point3, fitting};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

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
    let tol = 0.05;

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

    let mut wired = 0usize;
    let mut consumed = vec![false; segment_curves.len()];

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
            let edge_h = he.edge;
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
            let Some((seg_idx, d, _reversed)) = best else {
                continue;
            };
            if d > tol {
                continue;
            }

            let curve = segment_curves[seg_idx].clone();
            let (t0, t1) = curve.parameter_range();

            // Capture face surfaces for pcurve refit (this side
            // and twin) before we mutate the edge.
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
                    geometry: popkern_core::geometry::CurveGeometry::Nurbs(curve.clone()),
                    t_start: t0,
                    t_end: t1,
                });
                consumed[seg_idx] = true;
                wired += 1;
            }

            // Re-fit both half-edge pcurves from the new 3-D curve.
            // Without this, the cap-outer-boundary metric (which
            // reads pcurves directly) would still see the stale
            // pcurves on each side and report the original drift.
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
    }
    wired
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
