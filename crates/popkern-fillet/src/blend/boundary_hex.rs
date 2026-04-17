//! §5.3 Round 8 C1 hexagonal-boundary helpers.
//!
//! Extracted from `boundary.rs` to stay under the 500-line cap.
//! These functions drive the face-edge-bridge emission path in
//! `compute_blend_boundary` — they walk the incident planar faces
//! at a multi-fillet vertex, pair consecutive arc endpoints, and
//! construct straight-line bridge segments on the shared cube
//! face to form a closed hexagonal boundary.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::FaceHandle;

use super::boundary::BoundarySegment;

/// Reorient fillet arcs so consecutive arcs share a cube face via
/// a connecting endpoint pair, and the bridge from `arc[k].end`
/// to `arc[k+1].start` lies on that shared face.
///
/// Returns `(order, flip, oriented_arcs)`:
/// * `order[k]` — which input fillet index goes at position `k`
/// * `flip[k]` — whether arc `order[k]` was reversed to align
/// * `oriented_arcs[k]` — the arc as emitted (reversed if flipped)
///
/// Algorithm: fix arc 0 in the input orientation, then greedily
/// attach the next arc whose start matches the current arc's end
/// on some shared planar face (reversing if needed to match).
/// Returns `None` if no valid chain can be built (e.g., arcs
/// don't share faces — fall back to FilletEnd-only boundary).
#[allow(clippy::type_complexity)]
pub(super) fn reorient_arcs_for_hex_closure(
    fillet_indices: &[usize],
    fillet_arcs: &[NurbsCurve],
    fillet_endpoints: &[[Point3; 2]],
    planar_faces: &[(FaceHandle, Point3, Vec3)],
    tol_linear: f64,
) -> Option<(Vec<usize>, Vec<bool>, Vec<NurbsCurve>)> {
    let n = fillet_indices.len();
    if n < 2 {
        return None;
    }
    let mut order = vec![0usize; n];
    let mut flip = vec![false; n];
    let mut used = vec![false; n];
    order[0] = 0;
    used[0] = true;

    for k in 1..n {
        let prev_idx = order[k - 1];
        let prev_endpoints = if flip[k - 1] {
            [fillet_endpoints[prev_idx][1], fillet_endpoints[prev_idx][0]]
        } else {
            fillet_endpoints[prev_idx]
        };
        let prev_end = prev_endpoints[1];

        let mut best: Option<(usize, bool)> = None;
        'search: for (cand_idx, ep) in fillet_endpoints.iter().enumerate() {
            if used[cand_idx] {
                continue;
            }
            for cand_flip in [false, true] {
                let cand_start = if cand_flip { ep[1] } else { ep[0] };
                for (_, origin, normal) in planar_faces {
                    let d1 = point_plane_distance(&prev_end, origin, normal);
                    let d2 = point_plane_distance(&cand_start, origin, normal);
                    if d1 <= tol_linear && d2 <= tol_linear {
                        best = Some((cand_idx, cand_flip));
                        break 'search;
                    }
                }
            }
        }
        let (next_idx, next_flip) = best?;
        order[k] = next_idx;
        flip[k] = next_flip;
        used[next_idx] = true;
    }

    let oriented_arcs: Vec<NurbsCurve> = (0..n)
        .map(|k| {
            let arc = &fillet_arcs[order[k]];
            if flip[k] {
                reverse_nurbs_curve(arc)
            } else {
                arc.clone()
            }
        })
        .collect();

    Some((order, flip, oriented_arcs))
}

/// Unsigned distance from `p` to the plane with `origin` and
/// unit-length `normal`.
pub(super) fn point_plane_distance(p: &Point3, origin: &Point3, normal: &Vec3) -> f64 {
    let d = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    d.dot(normal).abs()
}

/// Reverse a NURBS curve's parameter direction by reversing the
/// control-point and weight order and reflecting the knot vector.
fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    let (t0, t1) = curve.parameter_range();
    let mut new_knots: Vec<f64> = curve.knots.iter().rev().map(|&k| t0 + (t1 - k)).collect();
    for k in &mut new_knots {
        *k = k.clamp(t0, t1);
    }
    let new_cps: Vec<Point3> = curve.control_points.iter().rev().copied().collect();
    let new_weights: Vec<f64> = curve.weights.iter().rev().copied().collect();
    new_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    NurbsCurve {
        degree: curve.degree,
        knots: new_knots,
        control_points: new_cps,
        weights: new_weights,
    }
}

/// Fit a cubic NURBS curve through four sample points at
/// `t ∈ {0, 1/3, 2/3, 1}`. Closed-form two-equation solve for
/// interior Bezier CPs given boundary samples — same math as
/// `blend::gregory::fit_cubic_bezier_through` but starting from
/// Point3 samples rather than a NurbsCurve.
pub(super) fn fit_cubic_nurbs_through_4_points(pts: [Point3; 4]) -> Result<NurbsCurve, String> {
    let p0 = pts[0];
    let p3 = pts[3];
    let q1 = pts[1];
    let q2 = pts[2];
    let solve = |q1c: f64, q2c: f64, p0c: f64, p3c: f64| -> (f64, f64) {
        let a = 27.0 * q1c - 8.0 * p0c - p3c;
        let b = 27.0 * q2c - p0c - 8.0 * p3c;
        ((2.0 * a - b) / 18.0, (2.0 * b - a) / 18.0)
    };
    let (p1x, p2x) = solve(q1.x, q2.x, p0.x, p3.x);
    let (p1y, p2y) = solve(q1.y, q2.y, p0.y, p3.y);
    let (p1z, p2z) = solve(q1.z, q2.z, p0.z, p3.z);
    Ok(NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            p0,
            Point3::new(p1x, p1y, p1z),
            Point3::new(p2x, p2y, p2z),
            p3,
        ],
        weights: vec![1.0; 4],
    })
}

/// Build a face-edge bridge curve between two 3D points, picking
/// whichever planar face contains both within `tol_linear`. The
/// bridge is a straight line sampled through a midpoint and fitted
/// as a degree-1 NURBS.
pub(super) fn build_face_edge_bridge_between(
    pa: &Point3,
    pb: &Point3,
    planar_faces: &[(FaceHandle, Point3, Vec3)],
    tol_linear: f64,
) -> Option<BoundarySegment> {
    let mut best: Option<(FaceHandle, f64)> = None;
    for (face, origin, normal) in planar_faces {
        let da = point_plane_distance(pa, origin, normal);
        let db = point_plane_distance(pb, origin, normal);
        let worst = da.max(db);
        if worst <= tol_linear {
            match &best {
                Some((_, cur)) if *cur <= worst => {}
                _ => best = Some((*face, worst)),
            }
        }
    }
    let (face, _) = best?;
    let mid = Point3::new(
        (pa.x + pb.x) * 0.5,
        (pa.y + pb.y) * 0.5,
        (pa.z + pb.z) * 0.5,
    );
    let curve = crate::center_curve::fit_nurbs_through_points(&[*pa, mid, *pb]).ok()?;
    Some(BoundarySegment::FaceEdge { face, curve })
}
