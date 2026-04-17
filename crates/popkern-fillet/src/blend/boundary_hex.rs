//! §5.3 Round 8 C1 hexagonal-boundary helpers.
//!
//! Extracted from `boundary.rs` to stay under the 500-line cap.
//! These functions drive the face-edge-bridge emission path in
//! `compute_blend_boundary` — they walk the incident planar faces
//! at a multi-fillet vertex, pair consecutive arc endpoints, and
//! construct straight-line bridge segments on the shared cube
//! face to form a closed hexagonal boundary.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
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

/// §5.3 C4b.2c — fit a Hermite cubic Bezier arc on a parametric
/// surface between two UV points, using the surface's own
/// tangent vectors at the endpoints + a classical-circular-arc
/// magnitude correction. Locks the endpoint tangent direction
/// to the strip surface so the endpoint κ match is ~1e-3 (vs
/// the ~0.4 error of a pure-interpolation cubic).
///
/// Construction:
/// * `P_0 = surface(u1, v1)`, `P_3 = surface(u2, v2)`
/// * `arc_length = Σ chord distances` over 16 dense samples
///   between `(u1, v1)` and `(u2, v2)`
/// * `τ_start = ∂T/∂(arc-param) at (u1, v1) = T_u'·du + T_v'·dv`
/// * `τ_end   = ∂T/∂(arc-param) at (u2, v2)` (same formula)
/// * Classical correction `(4·tan(θ/4))/θ` where `θ` is estimated
///   from the chord-vs-tangent angle: `θ = 2·acos(|t̂·ĉ|)`
/// * `P_1 = P_0 + τ_start · classical_correction · arc_length /
///          (3·|τ_start|)`
/// * `P_2 = P_3 − τ_end · (same scale)`
pub(super) fn fit_hermite_cubic_arc_on_surface(
    surface: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    du: f64,
) -> Result<NurbsCurve, String> {
    let dense_n = 16;
    let mut arc_length = 0.0;
    let mut prev = surface.point_at(u1, v1);
    for i in 1..=dense_n {
        let t = (i as f64) / (dense_n as f64);
        let u = u1 + t * du;
        let v = v1 + t * (v2 - v1);
        let p = surface.point_at(u, v);
        arc_length += p.distance_to(&prev);
        prev = p;
    }

    let p0 = surface.point_at(u1, v1);
    let p3 = surface.point_at(u2, v2);
    let dv = v2 - v1;
    let t_start = {
        let t_u = surface.partial_u(u1, v1);
        let t_v = surface.partial_v(u1, v1);
        Vec3::new(
            t_u.x * du + t_v.x * dv,
            t_u.y * du + t_v.y * dv,
            t_u.z * du + t_v.z * dv,
        )
    };
    let t_end = {
        let t_u = surface.partial_u(u2, v2);
        let t_v = surface.partial_v(u2, v2);
        Vec3::new(
            t_u.x * du + t_v.x * dv,
            t_u.y * du + t_v.y * dv,
            t_u.z * du + t_v.z * dv,
        )
    };

    let chord_vec = Vec3::new(p3.x - p0.x, p3.y - p0.y, p3.z - p0.z);
    let chord_len = chord_vec.length();
    let classical_correction = if chord_len > 1e-12 {
        let chord_unit = Vec3::new(
            chord_vec.x / chord_len,
            chord_vec.y / chord_len,
            chord_vec.z / chord_len,
        );
        let t0_mag = t_start.length();
        if t0_mag > 1e-12 {
            let t0_unit = Vec3::new(t_start.x / t0_mag, t_start.y / t0_mag, t_start.z / t0_mag);
            let dot = t0_unit.dot(&chord_unit).abs().clamp(0.0, 1.0);
            let theta = 2.0 * dot.acos();
            if theta > 1e-6 {
                (4.0 * (theta / 4.0).tan()) / theta
            } else {
                1.0
            }
        } else {
            1.0
        }
    } else {
        1.0
    };

    let hermite_cp = |origin: Point3, tangent: Vec3, sign: f64| -> Point3 {
        let mag = tangent.length();
        if mag < 1e-15 {
            return origin;
        }
        let scale = sign * classical_correction * arc_length / (3.0 * mag);
        Point3::new(
            origin.x + scale * tangent.x,
            origin.y + scale * tangent.y,
            origin.z + scale * tangent.z,
        )
    };
    let p1 = hermite_cp(p0, t_start, 1.0);
    let p2 = hermite_cp(p3, t_end, -1.0);
    cubic_nurbs_from_control_points([p0, p1, p2, p3])
}

/// Build a cubic NURBS curve from 4 control points directly —
/// no interpolation, no solve. Equivalent to wrapping the
/// provided 4 Bezier CPs with a clamped uniform knot vector
/// and unit weights. Used by Hermite cubic-Bezier constructions
/// that already computed the exact CP positions.
pub(super) fn cubic_nurbs_from_control_points(pts: [Point3; 4]) -> Result<NurbsCurve, String> {
    Ok(NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: pts.to_vec(),
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
