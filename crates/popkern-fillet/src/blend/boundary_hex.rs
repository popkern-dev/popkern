//! §5.3 Round 8 C1 hexagonal-boundary helpers + §5.3 item #6
//! curved-support generalization.
//!
//! Extracted from `boundary.rs` to stay under the 500-line cap.
//! These functions drive the face-edge-bridge emission path in
//! `compute_blend_boundary` — they walk the incident faces at a
//! multi-fillet vertex, pair consecutive arc endpoints, and
//! construct bridge segments on the shared face (via
//! [`surface_bridge`](super::surface_bridge)) to form a closed
//! hexagonal / N-gon boundary.
//!
//! **Item #6 upgrade:** the planar-only restriction is lifted —
//! any face whose surface is well-defined (Plane, Cylinder,
//! Sphere, Cone, Torus, NURBS, Gregory) can host a bridge. The
//! point-on-face distance check generalises from plane-normal
//! projection to `point_on_surface_distance`.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::FaceHandle;

use super::boundary::BoundarySegment;
use super::surface_bridge::{
    build_straight_line_bridge, build_surface_bridge, point_on_surface_distance, surface_kind_name,
};

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
    face_supports: &[(FaceHandle, SurfaceGeometry)],
    tol_linear: f64,
) -> Option<(Vec<usize>, Vec<bool>, Vec<NurbsCurve>)> {
    let n = fillet_indices.len();
    if n < 2 {
        return None;
    }
    // §5.3 item #5 (follow-up): slack tolerance for mixed-
    // convexity saddles. On a notched cube the convex fillet's
    // arc endpoints sit on cube faces while the concave fillets'
    // endpoints sit on notch faces; boolean subtraction can
    // leave endpoints a few multiples of `tol_linear` off their
    // logical support face due to trim-curve endpoint resampling.
    // Walking with a relaxed 10× band closes these saddle hexes
    // while staying well below the fillet radius so we never
    // accidentally match an endpoint against an unrelated
    // distant face.
    let tol_slack = tol_linear.max(1e-5) * 10.0;
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

        // Pass 1: strict tolerance — well-formed BReps where both
        // endpoints lie exactly on a shared support face.
        // Pass 2: slack tolerance — saddle / mixed-convexity cases
        // where the boolean result leaves trim endpoints a few
        // multiples of `tol_linear` off their logical face.
        // Pass 3: nearest-endpoint fallback — §5.3 item #5 saddle
        // path. When no face-support shares both endpoints within
        // tolerance (a notched cube's trim endpoints may sit on
        // the fillet surface, not the support face), pick the
        // candidate whose start (or end, post-flip) minimises the
        // straight-line distance to `prev_end`. Hex-adjacency is
        // a topological property; geometrically-close endpoints
        // are correct up to a ridge discontinuity the bridge
        // builder handles via its straight-line fallback.
        let mut best: Option<(usize, bool)> = None;
        for &tol_try in &[tol_linear, tol_slack] {
            if best.is_some() {
                break;
            }
            'search: for (cand_idx, ep) in fillet_endpoints.iter().enumerate() {
                if used[cand_idx] {
                    continue;
                }
                for cand_flip in [false, true] {
                    let cand_start = if cand_flip { ep[1] } else { ep[0] };
                    for (_, surface) in face_supports {
                        let d1 = point_on_surface_distance(&prev_end, surface);
                        let d2 = point_on_surface_distance(&cand_start, surface);
                        if d1 <= tol_try && d2 <= tol_try {
                            best = Some((cand_idx, cand_flip));
                            break 'search;
                        }
                    }
                }
            }
        }
        if best.is_none() {
            // §5.3 item #5 saddle fallback: nearest-endpoint
            // matching. On mixed-convexity BReps the trim endpoints
            // can sit on the fillet surface with ~10⁻² offsets from
            // any support face, so the face-adjacency test above
            // rejects every pair. Hex-adjacency is ultimately a
            // topological property — picking the candidate whose
            // start is geometrically closest to `prev_end`
            // preserves the correct cycle order and lets the
            // bridge builder fill in a straight-line segment on
            // the best-effort face.
            let mut nearest: Option<(f64, usize, bool)> = None;
            for (cand_idx, ep) in fillet_endpoints.iter().enumerate() {
                if used[cand_idx] {
                    continue;
                }
                for cand_flip in [false, true] {
                    let cand_start = if cand_flip { ep[1] } else { ep[0] };
                    let d = prev_end.distance_to(&cand_start);
                    if nearest.is_none_or(|(best_d, _, _)| d < best_d) {
                        nearest = Some((d, cand_idx, cand_flip));
                    }
                }
            }
            if let Some((_, idx, fl)) = nearest {
                best = Some((idx, fl));
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
/// whichever supporting face contains both within `tol_linear`.
/// The bridge follows the face's underlying
/// [`SurfaceGeometry`] — straight line on a plane, surface arc on
/// a cylinder / sphere / cone / torus, UV-linear on a NURBS /
/// Gregory. Falls back to a straight-line cubic Bezier if the
/// surface bridge construction fails (e.g., degenerate surface
/// parameterisation).
///
/// **§5.3 item #6:** generalised from planes-only (where bridges
/// were straight-line degree-1 NURBS) to every supported surface
/// variant. The downstream `build_n_quad_decomposition` path
/// re-samples the bridge via `fit_cubic_bezier_through`, so the
/// key requirement is that the bridge (a) lies on the face
/// surface and (b) interpolates `pa` / `pb`.
pub(super) fn build_face_edge_bridge_between(
    pa: &Point3,
    pb: &Point3,
    face_supports: &[(FaceHandle, SurfaceGeometry)],
    tol_linear: f64,
) -> Option<(BoundarySegment, Option<&'static str>)> {
    // Pick the face minimising `max(d(pa), d(pb))`. Strict
    // tolerance gets the surface-hugging cubic bridge; a looser
    // "best effort" pick triggers the straight-line fallback
    // (§5.3 item #5 saddle path where trim endpoints can sit on
    // the fillet surface instead of the adjacent support face,
    // leaving residual ~10⁻² offsets that strict tol would
    // reject entirely).
    let mut best_strict: Option<(FaceHandle, &SurfaceGeometry, f64)> = None;
    let mut best_any: Option<(FaceHandle, &SurfaceGeometry, f64)> = None;
    for (face, surface) in face_supports {
        let da = point_on_surface_distance(pa, surface);
        let db = point_on_surface_distance(pb, surface);
        let worst = da.max(db);
        if worst <= tol_linear && best_strict.is_none_or(|(_, _, cur)| worst < cur) {
            best_strict = Some((*face, surface, worst));
        }
        if best_any.is_none_or(|(_, _, cur)| worst < cur) {
            best_any = Some((*face, surface, worst));
        }
    }

    if let Some((face, surface, _)) = best_strict {
        if let Some(curve) = build_surface_bridge(pa, pb, surface, tol_linear) {
            return Some((BoundarySegment::FaceEdge { face, curve }, None));
        }
        let curve = build_straight_line_bridge(pa, pb)?;
        return Some((
            BoundarySegment::FaceEdge { face, curve },
            Some(surface_kind_name(surface)),
        ));
    }

    // Saddle / mixed-convexity fallback: neither endpoint sits
    // cleanly on any face support within `tol_linear`. Emit a
    // straight-line bridge on the best-effort face so the hex
    // closes with `2N` segments (downstream Gregory fit samples
    // the curve + re-fits cubic, so bridge-off-surface is
    // acceptable here).
    let (face, surface, _) = best_any?;
    let curve = build_straight_line_bridge(pa, pb)?;
    Some((
        BoundarySegment::FaceEdge { face, curve },
        Some(surface_kind_name(surface)),
    ))
}
