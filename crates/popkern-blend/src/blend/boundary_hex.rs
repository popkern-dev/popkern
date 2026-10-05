//! §5.3 hexagonal-boundary helpers + §5.3 item #6
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

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
use popkern_math::PortableMath;

use super::boundary::BoundarySegment;
#[path = "boundary_hex_ordering.rs"]
mod boundary_hex_ordering;
use boundary_hex_ordering::{angular_arc_ordering, ring_bridge_score, successor_walk_ordering};

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

    // Valence >= 4: NEITHER candidate ordering is universally
    // trustworthy. The successor walk mispairs arcs whose endpoints
    // crowd around the apex (measured on the pentagonal pyramid: a
    // star polygon whose bridges cross the apex axis), while the
    // angular sort folds when the arc layout is asymmetric (measured
    // on pyramid4_mixed: order [1,0,2,3] paired DIAGONAL fillets, the
    // FaceEdge bridges grabbed a diagonal face, every sphere-fit read
    // amplitude 0.56 where adjacent faces give 1.0, and the whole
    // apex silently fell to the Open rung — no cap, no seams).
    //
    // Adjacency has a measurable definition: consecutive arc ends
    // must be bridgeable ON A COMMON SUPPORT FACE. Build both
    // candidates and keep the one with fewer unsupported bridges
    // (ties: smaller total bridge length). The walk stays the sole
    // path for n <= 3, where the cube hex-closure fixtures validate
    // it and wrong-end pairing is benign.
    if n >= 4 {
        let walk =
            successor_walk_ordering(fillet_arcs, fillet_endpoints, face_supports, tol_linear);
        let angular = angular_arc_ordering(fillet_arcs, fillet_endpoints);
        return match (walk, angular) {
            (Some(w), Some(a)) => {
                let score_w = ring_bridge_score(&w, fillet_endpoints, face_supports, tol_slack);
                let score_a = ring_bridge_score(&a, fillet_endpoints, face_supports, tol_slack);
                if popkern_core::trace::trace_enabled("boundary-hex") {
                    eprintln!(
                        "[boundary-hex] n={n} walk order={:?} score={score_w:?}; \
                         angular order={:?} score={score_a:?}",
                        w.0, a.0
                    );
                }
                if score_a < score_w { Some(a) } else { Some(w) }
            }
            (walk, angular) => walk.or(angular),
        };
    }

    successor_walk_ordering(fillet_arcs, fillet_endpoints, face_supports, tol_linear)
}

/// Reverse an arc's parameter direction (exact: control points,
/// weights, and reflected knots).
pub(crate) fn reverse_arc(curve: &NurbsCurve) -> NurbsCurve {
    boundary_hex_ordering::reverse_nurbs_curve(curve)
}

/// Exact chart-boundary iso between two surface points (PK-189).
///
/// A swept strip's terminal section is the loft's boundary row —
/// exact for BOTH circular and conic profiles — while the Hermite
/// cubic fit below is an approximation that misses a conic section
/// by the shoulder gap (measured 1.1e-2 on the rho = 0.35 cube
/// corner). When the two projected end points ARE the endpoints of
/// one of the four chart-boundary isos of a NURBS surface, return
/// that iso verbatim (oriented start → end); anything else falls
/// through to the Hermite fit.
pub(super) fn exact_boundary_iso_between(
    surface: &popkern_core::geometry::SurfaceGeometry,
    start: &Point3,
    end: &Point3,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::CurveEval as _;
    let popkern_core::geometry::SurfaceGeometry::Nurbs(n) = surface else {
        return None;
    };
    let strided = |offset: usize, stride: usize, count: usize, degree: usize, knots: &[f64]| {
        let mut control_points = Vec::with_capacity(count);
        let mut weights = Vec::with_capacity(count);
        for i in 0..count {
            let idx = offset + i * stride;
            if idx >= n.control_points.len() {
                return None;
            }
            control_points.push(n.control_points[idx]);
            weights.push(n.weights[idx]);
        }
        Some(NurbsCurve {
            degree,
            knots: knots.to_vec(),
            control_points,
            weights,
        })
    };
    // Control-point layout: index = u * num_v + v.
    let candidates = [
        // v = first / last (rows across u)
        strided(0, n.num_v, n.num_u, n.degree_u, &n.knots_u),
        strided(n.num_v - 1, n.num_v, n.num_u, n.degree_u, &n.knots_u),
        // u = first / last (columns across v)
        strided(0, 1, n.num_v, n.degree_v, &n.knots_v),
        strided((n.num_u - 1) * n.num_v, 1, n.num_v, n.degree_v, &n.knots_v),
    ];
    let scale = start.distance_to(end).max(1e-9);
    let gate = scale * 1.0e-6;
    for candidate in candidates.into_iter().flatten() {
        let (t0, t1) = candidate.parameter_range();
        let (c0, c1) = (candidate.point_at(t0), candidate.point_at(t1));
        if c0.distance_to(start) <= gate && c1.distance_to(end) <= gate {
            return Some(candidate);
        }
        if c0.distance_to(end) <= gate && c1.distance_to(start) <= gate {
            return Some(NurbsCurve {
                degree: candidate.degree,
                knots: reversed_clamped_knots(&candidate.knots),
                control_points: candidate.control_points.iter().rev().copied().collect(),
                weights: candidate.weights.iter().rev().copied().collect(),
            });
        }
    }
    None
}

/// Reverse a clamped knot vector onto the same parameter range.
pub(super) fn reversed_clamped_knots(knots: &[f64]) -> Vec<f64> {
    let (lo, hi) = (knots[0], knots[knots.len() - 1]);
    knots.iter().rev().map(|k| lo + hi - k).collect()
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
            let theta = 2.0 * dot.pacos();
            if theta > 1e-6 {
                (4.0 * (theta / 4.0).ptan()) / theta
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

/// Collect the set of face handles incident to `vertex` by walking
/// all half-edges originating from it.
pub(super) fn faces_at_vertex_set(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut seen: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let Some(loop_data) = brep.loops.get(he.loop_) else {
            continue;
        };
        let face = loop_data.face;
        if !seen.contains(&face) {
            seen.push(face);
        }
    }
    seen
}

/// Golden-section iterations of the closest-point polish: the bracket
/// shrinks by 0.618 per step, so 80 steps reach roundoff from any
/// sample bracket.
const CLOSEST_POINT_POLISH_ITERATIONS: usize = 80;

/// The parameter of the closest point on `curve` to `target` (PK-548).
/// It takes the best of `n_samples + 1` uniform samples, narrows it by
/// golden-section search over the two sample intervals around it, and
/// finishes with Newton on the stationarity condition, to roundoff.
///
/// The ring station push-back measures arclength from this station on
/// every rail. Unpolished, each station was quantized to its own sample
/// spacing: 128 or 512 samples over trims about 3 units long. On the
/// pyramid4 apex, the one rail stored the other way round landed its
/// pushed station 0.0066 off its three symmetric neighbours'. That made
/// the blend boundary concave at a ring corner (186°), where no quad
/// chart is regular.
///
/// The default end arcs (`extract_fillet_end_near_vertex`) keep their
/// coarse search on purpose. On an orthogonal corner the exact closest
/// station puts the arc in the neighbouring face's plane, a degenerate
/// configuration the sampled station stays clear of. The push-back adds
/// a clearance margin, so only it takes the exact station.
///
/// The polish assumes the distance is unimodal within the bracket, which
/// holds for rails that pass the vertex once at sample spacing.
pub(crate) fn closest_param_on_nurbs_curve(
    curve: &NurbsCurve,
    target: &Point3,
    n_samples: usize,
) -> f64 {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let samples = n_samples.max(1);
    let at = |i: usize| t0 + (t1 - t0) * (i as f64) / (samples as f64);
    let mut best = (f64::INFINITY, 0usize);
    for i in 0..=samples {
        let d = curve.point_at(at(i)).distance_to(target);
        if d < best.0 {
            best = (d, i);
        }
    }
    let i = best.1;
    let (mut lo, mut hi) = (at(i.saturating_sub(1)), at((i + 1).min(samples)));
    let (bracket_lo, bracket_hi) = (lo.min(hi), lo.max(hi));
    let distance = |t: f64| curve.point_at(t).distance_to(target);
    let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut a, mut b) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
    let (mut fa, mut fb) = (distance(a), distance(b));
    for _ in 0..CLOSEST_POINT_POLISH_ITERATIONS {
        if fa <= fb {
            hi = b;
            (b, fb) = (a, fa);
            a = hi - ratio * (hi - lo);
            fa = distance(a);
        } else {
            lo = a;
            (a, fa) = (b, fb);
            b = lo + ratio * (hi - lo);
            fb = distance(b);
        }
    }
    // The distance is flat at its minimum, so the golden section stops
    // near √ε. Newton on the stationarity condition
    // `g(t) = (C(t) − P)·C′(t) = 0` finishes to roundoff, clamped to the
    // sample bracket.
    let mut t = 0.5 * (lo + hi);
    for _ in 0..NEWTON_POLISH_ITERATIONS {
        let (p, d1, d2) = (
            curve.point_at(t),
            curve.derivative_at(t, 1),
            curve.derivative_at(t, 2),
        );
        let offset = p - *target;
        let g = offset.dot(&d1);
        let dg = d1.dot(&d1) + offset.dot(&d2);
        if dg.is_nan() || dg <= 0.0 {
            break;
        }
        let next = (t - g / dg).clamp(bracket_lo, bracket_hi);
        if (next - t).abs() <= f64::EPSILON * (1.0 + t.abs()) {
            t = next;
            break;
        }
        t = next;
    }
    if distance(t) <= best.0 { t } else { at(i) }
}

/// Newton iterations finishing the closest-point polish.
const NEWTON_POLISH_ITERATIONS: usize = 8;
