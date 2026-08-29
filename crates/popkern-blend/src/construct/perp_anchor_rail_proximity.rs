//! Trim-rail proximity helpers for the anchor-split stage anchor placement
//! (PK-185/PK-186 reflex hole-corner gate). Split from
//! `perp_anchor_positions.rs` (file-size limit).

use popkern_core::geometry::{NurbsCurve, Point3};

/// Minimum distance from `point` to the sampled trim curve.
pub(super) fn distance_to_curve(curve: &NurbsCurve, point: &Point3) -> f64 {
    use popkern_core::geometry::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let mut best = f64::INFINITY;
    for k in 0..=64 {
        let t = t0 + (t1 - t0) * f64::from(k) / 64.0;
        best = best.min(curve.point_at(t).distance_to(point));
    }
    best
}

/// Closest point on the trim curve to `point` (coarse sample plus
/// ternary refinement — the rail is smooth and unimodal near the
/// foot for the class this serves).
pub(super) fn closest_point_on_curve(curve: &NurbsCurve, point: &Point3) -> Point3 {
    use popkern_core::geometry::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let mut bt = t0;
    let mut bd = f64::INFINITY;
    for k in 0..=64 {
        let t = t0 + (t1 - t0) * f64::from(k) / 64.0;
        let d = curve.point_at(t).distance_to(point);
        if d < bd {
            bd = d;
            bt = t;
        }
    }
    let step = (t1 - t0) / 64.0;
    let (mut lo, mut hi) = ((bt - step).max(t0), (bt + step).min(t1));
    for _ in 0..40 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(point) < curve.point_at(m2).distance_to(point) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    curve.point_at(0.5 * (lo + hi))
}
