//! Station queries on center curves (PK-139): closest-point lookup
//! used by the per-station variable-radius center march.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};

/// Closest point on a NURBS curve to a target: dense sampling with
/// two rounds of local ternary refinement (analytic offset
/// intersections here are lines/arcs/low-degree fits — 64 samples
/// bracket the minimum comfortably).
pub(crate) fn closest_point_on_nurbs(curve: &NurbsCurve, target: &Point3) -> Point3 {
    let (t0, t1) = curve.parameter_range();
    let mut best_t = 0.5 * (t0 + t1);
    let mut best_d = f64::MAX;
    let n = 64usize;
    for i in 0..=n {
        let t = t0 + (t1 - t0) * i as f64 / n as f64;
        let d = curve.point_at(t).distance_to(target);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    let mut lo = (best_t - (t1 - t0) / n as f64).max(t0);
    let mut hi = (best_t + (t1 - t0) / n as f64).min(t1);
    for _ in 0..40 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(target) < curve.point_at(m2).distance_to(target) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    curve.point_at(0.5 * (lo + hi))
}

/// Select the intersection curve whose midpoint is nearest to the given point.
pub(crate) fn select_nearest_curve(
    intersections: &[popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection],
    target: &Point3,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;

    let mut best: Option<(f64, NurbsCurve)> = None;

    for isect in intersections {
        if let SurfaceSurfaceIntersection::Curve(cos) = isect {
            let curve = &cos.curve_3d;
            let (t0, t1) = curve.parameter_range();
            let mid_pt = curve.point_at((t0 + t1) * 0.5);
            let dist = mid_pt.distance_to(target);
            let is_closer = best.as_ref().is_none_or(|(d, _)| dist < *d);
            if is_closer {
                best = Some((dist, curve.clone()));
            }
        }
    }

    best.map(|(_, c)| c)
        .ok_or_else(|| "no curve intersections found in offset result".into())
}
