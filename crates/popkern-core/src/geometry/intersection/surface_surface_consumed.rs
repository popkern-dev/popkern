//! Seed-consumption test for the surface–surface marching driver:
//! distance from a candidate seed to the polylines already traced.

use crate::geometry::Point3;

/// True when `p` lies within `radius` of any segment of any traced
/// polyline. The point-to-segment metric (rather than point-to-sample)
/// is what makes consumption honest for adaptively-stepped traces —
/// sample spacing may exceed the consumption radius by orders of
/// magnitude while the chord stays within tolerance of the curve.
pub(super) fn point_near_polylines(p: &Point3, polylines: &[Vec<Point3>], radius: f64) -> bool {
    let radius_sq = radius * radius;
    polylines.iter().any(|line| {
        if line.len() == 1 {
            return line[0].distance_to(p) < radius;
        }
        line.windows(2).any(|seg| {
            let a = seg[0];
            let b = seg[1];
            let ab = b - a;
            let len_sq = ab.dot(&ab);
            let t = if len_sq > 0.0 {
                ((*p - a).dot(&ab) / len_sq).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let closest = a + ab * t;
            let d = *p - closest;
            d.dot(&d) < radius_sq
        })
    })
}
