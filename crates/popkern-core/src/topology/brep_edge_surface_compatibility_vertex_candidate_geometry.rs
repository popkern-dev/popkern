//! Geometry helpers for edge/surface vertex candidate classification.

use crate::geometry::curve::CurveEval;
use crate::geometry::{NurbsCurve, Point3};

pub(super) fn closest_boundary_parameter(
    curve: &NurbsCurve,
    target: Point3,
    t_start: f64,
    t_end: f64,
    samples: usize,
) -> Option<f64> {
    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for index in 0..samples {
        let f = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * f;
        let distance = curve.point_at(t).distance_to(&target);
        if distance < best_distance {
            best_distance = distance;
            best_index = index;
        }
    }
    if !best_distance.is_finite() {
        return None;
    }
    let step = (t_end - t_start) / (samples - 1) as f64;
    let mut lo = t_start + step * best_index.saturating_sub(1) as f64;
    let mut hi = t_start + step * (best_index + 1).min(samples - 1) as f64;
    for _ in 0..32 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(&target) <= curve.point_at(m2).distance_to(&target) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some(0.5 * (lo + hi))
}

pub(super) fn point_centroid(points: &[Point3]) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for point in points {
        x += point.x;
        y += point.y;
        z += point.z;
    }
    let scale = 1.0 / points.len() as f64;
    Point3::new(x * scale, y * scale, z * scale)
}
