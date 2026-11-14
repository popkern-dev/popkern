//! Exact boxes of the trigonometric analytic curves: circles, arcs and
//! ellipses.
//!
//! Each coordinate of `c + a cos t + b sin t` is `c_i + A_i cos(t - φ_i)`
//! (`A_i = |(a_i, b_i)|`, `φ_i = atan2(b_i, a_i)`), extreme at the
//! interval's ends and at `t = φ_i + kπ` inside it. A box sampled at
//! seventeen parameters misses an extreme between samples by up to the
//! chord's sagitta: a rim circle of radius √3 starting 0.9 rad from its
//! x extreme read 0.0114 short of it, and a boolean's hole filter took
//! the box as containing and left a record piece inside the hole
//! (PK-1182).

use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// The box of `center + a cos t + b sin t` over `t_range`.
pub(super) fn trig_curve_box(center: Point3, a: Vec3, b: Vec3, t_range: (f64, f64)) -> Aabb {
    use std::f64::consts::PI;
    let (lo, hi) = (t_range.0.min(t_range.1), t_range.0.max(t_range.1));
    let at = |t: f64| center + a * t.pcos() + b * t.psin();
    let mut points = vec![at(lo), at(hi)];
    for (ai, bi) in [(a.x, b.x), (a.y, b.y), (a.z, b.z)] {
        if ai == 0.0 && bi == 0.0 {
            continue;
        }
        let phi = bi.patan2(ai);
        let first = ((lo - phi) / PI).ceil();
        let last = ((hi - phi) / PI).floor();
        if !(first.is_finite() && last.is_finite()) {
            continue;
        }
        // Both extremes are reached once the interval spans a period.
        let count = ((last - first) as i64 + 1).clamp(0, 2);
        for k in 0..count {
            points.push(at(phi + (first + k as f64) * PI));
        }
    }
    Aabb::from_points(&points).unwrap_or(Aabb {
        min: points[0],
        max: points[0],
    })
}

#[cfg(test)]
#[path = "analytic_curve_box_tests.rs"]
mod tests;
