//! Footprint projection for centre-curve stations: the closest point
//! of a support surface to a ball centre, inside the face's parameter
//! window.
//!
//! Planes project exactly. Every other surface asks the surface's own
//! `closest_point` (exact on the analytic classes, multi-start
//! Gauss-Newton on NURBS), carries periodic coordinates into the
//! face window by whole periods and clamps aperiodic ones; a 20 × 20
//! grid search over the window is the fallback when the projection
//! declines (a point on a cylinder's axis). The grid used to be the
//! only arm: its nodes sit 0.33 rad apart on a full cylinder chart,
//! and a 0.25-rad arc piece of a seam-straddling boss chain had every
//! footprint snapped to the seam node — the offset-intersect sweep
//! twisted (PK-421).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Closest point of `surface` to `point` within the face window
/// `u_range × v_range`.
pub fn project_onto_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    point: &Point3,
) -> Result<Point3, String> {
    if let SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    } = surface
    {
        let d = *point - *origin;
        let projected = *point - *normal * d.dot(normal);
        let dp = projected - *origin;
        let u = dp.dot(u_axis).clamp(u_range.0, u_range.1);
        let v = dp.dot(v_axis).clamp(v_range.0, v_range.1);
        return Ok(*origin + *u_axis * u + *v_axis * v);
    }
    if let Some((u, v)) = surface.closest_point(point) {
        let (period_u, period_v) = surface.parameter_periodicity();
        let u = wrap_then_clamp(u, u_range, period_u);
        let v = wrap_then_clamp(v, v_range, period_v);
        return Ok(surface.point_at(u, v));
    }
    Ok(grid_nearest(surface, u_range, v_range, point))
}

/// Coarse 20 × 20 grid search over the window.
fn grid_nearest(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    point: &Point3,
) -> Point3 {
    let n = 20usize;
    let mut best = surface.point_at((u_range.0 + u_range.1) * 0.5, (v_range.0 + v_range.1) * 0.5);
    let mut best_dist = best.distance_to(point);
    for iu in 0..n {
        let u = u_range.0 + (u_range.1 - u_range.0) * iu as f64 / (n - 1) as f64;
        for iv in 0..n {
            let v = v_range.0 + (v_range.1 - v_range.0) * iv as f64 / (n - 1) as f64;
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best = sp;
            }
        }
    }
    best
}

/// `x` carried into `[lo, hi]` by whole periods of a periodic chart
/// axis, then clamped; aperiodic axes clamp directly. A projection
/// answers in the principal branch (u = −0.25 for a station just
/// behind the seam of a [0, 2π] chart); clamping that branch parks
/// the station ON the seam, a quarter radian away.
pub(crate) fn wrap_then_clamp(x: f64, (lo, hi): (f64, f64), period: Option<f64>) -> f64 {
    let mut x = x;
    if let Some(p) = period
        && p > 0.0
        && lo.is_finite()
        && hi.is_finite()
    {
        if x < lo {
            x += p * ((lo - x) / p).ceil();
        }
        if x > hi {
            x -= p * ((x - hi) / p).ceil();
        }
    }
    x.clamp(lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use std::f64::consts::TAU;

    #[test]
    fn a_station_behind_the_seam_projects_onto_its_own_azimuth() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        // Ball centre at azimuth −0.115 rad, radius 1.05: the
        // footprint is at the same azimuth on the unit cylinder.
        let c = Point3::new(1.05 * (-0.115f64).cos(), 1.05 * (-0.115f64).sin(), 0.7);
        let p = project_onto_surface(&cyl, (0.0, TAU), (-1.0, 3.0), &c).expect("projects");
        let expected = Point3::new((-0.115f64).cos(), (-0.115f64).sin(), 0.7);
        assert!(p.distance_to(&expected) < 1e-8, "{p:?} vs {expected:?}");
        // A window that excludes the azimuth clamps to its nearer end
        // in the periodic sense: [0.5, 1.0] keeps u = 0.5 (0.615 rad
        // away; 1.0 is 1.115 rad away the other way round).
        let p = project_onto_surface(&cyl, (0.5, 1.0), (-1.0, 3.0), &c).expect("projects");
        assert!(p.distance_to(&Point3::new(0.5f64.cos(), 0.5f64.sin(), 0.7)) < 1e-8);
    }

    #[test]
    fn wrap_then_clamp_carries_by_whole_periods() {
        assert!((wrap_then_clamp(-0.25, (0.0, TAU), Some(TAU)) - (TAU - 0.25)).abs() < 1e-12);
        assert!((wrap_then_clamp(TAU + 0.3, (0.0, TAU), Some(TAU)) - 0.3).abs() < 1e-12);
        assert_eq!(wrap_then_clamp(-0.25, (0.0, 1.0), None), 0.0);
        assert_eq!(wrap_then_clamp(0.4, (0.0, 1.0), Some(TAU)), 0.4);
    }
}
