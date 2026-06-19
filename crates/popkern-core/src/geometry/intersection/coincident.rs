//! Coincident region detection for surface-surface intersection.

use crate::geometry::Point3;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Detect coincident regions where two surfaces overlap.
///
/// Samples surface 1 on an 8×8 grid within its parameter range. For each
/// sample point, projects it onto surface 2 via Gauss-Newton closest-point
/// iteration and records the minimum distance. If more than 80% of samples
/// are within linear tolerance the full parameter ranges are reported as a
/// `CoincidentRegion`.
#[allow(dead_code)]
pub(crate) fn detect_coincident_regions(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceSurfaceIntersection> {
    if !parameter_box_is_valid(u1_range, v1_range) || !parameter_box_is_valid(u2_range, v2_range) {
        return Vec::new();
    }

    const GRID: usize = 8;
    let total = GRID * GRID;

    let mut within_count: usize = 0;

    for i in 0..GRID {
        for j in 0..GRID {
            let u1 = u1_range.0 + (u1_range.1 - u1_range.0) * (i as f64 + 0.5) / (GRID as f64);
            let v1 = v1_range.0 + (v1_range.1 - v1_range.0) * (j as f64 + 0.5) / (GRID as f64);

            let p1 = surface1.point_at(u1, v1);

            let dist = closest_point_distance(p1, surface2, u2_range, v2_range);

            if dist < tolerance.linear() {
                within_count += 1;
            }
        }
    }

    // Require > 80% of samples to be within tolerance.
    if within_count * 5 > total * 4 {
        vec![SurfaceSurfaceIntersection::CoincidentRegion {
            u1_range,
            v1_range,
            u2_range,
            v2_range,
        }]
    } else {
        vec![]
    }
}

/// Closest-point distance from a point to a surface, used by the
/// coincident fast-path for quick rejection sampling.
pub(crate) fn closest_point_distance_pub(
    p1: Point3,
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> f64 {
    closest_point_distance(p1, surface2, u2_range, v2_range)
}

/// Find the minimum distance from `p1` to `surface2` by projecting `p1` onto
/// `surface2` via Gauss-Newton iteration, starting from the centre of the
/// parameter domain. Returns the distance at the best found projection.
#[allow(dead_code)]
fn closest_point_distance(
    p1: Point3,
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> f64 {
    if !parameter_box_is_valid(u2_range, v2_range) {
        return f64::INFINITY;
    }

    // Try multiple starting points for robustness: centre plus four corners.
    let u_mid = 0.5 * (u2_range.0 + u2_range.1);
    let v_mid = 0.5 * (v2_range.0 + v2_range.1);

    let starts = [
        (u_mid, v_mid),
        (u2_range.0, v2_range.0),
        (u2_range.0, v2_range.1),
        (u2_range.1, v2_range.0),
        (u2_range.1, v2_range.1),
    ];

    let mut best = f64::INFINITY;

    for (u0, v0) in starts {
        let dist = gauss_newton_closest_point(p1, surface2, u0, v0, u2_range, v2_range);
        if dist < best {
            best = dist;
        }
    }

    best
}

/// Gauss-Newton closest-point projection of `p1` onto `surface2`.
///
/// Minimises `f(u2,v2) = |P1 - S2(u2,v2)|²` using at most 10 iterations.
/// Parameters are clamped to `[u2_range, v2_range]` at each step.
/// Returns the Euclidean distance at the converged (or best) solution.
#[allow(dead_code)]
fn gauss_newton_closest_point(
    p1: Point3,
    surface2: &SurfaceGeometry,
    u0: f64,
    v0: f64,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> f64 {
    if !parameter_box_is_valid(u2_range, v2_range) || !u0.is_finite() || !v0.is_finite() {
        return f64::INFINITY;
    }

    const MAX_ITER: usize = 10;
    const STEP_TOL: f64 = 1e-12;

    let mut u = clamp_to_range(u0, u2_range);
    let mut v = clamp_to_range(v0, v2_range);

    for _ in 0..MAX_ITER {
        let s = surface2.point_at(u, v);
        let diff = p1 - s; // Vec3: (P1 - S2)

        let su = surface2.partial_u(u, v);
        let sv = surface2.partial_v(u, v);

        // Gradient of f = -2 * (P1-S2)·S2_u  and  -2 * (P1-S2)·S2_v
        // For the Newton step we solve the 2×2 Gauss-Newton normal equations:
        //   [su·su  su·sv] [Δu]   [(P1-S2)·su]
        //   [sv·su  sv·sv] [Δv] = [(P1-S2)·sv]
        let rhs_u = diff.dot(&su);
        let rhs_v = diff.dot(&sv);

        let a00 = su.dot(&su);
        let a01 = su.dot(&sv);
        let a11 = sv.dot(&sv);
        let det = a00 * a11 - a01 * a01;

        if det.abs() < f64::EPSILON {
            // Degenerate Jacobian — cannot improve further.
            break;
        }

        let du = (a11 * rhs_u - a01 * rhs_v) / det;
        let dv = (a00 * rhs_v - a01 * rhs_u) / det;

        u = clamp_to_range(u + du, u2_range);
        v = clamp_to_range(v + dv, v2_range);

        if du * du + dv * dv < STEP_TOL {
            break;
        }
    }

    let s_final = surface2.point_at(u, v);
    p1.distance_to(&s_final)
}

fn parameter_box_is_valid(u_range: (f64, f64), v_range: (f64, f64)) -> bool {
    parameter_range_is_valid(u_range) && parameter_range_is_valid(v_range)
}

fn parameter_range_is_valid(range: (f64, f64)) -> bool {
    range.0.is_finite() && range.1.is_finite() && range.0 < range.1
}

fn clamp_to_range(value: f64, range: (f64, f64)) -> f64 {
    debug_assert!(parameter_range_is_valid(range));
    value.clamp(range.0, range.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Vec3;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn identical_planes_coincident() {
        let s = xy_plane();
        let regions = detect_coincident_regions(
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(!regions.is_empty());
        assert!(matches!(
            &regions[0],
            SurfaceSurfaceIntersection::CoincidentRegion { .. }
        ));
    }

    #[test]
    fn perpendicular_planes_not_coincident() {
        // XY plane (z = 0)
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // XZ plane (y = 0)
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let regions = detect_coincident_regions(
            &s1,
            (0.0, 1.0),
            (0.0, 1.0),
            &s2,
            (0.0, 1.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(regions.is_empty());
    }

    #[test]
    fn parallel_offset_planes_not_coincident() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // Same orientation but offset by 1 along Z.
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let regions = detect_coincident_regions(
            &s1,
            (0.0, 1.0),
            (0.0, 1.0),
            &s2,
            (0.0, 1.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(regions.is_empty());
    }

    #[test]
    fn reversed_projection_range_declines_without_panic() {
        let s = xy_plane();
        let distance =
            closest_point_distance(Point3::new(0.25, 0.25, 0.0), &s, (0.0, -1.0), (0.0, 1.0));

        assert!(distance.is_infinite());
    }

    #[test]
    fn invalid_coincident_region_range_declines_without_panic() {
        let s = xy_plane();
        let regions = detect_coincident_regions(
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &s,
            (0.0, -1.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );

        assert!(regions.is_empty());
    }
}
