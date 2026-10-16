//! Line-sphere intersection algorithm.

use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Intersect a line `P + t*D` with a sphere of given center `C` and radius `R`.
///
/// Solves `|P + t*D - C|² = R²`, a quadratic in `t`.
///
/// Surface parameters follow the Z-axis spherical convention used in `analytic_surfaces`:
/// - `u` = longitude = `atan2(Q.y, Q.x)`, normalised to `[0, 2π)`
/// - `v` = colatitude = `acos(Q.z / R)` in `[0, π]`
///
/// where `Q = intersection_point - center`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_sphere(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    sphere_center: &Point3,
    sphere_radius: f64,
    sphere_axis: &Vec3,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // W = line_origin - sphere_center
    let w: Vec3 = *line_origin - *sphere_center;
    let d = *line_dir;
    let r = sphere_radius;

    // |W + t*D|² = R²
    // t²*(D·D) + 2t*(W·D) + (W·W - R²) = 0
    let a = d.dot(&d);
    let b = 2.0 * w.dot(&d);
    let c = w.dot(&w) - r * r;

    let disc = b * b - 4.0 * a * c;
    let eps2 = tolerance.linear() * tolerance.linear();

    let candidates: Vec<f64> = if a < eps2 {
        // Degenerate: zero-length direction — no parameterised intersection
        return Vec::new();
    } else if disc < -eps2 {
        return Vec::new();
    } else if disc <= eps2 {
        // Tangent
        vec![-b / (2.0 * a)]
    } else {
        let sqrt_disc = disc.sqrt();
        vec![(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)]
    };

    let tau = std::f64::consts::TAU;

    let mut results = Vec::new();

    for t in candidates {
        if t < t_range.0 - tolerance.linear() || t > t_range.1 + tolerance.linear() {
            continue;
        }

        let point = *line_origin + d * t;
        // In the sphere's chart frame (identity for a +Z pole).
        let q: Vec3 = crate::geometry::analytic_surfaces::sphere_to_local(
            sphere_axis,
            point - *sphere_center,
        );

        // Colatitude v = acos(q.z / R), clamped to handle floating-point noise
        let cos_v = (q.z / r).clamp(-1.0, 1.0);
        let v = cos_v.pacos(); // in [0, π]

        // Longitude u = atan2(q.y, q.x), normalised to [0, 2π)
        let u_raw = q.y.patan2(q.x);
        let u = ((u_raw % tau) + tau) % tau;

        // Filter by surface v_range
        if v < v_range.0 - tolerance.angular() || v > v_range.1 + tolerance.angular() {
            continue;
        }

        // Filter by surface u_range (handle full-circle)
        if !u_in_range(u, u_range, tolerance.angular()) {
            continue;
        }

        results.push(SurfaceCurveIntersection::Point { point, t, u, v });
    }

    // Sort by t parameter
    results.sort_by(|a, b| {
        a.t_val()
            .partial_cmp(&b.t_val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results
}

/// Returns true if `u` is within the given u_range (handles full-circle wrap-around).
fn u_in_range(u: f64, u_range: (f64, f64), angular_tol: f64) -> bool {
    let tau = std::f64::consts::TAU;
    let (lo, hi) = u_range;
    if hi - lo >= tau - angular_tol {
        return true;
    }
    u >= lo - angular_tol && u <= hi + angular_tol
}

/// Helper to extract `t` from a `SurfaceCurveIntersection::Point`.
trait TVal {
    fn t_val(&self) -> f64;
}

impl TVal for SurfaceCurveIntersection {
    fn t_val(&self) -> f64 {
        match self {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range, .. } => t_range.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, TAU};

    fn default_tol() -> Tolerance {
        Tolerance::default()
    }

    /// Unit sphere at origin.
    fn unit_sphere() -> (Point3, f64) {
        (Point3::new(0.0, 0.0, 0.0), 1.0)
    }

    #[test]
    fn through_center_two_intersections() {
        // Line along X through centre of unit sphere: hits at x=−1 and x=+1.
        let (center, radius) = unit_sphere();
        let origin = Point3::new(-3.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            2,
            "expected 2 intersections, got {:?}",
            results
        );
        // t values should be 2 and 4 (origin at x=−3, hits x=−1 at t=2, x=+1 at t=4)
        let ts: Vec<f64> = results
            .iter()
            .map(|r| match r {
                SurfaceCurveIntersection::Point { t, .. } => *t,
                _ => panic!("unexpected overlap"),
            })
            .collect();
        assert!((ts[0] - 2.0).abs() < 1e-6, "t0 expected 2.0, got {}", ts[0]);
        assert!((ts[1] - 4.0).abs() < 1e-6, "t1 expected 4.0, got {}", ts[1]);
    }

    #[test]
    fn tangent_one_intersection() {
        // Line tangent to unit sphere: origin=(−2, 1, 0), dir=(1,0,0) — touches at (0,1,0).
        let (center, radius) = unit_sphere();
        let origin = Point3::new(-2.0, 1.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            1,
            "expected 1 tangent intersection, got {:?}",
            results
        );
        if let SurfaceCurveIntersection::Point { point, .. } = &results[0] {
            assert!(
                (point.distance_to(&Point3::new(0.0, 1.0, 0.0))) < 1e-5,
                "tangent point wrong: {:?}",
                point
            );
        }
    }

    #[test]
    fn miss_no_intersection() {
        // Line that misses the sphere: origin=(−2, 2, 0), dir=(1,0,0) — y=2 > radius=1.
        let (center, radius) = unit_sphere();
        let origin = Point3::new(-2.0, 2.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert!(
            results.is_empty(),
            "expected no intersection, got {:?}",
            results
        );
    }

    #[test]
    fn uv_params_match_point_at() {
        use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

        let (center, radius) = unit_sphere();
        let origin = Point3::new(-3.0, 0.0, 0.5);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        let sphere = SurfaceGeometry::Sphere {
            center,
            radius,
            axis: crate::geometry::Vec3::Z,
        };
        assert!(!results.is_empty(), "expected intersections");
        for r in &results {
            if let SurfaceCurveIntersection::Point { point, u, v, .. } = r {
                let surface_pt = sphere.point_at(*u, *v);
                assert!(
                    point.distance_to(&surface_pt) < 1e-5,
                    "surface params don't reproduce point: pt={:?}, surface_pt={:?}, u={}, v={}",
                    point,
                    surface_pt,
                    u,
                    v
                );
            }
        }
    }

    #[test]
    fn north_pole_v_near_zero() {
        // Line along Z through north pole: should hit at z=+1 (v≈0) and z=−1 (v≈π).
        let (center, radius) = unit_sphere();
        let origin = Point3::new(0.0, 0.0, -3.0);
        let dir = Vec3::new(0.0, 0.0, 1.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            2,
            "expected 2 intersections, got {:?}",
            results
        );
        let vs: Vec<f64> = results
            .iter()
            .map(|r| match r {
                SurfaceCurveIntersection::Point { v, .. } => *v,
                _ => panic!("overlap unexpected"),
            })
            .collect();
        // First hit: south pole v=π, second hit: north pole v=0
        assert!(
            vs[0] > PI - 0.01,
            "first hit should be south pole, v={}",
            vs[0]
        );
        assert!(vs[1] < 0.01, "second hit should be north pole, v={}", vs[1]);
    }

    #[test]
    fn offset_sphere_center() {
        // Sphere at (5,0,0) radius 2, line along X.
        let center = Point3::new(5.0, 0.0, 0.0);
        let radius = 2.0;
        let origin = Point3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_sphere(
            &origin,
            &dir,
            (0.0, 20.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            2,
            "expected 2 intersections, got {:?}",
            results
        );
        let ts: Vec<f64> = results
            .iter()
            .map(|r| match r {
                SurfaceCurveIntersection::Point { t, .. } => *t,
                _ => panic!(),
            })
            .collect();
        assert!((ts[0] - 3.0).abs() < 1e-6, "t0 expected 3.0, got {}", ts[0]);
        assert!((ts[1] - 7.0).abs() < 1e-6, "t1 expected 7.0, got {}", ts[1]);
    }

    #[test]
    fn zero_direction_returns_empty() {
        // Zero-length direction vector should return no intersections (line 43)
        let (center, radius) = unit_sphere();
        let results = intersect_line_sphere(
            &Point3::new(1.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 0.0),
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert!(results.is_empty());
    }

    #[test]
    fn t_range_filters_out() {
        // Line through sphere but t_range excludes the intersections (line 63)
        let (center, radius) = unit_sphere();
        let results = intersect_line_sphere(
            &Point3::new(-3.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (5.0, 10.0), // intersections at t=2 and t=4, but range starts at 5
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert!(results.is_empty());
    }

    #[test]
    fn v_range_filters_out() {
        // Restrict v_range to exclude south pole (line 79)
        let (center, radius) = unit_sphere();
        // Line along Z axis hits north pole (v=0) and south pole (v=PI)
        let results = intersect_line_sphere(
            &Point3::new(0.0, 0.0, -3.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, 0.1), // v_range excludes most of sphere
            &default_tol(),
        );
        // Only the north pole (v=0) should pass the v_range filter
        assert!(results.len() <= 1);
    }

    #[test]
    fn u_range_filters_out() {
        // Line through sphere with restricted u_range (line 84/102)
        let (center, radius) = unit_sphere();
        // Line along X from -3 through sphere. Hits at (-1,0,0) u=PI and (1,0,0) u=0
        let results = intersect_line_sphere(
            &Point3::new(-3.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &center,
            radius,
            &Vec3::Z,
            (0.5, 2.5), // restricted u_range, excludes u=0 and u=PI
            (0.0, PI),
            &default_tol(),
        );
        assert!(results.is_empty());
    }
}
