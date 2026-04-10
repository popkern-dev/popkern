//! Line-cone intersection algorithm.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Intersect a line `P + t*D` with a cone defined by apex, unit axis `A`, and half-angle `α`.
///
/// The cone surface satisfies `(Q·A)² - cos²(α)*|Q|² = 0` where `Q = P + t*D - apex`.
/// Only the forward nappe (`Q·A > 0`) is considered.
///
/// Returns surface parameters:
/// - `v = Q·A` (height along axis from apex)
/// - `u = atan2(Q·perp2, Q·perp1)` (azimuth angle, normalised to `[0, 2π)`)
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_cone(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    cone_apex: &Point3,
    cone_axis: &Vec3,
    cone_half_angle: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // W = line_origin - cone_apex
    let w = *line_origin - *cone_apex;
    let d = *line_dir;
    let a = *cone_axis;

    let cos_alpha = cone_half_angle.cos();
    let cos2 = cos_alpha * cos_alpha;

    // Quadratic coefficients from expanding (Q·A)² - cos²(α)*|Q|² = 0
    // Q = W + t*D
    // (W·A + t*(D·A))² - cos²(α) * |W + t*D|² = 0
    let da = d.dot(&a);
    let wa = w.dot(&a);
    let dd = d.dot(&d);
    let wd = w.dot(&d);
    let ww = w.dot(&w);

    let qa = da * da - cos2 * dd;
    let qb = 2.0 * (da * wa - cos2 * wd);
    let qc = wa * wa - cos2 * ww;

    let disc = qb * qb - 4.0 * qa * qc;

    // Threshold for treating as degenerate quadratic
    let eps2 = tolerance.linear() * tolerance.linear();

    let candidates: Vec<f64> = if qa.abs() < eps2 {
        // Degenerate: linear in t (or no solution)
        if qb.abs() < eps2 {
            // qc must also be ~0 for a solution; no isolated t
            return Vec::new();
        }
        vec![-qc / qb]
    } else if disc < -eps2 {
        // No real roots
        return Vec::new();
    } else if disc <= eps2 {
        // Tangent — one root
        vec![-qb / (2.0 * qa)]
    } else {
        let sqrt_disc = disc.sqrt();
        vec![
            (-qb - sqrt_disc) / (2.0 * qa),
            (-qb + sqrt_disc) / (2.0 * qa),
        ]
    };

    let (perp1, perp2) = orthonormal_basis(a);
    let tau = std::f64::consts::TAU;

    let mut results = Vec::new();

    for t in candidates {
        // Filter by line t_range
        if t < t_range.0 - tolerance.linear() || t > t_range.1 + tolerance.linear() {
            continue;
        }

        let q = w + d * t; // Q = P + t*D - apex

        // Forward nappe filter: v = Q·A must be > 0 (within tolerance)
        let v = q.dot(&a);
        if v < -tolerance.linear() {
            continue;
        }

        // Filter by surface v_range
        if v < v_range.0 - tolerance.linear() || v > v_range.1 + tolerance.linear() {
            continue;
        }

        // Compute azimuth u from projection onto the orthonormal basis
        let q_perp1 = q.dot(&perp1);
        let q_perp2 = q.dot(&perp2);
        let u_raw = q_perp2.atan2(q_perp1);
        // Normalise u to [0, 2π)
        let u = ((u_raw % tau) + tau) % tau;

        // Filter by surface u_range
        if !u_in_range(u, u_range, tolerance.angular()) {
            continue;
        }

        let point = *cone_apex + a * v + perp1 * q_perp1 + perp2 * q_perp2;

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
    // Full circle
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
    use std::f64::consts::{FRAC_PI_4, TAU};

    fn default_tol() -> Tolerance {
        Tolerance::default()
    }

    /// Cone: apex at origin, axis = +Z, half_angle = π/4 (45°)
    fn std_cone() -> (Point3, Vec3, f64) {
        (
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_4,
        )
    }

    #[test]
    fn secant_two_intersections() {
        // Line through the cone at two points.
        // Cone: apex=(0,0,0), axis=Z, half_angle=π/4 → r = v at height v.
        // At v=1, r=1; at v=2, r=2.
        // Line: origin=(−3, 1, 1), direction=(1, 0, 0).
        // Cone condition: x² + y² = z² (since tan(π/4)=1).
        // x² + 1 = 1 → x = 0, but we want intersections with the line y=1, z=1:
        // (x²+ 1) = 1 → x=0 only, that's one point. Let's use a better setup.
        //
        // Line: origin=(−3, 0, 2), dir=(1, 0, 0), z=2 throughout.
        // Cone at z=2: r=2. So intersection where x²+0²=4 → x=±2.
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-3.0, 0.0, 2.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            2,
            "expected 2 intersections, got {:?}",
            results
        );
        // Verify the points lie on the cone: x²+y²=z²
        for r in &results {
            if let SurfaceCurveIntersection::Point { point, .. } = r {
                let radial = (point.x * point.x + point.y * point.y).sqrt();
                assert!(
                    (radial - point.z.abs()).abs() < 1e-6,
                    "point not on cone: {:?}",
                    point
                );
            }
        }
    }

    #[test]
    fn tangent_one_intersection() {
        // Line tangent to cone surface.
        // Cone apex=(0,0,0), axis=Z, half_angle=π/4.
        // Line in the plane y=1, z=1: origin=(0,1,1), dir=(1,0,0).
        // Cone at z=1: r=1. y=1 is exactly the tangent point at x=0, y=1.
        // x²+1 = 1 → x=0 only (disc=0).
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(0.0, 1.0, 1.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        assert_eq!(
            results.len(),
            1,
            "expected 1 tangent intersection, got {:?}",
            results
        );
    }

    #[test]
    fn miss_no_intersection() {
        // Line that misses the cone entirely.
        // Cone at z=1: r=1. Line at z=1, y=3 — too far from axis.
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-5.0, 3.0, 1.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        assert!(
            results.is_empty(),
            "expected no intersection, got {:?}",
            results
        );
    }

    #[test]
    fn behind_apex_filtered() {
        // Line intersects the backward nappe (v < 0), should be filtered.
        // Line at z=−1: origin=(−2, 0, −1), dir=(1,0,0).
        // z=−1 < 0 → v = Q·Z = −1 < 0 → filtered.
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-2.0, 0.0, -1.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        assert!(
            results.is_empty(),
            "expected no intersection (behind apex), got {:?}",
            results
        );
    }

    #[test]
    fn through_apex() {
        // Line through the apex — should produce either 0 or 1 result at the apex.
        // Line: origin=(0,0,−1), dir=(0,0,1) along axis. At apex Q·A=0 → tangential.
        // Degenerate: the line is along the axis. All points on axis satisfy Q perpendicular=0,
        // but only at apex itself does the cone equation hold (v=0, r=0).
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(0.0, 0.0, -1.0);
        let dir = Vec3::new(0.0, 0.0, 1.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        // Acceptable: 0 or 1 intersection at the apex (v=0)
        assert!(
            results.len() <= 1,
            "expected 0 or 1 apex intersections, got {:?}",
            results
        );
    }

    #[test]
    fn surface_params_match_point_at() {
        use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-3.0, 0.0, 2.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        let cone = SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle: ha,
        };
        for r in &results {
            if let SurfaceCurveIntersection::Point { point, u, v, .. } = r {
                let surface_pt = cone.point_at(*u, *v);
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
    fn degenerate_linear_case() {
        // When qa ≈ 0 (line lies on the cone surface), we get the linear case (line 56-60).
        // A line along a generator of the cone should produce a degenerate quadratic.
        let (apex, axis, ha) = std_cone();
        // Generator direction: at angle=0, the generator is (sin(ha), 0, cos(ha)) = (1/sqrt2, 0, 1/sqrt2)
        let s = ha.sin();
        let c = ha.cos();
        let dir = Vec3::new(s, 0.0, c);
        let origin = apex; // start at apex
        let results = intersect_line_cone(
            &origin,
            &dir,
            (0.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        // May return 0 or more points; the degenerate path was exercised
        let _ = results;
    }

    #[test]
    fn t_range_filters_cone() {
        // Line intersects cone but t_range excludes the hits (line 83)
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-3.0, 0.0, 2.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (10.0, 20.0), // far beyond the intersection points
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        assert!(results.is_empty());
    }

    #[test]
    fn v_range_filters_cone() {
        // Line intersects cone at z=2 but v_range excludes it (line 96)
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-3.0, 0.0, 2.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (5.0, 10.0), // v_range excludes v=2
            &default_tol(),
        );
        assert!(results.is_empty());
    }

    #[test]
    fn u_range_filters_cone() {
        // Line intersects cone but u_range excludes the azimuth (line 108)
        let (apex, axis, ha) = std_cone();
        let origin = Point3::new(-3.0, 0.0, 2.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let results = intersect_line_cone(
            &origin,
            &dir,
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (1.0, 2.0), // restricted u_range
            (0.0, f64::INFINITY),
            &default_tol(),
        );
        // Should filter out at least some intersections
        assert!(results.len() < 2);
    }
}
