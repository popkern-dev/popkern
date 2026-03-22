//! Line-torus intersection algorithm.
//!
//! Produces a degree-4 polynomial from the torus implicit equation and finds
//! real roots via sign-change bracketing + Newton-bisection refinement.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Evaluate a quartic polynomial `a4*t^4 + a3*t^3 + a2*t^2 + a1*t + a0` using
/// Horner's method for numerical stability.
fn eval_quartic(coeffs: &[f64; 5], t: f64) -> f64 {
    // coeffs = [a0, a1, a2, a3, a4]
    let (a0, a1, a2, a3, a4) = (coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]);
    ((((a4 * t + a3) * t) + a2) * t + a1) * t + a0
}

/// Evaluate the derivative of the quartic: `4*a4*t^3 + 3*a3*t^2 + 2*a2*t + a1`.
fn eval_quartic_deriv(coeffs: &[f64; 5], t: f64) -> f64 {
    let (a1, a2, a3, a4) = (coeffs[1], coeffs[2], coeffs[3], coeffs[4]);
    (((4.0 * a4 * t) + 3.0 * a3) * t + 2.0 * a2) * t + a1
}

/// Find real roots of a quartic polynomial within `[lo, hi]` using sign-change
/// bracketing followed by Newton-bisection hybrid refinement.
fn find_quartic_roots(coeffs: &[f64; 5], lo: f64, hi: f64) -> Vec<f64> {
    // Number of sample points for bracketing
    let n_samples = 2048;
    let dt = (hi - lo) / n_samples as f64;

    let mut roots = Vec::with_capacity(4);
    let mut prev_t = lo;
    let mut prev_val = eval_quartic(coeffs, lo);

    for i in 1..=n_samples {
        let t = lo + i as f64 * dt;
        let val = eval_quartic(coeffs, t);

        if prev_val.is_sign_positive() != val.is_sign_positive() && prev_val != 0.0 && val != 0.0 {
            // Sign change — bracket contains a root. Refine with Newton-bisection.
            if let Some(root) = refine_root(coeffs, prev_t, t) {
                // Avoid duplicate roots (e.g., a double root straddling two brackets)
                if roots.last().map_or(true, |&last: &f64| (last - root).abs() > 1e-10) {
                    roots.push(root);
                }
            }
        }

        // Check for exact zero at sample point
        if val.abs() < 1e-15 {
            if roots.last().map_or(true, |&last: &f64| (last - t).abs() > 1e-10) {
                roots.push(t);
            }
        }

        prev_t = t;
        prev_val = val;

        if roots.len() >= 4 {
            break;
        }
    }

    roots
}

/// Refine a root bracketed in `[a, b]` using Newton-bisection hybrid.
fn refine_root(coeffs: &[f64; 5], mut a: f64, mut b: f64) -> Option<f64> {
    let mut fa = eval_quartic(coeffs, a);
    let _fb = eval_quartic(coeffs, b);

    let mut root = 0.5 * (a + b);
    let max_iter = 64;

    for _ in 0..max_iter {
        let f_root = eval_quartic(coeffs, root);

        if f_root.abs() < 1e-14 {
            return Some(root);
        }

        // Try Newton step
        let df = eval_quartic_deriv(coeffs, root);
        let mut used_newton = false;
        if df.abs() > 1e-30 {
            let newton = root - f_root / df;
            if newton > a && newton < b {
                root = newton;
                used_newton = true;
            }
        }

        if !used_newton {
            // Fall back to bisection
            root = 0.5 * (a + b);
        }

        let f_mid = eval_quartic(coeffs, root);
        if fa.is_sign_positive() != f_mid.is_sign_positive() {
            b = root;
        } else {
            a = root;
            fa = f_mid;
        }

        if (b - a).abs() < 1e-14 {
            return Some(0.5 * (a + b));
        }
    }

    Some(0.5 * (a + b))
}

/// Normalize an angle into `[0, 2*pi)`.
fn normalize_angle(a: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let r = a % tau;
    if r < 0.0 { r + tau } else { r }
}

/// Check if an angle is within `[range_start, range_end]`, handling wrapping.
fn angle_in_range(angle: f64, start: f64, end: f64) -> bool {
    let tau = std::f64::consts::TAU;

    // If the raw range spans a full circle (or more), accept everything.
    if (end - start).abs() >= tau - 1e-12 {
        return true;
    }

    let a = normalize_angle(angle);
    let s = normalize_angle(start);
    let e = normalize_angle(end);

    if s <= e {
        a >= s - 1e-10 && a <= e + 1e-10
    } else {
        a >= s - 1e-10 || a <= e + 1e-10
    }
}

/// Intersect a line `P + t*D` with a torus defined by center, axis, major
/// radius `R`, and minor (tube) radius `r`.
///
/// The torus surface is parameterised as:
///   `u` = angle around major circle (longitude), `v` = angle around tube (latitude).
///
/// Returns intersection points filtered to the given parameter ranges.
#[allow(dead_code)]
pub(crate) fn intersect_line_torus(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    torus_center: &Point3,
    torus_axis: &Vec3,
    torus_major_radius: f64,
    torus_minor_radius: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    let big_r = torus_major_radius;
    let small_r = torus_minor_radius;

    // Build the local coordinate frame: axis along Z, perp1 along X, perp2 along Y.
    let axis = match torus_axis.normalized() {
        Some(a) => a,
        None => return Vec::new(),
    };
    let (perp1, perp2) = orthonormal_basis(axis);

    // Transform line origin into the torus local frame.
    let offset = *line_origin - *torus_center;
    let p_local = Vec3::new(
        offset.dot(&perp1),
        offset.dot(&perp2),
        offset.dot(&axis),
    );

    // Transform line direction into the torus local frame.
    let d_local = Vec3::new(
        line_dir.dot(&perp1),
        line_dir.dot(&perp2),
        line_dir.dot(&axis),
    );

    // Torus implicit equation (axis along Z, center at origin):
    //   (x² + y² + z² + R² - r²)² = 4R²(x² + y²)
    //
    // Substitute ray: x = px + t*dx, y = py + t*dy, z = pz + t*dz
    //
    // Let:
    //   alpha = D·D (all 3 components)
    //   beta  = 2*(P·D)
    //   gamma = P·P + R² - r²
    //   delta = dx² + dy²  (XY projection of D)
    //   eps   = 2*(px*dx + py*dy) (XY projection dot)
    //   zeta  = px² + py²  (XY projection of P)
    //
    //   sum(t) = alpha*t² + beta*t + gamma  = |P+tD|² + R² - r²
    //   The equation is: sum(t)² = 4*R²*(delta*t² + eps*t + zeta)
    //
    // Expand sum² and subtract RHS to get quartic coefficients.

    let alpha = d_local.dot(&d_local);
    let beta = 2.0 * p_local.dot(&d_local);
    let gamma = p_local.dot(&p_local) + big_r * big_r - small_r * small_r;

    let delta = d_local.x * d_local.x + d_local.y * d_local.y;
    let eps = 2.0 * (p_local.x * d_local.x + p_local.y * d_local.y);
    let zeta = p_local.x * p_local.x + p_local.y * p_local.y;

    let four_r2 = 4.0 * big_r * big_r;

    // sum(t)² = (alpha*t² + beta*t + gamma)²
    // Expand:
    //   a4 = alpha²
    //   a3 = 2*alpha*beta
    //   a2 = beta² + 2*alpha*gamma
    //   a1 = 2*beta*gamma
    //   a0 = gamma²
    //
    // Subtract 4*R²*(delta*t² + eps*t + zeta):
    //   a4 = alpha²
    //   a3 = 2*alpha*beta
    //   a2 = beta² + 2*alpha*gamma - 4*R²*delta
    //   a1 = 2*beta*gamma - 4*R²*eps
    //   a0 = gamma² - 4*R²*zeta

    let a4 = alpha * alpha;
    let a3 = 2.0 * alpha * beta;
    let a2 = beta * beta + 2.0 * alpha * gamma - four_r2 * delta;
    let a1 = 2.0 * beta * gamma - four_r2 * eps;
    let a0 = gamma * gamma - four_r2 * zeta;

    // Degenerate: if a4 is essentially zero, the line direction is zero-length
    if a4.abs() < 1e-30 {
        return Vec::new();
    }

    // Normalize coefficients by a4 for better numerical conditioning
    let coeffs = [a0 / a4, a1 / a4, a2 / a4, a3 / a4, 1.0];

    // Search range: use t_range with generous padding to catch roots near boundaries
    let search_pad = (t_range.1 - t_range.0).max(10.0);
    let search_lo = t_range.0 - search_pad;
    let search_hi = t_range.1 + search_pad;

    let raw_roots = find_quartic_roots(&coeffs, search_lo, search_hi);

    let mut results = Vec::new();

    for root in raw_roots {
        // Check t is within the line's parameter range (with small tolerance)
        if root < t_range.0 - tolerance.parametric || root > t_range.1 + tolerance.parametric {
            continue;
        }

        // Clamp to exact range boundary
        let t = root.clamp(t_range.0, t_range.1);

        // Compute intersection point in local frame
        let local_pt = Vec3::new(
            p_local.x + t * d_local.x,
            p_local.y + t * d_local.y,
            p_local.z + t * d_local.z,
        );

        // Verify the point actually lies on the torus (within tolerance)
        let dist_xy = (local_pt.x * local_pt.x + local_pt.y * local_pt.y).sqrt();
        let tube_dist = ((dist_xy - big_r).powi(2) + local_pt.z * local_pt.z).sqrt();
        if (tube_dist - small_r).abs() > tolerance.linear.max(1e-6) {
            continue;
        }

        // Compute world-space intersection point
        let world_pt = *torus_center
            + perp1 * local_pt.x
            + perp2 * local_pt.y
            + axis * local_pt.z;

        // Compute torus surface parameters
        // u = angle around major circle = atan2(y_local, x_local)
        let u = normalize_angle(local_pt.y.atan2(local_pt.x));

        // v = angle around tube cross-section
        // v = atan2(z_local, dist_to_axis - R)
        let v = normalize_angle(local_pt.z.atan2(dist_xy - big_r));

        // Filter to u and v ranges
        if !angle_in_range(u, u_range.0, u_range.1) {
            continue;
        }
        if !angle_in_range(v, v_range.0, v_range.1) {
            continue;
        }

        results.push(SurfaceCurveIntersection::Point {
            point: world_pt,
            t,
            u,
            v,
        });
    }

    // Sort by t parameter
    results.sort_by(|a, b| {
        let ta = match a {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        let tb = match b {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn line_through_torus_four_points() {
        // Line through the center of a torus, perpendicular to axis
        // Should hit the torus 4 times (enter outer tube, exit inner, enter inner, exit outer)
        let r = intersect_line_torus(
            &Point3::new(-5.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-5.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            3.0,
            1.0, // major=3, minor=1
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 4, "expected 4 intersections, got {}: {:?}", r.len(), r);

        // Verify the t values correspond to the expected hit points:
        // Outer ring at x = -(R+r) = -4, inner at x = -(R-r) = -2,
        // then inner at x = R-r = 2, outer at x = R+r = 4.
        let t_values: Vec<f64> = r.iter().map(|i| match i {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            _ => panic!("expected Point"),
        }).collect();

        // Line origin is at x=-5, direction is (1,0,0), so t = x - (-5) = x + 5
        let expected_t = [1.0, 3.0, 7.0, 9.0]; // x = -4, -2, 2, 4
        for (actual, expected) in t_values.iter().zip(expected_t.iter()) {
            assert!(
                (actual - expected).abs() < 1e-6,
                "expected t={}, got t={}",
                expected,
                actual
            );
        }
    }

    #[test]
    fn line_missing_torus() {
        // Line far above the torus
        let r = intersect_line_torus(
            &Point3::new(-5.0, 0.0, 10.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-5.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            3.0,
            1.0,
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty(), "expected no intersections, got {}", r.len());
    }

    #[test]
    fn line_tangent_to_torus() {
        // Line tangent to the outer edge of the torus at y = R+r = 4
        let r = intersect_line_torus(
            &Point3::new(-5.0, 4.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-5.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            3.0,
            1.0,
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        // Tangent to outer edge at y=4 = R+r = 3+1
        assert!(
            r.len() <= 2,
            "tangent should give at most 2 intersections, got {}",
            r.len()
        );
    }

    #[test]
    fn line_two_points_through_tube() {
        // Line parallel to axis, passing through the tube at x=R=3, y=0
        // Should hit twice: top and bottom of the tube.
        let r = intersect_line_torus(
            &Point3::new(3.0, 0.0, -5.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (-5.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            3.0,
            1.0,
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 2, "expected 2 intersections, got {}: {:?}", r.len(), r);

        // Should hit at z = -1 and z = 1 (tube centered at distance R from axis)
        let t_values: Vec<f64> = r.iter().map(|i| match i {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            _ => panic!("expected Point"),
        }).collect();
        // Line origin z=-5, dir (0,0,1), so t = z + 5
        assert!((t_values[0] - 4.0).abs() < 1e-6, "expected t=4, got {}", t_values[0]);
        assert!((t_values[1] - 6.0).abs() < 1e-6, "expected t=6, got {}", t_values[1]);
    }

    #[test]
    fn tilted_torus_axis() {
        // Torus with axis along Y instead of Z — verify the rotation works
        let r = intersect_line_torus(
            &Point3::new(-5.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-5.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 1.0, 0.0), // axis along Y
            3.0,
            1.0,
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        // Line in XZ plane through center, perpendicular to Y axis
        // Same geometry as the 4-point test, just rotated
        assert_eq!(r.len(), 4, "expected 4 intersections with Y-axis torus, got {}", r.len());
    }

    #[test]
    fn offset_torus_center() {
        // Torus centered at (10, 0, 0), line through its center along X
        let r = intersect_line_torus(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 20.0),
            &Point3::new(10.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            3.0,
            1.0,
            (0.0, TAU),
            (0.0, TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 4, "expected 4 intersections, got {}", r.len());

        let t_values: Vec<f64> = r.iter().map(|i| match i {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            _ => panic!("expected Point"),
        }).collect();
        // Hits at x = 10-4, 10-2, 10+2, 10+4 → t = 6, 8, 12, 14
        let expected = [6.0, 8.0, 12.0, 14.0];
        for (actual, exp) in t_values.iter().zip(expected.iter()) {
            assert!(
                (actual - exp).abs() < 1e-6,
                "expected t={}, got t={}",
                exp,
                actual
            );
        }
    }
}
