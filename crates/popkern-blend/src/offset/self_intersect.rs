//! Offset self-intersection detection.
//!
//! When an inward offset exceeds the local curvature radius the offset
//! surface folds back on itself. This module checks whether a proposed
//! offset distance is safe, and if not, returns the maximum safe radius.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Check whether offsetting `surface` by `distance` would cause
/// self-intersection.
///
/// Returns `Some(max_safe_radius)` if the offset *would* self-intersect,
/// indicating the largest safe offset distance. Returns `None` when the
/// offset is safe.
///
/// For analytic primitives the curvature radius is known exactly. For
/// NURBS and cone surfaces the curvature is estimated by sampling second
/// derivatives over the given UV range.
pub fn check_offset_self_intersection(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Option<f64> {
    // The veto applies when the offset moves TOWARD the center of
    // curvature — for radius-form surfaces that is the radius-shrinking
    // direction, which depends on the stored normal's radial sign
    // (contract: `inward` means "along -n"). Radius-growing offsets on
    // analytic surfaces never self-intersect.
    let sign = if inward { -1.0 } else { 1.0 };
    match surface {
        SurfaceGeometry::Plane { .. } => return None,
        SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Sphere { .. }
        | SurfaceGeometry::Torus { .. } => {
            let radial = super::analytic::radial_normal_sign(surface)?;
            if sign * radial > 0.0 {
                return None; // radius grows: safe
            }
        }
        _ => {} // NURBS/cone: sampled curvature check both ways
    }

    let min_curvature_radius = min_curvature_radius(surface, u_range, v_range);

    if distance < min_curvature_radius {
        None
    } else {
        // Return the max safe radius, slightly inside the curvature radius.
        let safe = min_curvature_radius * 0.99;
        Some(if safe > 0.0 { safe } else { 0.0 })
    }
}

/// Compute the minimum curvature radius of a surface.
///
/// For analytic surfaces this is exact. For NURBS/Cone it is estimated by
/// sampling.
fn min_curvature_radius(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> f64 {
    match surface {
        SurfaceGeometry::Plane { .. } => f64::INFINITY,
        SurfaceGeometry::Cylinder { radius, .. } => *radius,
        SurfaceGeometry::Sphere { radius, .. } => *radius,
        SurfaceGeometry::Torus { minor_radius, .. } => *minor_radius,
        _ => sample_min_curvature_radius(surface, u_range, v_range),
    }
}

/// Estimate minimum curvature radius by sampling second derivatives.
///
/// Uses a 16x16 grid over the UV domain. At each sample the principal
/// curvature is approximated from the second fundamental form.
fn sample_min_curvature_radius(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> f64 {
    let n = 16_usize;
    let mut min_r = f64::INFINITY;
    let h = 1e-5;

    for i in 0..n {
        for j in 0..n {
            let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n - 1) as f64);
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n - 1) as f64);

            // First derivatives (already available via SurfaceEval).
            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);

            // Second derivatives via finite differences on first derivatives.
            let suu = (surface.partial_u(u + h, v) - surface.partial_u(u - h, v)) * (0.5 / h);
            let svv = (surface.partial_v(u, v + h) - surface.partial_v(u, v - h)) * (0.5 / h);
            let suv = (surface.partial_u(u, v + h) - surface.partial_u(u, v - h)) * (0.5 / h);

            // Surface normal.
            let normal = su.cross(&sv);
            let normal_len = normal.length();
            if normal_len < 1e-12 {
                continue; // degenerate point
            }
            let n_hat = normal * (1.0 / normal_len);

            // Coefficients of the second fundamental form.
            let l = suu.dot(&n_hat);
            let m = suv.dot(&n_hat);
            let n_coeff = svv.dot(&n_hat);

            // Coefficients of the first fundamental form.
            let e = su.dot(&su);
            let f = su.dot(&sv);
            let g = sv.dot(&sv);

            // Mean curvature H and Gaussian curvature K from the shape operator.
            let denom = e * g - f * f;
            if denom.abs() < 1e-20 {
                continue;
            }
            let h_curv = (e * n_coeff - 2.0 * f * m + g * l) / (2.0 * denom);
            let k_curv = (l * n_coeff - m * m) / denom;

            // Principal curvatures from H and K: kappa = H ± sqrt(H² - K).
            let disc = h_curv * h_curv - k_curv;
            let sqrt_disc = if disc > 0.0 { disc.sqrt() } else { 0.0 };
            let kappa1 = (h_curv + sqrt_disc).abs();
            let kappa2 = (h_curv - sqrt_disc).abs();
            let max_kappa = kappa1.max(kappa2);

            if max_kappa > 1e-12 {
                let r = 1.0 / max_kappa;
                if r < min_r {
                    min_r = r;
                }
            }
        }
    }

    min_r
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};

    #[test]
    fn tight_cylinder_detects_self_intersection() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        // Inward offset of 3.0 on a radius-2 cylinder should detect.
        let result = check_offset_self_intersection(
            &cyl,
            3.0,
            true,
            (0.0, std::f64::consts::TAU),
            (-1.0, 1.0),
        );
        assert!(result.is_some());
        let max_safe = result.unwrap();
        assert!(max_safe < 2.0, "max_safe={max_safe} should be < 2.0");
    }

    #[test]
    fn large_cylinder_is_safe() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 10.0,
        };
        let result = check_offset_self_intersection(
            &cyl,
            1.0,
            true,
            (0.0, std::f64::consts::TAU),
            (-1.0, 1.0),
        );
        assert!(result.is_none());
    }

    #[test]
    fn plane_always_safe() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = check_offset_self_intersection(&plane, 1000.0, true, (-1.0, 1.0), (-1.0, 1.0));
        assert!(result.is_none());
    }

    /// Cone falls through to the sampling path. A narrow cone (small half-angle)
    /// has high curvature near the apex — a large inward offset should detect.
    #[test]
    fn cone_inward_offset_detects_self_intersection() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: std::f64::consts::FRAC_PI_6, // 30 degrees
        };
        // Sample away from the apex where curvature is finite but non-trivial.
        let result = check_offset_self_intersection(
            &cone,
            50.0,
            true,
            (0.0, std::f64::consts::TAU),
            (0.5, 2.0),
        );
        // A very large offset on a cone should self-intersect.
        assert!(
            result.is_some(),
            "large inward offset on cone should detect self-intersection"
        );
        let max_safe = result.unwrap();
        assert!(max_safe > 0.0, "max safe radius should be positive");
    }

    /// Cone with a safe (small) inward offset should report None.
    #[test]
    fn cone_safe_offset() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: std::f64::consts::FRAC_PI_6,
        };
        // Far from the apex (v in [5, 10]) the curvature radius is large.
        let result = check_offset_self_intersection(
            &cone,
            0.001,
            true,
            (0.0, std::f64::consts::TAU),
            (5.0, 10.0),
        );
        assert!(
            result.is_none(),
            "tiny offset on cone far from apex should be safe"
        );
    }

    /// NURBS surface (bilinear flat patch) should behave like a plane —
    /// infinite curvature radius, always safe for inward offset.
    #[test]
    fn nurbs_flat_patch_safe() {
        use popkern_core::geometry::NurbsSurface;
        let nurbs = NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0; 4],
        };
        let surface = SurfaceGeometry::Nurbs(nurbs);
        let result = check_offset_self_intersection(&surface, 100.0, true, (0.0, 1.0), (0.0, 1.0));
        assert!(
            result.is_none(),
            "flat NURBS patch should be safe for any offset"
        );
    }

    /// NURBS surface with high curvature (curved inward) should detect
    /// self-intersection for a large offset.
    #[test]
    fn nurbs_curved_detects_self_intersection() {
        use popkern_core::geometry::NurbsSurface;
        // Create a patch that bulges upward, giving it noticeable curvature.
        let nurbs = NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                // Row 0
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                // Row 1 — bulges upward
                Point3::new(0.0, 0.5, 1.0),
                Point3::new(0.5, 0.5, 1.0),
                Point3::new(1.0, 0.5, 1.0),
                // Row 2
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        };
        let surface = SurfaceGeometry::Nurbs(nurbs);
        // Large inward offset on a surface with curvature ~1 should detect.
        let result = check_offset_self_intersection(&surface, 100.0, true, (0.0, 1.0), (0.0, 1.0));
        assert!(
            result.is_some(),
            "curved NURBS with large offset should detect self-intersection"
        );
    }

    /// Outward offset on analytic surfaces always returns None.
    #[test]
    fn outward_offset_analytic_always_safe() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: popkern_core::geometry::Vec3::Z,
        };
        // The standalone Sphere parametrization's normal points toward
        // the center, so the radius-GROWING (always-safe) direction is
        // the along--n flag (inward=true); along +n shrinks and must
        // veto past the radius.
        let grow = check_offset_self_intersection(
            &sphere,
            1000.0,
            true,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
        );
        assert!(
            grow.is_none(),
            "radius-growing offset on sphere is always safe"
        );
        let shrink = check_offset_self_intersection(
            &sphere,
            1000.0,
            false,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
        );
        assert!(
            shrink.is_some(),
            "radius-shrinking offset past the sphere radius must veto"
        );
    }

    /// Outward offset on a cone still goes through sampling (non-analytic path).
    #[test]
    fn outward_offset_cone_sampling_path() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: std::f64::consts::FRAC_PI_6,
        };
        // Outward offset on a cone — should go through sampling path but still be safe
        // for a small offset far from apex.
        let result = check_offset_self_intersection(
            &cone,
            0.001,
            false,
            (0.0, std::f64::consts::TAU),
            (5.0, 10.0),
        );
        assert!(
            result.is_none(),
            "small outward offset on cone should be safe"
        );
    }

    /// The safe radius returned should be slightly below the curvature radius.
    #[test]
    fn safe_radius_is_below_curvature_radius() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 5.0,
            axis: popkern_core::geometry::Vec3::Z,
        };
        // Radius-shrinking direction for the sphere parametrization is
        // along +n (inward=false; its stored normal points center-ward).
        let result = check_offset_self_intersection(
            &sphere,
            10.0,
            false,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
        );
        assert!(result.is_some());
        let max_safe = result.unwrap();
        // Should be 5.0 * 0.99 = 4.95
        assert!((max_safe - 4.95).abs() < 0.01, "max_safe={max_safe}");
    }
}
