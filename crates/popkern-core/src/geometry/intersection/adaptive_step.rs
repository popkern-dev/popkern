//! Curvature-based adaptive step size for intersection curve marching.
//!
//! Uses second fundamental form information (via `second_partial_*`) to
//! estimate the maximum curvature of the intersection curve at a point.
//! The step size is inversely proportional to the curvature: smaller steps
//! in high-curvature regions, larger in flat regions.

use crate::geometry::Vec3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Estimate the curvature of the intersection curve from surface curvatures.
///
/// Uses the principal curvatures of both surfaces at the intersection point
/// to bound the intersection curve curvature. The intersection curve curvature
/// is bounded by the sum of the maximum principal curvatures of both surfaces.
///
/// Returns the estimated curvature magnitude (always non-negative).
pub(crate) fn estimate_intersection_curvature(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
) -> f64 {
    let kmax1 = max_normal_curvature(surface1, u1, v1);
    let kmax2 = max_normal_curvature(surface2, u2, v2);

    // The intersection curve curvature is bounded by the sum of surface curvatures
    // (Meusnier's theorem gives the relationship between normal and geodesic curvature)
    kmax1 + kmax2
}

/// Compute the maximum normal curvature of a surface at (u, v).
///
/// Uses the eigenvalues of the shape operator (second fundamental form / first
/// fundamental form). For numerical stability, uses a simplified estimate:
/// max(|Suu|, |Svv|, |Suv|) / min(|Su|, |Sv|)^2
fn max_normal_curvature(surface: &SurfaceGeometry, u: f64, v: f64) -> f64 {
    let su = surface.partial_u(u, v);
    let sv = surface.partial_v(u, v);

    let su_len = su.length();
    let sv_len = sv.length();

    // Degenerate surface point — can't compute curvature
    if su_len < 1e-12 || sv_len < 1e-12 {
        return 0.0;
    }

    let suu = surface.second_partial_uu(u, v);
    let suv = surface.second_partial_uv(u, v);
    let svv = surface.second_partial_vv(u, v);

    let n = su.cross(&sv);
    let n_len = n.length();
    if n_len < 1e-12 {
        return 0.0;
    }
    let n_unit = n * (1.0 / n_len);

    // Second fundamental form coefficients: L, M, N
    let big_l = suu.dot(&n_unit);
    let big_m = suv.dot(&n_unit);
    let big_n = svv.dot(&n_unit);

    // First fundamental form coefficients: E, F, G
    let big_e = su.dot(&su);
    let big_f = su.dot(&sv);
    let big_g = sv.dot(&sv);

    // Principal curvatures are eigenvalues of the shape operator:
    // det(L-κE, M-κF; M-κF, N-κG) = 0
    // => (EG-F²)κ² - (EN+GL-2FM)κ + (LN-M²) = 0
    let denom = big_e * big_g - big_f * big_f;
    if denom.abs() < 1e-20 {
        return 0.0;
    }

    let trace = (big_e * big_n + big_g * big_l - 2.0 * big_f * big_m) / denom;
    let det = (big_l * big_n - big_m * big_m) / denom;

    // κ = (trace ± sqrt(trace² - 4*det)) / 2
    let discriminant = trace * trace - 4.0 * det;
    if discriminant < 0.0 {
        // Complex eigenvalues shouldn't happen for real surfaces, but be safe
        return trace.abs() * 0.5;
    }

    let sqrt_d = discriminant.sqrt();
    let k1 = (trace + sqrt_d) * 0.5;
    let k2 = (trace - sqrt_d) * 0.5;

    k1.abs().max(k2.abs())
}

/// Compute the adaptive step size given the intersection curvature.
///
/// The step is chosen so that the chord deviation (sag) from the true curve
/// is approximately `tolerance / 2`. For a circular arc of radius R = 1/κ,
/// the sag for a chord of length h is approximately h²κ/8.
/// Setting sag = tol/2: h = 2*sqrt(tol/κ).
///
/// The result is clamped to `[min_step, max_step]`.
pub(crate) fn curvature_based_step(
    curvature: f64,
    tolerance: f64,
    min_step: f64,
    max_step: f64,
) -> f64 {
    if curvature < 1e-12 {
        // Flat region — use maximum step
        return max_step;
    }

    // h = 2 * sqrt(tol / κ)
    let ideal = 2.0 * (tolerance / curvature).sqrt();
    ideal.clamp(min_step, max_step)
}

/// Compute the Jacobian condition number proxy for near-tangent detection.
///
/// Returns `|N1 × N2| / (|N1| * |N2|)` — the sine of the angle between
/// the surface normals. Values near zero indicate near-tangent configuration.
#[allow(dead_code)]
pub(crate) fn normal_cross_ratio(n1: &Vec3, n2: &Vec3) -> f64 {
    let cross_len = n1.cross(n2).length();
    let denom = n1.length() * n2.length();
    if denom < 1e-30 {
        return 0.0;
    }
    cross_len / denom
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    #[test]
    fn plane_has_zero_curvature() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let k = max_normal_curvature(&plane, 0.0, 0.0);
        assert!(k < 1e-10, "plane curvature should be zero, got {k}");
    }

    #[test]
    fn sphere_curvature_is_inverse_radius() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 2.0,
            axis: crate::geometry::Vec3::Z,
        };
        // At equator, away from poles
        let k = max_normal_curvature(&sphere, 0.5, std::f64::consts::FRAC_PI_2);
        // Principal curvature of a sphere is 1/r = 0.5
        assert!(
            (k - 0.5).abs() < 0.05,
            "sphere curvature expected ~0.5, got {k}"
        );
    }

    #[test]
    fn cylinder_curvature_is_inverse_radius() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        // Cylinder has κ1 = 1/r, κ2 = 0
        let k = max_normal_curvature(&cyl, 0.0, 0.0);
        assert!(
            (k - 1.0).abs() < 0.1,
            "cylinder curvature expected ~1.0, got {k}"
        );
    }

    #[test]
    fn flat_region_gives_max_step() {
        let step = curvature_based_step(0.0, 1e-6, 1e-6, 1e-3);
        assert!((step - 1e-3).abs() < 1e-10);
    }

    #[test]
    fn high_curvature_gives_small_step() {
        let step = curvature_based_step(100.0, 1e-6, 1e-8, 1e-3);
        // h = 2*sqrt(1e-6 / 100) = 2*sqrt(1e-8) = 2e-4
        assert!(step < 1e-3, "expected step < max, got {step}");
        assert!(step > 1e-8, "expected step > min, got {step}");
    }

    #[test]
    fn normal_cross_ratio_perpendicular() {
        let n1 = Vec3::new(1.0, 0.0, 0.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let ratio = normal_cross_ratio(&n1, &n2);
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn normal_cross_ratio_parallel() {
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let n2 = Vec3::new(0.0, 0.0, 1.0);
        let ratio = normal_cross_ratio(&n1, &n2);
        assert!(ratio < 1e-10);
    }
}
