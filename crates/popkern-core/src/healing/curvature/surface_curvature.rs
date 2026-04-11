//! Surface curvature computation from the shape operator.
//!
//! All functions accept any surface implementing [`SurfaceEval`] and
//! return `None` at degenerate points (poles, collapsed parameter lines)
//! where the shape operator is undefined.

use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::fundamental_forms::shape_operator;

/// Aggregate surface curvature at a single point.
///
/// Contains principal curvatures, Gaussian and mean curvature,
/// and the minimum radius of curvature.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceCurvature {
    /// First principal curvature (larger in absolute value).
    pub kappa_1: f64,
    /// Second principal curvature.
    pub kappa_2: f64,
    /// Gaussian curvature: K = kappa_1 * kappa_2.
    pub gaussian: f64,
    /// Mean curvature: H = (kappa_1 + kappa_2) / 2.
    pub mean: f64,
    /// Minimum radius of curvature: 1 / max(|kappa_1|, |kappa_2|).
    /// Infinite (f64::INFINITY) when both principal curvatures are zero.
    pub min_radius: f64,
}

/// Compute all surface curvature quantities at a point.
///
/// Returns `None` at degenerate points where the shape operator
/// cannot be computed (singular normal, collapsed parameter lines).
pub fn surface_curvature_at(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<SurfaceCurvature> {
    let so = shape_operator(surface, u, v)?;
    let k1 = so.kappa_1;
    let k2 = so.kappa_2;
    let max_abs = k1.abs().max(k2.abs());
    let min_r = if max_abs > 1e-15 {
        1.0 / max_abs
    } else {
        f64::INFINITY
    };
    Some(SurfaceCurvature {
        kappa_1: k1,
        kappa_2: k2,
        gaussian: k1 * k2,
        mean: (k1 + k2) / 2.0,
        min_radius: min_r,
    })
}

/// Principal curvatures at a surface point.
///
/// Returns `(kappa_1, kappa_2)` where `kappa_1 >= kappa_2` in absolute
/// value. Returns `None` at degenerate points.
pub fn principal_curvatures(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<(f64, f64)> {
    let so = shape_operator(surface, u, v)?;
    Some((so.kappa_1, so.kappa_2))
}

/// Gaussian curvature at a surface point.
///
/// K = kappa_1 * kappa_2 = det(II) / det(I).
///
/// - Positive on elliptic points (sphere-like).
/// - Zero on parabolic points (cylinder-like).
/// - Negative on hyperbolic points (saddle-like).
///
/// Returns `None` at degenerate points.
pub fn gaussian_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    Some(so.kappa_1 * so.kappa_2)
}

/// Mean curvature at a surface point.
///
/// H = (kappa_1 + kappa_2) / 2.
///
/// Zero on minimal surfaces (soap films). For a sphere of radius r,
/// H = 1/r (up to sign convention).
///
/// Returns `None` at degenerate points.
pub fn mean_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    Some((so.kappa_1 + so.kappa_2) / 2.0)
}

/// Minimum radius of curvature at a surface point.
///
/// Defined as `1 / max(|kappa_1|, |kappa_2|)`. Returns
/// `f64::INFINITY` when both principal curvatures are zero (flat surface).
///
/// Returns `None` at degenerate points.
pub fn min_radius_of_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    let max_abs = so.kappa_1.abs().max(so.kappa_2.abs());
    if max_abs > 1e-15 {
        Some(1.0 / max_abs)
    } else {
        Some(f64::INFINITY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    #[test]
    fn sphere_gaussian_curvature() {
        let r = 3.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
        };
        // Evaluate away from poles: u=1.0, v=pi/2 (equator)
        let k = gaussian_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
        let expected = 1.0 / (r * r);
        let got = k.expect("sphere should have well-defined curvature at equator");
        assert!(
            (got - expected).abs() < 0.01,
            "Gaussian curvature: got {got}, expected {expected}"
        );
    }

    #[test]
    fn sphere_mean_curvature() {
        let r = 2.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
        };
        let h = mean_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
        let expected = 1.0 / r;
        let got = h.expect("sphere should have well-defined curvature at equator");
        // Mean curvature may be negative depending on normal orientation; check absolute value
        assert!(
            (got.abs() - expected).abs() < 0.01,
            "Mean curvature: got {got}, expected +/-{expected}"
        );
    }

    #[test]
    fn sphere_surface_curvature_at() {
        let r = 4.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(1.0, 2.0, 3.0),
            radius: r,
        };
        let sc = surface_curvature_at(&sphere, 0.5, std::f64::consts::FRAC_PI_2);
        let c = sc.expect("sphere curvature should be defined at equator");
        let expected_k = 1.0 / (r * r);
        let expected_h = 1.0 / r;
        assert!(
            (c.gaussian - expected_k).abs() < 0.01,
            "K: got {}, expected {expected_k}",
            c.gaussian
        );
        assert!(
            (c.mean.abs() - expected_h).abs() < 0.01,
            "H: got {}, expected +/-{expected_h}",
            c.mean
        );
        assert!(
            (c.min_radius - r).abs() < 0.1,
            "min_radius: got {}, expected {r}",
            c.min_radius
        );
    }

    #[test]
    fn plane_zero_curvature() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let k = gaussian_curvature(&plane, 0.5, 0.5);
        assert!(
            k.expect("plane curvature should be defined").abs() < 1e-10,
            "Plane should have K=0"
        );
        let h = mean_curvature(&plane, 0.5, 0.5);
        assert!(
            h.expect("plane curvature should be defined").abs() < 1e-10,
            "Plane should have H=0"
        );
        let mr = min_radius_of_curvature(&plane, 0.5, 0.5);
        assert_eq!(
            mr.expect("plane min_radius should be defined"),
            f64::INFINITY,
            "Plane should have infinite min radius"
        );
    }

    #[test]
    fn cylinder_curvature() {
        let r = 5.0;
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: r,
        };
        // Evaluate at arbitrary point on cylinder (u=1.0, v=0.0)
        let k = gaussian_curvature(&cyl, 1.0, 0.0);
        let got_k = k.expect("cylinder curvature should be defined");
        assert!(got_k.abs() < 0.01, "Cylinder should have K=0, got {got_k}");

        let h = mean_curvature(&cyl, 1.0, 0.0);
        let got_h = h.expect("cylinder curvature should be defined");
        let expected_h = 1.0 / (2.0 * r);
        assert!(
            (got_h.abs() - expected_h).abs() < 0.01,
            "Cylinder mean curvature: got {got_h}, expected +/-{expected_h}"
        );

        let mr = min_radius_of_curvature(&cyl, 1.0, 0.0);
        let got_mr = mr.expect("cylinder min_radius should be defined");
        assert!(
            (got_mr - r).abs() < 0.1,
            "Cylinder min_radius: got {got_mr}, expected {r}"
        );
    }

    #[test]
    fn principal_curvatures_sphere() {
        let r = 2.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
        };
        let (k1, k2) = principal_curvatures(&sphere, 1.0, std::f64::consts::FRAC_PI_2)
            .expect("should be defined at equator");
        let expected = 1.0 / r;
        assert!(
            (k1.abs() - expected).abs() < 0.01,
            "kappa_1: got {k1}, expected +/-{expected}"
        );
        assert!(
            (k2.abs() - expected).abs() < 0.01,
            "kappa_2: got {k2}, expected +/-{expected}"
        );
    }
}
