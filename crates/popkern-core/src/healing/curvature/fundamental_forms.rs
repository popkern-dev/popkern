//! First and second fundamental forms of differential geometry.
//!
//! These are the mathematical building blocks for all curvature analysis.
//! They operate on any surface that implements [`SurfaceEval`], making them
//! representation-agnostic.

use crate::geometry::Vec3;
use crate::geometry::surface::SurfaceEval;

/// First fundamental form coefficients at a surface point.
///
/// Encodes the intrinsic metric: distances, angles, and area element
/// on the surface. Computed from the first partial derivatives.
///
/// - `E = S_u . S_u`
/// - `F = S_u . S_v`
/// - `G = S_v . S_v`
#[derive(Debug, Clone, Copy)]
pub struct FirstFundamentalForm {
    /// E = S_u . S_u
    pub e: f64,
    /// F = S_u . S_v
    pub f: f64,
    /// G = S_v . S_v
    pub g: f64,
}

impl FirstFundamentalForm {
    /// Area element: sqrt(EG - F^2). Zero at degenerate points.
    pub fn area_element(&self) -> f64 {
        let disc = self.e * self.g - self.f * self.f;
        if disc > 0.0 { disc.sqrt() } else { 0.0 }
    }

    /// Determinant: EG - F^2.
    pub fn determinant(&self) -> f64 {
        self.e * self.g - self.f * self.f
    }
}

/// Second fundamental form coefficients at a surface point.
///
/// Encodes the extrinsic curvature — how the surface curves in 3D space.
/// Computed from the second partial derivatives projected onto the normal.
///
/// - `L = S_uu . n`
/// - `M = S_uv . n`
/// - `N = S_vv . n`
#[derive(Debug, Clone, Copy)]
pub struct SecondFundamentalForm {
    /// L = S_uu . n
    pub l: f64,
    /// M = S_uv . n
    pub m: f64,
    /// N = S_vv . n (named `n_val` to avoid conflict with normal vector `n`)
    pub n_val: f64,
}

impl SecondFundamentalForm {
    /// Determinant: LN - M^2.
    pub fn determinant(&self) -> f64 {
        self.l * self.n_val - self.m * self.m
    }
}

/// Shape operator eigenvalue decomposition.
///
/// The shape operator W = I^{-1} * II has eigenvalues equal to the
/// principal curvatures and eigenvectors equal to the principal directions.
#[derive(Debug, Clone, Copy)]
pub struct ShapeOperator {
    /// First principal curvature (larger in absolute value).
    pub kappa_1: f64,
    /// Second principal curvature.
    pub kappa_2: f64,
    /// First principal direction (in 3D, along the surface).
    pub dir_1: Vec3,
    /// Second principal direction (in 3D, perpendicular to dir_1 on the surface).
    pub dir_2: Vec3,
}

/// Compute the first fundamental form at a surface point.
pub fn first_fundamental_form(surface: &dyn SurfaceEval, u: f64, v: f64) -> FirstFundamentalForm {
    let su = surface.partial_u(u, v);
    let sv = surface.partial_v(u, v);
    FirstFundamentalForm {
        e: su.dot(&su),
        f: su.dot(&sv),
        g: sv.dot(&sv),
    }
}

/// Compute the second fundamental form at a surface point.
///
/// Returns `None` at singular points where the surface normal is
/// undefined (pole, collapsed parameter line).
pub fn second_fundamental_form(
    surface: &dyn SurfaceEval,
    u: f64,
    v: f64,
) -> Option<SecondFundamentalForm> {
    let n = surface.normal_at(u, v);
    let n_len = n.length();
    if n_len < 1e-15 {
        return None;
    }
    let n_unit = n * (1.0 / n_len);

    let suu = surface.second_partial_uu(u, v);
    let suv = surface.second_partial_uv(u, v);
    let svv = surface.second_partial_vv(u, v);

    Some(SecondFundamentalForm {
        l: suu.dot(&n_unit),
        m: suv.dot(&n_unit),
        n_val: svv.dot(&n_unit),
    })
}

/// Compute the shape operator (Weingarten map) at a surface point.
///
/// Returns principal curvatures and directions. Returns `None` at
/// singular points or where the first fundamental form is degenerate.
pub fn shape_operator(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<ShapeOperator> {
    let ff1 = first_fundamental_form(surface, u, v);
    let ff2 = second_fundamental_form(surface, u, v)?;

    let det_i = ff1.determinant();
    if det_i.abs() < 1e-30 {
        return None;
    }

    // Shape operator W = I^{-1} * II
    // W = (1/det_I) * [[G, -F], [-F, E]] * [[L, M], [M, N]]
    let inv_det = 1.0 / det_i;
    let w11 = inv_det * (ff1.g * ff2.l - ff1.f * ff2.m);
    let w12 = inv_det * (ff1.g * ff2.m - ff1.f * ff2.n_val);
    let w21 = inv_det * (ff1.e * ff2.m - ff1.f * ff2.l);
    let w22 = inv_det * (ff1.e * ff2.n_val - ff1.f * ff2.m);

    // Eigenvalues of 2x2 matrix [[w11, w12], [w21, w22]]
    // Using quadratic formula on characteristic polynomial:
    // lambda^2 - trace*lambda + det = 0
    let trace = w11 + w22;
    let det_w = w11 * w22 - w12 * w21;
    let disc = trace * trace - 4.0 * det_w;
    let disc_sqrt = if disc > 0.0 { disc.sqrt() } else { 0.0 };

    let kappa_1 = 0.5 * (trace + disc_sqrt);
    let kappa_2 = 0.5 * (trace - disc_sqrt);

    // Principal directions in parameter space, mapped to 3D
    let su = surface.partial_u(u, v);
    let sv = surface.partial_v(u, v);

    let (dir_1, dir_2) = if w12.abs() > 1e-15 {
        // Eigenvector for kappa_1: (kappa_1 - w22, w21)
        let d1_u = kappa_1 - w22;
        let d1_v = w21;
        let d1 = su * d1_u + sv * d1_v;

        let d2_u = kappa_2 - w22;
        let d2_v = w21;
        let d2 = su * d2_u + sv * d2_v;
        (d1.normalized().unwrap_or(su), d2.normalized().unwrap_or(sv))
    } else if w21.abs() > 1e-15 {
        // Use rows from (W - lambda*I) = 0
        let d1_u = w12;
        let d1_v = kappa_1 - w11;
        let d1 = su * d1_u + sv * d1_v;

        let d2_u = w12;
        let d2_v = kappa_2 - w11;
        let d2 = su * d2_u + sv * d2_v;
        (d1.normalized().unwrap_or(su), d2.normalized().unwrap_or(sv))
    } else {
        // Diagonal matrix — principal directions align with parameter directions
        (
            su.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0)),
            sv.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0)),
        )
    };

    Some(ShapeOperator {
        kappa_1,
        kappa_2,
        dir_1,
        dir_2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::geometry::surface::SurfaceGeometry;

    #[test]
    fn plane_has_zero_second_form() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let ff2 = second_fundamental_form(&plane, 0.5, 0.5);
        let form = ff2.expect("plane normal should be defined");
        assert!(form.l.abs() < 1e-10);
        assert!(form.m.abs() < 1e-10);
        assert!(form.n_val.abs() < 1e-10);
    }

    #[test]
    fn plane_first_form_is_identity() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let ff1 = first_fundamental_form(&plane, 0.0, 0.0);
        assert!((ff1.e - 1.0).abs() < 1e-10);
        assert!(ff1.f.abs() < 1e-10);
        assert!((ff1.g - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_principal_curvatures_are_1_over_r() {
        let r = 2.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
            axis: crate::geometry::Vec3::Z,
        };
        // Evaluate away from poles (u=1.0, v=pi/2 is equator)
        let so = shape_operator(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
        let op = so.expect("sphere should have well-defined shape operator at equator");
        // Principal curvatures of a sphere are both 1/r (or -1/r depending on normal convention)
        let expected = 1.0 / r;
        assert!(
            (op.kappa_1.abs() - expected).abs() < 0.01,
            "kappa_1={}, expected +/-{}",
            op.kappa_1,
            expected
        );
        assert!(
            (op.kappa_2.abs() - expected).abs() < 0.01,
            "kappa_2={}, expected +/-{}",
            op.kappa_2,
            expected
        );
    }

    #[test]
    fn area_element_positive_for_non_degenerate() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let ff1 = first_fundamental_form(&plane, 0.0, 0.0);
        assert!(ff1.area_element() > 0.0);
        assert!((ff1.area_element() - 1.0).abs() < 1e-10);
    }
}
