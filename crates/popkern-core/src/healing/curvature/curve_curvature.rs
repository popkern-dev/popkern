//! Curve curvature and torsion computation.
//!
//! Uses the Frenet-Serret formulas for space curves:
//! - Curvature kappa = |C' x C''| / |C'|^3
//! - Torsion tau = (C' x C'') . C''' / |C' x C''|^2

use crate::geometry::curve::CurveEval;

/// Curvature of a parametric curve at parameter `t`.
///
/// Computed as kappa = |C' x C''| / |C'|^3 where C' and C'' are the
/// first and second derivatives. Returns `None` if the tangent vector
/// has near-zero length (degenerate point or zero-speed parameterization).
pub fn curve_curvature_at(curve: &dyn CurveEval, t: f64) -> Option<f64> {
    let d1 = curve.derivative_at(t, 1);
    let d2 = curve.derivative_at(t, 2);

    let d1_len = d1.length();
    if d1_len < 1e-15 {
        return None;
    }

    let cross = d1.cross(&d2);
    let kappa = cross.length() / (d1_len * d1_len * d1_len);
    Some(kappa)
}

/// Torsion of a parametric curve at parameter `t`.
///
/// Computed as tau = (C' x C'') . C''' / |C' x C''|^2.
/// Measures how the curve twists out of its osculating plane.
///
/// Returns `None` if the cross product C' x C'' has near-zero length
/// (straight line, inflection point, or degenerate parameterization).
pub fn curve_torsion_at(curve: &dyn CurveEval, t: f64) -> Option<f64> {
    let d1 = curve.derivative_at(t, 1);
    let d2 = curve.derivative_at(t, 2);
    let d3 = curve.derivative_at(t, 3);

    let cross = d1.cross(&d2);
    let cross_len_sq = cross.length_squared();
    if cross_len_sq < 1e-30 {
        return None;
    }

    let tau = cross.dot(&d3) / cross_len_sq;
    Some(tau)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::{Point3, Vec3};

    #[test]
    fn line_has_zero_curvature() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 2.0, 3.0),
        };
        let kappa = curve_curvature_at(&line, 0.5);
        let got = kappa.expect("line curvature should be defined");
        assert!(
            got.abs() < 1e-10,
            "Line should have zero curvature, got {got}"
        );
    }

    #[test]
    fn line_torsion_is_none() {
        // Torsion is undefined for a straight line (cross product is zero)
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let tau = curve_torsion_at(&line, 0.5);
        assert!(tau.is_none(), "Torsion should be None for a straight line");
    }

    #[test]
    fn circle_curvature_is_1_over_r() {
        let r = 3.0;
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: r,
        };
        let kappa = curve_curvature_at(&circle, 1.0);
        let expected = 1.0 / r;
        let got = kappa.expect("circle curvature should be defined");
        assert!(
            (got - expected).abs() < 1e-10,
            "Circle curvature: got {got}, expected {expected}"
        );
    }

    #[test]
    fn circle_torsion_is_zero() {
        // A planar circle has zero torsion
        let r = 2.0;
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: r,
        };
        let tau = curve_torsion_at(&circle, 1.0);
        let got = tau.expect("circle torsion should be defined");
        assert!(
            got.abs() < 1e-10,
            "Circle should have zero torsion, got {got}"
        );
    }

    #[test]
    fn circle_curvature_different_radius() {
        let r = 0.5;
        let circle = CurveGeometry::Circle {
            center: Point3::new(1.0, 2.0, 3.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            radius: r,
        };
        let kappa = curve_curvature_at(&circle, 2.5);
        let expected = 1.0 / r;
        let got = kappa.expect("circle curvature should be defined");
        assert!(
            (got - expected).abs() < 1e-10,
            "Circle curvature: got {got}, expected {expected}"
        );
    }
}
