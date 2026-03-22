//! Revolution geometry helpers — point rotation, revolution curves.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::TAU;

/// Rotate a point around an axis using Rodrigues' rotation formula.
pub fn rotate_point(point: &Point3, axis_origin: &Point3, axis_dir: &Vec3, angle: f64) -> Point3 {
    let v = *point - *axis_origin;
    let k = axis_dir.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let rotated = v * angle.cos() + k.cross(&v) * angle.sin() + k * k.dot(&v) * (1.0 - angle.cos());
    *axis_origin + rotated
}

/// Create a revolution curve for a point revolving around an axis.
/// Returns Circle for full 360°, Arc for partial.
pub fn revolution_curve(
    point: &Point3,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    angle: f64,
) -> CurveGeometry {
    let axis = axis_dir.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let v = *point - *axis_origin;
    let proj = axis * v.dot(&axis);
    let center = *axis_origin + proj;
    let radial = v - proj;
    let radius = radial.length();

    if (angle - TAU).abs() < 1e-10 {
        CurveGeometry::Circle { center, normal: axis, radius }
    } else {
        CurveGeometry::Arc {
            center,
            normal: axis,
            radius,
            start_angle: 0.0,
            end_angle: angle,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    #[test]
    fn rotate_90_degrees() {
        let p = Point3::new(1.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(0.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), FRAC_PI_2);
        assert!((r.x).abs() < 1e-10);
        assert!((r.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rotate_360_returns_to_start() {
        let p = Point3::new(1.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(0.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), TAU);
        assert!((r.x - 1.0).abs() < 1e-10);
        assert!((r.y).abs() < 1e-10);
    }

    #[test]
    fn rotate_around_offset_axis() {
        let p = Point3::new(2.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(1.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), PI);
        assert!((r.x).abs() < 1e-10);
        assert!((r.y).abs() < 1e-10);
    }

    #[test]
    fn revolution_curve_full_is_circle() {
        let c = revolution_curve(
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        assert!(matches!(c, CurveGeometry::Circle { .. }));
    }

    #[test]
    fn revolution_curve_partial_is_arc() {
        let c = revolution_curve(
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        );
        assert!(matches!(c, CurveGeometry::Arc { .. }));
    }

    #[test]
    fn revolution_curve_radius_correct() {
        let c = revolution_curve(
            &Point3::new(3.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        if let CurveGeometry::Circle { radius, .. } = c {
            assert!((radius - 3.0).abs() < 1e-10);
        } else {
            panic!("Expected Circle");
        }
    }

    #[test]
    fn revolution_curve_offset_point() {
        // Point at (2, 0, 3), axis is Z through origin. Should orbit at radius 2, center at (0,0,3)
        let c = revolution_curve(
            &Point3::new(2.0, 0.0, 3.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        if let CurveGeometry::Circle { center, radius, .. } = c {
            assert!((center.z - 3.0).abs() < 1e-10);
            assert!((radius - 2.0).abs() < 1e-10);
        } else {
            panic!("Expected Circle");
        }
    }
}
