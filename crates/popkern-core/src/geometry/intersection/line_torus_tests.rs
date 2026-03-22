//! Tests for line-torus intersection.

#[cfg(test)]
mod tests {
    use crate::geometry::intersection::SurfaceCurveIntersection;
    use crate::geometry::intersection::line_torus::intersect_line_torus;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
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
        assert_eq!(
            r.len(),
            4,
            "expected 4 intersections, got {}: {:?}",
            r.len(),
            r
        );

        // Verify the t values correspond to the expected hit points:
        // Outer ring at x = -(R+r) = -4, inner at x = -(R-r) = -2,
        // then inner at x = R-r = 2, outer at x = R+r = 4.
        let t_values: Vec<f64> = r
            .iter()
            .map(|i| match i {
                SurfaceCurveIntersection::Point { t, .. } => *t,
                _ => panic!("expected Point"),
            })
            .collect();

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
        assert_eq!(
            r.len(),
            2,
            "expected 2 intersections, got {}: {:?}",
            r.len(),
            r
        );

        // Should hit at z = -1 and z = 1 (tube centered at distance R from axis)
        let t_values: Vec<f64> = r
            .iter()
            .map(|i| match i {
                SurfaceCurveIntersection::Point { t, .. } => *t,
                _ => panic!("expected Point"),
            })
            .collect();
        // Line origin z=-5, dir (0,0,1), so t = z + 5
        assert!(
            (t_values[0] - 4.0).abs() < 1e-6,
            "expected t=4, got {}",
            t_values[0]
        );
        assert!(
            (t_values[1] - 6.0).abs() < 1e-6,
            "expected t=6, got {}",
            t_values[1]
        );
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
        assert_eq!(
            r.len(),
            4,
            "expected 4 intersections with Y-axis torus, got {}",
            r.len()
        );
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

        let t_values: Vec<f64> = r
            .iter()
            .map(|i| match i {
                SurfaceCurveIntersection::Point { t, .. } => *t,
                _ => panic!("expected Point"),
            })
            .collect();
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
