//! Tests for the analytic plane-plane intersection solver.

#[cfg(test)]
mod tests {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::intersection::plane_plane::{solve_2x2, try_plane_plane};
    use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    fn xy_plane_params() -> (Point3, Vec3, Vec3, Vec3) {
        (
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
    }

    fn xz_plane_params() -> (Point3, Vec3, Vec3, Vec3) {
        (
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
    }

    fn default_tol() -> Tolerance {
        Tolerance::default()
    }

    /// Two perpendicular planes (XY and XZ) should produce a line along the X-axis.
    #[test]
    fn perpendicular_planes_give_line() {
        let (o1, n1, u1, v1) = xy_plane_params();
        let (o2, n2, u2, v2) = xz_plane_params();
        let range = (-2.0_f64, 2.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            o2,
            n2,
            u2,
            v2,
            range,
            range,
            &default_tol(),
        )
        .unwrap();
        assert_eq!(result.len(), 1, "Expected exactly one intersection curve");
        assert!(
            matches!(&result[0], SurfaceSurfaceIntersection::Curve(_)),
            "Expected Curve variant"
        );
    }

    /// Parallel planes separated by distance 10 should produce no intersection.
    #[test]
    fn parallel_separated_planes_give_empty() {
        let (o1, n1, u1, v1) = xy_plane_params();
        // Same orientation but shifted by 10 in Z — well beyond linear tolerance.
        let o2 = Point3::new(0.0, 0.0, 10.0);
        let range = (0.0_f64, 1.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            o2,
            n1,
            u1,
            v1,
            range,
            range,
            &default_tol(),
        )
        .unwrap();
        assert!(
            result.is_empty(),
            "Parallel separated planes should not intersect"
        );
    }

    /// Coincident planes (same plane) should produce a CoincidentRegion.
    #[test]
    fn coincident_planes_give_region() {
        let (o1, n1, u1, v1) = xy_plane_params();
        let range = (0.0_f64, 1.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            &default_tol(),
        )
        .unwrap();
        assert_eq!(result.len(), 1);
        assert!(matches!(
            &result[0],
            SurfaceSurfaceIntersection::CoincidentRegion { .. }
        ));
    }

    /// Two planes at 45° should intersect with a line whose direction is along the X-axis.
    #[test]
    fn planes_at_45_degrees() {
        // Plane 1: XY plane (z = 0), normal = (0, 0, 1)
        let o1 = Point3::new(0.0, 0.0, 0.0);
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let u1 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(0.0, 1.0, 0.0);

        // Plane 2: tilted 45° around the X-axis so normal = (0, -sin45, cos45).
        let s = 1.0_f64 / 2.0_f64.sqrt();
        let n2 = Vec3::new(0.0, -s, s);
        let u2 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, s, s); // v2 ⊥ n2, in YZ plane
        let o2 = Point3::new(0.0, 0.0, 0.0);

        let range = (-3.0_f64, 3.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            o2,
            n2,
            u2,
            v2,
            range,
            range,
            &default_tol(),
        )
        .unwrap();

        assert_eq!(result.len(), 1, "Expected one intersection curve");
        if let SurfaceSurfaceIntersection::Curve(cos) = &result[0] {
            // The intersection direction should be (±1, 0, 0) (along X-axis).
            let p_start = cos.curve_3d.point_at(cos.curve_3d.parameter_range().0);
            let p_end = cos.curve_3d.point_at(cos.curve_3d.parameter_range().1);
            let tangent = p_end - p_start;
            let tangent_norm = tangent.normalized().unwrap();
            // Direction should be along X.
            assert!(
                tangent_norm.y.abs() < 1e-8 && tangent_norm.z.abs() < 1e-8,
                "Tangent should be along X, got {:?}",
                tangent_norm
            );
        } else {
            panic!("Expected Curve variant");
        }
    }

    /// Verify that UV curves produce correct projections of known 3D points.
    #[test]
    fn uv_curves_are_correct() {
        // XY plane (surface 1): origin at (0,0,0), u=X, v=Y
        let o1 = Point3::new(0.0, 0.0, 0.0);
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let u1 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(0.0, 1.0, 0.0);

        // XZ plane (surface 2): origin at (0,0,0), u=X, v=Z
        let o2 = Point3::new(0.0, 0.0, 0.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let u2 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);

        let range = (-1.0_f64, 1.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            range,
            range,
            o2,
            n2,
            u2,
            v2,
            range,
            range,
            &default_tol(),
        )
        .unwrap();

        assert_eq!(result.len(), 1);
        if let SurfaceSurfaceIntersection::Curve(cos) = &result[0] {
            let (t_lo, t_hi) = cos.curve_3d.parameter_range();

            // Sample at t_lo, midpoint, and t_hi.
            for frac in [0.0_f64, 0.5, 1.0] {
                let t = t_lo + frac * (t_hi - t_lo);

                let p3d = cos.curve_3d.point_at(t);
                let uv1 = cos.curve_uv1.point_at(t);
                let uv2 = cos.curve_uv2.point_at(t);

                // p3d should lie on surface 1 (z ≈ 0).
                assert!(p3d.z.abs() < 1e-8, "3D point z should be 0: got {}", p3d.z);

                // UV1: u = p3d.x, v = p3d.y (since u1=X, v1=Y, origin at 0).
                assert!(
                    (uv1.x - p3d.x).abs() < 1e-8,
                    "UV1.u mismatch: {} vs {}",
                    uv1.x,
                    p3d.x
                );
                assert!(
                    (uv1.y - p3d.y).abs() < 1e-8,
                    "UV1.v mismatch: {} vs {}",
                    uv1.y,
                    p3d.y
                );

                // p3d should lie on surface 2 (y ≈ 0).
                assert!(p3d.y.abs() < 1e-8, "3D point y should be 0: got {}", p3d.y);

                // UV2: u = p3d.x, v = p3d.z (since u2=X, v2=Z, origin at 0).
                assert!(
                    (uv2.x - p3d.x).abs() < 1e-8,
                    "UV2.u mismatch: {} vs {}",
                    uv2.x,
                    p3d.x
                );
                assert!(
                    (uv2.y - p3d.z).abs() < 1e-8,
                    "UV2.v mismatch: {} vs {}",
                    uv2.y,
                    p3d.z
                );
            }
        } else {
            panic!("Expected Curve variant");
        }
    }

    /// Planes that intersect within a UV box offset from the origin in the u direction.
    ///
    /// The XY ∩ XZ intersection is the X-axis (y=0, z=0).
    /// Surface 1 (XY): u=x, v=y → v is always 0 on the line.
    /// Surface 2 (XZ): u=x, v=z → v is always 0 on the line.
    /// Use u_range [10, 11] and v_range [-1, 1] so the line is within the box.
    #[test]
    fn line_within_offset_uv_box() {
        let o1 = Point3::new(0.0, 0.0, 0.0);
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let u1 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(0.0, 1.0, 0.0);

        let o2 = Point3::new(0.0, 0.0, 0.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let u2 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);

        // u offset to [10, 11]; v contains 0 so the intersection line passes through.
        let u_range = (10.0_f64, 11.0_f64);
        let v_range = (-1.0_f64, 1.0_f64);
        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            u_range,
            v_range,
            o2,
            n2,
            u2,
            v2,
            u_range,
            v_range,
            &default_tol(),
        )
        .unwrap();
        // The intersection line (X-axis, x ∈ [10, 11], y=0, z=0) is within both UV boxes.
        assert_eq!(
            result.len(),
            1,
            "Should find intersection within offset UV box"
        );
    }

    /// Intersection line that does not reach the UV box should give empty result.
    #[test]
    fn line_outside_uv_box_gives_empty() {
        // XY plane, UV range [100, 101] x [0, 1]
        let o1 = Point3::new(0.0, 0.0, 0.0);
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let u1 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(0.0, 1.0, 0.0);

        // YZ plane (x=0), UV range [0, 1] x [0, 1].
        // Intersection is the Y-axis (x=0, z=0), which does NOT overlap with u1 range [100,101].
        let o2 = Point3::new(0.0, 0.0, 0.0);
        let n2 = Vec3::new(1.0, 0.0, 0.0);
        let u2 = Vec3::new(0.0, 1.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);

        let result = try_plane_plane(
            o1,
            n1,
            u1,
            v1,
            (100.0, 101.0),
            (0.0, 1.0),
            o2,
            n2,
            u2,
            v2,
            (0.0, 1.0),
            (0.0, 1.0),
            &default_tol(),
        )
        .unwrap();
        assert!(
            result.is_empty(),
            "Intersection line outside UV box should give empty result"
        );
    }

    /// Solve 2×2 helper: verify basic correctness.
    #[test]
    fn solve_2x2_basic() {
        // 1*x + 0*y = 3  => x=3
        // 0*x + 1*y = 7  => y=7
        let (x, y) = solve_2x2(1.0, 0.0, 3.0, 0.0, 1.0, 7.0).unwrap();
        assert!((x - 3.0).abs() < 1e-10);
        assert!((y - 7.0).abs() < 1e-10);
    }

    /// Degenerate 2×2 system should return an error.
    #[test]
    fn solve_2x2_degenerate() {
        // Both rows identical: 1*x + 1*y = 1 (singular)
        let result = solve_2x2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(result.is_err());
    }

    /// Verify the fast path is exercised from the public API (intersect_surfaces).
    #[test]
    fn public_api_two_perpendicular_planes() {
        use crate::geometry::intersection::surface_surface::{
            SurfaceSurfaceIntersection, intersect_surfaces,
        };
        use crate::geometry::surface::SurfaceGeometry;

        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let result = intersect_surfaces(
            &s1,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &s2,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &default_tol(),
        )
        .unwrap();
        assert_eq!(result.intersections.len(), 1);
        assert!(matches!(
            &result.intersections[0],
            SurfaceSurfaceIntersection::Curve(_)
        ));
    }

    /// Fast path: parallel separated planes via public API.
    #[test]
    fn public_api_parallel_separated_planes() {
        use crate::geometry::intersection::surface_surface::intersect_surfaces;
        use crate::geometry::surface::SurfaceGeometry;

        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = intersect_surfaces(
            &s1,
            (0.0, 1.0),
            (0.0, 1.0),
            &s2,
            (0.0, 1.0),
            (0.0, 1.0),
            &default_tol(),
        )
        .unwrap();
        assert!(result.intersections.is_empty());
    }
}
