//! Tests for the analytic plane-sphere intersection solver.

use super::*;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::{PI, TAU};

fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-6)
        .with_angular(1e-6)
        .with_parametric(1e-9)
}

fn unit_sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    }
}

fn xy_plane(z: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// Evaluate a NURBS curve at parameter t and return the 3D point.
fn eval_nurbs(curve: &NurbsCurve, t: f64) -> Point3 {
    use crate::geometry::curve::CurveEval;
    curve.point_at(t)
}

/// Unwrap an `Option` in tests without using `.unwrap()` or `.expect()`.
fn must_some<T>(opt: Option<T>, msg: &str) -> T {
    match opt {
        Some(v) => v,
        None => panic!("{msg}"),
    }
}

#[test]
fn plane_sphere_equator() {
    // z=0 plane through unit sphere -> circle of radius 1
    let plane = xy_plane(0.0);
    let sphere = unit_sphere();
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &tol,
    );

    let intersections = must_some(result, "should match plane-sphere pair");
    assert_eq!(intersections.len(), 1, "should find exactly one circle");

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed, "intersection circle should be closed");

        // Sample the 3D curve and verify points lie on both surfaces.
        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..8 {
            let t = t_min + (t_max - t_min) * (i as f64) / 8.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // Should lie on sphere (distance from origin == 1)
            let dist = pt.distance_to(&Point3::new(0.0, 0.0, 0.0));
            assert!(
                (dist - 1.0).abs() < 1e-4,
                "point at t={t} not on sphere: dist={dist}"
            );

            // Should lie on z=0 plane
            assert!(
                pt.z.abs() < 1e-4,
                "point at t={t} not on z=0 plane: z={}",
                pt.z
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn plane_sphere_no_intersection() {
    // z=2 plane, unit sphere -> no intersection
    let plane = xy_plane(2.0);
    let sphere = unit_sphere();
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &tol,
    );

    let intersections = must_some(result, "should match plane-sphere pair");
    assert!(
        intersections.is_empty(),
        "plane at z=2 should not intersect unit sphere"
    );
}

#[test]
fn plane_sphere_tangent() {
    // z=1 plane, unit sphere -> tangent (degenerate, treated as empty)
    let plane = xy_plane(1.0);
    let sphere = unit_sphere();
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &tol,
    );

    let intersections = must_some(result, "should match plane-sphere pair");
    assert!(
        intersections.is_empty(),
        "tangent plane should yield empty (degenerate)"
    );
}

#[test]
fn plane_sphere_offset() {
    // z=0.5 plane, unit sphere -> circle of radius sqrt(0.75)
    let plane = xy_plane(0.5);
    let sphere = unit_sphere();
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &tol,
    );

    let intersections = must_some(result, "should match plane-sphere pair");
    assert_eq!(intersections.len(), 1);

    let expected_radius = (0.75_f64).sqrt();

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        // Sample the 3D curve and check radius from z-axis and z-value.
        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..8 {
            let t = t_min + (t_max - t_min) * (i as f64) / 8.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // Distance from z-axis should equal expected_radius
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - expected_radius).abs() < 1e-4,
                "radial distance {radial} != expected {expected_radius}"
            );

            // z-coordinate should be 0.5
            assert!((pt.z - 0.5).abs() < 1e-4, "z={} should be 0.5", pt.z);

            // Point should lie on unit sphere
            let dist = pt.distance_to(&Point3::new(0.0, 0.0, 0.0));
            assert!(
                (dist - 1.0).abs() < 1e-4,
                "point not on sphere: dist={dist}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn plane_sphere_swapped_order() {
    // Sphere first, plane second -- should still work.
    let plane = xy_plane(0.0);
    let sphere = unit_sphere();
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &tol,
    );

    let intersections = must_some(result, "should match sphere-plane pair");
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed);

        // UV1 should be sphere UV (linear u sweep at constant v)
        assert_eq!(cos.curve_uv1.degree, 1);
        // UV2 should be plane UV (circle)
        assert_eq!(cos.curve_uv2.degree, 2);
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn plane_sphere_samples_on_both_surfaces() {
    // Verify samples on the result curve lie on both surfaces within tolerance.
    let plane = xy_plane(0.3);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &tol,
    );

    let intersections = must_some(result, "should match");
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..16 {
            let t = t_min + (t_max - t_min) * (i as f64) / 16.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // On plane: z should be 0.3
            assert!((pt.z - 0.3).abs() < 1e-3, "point not on plane: z={}", pt.z);

            // On sphere: distance from origin should be 2.0
            let dist = pt.distance_to(&Point3::new(0.0, 0.0, 0.0));
            assert!(
                (dist - 2.0).abs() < 1e-3,
                "point not on sphere: dist={dist}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn returns_none_for_non_plane_sphere() {
    // Two planes should return None (not handled by this solver).
    let p1 = xy_plane(0.0);
    let p2 = xy_plane(1.0);
    let tol = test_tolerance();

    let result = try_plane_sphere(
        &p1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &p2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &tol,
    );
    assert!(result.is_none(), "two planes should not match");
}
