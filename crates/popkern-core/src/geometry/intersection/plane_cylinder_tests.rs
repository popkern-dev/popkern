//! Tests for the analytic plane-cylinder intersection solver.

use super::*;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-6)
        .with_angular(1e-6)
        .with_parametric(1e-9)
}

fn z_axis_cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
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
fn perpendicular_plane_cylinder_circle() {
    // z=3 plane through a unit cylinder -> circle of radius 1 at height 3.
    let plane = xy_plane(3.0);
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &tol,
    );

    let intersections = must_some(result, "should match plane-cylinder pair");
    assert_eq!(intersections.len(), 1, "should find exactly one circle");

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed, "circle should be closed");

        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..8 {
            let t = t_min + (t_max - t_min) * (i as f64) / 8.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // Should lie on cylinder (radial distance == 1)
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 1.0).abs() < 1e-4,
                "point at t={t} not on cylinder: radial={radial}"
            );

            // Should lie on z=3 plane.
            assert!(
                (pt.z - 3.0).abs() < 1e-4,
                "point at t={t} not on plane: z={}",
                pt.z
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn parallel_non_intersecting_empty() {
    // YZ plane at x=5, unit cylinder along Z axis -> no intersection.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(5.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &tol,
    );

    let intersections = must_some(result, "should match plane-cylinder pair");
    assert!(
        intersections.is_empty(),
        "plane at x=5 should not intersect unit cylinder"
    );
}

#[test]
fn parallel_intersecting_two_lines() {
    // YZ plane at x=0 (through the cylinder axis), unit cylinder -> two lines.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &tol,
    );

    let intersections = must_some(result, "should match plane-cylinder pair");
    assert_eq!(
        intersections.len(),
        2,
        "should find two lines for x=0 plane through unit cylinder"
    );

    for (idx, isect) in intersections.iter().enumerate() {
        if let SurfaceSurfaceIntersection::Curve(cos) = isect {
            assert!(!cos.is_closed, "lines should not be closed");

            // Sample and verify points lie on both surfaces.
            let (t_min, t_max) = (
                cos.curve_3d.knots[0],
                cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
            );
            for i in 0..5 {
                let t = t_min + (t_max - t_min) * (i as f64) / 4.0;
                let pt = eval_nurbs(&cos.curve_3d, t);

                // On plane: x should be 0.
                assert!(
                    pt.x.abs() < 1e-3,
                    "line {idx} point at t={t} not on plane: x={}",
                    pt.x
                );

                // On cylinder: radial distance should be 1.
                let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
                assert!(
                    (radial - 1.0).abs() < 1e-3,
                    "line {idx} point at t={t} not on cylinder: radial={radial}"
                );
            }
        } else {
            panic!("expected Curve intersection for line {idx}");
        }
    }
}

#[test]
fn oblique_plane_cylinder_ellipse() {
    // A 45-degree tilted plane through the origin cutting a unit cylinder.
    // Normal = (0, 0.707, 0.707), so cos_alpha = 0.707, sin_alpha = 0.707.
    // Semi-minor = 1, semi-major = 1/sin(45) = sqrt(2).
    let n = Vec3::new(0.0, 1.0, 1.0)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let u_axis = Vec3::new(1.0, 0.0, 0.0);
    let v_raw = n.cross(&u_axis);
    let v_axis = v_raw.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));

    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: n,
        u_axis,
        v_axis,
    };
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-5.0, 5.0),
        &tol,
    );

    let intersections = must_some(result, "should match plane-cylinder pair");
    assert_eq!(intersections.len(), 1, "should find exactly one ellipse");

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed, "ellipse should be closed");

        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..16 {
            let t = t_min + (t_max - t_min) * (i as f64) / 16.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // On cylinder: radial distance should be 1.
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 1.0).abs() < 1e-3,
                "point at t={t} not on cylinder: radial={radial}"
            );

            // On plane: (pt - origin) . normal should be ~0.
            let d = (pt - Point3::new(0.0, 0.0, 0.0)).dot(&n);
            assert!(
                d.abs() < 1e-3,
                "point at t={t} not on plane: signed_dist={d}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn samples_lie_on_both_surfaces() {
    // Perpendicular case: z=0 plane, radius=2 cylinder.
    let plane = xy_plane(0.0);
    let cyl = z_axis_cylinder(2.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-5.0, 5.0),
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

            // On plane: z should be 0.
            assert!(pt.z.abs() < 1e-3, "not on plane: z={}", pt.z);

            // On cylinder: radial distance should be 2.
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 2.0).abs() < 1e-3,
                "not on cylinder: radial={radial}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn swapped_order_works() {
    // Cylinder first, plane second.
    let plane = xy_plane(1.0);
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &tol,
    );

    let intersections = must_some(result, "should match cylinder-plane pair");
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed);
        // UV1 should be cylinder UV (linear u sweep at constant v)
        assert_eq!(cos.curve_uv1.degree, 1);
        // UV2 should be plane UV (circle)
        assert_eq!(cos.curve_uv2.degree, 2);
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn returns_none_for_non_plane_cylinder() {
    // Two planes should return None.
    let p1 = xy_plane(0.0);
    let p2 = xy_plane(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
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

#[test]
fn parallel_tangent_is_empty() {
    // YZ plane at x=1, unit cylinder -> tangent (degenerate, empty).
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = z_axis_cylinder(1.0);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &tol,
    );

    let intersections = must_some(result, "should match plane-cylinder pair");
    assert!(
        intersections.is_empty(),
        "tangent plane should yield empty (degenerate)"
    );
}

#[test]
fn oblique_samples_on_both_surfaces_general() {
    // Oblique at 30 degrees: normal = (0, sin30, cos30) = (0, 0.5, 0.866).
    let angle = std::f64::consts::FRAC_PI_6; // 30 degrees
    let n = Vec3::new(0.0, angle.sin(), angle.cos());
    let u_axis = Vec3::new(1.0, 0.0, 0.0);
    let v_raw = n.cross(&u_axis);
    let v_axis = v_raw.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));

    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: n,
        u_axis,
        v_axis,
    };
    let cyl = z_axis_cylinder(1.5);
    let tol = test_tolerance();

    let result = try_plane_cylinder(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &tol,
    );

    let intersections = must_some(result, "should match");
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed);

        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..16 {
            let t = t_min + (t_max - t_min) * (i as f64) / 16.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // On cylinder.
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 1.5).abs() < 1e-2,
                "not on cylinder: radial={radial}"
            );

            // On plane.
            let d = (pt - Point3::new(0.0, 0.0, 0.0)).dot(&n);
            assert!(d.abs() < 1e-2, "not on plane: d={d}");
        }
    } else {
        panic!("expected Curve intersection");
    }
}
