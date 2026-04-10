//! Tests for the analytic plane-cone intersection solver.

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

/// A cone along +Z with apex at origin, half_angle = 45 degrees.
fn z_axis_cone() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    }
}

/// A cone matching the constructor pattern: apex at (0,0,height), axis=(0,0,-1).
fn downward_cone(height: f64, radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, height),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: (radius / height).atan(),
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

/// Check whether a 3D point lies on a cone surface.
fn point_on_cone(pt: Point3, apex: Point3, axis: Vec3, half_angle: f64, tol: f64) -> bool {
    let d = pt - apex;
    let along_axis = d.dot(&axis);
    if along_axis < -tol {
        return false; // behind the apex
    }
    let radial_sq = d.dot(&d) - along_axis * along_axis;
    let expected_r = along_axis * half_angle.tan();
    (radial_sq.sqrt() - expected_r.abs()).abs() < tol
}

#[test]
fn perpendicular_plane_cone_circle() {
    // z=2 plane through a 45-degree cone with apex at origin.
    // At v=2, radius = 2*tan(45) = 2.
    let plane = xy_plane(2.0);
    let cone = z_axis_cone();
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone,
        (0.0, TAU),
        (0.0, 10.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should match plane-cone pair"),
    };
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

            // Should lie on cone (radial distance == 2 at z=2).
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - 2.0).abs() < 1e-4,
                "point at t={t} not on cone: radial={radial}"
            );

            // Should lie on z=2 plane.
            assert!(
                (pt.z - 2.0).abs() < 1e-4,
                "point at t={t} not on plane: z={}",
                pt.z
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn perpendicular_downward_cone_circle() {
    // Downward cone: apex=(0,0,5), axis=(0,0,-1), half_angle=atan(1/5).
    // Plane at z=3 intersects at v=2 (distance from apex along axis).
    // Radius at v=2: 2*tan(atan(1/5)) = 2/5 = 0.4.
    let cone = downward_cone(5.0, 1.0);
    let plane = xy_plane(3.0);
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone,
        (0.0, TAU),
        (0.0, 10.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should match plane-cone pair"),
    };
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed);

        let expected_radius = 2.0 * (1.0_f64 / 5.0).atan().tan();
        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..8 {
            let t = t_min + (t_max - t_min) * (i as f64) / 8.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - expected_radius).abs() < 1e-4,
                "radial={radial}, expected={expected_radius}"
            );
            assert!((pt.z - 3.0).abs() < 1e-4, "z={}", pt.z);
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn oblique_plane_cone_ellipse() {
    // A tilted plane cutting the cone obliquely.
    // Cone: apex at origin, axis=+Z, half_angle=30 degrees.
    // Plane: normal tilted 20 degrees from Z, at z=3.
    let ha = std::f64::consts::FRAC_PI_6; // 30 degrees
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: ha,
    };

    let tilt = 20.0_f64.to_radians();
    let n = Vec3::new(0.0, tilt.sin(), tilt.cos());
    let u_axis = Vec3::new(1.0, 0.0, 0.0);
    let v_raw = n.cross(&u_axis);
    let v_axis = v_raw.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));

    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 3.0),
        normal: n,
        u_axis,
        v_axis,
    };
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &cone,
        (0.0, TAU),
        (0.0, 20.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should find elliptic intersection"),
    };
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

            // On cone check.
            assert!(
                point_on_cone(
                    pt,
                    Point3::new(0.0, 0.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    ha,
                    0.1
                ),
                "point at t={t} not on cone: {:?}",
                pt
            );

            // On plane check.
            let d = (pt - Point3::new(0.0, 0.0, 3.0)).dot(&n);
            assert!(
                d.abs() < 0.1,
                "point at t={t} not on plane: signed_dist={d}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn no_intersection_plane_below_apex() {
    // Plane at z=-1, cone apex at origin with axis +Z.
    // The cone only exists for v > 0 (z > 0), so no intersection.
    let plane = xy_plane(-1.0);
    let cone = z_axis_cone();
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone,
        (0.0, TAU),
        (0.0, 10.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should match plane-cone pair"),
    };
    assert!(
        intersections.is_empty(),
        "plane below apex should yield no intersection"
    );
}

#[test]
fn samples_lie_on_both_surfaces() {
    // Perpendicular case: z=5 plane, cone apex at origin, ha=30deg.
    let ha = std::f64::consts::FRAC_PI_6;
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: ha,
    };
    let plane = xy_plane(5.0);
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &cone,
        (0.0, TAU),
        (0.0, 20.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should match"),
    };
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        let expected_r = 5.0 * ha.tan();
        let (t_min, t_max) = (
            cos.curve_3d.knots[0],
            cos.curve_3d.knots[cos.curve_3d.knots.len() - 1],
        );
        for i in 0..16 {
            let t = t_min + (t_max - t_min) * (i as f64) / 16.0;
            let pt = eval_nurbs(&cos.curve_3d, t);

            // On plane.
            assert!((pt.z - 5.0).abs() < 1e-3, "not on plane: z={}", pt.z);

            // On cone.
            let radial = (pt.x * pt.x + pt.y * pt.y).sqrt();
            assert!(
                (radial - expected_r).abs() < 1e-3,
                "not on cone: radial={radial}, expected={expected_r}"
            );
        }
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn swapped_ordering_works() {
    // Cone first, plane second.
    let plane = xy_plane(2.0);
    let cone = z_axis_cone();
    let tol = test_tolerance();

    let result = try_plane_cone(
        &cone,
        (0.0, TAU),
        (0.0, 10.0),
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &tol,
    );

    let intersections = match result {
        Some(v) => v,
        None => panic!("should match cone-plane pair"),
    };
    assert_eq!(intersections.len(), 1);

    if let SurfaceSurfaceIntersection::Curve(ref cos) = intersections[0] {
        assert!(cos.is_closed);
        // UV1 should be cone UV, UV2 should be plane UV.
        // Cone UV for perpendicular circle: degree-2 rational (projected).
        // Plane UV for perpendicular circle: degree-2 rational.
        assert_eq!(cos.curve_uv1.degree, 2);
        assert_eq!(cos.curve_uv2.degree, 2);
    } else {
        panic!("expected Curve intersection");
    }
}

#[test]
fn degenerate_plane_through_apex_returns_none() {
    // Plane at z=0 passes through the apex at origin.
    let plane = xy_plane(0.0);
    let cone = z_axis_cone();
    let tol = test_tolerance();

    let result = try_plane_cone(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone,
        (0.0, TAU),
        (0.0, 10.0),
        &tol,
    );

    // Should return None (degenerate, falls through to marching).
    assert!(
        result.is_none(),
        "plane through apex should return None for marching fallback"
    );
}

#[test]
fn returns_none_for_non_plane_cone() {
    // Two planes should return None.
    let p1 = xy_plane(0.0);
    let p2 = xy_plane(1.0);
    let tol = test_tolerance();

    let result = try_plane_cone(
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
