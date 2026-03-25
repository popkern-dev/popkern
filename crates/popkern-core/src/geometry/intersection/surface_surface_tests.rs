use super::*;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{PI, TAU};

/// A coarser tolerance so marching takes fewer, larger steps in tests.
fn test_tolerance() -> Tolerance {
    Tolerance {
        linear: 1e-3,
        angular: 1e-3,
        parametric: 1e-6,
    }
}

#[test]
fn two_perpendicular_planes_one_curve() {
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
        &test_tolerance(),
    )
    .unwrap();
    assert!(!result.intersections.is_empty(), "Should find intersection");
    assert!(matches!(
        &result.intersections[0],
        SurfaceSurfaceIntersection::Curve(_)
    ));
}

#[test]
fn coincident_planes_region() {
    let s = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = intersect_surfaces(
        &s,
        (0.0, 1.0),
        (0.0, 1.0),
        &s,
        (0.0, 1.0),
        (0.0, 1.0),
        &test_tolerance(),
    )
    .unwrap();
    assert!(
        result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. }))
    );
}

#[test]
fn non_intersecting_empty() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 10.0),
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
        &test_tolerance(),
    )
    .unwrap();
    assert!(result.intersections.is_empty());
}

#[test]
fn coincident_cylinders_return_region_not_hang() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.5, 0.5, -0.5),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.15,
    };
    let tol = Tolerance {
        linear: 1e-3,
        angular: 1e-3,
        parametric: 1e-6,
    };
    let result = intersect_surfaces(
        &cyl,
        (0.0, TAU),
        (-1.0, 2.0),
        &cyl,
        (0.0, TAU),
        (-1.0, 2.0),
        &tol,
    )
    .unwrap();
    assert!(
        result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. })),
        "identical cylinders should produce a CoincidentRegion, got {:?}",
        result.intersections
    );
}

#[test]
fn coincident_planes_return_region() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = test_tolerance();
    let result = intersect_surfaces(
        &plane,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &tol,
    )
    .unwrap();
    assert!(
        result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. })),
        "identical planes should produce a CoincidentRegion"
    );
}

#[test]
fn plane_and_sphere() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let result = intersect_surfaces(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &test_tolerance(),
    )
    .unwrap();
    // Should find an intersection curve (circle of radius 1 in XY plane)
    assert!(
        !result.intersections.is_empty(),
        "Should find plane-sphere intersection"
    );
}
