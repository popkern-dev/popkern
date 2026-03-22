#[cfg(test)]
use super::*;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[test]
fn invalid_t_range_fails() {
    let s = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let c = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(
        intersect_surface_curve(
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &c,
            (1.0, 0.0),
            &Tolerance::default()
        )
        .is_err()
    );
}

// -----------------------------------------------------------------------
// dispatch_line_cylinder: verify analytic path
// -----------------------------------------------------------------------
#[test]
fn dispatch_line_cylinder() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let line = CurveGeometry::Line {
        origin: Point3::new(2.0, 0.0, 0.5),
        direction: Vec3::new(-1.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();
    let result = intersect_surface_curve(
        &cyl,
        (0.0, std::f64::consts::TAU),
        (-5.0, 5.0),
        &line,
        (0.0, 4.0),
        &tol,
    )
    .unwrap();
    let pts: Vec<_> = result
        .intersections
        .iter()
        .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
        .collect();
    assert!(
        pts.len() >= 2,
        "Expected 2 intersections for line through cylinder, got {}",
        pts.len()
    );
}

// -----------------------------------------------------------------------
// dispatch_line_cone: verify analytic path
// -----------------------------------------------------------------------
#[test]
fn dispatch_line_cone() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let line = CurveGeometry::Line {
        origin: Point3::new(2.0, 0.0, 1.0),
        direction: Vec3::new(-1.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();
    let result = intersect_surface_curve(
        &cone,
        (0.0, std::f64::consts::TAU),
        (0.0, 10.0),
        &line,
        (0.0, 4.0),
        &tol,
    )
    .unwrap();
    let pts: Vec<_> = result
        .intersections
        .iter()
        .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
        .collect();
    assert!(
        pts.len() >= 1,
        "Expected at least 1 intersection for line through cone, got {}",
        pts.len()
    );
}

// -----------------------------------------------------------------------
// dispatch_line_torus: verify analytic path
// -----------------------------------------------------------------------
#[test]
fn dispatch_line_torus() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    let line = CurveGeometry::Line {
        origin: Point3::new(5.0, 0.0, 0.0),
        direction: Vec3::new(-1.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();
    let result = intersect_surface_curve(
        &torus,
        (0.0, std::f64::consts::TAU),
        (0.0, std::f64::consts::TAU),
        &line,
        (0.0, 10.0),
        &tol,
    )
    .unwrap();
    let pts: Vec<_> = result
        .intersections
        .iter()
        .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
        .collect();
    assert!(
        pts.len() >= 2,
        "Expected at least 2 intersections for line through torus, got {}",
        pts.len()
    );
}
