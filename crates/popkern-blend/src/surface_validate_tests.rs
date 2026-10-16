//! Unit tests for [`crate::surface_validate::validate_fillet_surface`]
//! — the §13.2 post-construction self-intersection detector.

use super::*;
use popkern_core::geometry::Vec3;
use std::f64::consts::{PI, TAU};

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// primitive-level tests
// ---------------------------------------------------------------------------

#[test]
fn valid_torus_passes_validation() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    assert!(validate_fillet_surface(&torus, (0.0, 1.0), (0.0, 1.0), &unit_tol()).is_ok());
}

#[test]
fn spindle_torus_is_rejected() {
    // minor > major → spindle torus, inner tube self-intersects.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 0.3,
        minor_radius: 0.5,
    };
    let result = validate_fillet_surface(&torus, (0.0, 1.0), (0.0, 1.0), &unit_tol());
    match result {
        Err(SurfaceValidationIssue::SpindleTorus {
            major_radius,
            minor_radius,
        }) => {
            assert!((major_radius - 0.3).abs() < 1e-12);
            assert!((minor_radius - 0.5).abs() < 1e-12);
        }
        other => panic!("expected SpindleTorus, got {other:?}"),
    }
}

#[test]
fn torus_with_minor_equals_major_passes() {
    // Edge case: minor == major is a horn torus (tangent to itself
    // at a single circle). We accept it as the boundary of the
    // valid region; reject only strict inequality.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 0.5,
        minor_radius: 0.5,
    };
    assert!(validate_fillet_surface(&torus, (0.0, 1.0), (0.0, 1.0), &unit_tol()).is_ok());
}

#[test]
fn plane_passes_validation() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(validate_fillet_surface(&plane, (-1.0, 1.0), (-1.0, 1.0), &unit_tol()).is_ok());
}

#[test]
fn cylinder_passes_validation() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert!(validate_fillet_surface(&cyl, (0.0, TAU), (-1.0, 1.0), &unit_tol()).is_ok());
}

#[test]
fn sphere_passes_validation() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    assert!(validate_fillet_surface(&sphere, (0.0, TAU), (0.0, PI), &unit_tol()).is_ok());
}

#[test]
fn cone_passes_validation() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.5_f64.atan(),
    };
    assert!(validate_fillet_surface(&cone, (0.0, TAU), (0.0, 2.0), &unit_tol()).is_ok());
}

// ---------------------------------------------------------------------------
// NURBS surface tests
// ---------------------------------------------------------------------------

/// Build a small flat NURBS sheet for sampling. Degree-1 in both
/// directions, 3×3 control points, knots [0, 0, 0.5, 1, 1] on each
/// axis — a 1×1 unit square in the z = 0 plane. This surface has no
/// twist and no self-intersection; validation should pass cleanly.
fn flat_nurbs_sheet() -> SurfaceGeometry {
    use popkern_core::geometry::NurbsSurface;
    let n = 3;
    let mut cps = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            cps.push(Point3::new(i as f64 * 0.5, j as f64 * 0.5, 0.0));
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: cps,
        num_u: n,
        num_v: n,
        weights: vec![1.0; n * n],
    })
}

#[test]
fn flat_nurbs_sheet_passes_validation() {
    let sheet = flat_nurbs_sheet();
    assert!(
        validate_fillet_surface(&sheet, (0.0, 1.0), (0.0, 1.0), &unit_tol()).is_ok(),
        "a flat NURBS sheet should pass §13.2 validation"
    );
}
