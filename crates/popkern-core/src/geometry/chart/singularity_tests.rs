use super::*;
use crate::geometry::{Point3, Vec3};

fn bilinear(corners: [Point3; 4]) -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: corners.to_vec(),
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    })
}

#[test]
fn poles_and_apexes_are_singular_rows() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 3.0,
        axis: Vec3::Z,
    };
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::Z,
        half_angle: 0.3,
    };
    assert!(has_singular_row(&sphere, 1e-6));
    assert!(has_singular_row(&cone, 1e-6));
}

#[test]
fn a_torus_is_singular_only_when_its_tube_reaches_the_axis() {
    let torus = |minor| SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::Z,
        major_radius: 2.0,
        minor_radius: minor,
    };
    assert!(!has_singular_row(&torus(0.5), 1e-6));
    assert!(has_singular_row(&torus(2.0), 1e-6));
    assert!(has_singular_row(&torus(2.5), 1e-6));
}

#[test]
fn planes_and_cylinders_have_none() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::Z,
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::Z,
        radius: 0.1,
    };
    assert!(!has_singular_row(&plane, 1e-6));
    assert!(!has_singular_row(&cylinder, 1e-6));
}

#[test]
fn a_nurbs_row_of_coincident_points_is_singular() {
    let p = |x, y| Point3::new(x, y, 0.0);
    // Grid (u, v): (0,0) (0,1) (1,0) (1,1), v fastest.
    assert!(!has_singular_row(
        &bilinear([p(0.0, 0.0), p(0.0, 1.0), p(1.0, 0.0), p(1.0, 1.0)]),
        1e-6
    ));
    // The u = 0 row collapsed to one point (a triangle).
    assert!(has_singular_row(
        &bilinear([p(0.0, 0.0), p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0)]),
        1e-6
    ));
    // The v = 1 column collapsed within the tolerance.
    assert!(has_singular_row(
        &bilinear([p(0.0, 0.0), p(1.0, 1.0), p(1.0, 0.0), p(1.0, 1.0 + 5e-7)]),
        1e-6
    ));
    // Apart by more than it, not.
    assert!(!has_singular_row(
        &bilinear([p(0.0, 0.0), p(1.0, 1.0), p(1.0, 0.0), p(1.0, 1.0 + 1e-3)]),
        1e-6
    ));
}
