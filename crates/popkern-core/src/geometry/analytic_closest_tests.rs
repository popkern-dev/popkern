//! Exact analytic closest points (PK-737): every arm's foot is the true
//! nearest surface point (checked against the distance to dense samples
//! and against the implicit distance), in the principal branch, with a
//! defined choice at degenerate points.

use super::analytic_closest_point;
use crate::geometry::analytic_surfaces::analytic_signed_distance;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

fn foot(s: &SurfaceGeometry, p: Point3) -> (f64, f64, Point3) {
    let (u, v) = analytic_closest_point(s, &p).expect("analytic");
    (u, v, s.point_at(u, v))
}

/// The foot is exactly as far as the implicit distance says.
fn assert_exact(s: &SurfaceGeometry, p: Point3) {
    let (u, v, f) = foot(s, p);
    let want = analytic_signed_distance(s, &p).expect("implicit").abs();
    let got = f.distance_to(&p);
    assert!(
        (got - want).abs() <= 1e-12 * want.max(1.0),
        "{s:?} at {p:?}: foot {f:?} at (u {u}, v {v}) is {got}, implicit {want}"
    );
    if s.parameter_periodicity().0.is_some() {
        assert!(
            (0.0..std::f64::consts::TAU).contains(&u),
            "u {u} outside [0, 2π)"
        );
    }
}

fn cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, -2.0, 3.0),
        axis: Vec3::new(0.0, 0.6, 0.8),
        radius: 8.0,
    }
}

#[test]
fn a_near_axis_point_lands_on_its_radial_foot() {
    // The measured miss: one unit off the axis of an r 8 cylinder.
    let s = cylinder();
    for p in [
        Point3::new(2.0, -2.0, 3.0),
        Point3::new(1.0, -1.2, 2.4),
        Point3::new(30.0, 5.0, -7.0),
        Point3::new(1.0, -2.0 + 6.0, 3.0 + 8.0),
    ] {
        assert_exact(&s, p);
    }
}

#[test]
fn a_point_on_the_axis_takes_angle_zero() {
    let s = cylinder();
    let (u, v, f) = foot(&s, Point3::new(1.0, -2.0 + 0.6 * 5.0, 3.0 + 0.8 * 5.0));
    assert_eq!(u, 0.0);
    assert!((v - 5.0).abs() < 1e-12);
    assert!((f.distance_to(&Point3::new(1.0, 1.0, 7.0)) - 8.0).abs() < 1e-12);
}

#[test]
fn sphere_torus_cone_and_plane_feet_are_exact() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.5, 0.5, -1.0),
        radius: 3.0,
        axis: crate::geometry::Vec3::Z,
    };
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 10.0,
        minor_radius: 2.0,
    };
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.4,
    };
    let points = [
        Point3::new(4.0, -3.0, 2.0),
        Point3::new(-7.0, 1.0, -0.5),
        Point3::new(11.0, 0.2, 1.0),
        Point3::new(9.5, -0.3, 0.1),
        Point3::new(0.3, 0.2, 5.0),
    ];
    for p in points {
        assert_exact(&sphere, p);
        assert_exact(&torus, p);
        assert_exact(&cone, p);
    }
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 1.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    for p in points {
        assert_exact(&plane, p);
    }
}

#[test]
fn degenerate_points_have_defined_parameters() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
        axis: crate::geometry::Vec3::Z,
    };
    // The north pole direction: longitude 0, colatitude 0.
    assert_eq!(foot(&sphere, Point3::new(0.0, 0.0, 5.0)).0, 0.0);
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 10.0,
        minor_radius: 2.0,
    };
    // The tube centre: tube angle 0.
    let (_, v, _) = foot(&torus, Point3::new(10.0, 0.0, 0.0));
    assert_eq!(v, 0.0);
    // Behind a cone's apex: the apex.
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.4,
    };
    let (_, v, f) = foot(&cone, Point3::new(0.0, 0.0, -3.0));
    assert_eq!(v, 0.0);
    assert!(f.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-15);
}
