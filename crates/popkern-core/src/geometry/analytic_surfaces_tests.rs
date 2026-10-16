use super::*;
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn plane_point_at_origin() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p = plane.point_at(0.0, 0.0);
    assert!((p.x).abs() < 1e-10);
    assert!((p.y).abs() < 1e-10);
}

#[test]
fn plane_point_at_offset() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p = plane.point_at(3.0, 4.0);
    assert!((p.x - 3.0).abs() < 1e-10);
    assert!((p.y - 4.0).abs() < 1e-10);
}

#[test]
fn plane_normal_is_constant() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let n = plane.normal_at(5.0, 10.0);
    assert!((n.x).abs() < 1e-10);
    assert!((n.y).abs() < 1e-10);
    assert!(n.z > 0.0);
}

#[test]
fn sphere_point_at_north_pole() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let p = sphere.point_at(0.0, 0.0);
    let dist = ((p.x).powi(2) + (p.y).powi(2) + (p.z).powi(2)).sqrt();
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn sphere_all_points_on_sphere() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: 5.0,
        axis: crate::geometry::Vec3::Z,
    };
    for u in [0.0, 1.0, PI, 5.0] {
        for v in [0.0, 0.5, FRAC_PI_2, PI] {
            let p = sphere.point_at(u, v);
            let dist = ((p.x - 1.0).powi(2) + (p.y - 2.0).powi(2) + (p.z - 3.0).powi(2)).sqrt();
            assert!(
                (dist - 5.0).abs() < 1e-10,
                "Point at u={u}, v={v} not on sphere: dist={dist}"
            );
        }
    }
}

#[test]
fn cylinder_point_on_surface() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let p = cyl.point_at(0.0, 5.0);
    let radial_dist = ((p.x).powi(2) + (p.y).powi(2)).sqrt();
    assert!((radial_dist - 2.0).abs() < 1e-10);
    assert!((p.z - 5.0).abs() < 1e-10);
}

#[test]
fn torus_point_on_surface() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    let p = torus.point_at(0.0, 0.0);
    assert!((p.x - 4.0).abs() < 1e-10);
    assert!((p.y).abs() < 1e-10);
    assert!((p.z).abs() < 1e-10);
}

#[test]
fn plane_partials_and_range() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!((plane.partial_u(0.0, 0.0).x - 1.0).abs() < 1e-10);
    assert!((plane.partial_v(0.0, 0.0).y - 1.0).abs() < 1e-10);
    let ((u0, u1), (v0, v1)) = plane.parameter_range();
    assert!(u0.is_infinite() && u0 < 0.0);
    assert!(u1.is_infinite() && u1 > 0.0);
    assert!(v0.is_infinite() && v0 < 0.0);
    assert!(v1.is_infinite() && v1 > 0.0);
}

#[test]
fn cylinder_partials_normal_range() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    assert!((cyl.partial_u(0.0, 0.0).length() - 2.0).abs() < 1e-10);
    assert!((cyl.partial_v(0.0, 0.0).z - 1.0).abs() < 1e-10);
    assert!(cyl.normal_at(0.0, 0.0).length() > 1e-10);
    let ((u0, u1), (v0, v1)) = cyl.parameter_range();
    assert!((u0).abs() < 1e-10 && (u1 - TAU).abs() < 1e-10);
    assert!(v0.is_infinite() && v1.is_infinite());
}

#[test]
fn cone_partials_normal_point_range() {
    let ha = FRAC_PI_2 / 4.0;
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: ha,
    };
    assert!(cone.partial_u(0.0, 1.0).length() > 1e-10);
    assert!(cone.partial_v(0.0, 1.0).length() > 1e-10);
    assert!(cone.normal_at(0.0, 1.0).length() > 1e-10);
    let p = cone.point_at(0.0, 1.0);
    let radial = ((p.x).powi(2) + (p.y).powi(2)).sqrt();
    assert!((radial - ha.tan()).abs() < 1e-10);
    assert!((p.z - 1.0).abs() < 1e-10);
    let ((u0, u1), (v0, v1)) = cone.parameter_range();
    assert!((u0).abs() < 1e-10 && (u1 - TAU).abs() < 1e-10);
    assert!((v0).abs() < 1e-10 && v1.is_infinite());
}

#[test]
fn sphere_partials_normal_range() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    assert!(sphere.partial_u(0.0, FRAC_PI_2).length() > 1e-10);
    assert!(sphere.partial_v(0.0, FRAC_PI_2).length() > 1e-10);
    assert!(sphere.normal_at(0.0, FRAC_PI_2).length() > 1e-10);
    let ((u0, u1), (v0, v1)) = sphere.parameter_range();
    assert!((u0).abs() < 1e-10 && (u1 - TAU).abs() < 1e-10);
    assert!((v0).abs() < 1e-10 && (v1 - PI).abs() < 1e-10);
}

#[test]
fn torus_partials_normal_range() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    assert!(torus.partial_u(0.0, 0.0).length() > 1e-10);
    assert!(torus.partial_v(0.0, 0.0).length() > 1e-10);
    assert!(torus.normal_at(0.0, 0.0).length() > 1e-10);
    let ((u0, u1), (v0, v1)) = torus.parameter_range();
    assert!((u0).abs() < 1e-10 && (u1 - TAU).abs() < 1e-10);
    assert!((v0).abs() < 1e-10 && (v1 - TAU).abs() < 1e-10);
}

#[test]
fn nurbs_dispatch_all_methods() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let nurbs = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    });
    assert!((nurbs.point_at(0.5, 0.5).x - 0.5).abs() < 1e-4);
    assert!(nurbs.partial_u(0.5, 0.5).length() > 0.5);
    assert!(nurbs.partial_v(0.5, 0.5).length() > 0.5);
    assert!(nurbs.normal_at(0.5, 0.5).length() > 1e-10);
    let ((u0, u1), (v0, v1)) = nurbs.parameter_range();
    assert!((u0).abs() < 1e-10 && (u1 - 1.0).abs() < 1e-10);
    assert!((v0).abs() < 1e-10 && (v1 - 1.0).abs() < 1e-10);
}

#[test]
fn orthonormal_basis_all_axes() {
    for (axis, desc) in [
        (Vec3::new(0.0, 1.0, 0.0), "x smallest"),
        (Vec3::new(1.0, 0.0, 0.0), "y smallest"),
        (Vec3::new(0.0, 0.0, 1.0), "z smallest"),
    ] {
        let (p1, p2) = orthonormal_basis(axis);
        assert!((p1.dot(&axis)).abs() < 1e-10, "{desc}: p1 not perp");
        assert!((p2.dot(&axis)).abs() < 1e-10, "{desc}: p2 not perp");
    }
}

#[test]
fn orthonormal_basis_z_reference_branch() {
    let axis = Vec3::new(1.0, 1.0, 0.0).normalized().unwrap();
    let (p1, p2) = orthonormal_basis(axis);
    assert!((p1.dot(&axis)).abs() < 1e-10, "p1 not perpendicular");
    assert!((p2.dot(&axis)).abs() < 1e-10, "p2 not perpendicular");
    assert!((p1.dot(&p2)).abs() < 1e-10, "p1 and p2 not perpendicular");
}
