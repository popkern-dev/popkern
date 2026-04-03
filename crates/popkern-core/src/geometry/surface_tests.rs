use super::*;
use crate::geometry::surface::SurfaceEval;

#[test]
fn plane_is_not_periodic() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(plane.is_periodic(), (false, false));
}

#[test]
fn cylinder_is_periodic_in_u_only() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert_eq!(cyl.is_periodic(), (true, false));
}

#[test]
fn cone_is_periodic_in_u_only() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    assert_eq!(cone.is_periodic(), (true, false));
}

#[test]
fn sphere_is_periodic_in_u_only() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    assert_eq!(sphere.is_periodic(), (true, false));
}

#[test]
fn torus_is_periodic_in_both() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    assert_eq!(torus.is_periodic(), (true, true));
}

#[test]
fn periodicity_returns_correct_periods() {
    use std::f64::consts::TAU;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(plane.periodicity(), (None, None));

    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert_eq!(cyl.periodicity(), (Some(TAU), None));

    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    assert_eq!(torus.periodicity(), (Some(TAU), Some(TAU)));
}

#[test]
fn nurbs_is_not_periodic() {
    let nurbs = crate::geometry::nurbs_surface::NurbsSurface {
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
    };
    let geom = SurfaceGeometry::Nurbs(nurbs);
    assert_eq!(geom.is_periodic(), (false, false));
}

#[test]
fn plane_to_nurbs_preserves_shape() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let nurbs = plane.to_nurbs_in_range((-1.0, 1.0), (-1.0, 1.0));
    for u in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
        for v in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
            let p1 = plane.point_at(u, v);
            let p2 = nurbs.point_at(u, v);
            assert!(
                p1.distance_to(&p2) < 1e-8,
                "Plane-to-NURBS mismatch at u={u}, v={v}: p1={p1:?}, p2={p2:?}"
            );
        }
    }
}

#[test]
fn nurbs_to_nurbs_is_clone() {
    let nurbs = NurbsSurface {
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
    };
    let geom = SurfaceGeometry::Nurbs(nurbs.clone());
    let promoted = geom.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0));
    assert_eq!(promoted.degree_u, nurbs.degree_u);
}

#[test]
fn sphere_to_nurbs_approximate() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let nurbs = sphere.to_nurbs_in_range((0.5, 2.5), (0.5, 2.5));
    for u in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
        for v in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
            let p = nurbs.point_at(u, v);
            let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (r - 1.0).abs() < 0.1,
                "Point not on unit sphere: r={r} at u={u}, v={v}"
            );
        }
    }
}

#[test]
fn sphere_to_nurbs_default_range() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let nurbs = sphere.to_nurbs();
    let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
    let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Point not on unit sphere: r={r}");
}

#[test]
fn cylinder_to_nurbs_in_range() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let nurbs = cyl.to_nurbs_in_range((0.0, std::f64::consts::TAU), (-1.0, 1.0));
    let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, 0.0);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Point not on cylinder: r={r}");
}
