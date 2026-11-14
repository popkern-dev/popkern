//! The typed intersection's curves are the untyped intersection's, bit
//! for bit: contacts are reported beside them and change nothing a
//! tangency-silent caller reads (PK-1126, the plan's S2 exit).

use super::*;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{PI, TAU};

/// A chart window.
type Window = ((f64, f64), (f64, f64));
/// A surface pair, each over its window.
type Case = (SurfaceGeometry, Window, SurfaceGeometry, Window);

#[test]
fn typed_curves_are_the_untyped_curves() {
    let tol = Tolerance::default();
    let sphere = |c: (f64, f64, f64), r: f64| SurfaceGeometry::Sphere {
        center: Point3::new(c.0, c.1, c.2),
        radius: r,
        axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let z_cyl = |x: f64, r: f64| SurfaceGeometry::Cylinder {
        origin: Point3::new(x, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(-20.0, -20.0, -5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 3.0,
    };
    let full_sphere = ((0.0, TAU), (0.0, PI));
    let rod = ((0.0, TAU), (0.0, 24.0));
    let sheet = ((0.0, 40.0), (0.0, 40.0));
    // Viviani (crossing), a ball in its band (tangent curve), a ball
    // on a face (touch), a torus crossing a plane, a plain pierce.
    let cases: Vec<Case> = vec![
        (
            z_cyl(5.0, 5.0),
            rod,
            sphere((0.0, 0.0, 0.0), 10.0),
            full_sphere,
        ),
        (
            z_cyl(0.0, 5.0),
            rod,
            sphere((0.0, 0.0, 0.0), 5.0),
            full_sphere,
        ),
        (
            plane.clone(),
            sheet,
            sphere((0.0, 0.0, 0.0), 5.0),
            full_sphere,
        ),
        (plane, sheet, torus, ((0.0, TAU), (0.0, TAU))),
        (
            z_cyl(2.0, 3.0),
            rod,
            sphere((0.0, 0.0, 0.0), 8.0),
            full_sphere,
        ),
    ];
    for (k, (s1, (u1, v1), s2, (u2, v2))) in cases.iter().enumerate() {
        let plain = intersect_surfaces(s1, *u1, *v1, s2, *u2, *v2, &tol).expect("plain");
        let typed =
            intersect_surfaces_typed(s1, *u1, *v1, s2, *u2, *v2, &tol, 20.0).expect("typed");
        assert_eq!(
            format!("{:?}", typed.curves),
            format!("{plain:?}"),
            "case {k}: the typed curves differ"
        );
    }
}
