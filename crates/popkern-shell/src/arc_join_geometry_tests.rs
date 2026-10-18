use super::*;
use popkern_core::geometry::curve::CurveOnEdge;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn great_arc_reproduces_endpoints_and_stays_on_sphere() {
    let c = Point3::new(1.0, 2.0, 3.0);
    let r = 0.25;
    let a = c + Vec3::new(r, 0.0, 0.0);
    let b = c + Vec3::new(0.0, 0.0, r);
    let arc = great_arc(c, r, a, b, &tol()).unwrap();
    assert!(arc.geometry.point_at(arc.t_start).distance_to(&a) < 1e-12);
    assert!(arc.geometry.point_at(arc.t_end).distance_to(&b) < 1e-12);
    // Quarter arc, minor by construction.
    assert!((arc.t_end - arc.t_start - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    let mid = arc.geometry.point_at(0.5 * (arc.t_start + arc.t_end));
    assert!((mid.distance_to(&c) - r).abs() < 1e-12);
}

#[test]
fn great_arc_declines_antipodal_endpoints() {
    let c = Point3::new(0.0, 0.0, 0.0);
    let r = 1.0;
    let a = Point3::new(r, 0.0, 0.0);
    let b = Point3::new(-r, 0.0, 0.0);
    assert!(great_arc(c, r, a, b, &tol()).is_none());
}

#[test]
fn offset_rail_translates_a_line_on_a_plane_offset() {
    let clone = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 1.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    // Offset of the z=1 plane moved up by 0.2.
    let offset = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.2),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let (rail, dev) = offset_rail(&clone, &offset, None, &tol()).unwrap();
    assert_eq!(dev, 0.0, "similarity arm is exact");
    let p = rail.geometry.point_at(0.5);
    assert!(p.distance_to(&Point3::new(0.5, 0.0, 1.2)) < 1e-12, "{p:?}");
}

#[test]
fn strip_surface_promotes_line_to_cylinder_and_circle_to_torus() {
    let line = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 1.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let dummy = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s = strip_surface(
        &line,
        0.1,
        &dummy,
        &dummy,
        &line,
        &line,
        &Tolerance::default(),
    )
    .unwrap();
    assert!(matches!(s, SurfaceGeometry::Cylinder { radius, .. } if (radius - 0.1).abs() < 1e-15));

    let circle = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let s = strip_surface(
        &circle,
        0.1,
        &dummy,
        &dummy,
        &circle,
        &circle,
        &Tolerance::default(),
    )
    .unwrap();
    assert!(matches!(
        s,
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } if (major_radius - 1.0).abs() < 1e-15 && (minor_radius - 0.1).abs() < 1e-15
    ));

    // Spindle guard: ball radius >= spine radius declines.
    assert!(
        strip_surface(
            &circle,
            1.5,
            &dummy,
            &dummy,
            &circle,
            &circle,
            &Tolerance::default()
        )
        .is_err()
    );
}
