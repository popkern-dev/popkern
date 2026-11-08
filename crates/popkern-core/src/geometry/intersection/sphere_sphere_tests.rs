use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;

fn ball(c: (f64, f64, f64), r: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(c.0, c.1, c.2),
        radius: r,
        axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

/// Balls of radius 5 at the origin and (4.8, 0, 3.6) meet in the circle
/// of radius 4 about (2.4, 0, 1.8), through the first's pole (0, 0, 5)
/// and the second's (4.8, 0, −1.4); each pcurve maps onto the 3D circle
/// at the same parameter, through the poles too.
#[test]
fn balls_meeting_through_a_pole_meet_in_one_circle() {
    let tol = Tolerance::default();
    let (s1, s2) = (ball((0.0, 0.0, 0.0), 5.0), ball((4.8, 0.0, 3.6), 5.0));
    let found = try_sphere_sphere(&s1, &s2, &tol).expect("two spheres");
    assert_eq!(found.len(), 1);
    let SurfaceSurfaceIntersection::Curve(c) = &found[0] else {
        panic!("a curve");
    };
    let centre = Point3::new(2.4, 0.0, 1.8);
    let (t0, t1) = c.curve_3d.parameter_range();
    let mut through_poles = [false, false];
    for k in 0..=400 {
        let t = t0 + (t1 - t0) * f64::from(k) / 400.0;
        let p = c.curve_3d.point_at(t);
        assert!((p.distance_to(&centre) - 4.0).abs() < 1e-9, "{p:?}");
        let q1 = c.curve_uv1.point_at(t);
        let q2 = c.curve_uv2.point_at(t);
        assert!(s1.point_at(q1.x, q1.y).distance_to(&p) < 1e-5, "uv1 at {t}");
        assert!(s2.point_at(q2.x, q2.y).distance_to(&p) < 1e-5, "uv2 at {t}");
        through_poles[0] |= p.distance_to(&Point3::new(0.0, 0.0, 5.0)) < 0.05;
        through_poles[1] |= p.distance_to(&Point3::new(4.8, 0.0, -1.4)) < 0.05;
    }
    assert_eq!(through_poles, [true, true]);
}

/// Touching, apart, nested, and concentric spheres emit no curve;
/// coincident ones fall through to the coincidence detector.
#[test]
fn touching_apart_and_coincident_spheres() {
    let tol = Tolerance::default();
    let a = ball((0.0, 0.0, 0.0), 5.0);
    for b in [
        ball((10.0, 0.0, 0.0), 5.0),
        ball((20.0, 0.0, 0.0), 5.0),
        ball((1.0, 0.0, 0.0), 2.0),
        ball((0.0, 0.0, 0.0), 3.0),
        ball((3.0, 0.0, 0.0), 2.0),
    ] {
        assert_eq!(try_sphere_sphere(&a, &b, &tol).map(|v| v.len()), Some(0));
    }
    assert!(try_sphere_sphere(&a, &ball((0.0, 0.0, 0.0), 5.0), &tol).is_none());
}
