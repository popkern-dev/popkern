use super::*;
use crate::geometry::revolution::revolution_surface;
use std::f64::consts::TAU;

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// The surface `profile` revolved a full turn about z, its `v` on
/// [0, 1] as the revolve builds it.
fn revolved(profile: &NurbsCurve) -> SurfaceGeometry {
    let mut s =
        revolution_surface(profile, &Point3::new(0.0, 0.0, 0.0), &Vec3::Z, TAU).expect("revolve");
    for k in &mut s.knots_v {
        *k /= TAU;
    }
    SurfaceGeometry::Nurbs(s)
}

fn ball(z: f64, r: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, z),
        radius: r,
        axis: Vec3::Z,
    }
}

fn plane(z: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::Z,
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

const WIDE: Window = ((-1e3, 1e3), (-1e3, 1e3));
const CHART: Window = ((0.0, 1.0), (0.0, 1.0));

/// Every curve is a circle of radius `r` at height `z` about the axis,
/// and both pcurves map onto it at the curve's own parameter.
fn assert_circle(
    found: &[SurfaceSurfaceIntersection],
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    z: f64,
    r: f64,
) {
    assert_eq!(found.len(), 1, "one circle");
    let SurfaceSurfaceIntersection::Curve(c) = &found[0] else {
        panic!("a curve");
    };
    assert!(c.is_closed);
    let (t0, t1) = c.curve_3d.parameter_range();
    for k in 0..=200 {
        let t = t0 + (t1 - t0) * f64::from(k) / 200.0;
        let p = c.curve_3d.point_at(t);
        assert!(
            (p.x.hypot(p.y) - r).abs() < 1e-12 * r.max(1.0) && (p.z - z).abs() < 1e-12,
            "{p:?}"
        );
        let (q1, q2) = (c.curve_uv1.point_at(t), c.curve_uv2.point_at(t));
        assert!(s1.point_at(q1.x, q1.y).distance_to(&p) < 1e-6, "uv1 at {t}");
        assert!(s2.point_at(q2.x, q2.y).distance_to(&p) < 1e-6, "uv2 at {t}");
    }
}

/// A cup's revolved top (the radial line from the axis to radius 7.07
/// at z = 5) meets a ball of radius 5 centred there in its equator, in
/// either argument order: exactly, as the surface's iso-circle.
#[test]
fn a_revolved_disc_meets_a_coaxial_ball_in_its_equator() {
    let tol = Tolerance::default();
    let disc = revolved(&line(
        Point3::new(0.0, 0.0, 5.0),
        Point3::new(50f64.sqrt(), 0.0, 5.0),
    ));
    let b = ball(5.0, 5.0);
    let found = try_revolved_nurbs(&disc, CHART, &b, WIDE, &tol).expect("coaxial");
    assert_circle(&found, &disc, &b, 5.0, 5.0);
    let found = try_revolved_nurbs(&b, WIDE, &disc, CHART, &tol).expect("either order");
    assert_circle(&found, &b, &disc, 5.0, 5.0);
}

/// The paraboloid x² + y² = 8z (a parabola, a polynomial quadratic,
/// revolved) cut 1e-4 above its tip: the circle of radius √(8e-4),
/// beside the chart's pole.
#[test]
fn a_paraboloid_cut_just_above_its_tip_is_one_exact_circle() {
    let tol = Tolerance::default();
    let parabola = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(8.0, 0.0, 8.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let cup = revolved(&parabola);
    let p = plane(1e-4);
    let found = try_revolved_nurbs(&cup, CHART, &p, WIDE, &tol).expect("coaxial");
    assert_circle(&found, &cup, &p, 1e-4, 8e-4f64.sqrt());
    // Its tip touches the plane through it, a point: nothing emitted.
    let found = try_revolved_nurbs(&cup, CHART, &plane(0.0), WIDE, &tol).expect("coaxial");
    assert!(found.is_empty(), "{}", found.len());
}

/// A ball off the axis, a tilted plane, and a surface that is not a
/// revolution keep the marching path; the disc's own plane is a
/// coincidence, left to the coincidence detector.
#[test]
fn off_axis_pairs_non_revolutions_and_coincidences_decline() {
    let tol = Tolerance::default();
    let disc = revolved(&line(
        Point3::new(0.0, 0.0, 5.0),
        Point3::new(7.0, 0.0, 5.0),
    ));
    let off = SurfaceGeometry::Sphere {
        center: Point3::new(0.5, 0.0, 5.0),
        radius: 5.0,
        axis: Vec3::Z,
    };
    assert!(try_revolved_nurbs(&disc, CHART, &off, WIDE, &tol).is_none());
    let tilted = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        normal: Vec3::new(0.0, 0.6, 0.8),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.8, -0.6),
    };
    assert!(try_revolved_nurbs(&disc, CHART, &tilted, WIDE, &tol).is_none());
    assert!(try_revolved_nurbs(&disc, CHART, &plane(5.0), WIDE, &tol).is_none());
    // A bilinear patch: no arcs in v.
    let patch = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    assert!(try_revolved_nurbs(&patch, CHART, &ball(0.0, 1.0), WIDE, &tol).is_none());
}

/// A revolve authors the angle in u (PK-455): the same paraboloid with
/// its parameter directions exchanged gives the same circle, its pcurve
/// the line at the root's v.
#[test]
fn a_surface_with_the_angle_in_u_is_read_through_its_transpose() {
    let tol = Tolerance::default();
    let parabola = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(8.0, 0.0, 8.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let SurfaceGeometry::Nurbs(s) = revolved(&parabola) else {
        panic!("nurbs");
    };
    let cup = SurfaceGeometry::Nurbs(transpose(&s));
    let p = plane(1e-4);
    let found = try_revolved_nurbs(&cup, CHART, &p, WIDE, &tol).expect("coaxial");
    assert_circle(&found, &cup, &p, 1e-4, 8e-4f64.sqrt());
    // And in the periodic form a revolve serves (PK-890: its knot
    // vector runs past the domain, a wrap-around control point, an
    // even count): the domain and the on-circle points are read from
    // the knots and weights, not the vector's ends and even indices.
    let periodic = transpose(&s)
        .periodic_u_from_closed_clamped(tol.linear())
        .expect("periodic form");
    assert!(periodic.num_u.is_multiple_of(2), "{}", periodic.num_u);
    let cup = SurfaceGeometry::Nurbs(periodic);
    let found = try_revolved_nurbs(&cup, CHART, &p, WIDE, &tol).expect("coaxial");
    assert_circle(&found, &cup, &p, 1e-4, 8e-4f64.sqrt());
}
