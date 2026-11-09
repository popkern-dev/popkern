use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;

const WIDE: ((f64, f64), (f64, f64)) = ((-100.0, 100.0), (-100.0, 100.0));

fn ball(c: Point3, r: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: c,
        radius: r,
        axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn cone(apex: Point3, half_angle: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex,
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle,
    }
}

/// Every curve's 3D points and both pcurves' images lie on a circle of
/// one of `radii` about the z axis through `o`, at the matching height.
fn assert_circles(
    found: &[SurfaceSurfaceIntersection],
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    o: Point3,
    rings: &[(f64, f64)],
) {
    assert_eq!(found.len(), rings.len());
    for (f, &(z, r)) in found.iter().zip(rings) {
        let SurfaceSurfaceIntersection::Curve(c) = f else {
            panic!("a curve");
        };
        let on_ring = |p: Point3| {
            let q = p - o;
            ((q.x * q.x + q.y * q.y).sqrt() - r)
                .abs()
                .max((q.z - z).abs())
        };
        for k in 0..=200 {
            let s = f64::from(k) / 200.0;
            let at = |crv: &crate::geometry::nurbs_curve::NurbsCurve| {
                let (t0, t1) = crv.parameter_range();
                crv.point_at(t0 + (t1 - t0) * s)
            };
            assert!(on_ring(at(&c.curve_3d)) < 1e-9);
            let q1 = at(&c.curve_uv1);
            let q2 = at(&c.curve_uv2);
            assert!(on_ring(s1.point_at(q1.x, q1.y)) < 1e-6, "uv1 at {s}");
            assert!(on_ring(s2.point_at(q2.x, q2.y)) < 1e-6, "uv2 at {s}");
        }
    }
}

/// A ball of radius 8 centred 10 up a 45° cone's axis meets it in two
/// circles, at v = (20 ± √112) / 4, in either argument order and a
/// thousand units from the origin alike.
#[test]
fn a_ball_on_the_axis_meets_the_cone_in_two_circles() {
    let tol = Tolerance::default();
    let q = std::f64::consts::FRAC_PI_4;
    for o in [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1234.5, -678.25, 910.125),
    ] {
        let (b, k) = (ball(o + Vec3::new(0.0, 0.0, 10.0), 8.0), cone(o, q));
        let rings: Vec<(f64, f64)> = [(20.0 - 112f64.sqrt()) / 4.0, (20.0 + 112f64.sqrt()) / 4.0]
            .iter()
            .map(|&v| (v, v))
            .collect();
        let found = try_sphere_cone(&b, &k, WIDE, &tol).expect("on the axis");
        assert_circles(&found, &b, &k, o, &rings);
        let found = try_sphere_cone(&k, &b, WIDE, &tol).expect("either order");
        assert_circles(&found, &k, &b, o, &rings);
    }
}

/// The cone's v window clips the circles: [0, 5] keeps the lower one.
#[test]
fn the_cone_window_clips_the_circles() {
    let tol = Tolerance::default();
    let o = Point3::new(0.0, 0.0, 0.0);
    let (b, k) = (
        ball(Point3::new(0.0, 0.0, 10.0), 8.0),
        cone(o, std::f64::consts::FRAC_PI_4),
    );
    let v = (20.0 - 112f64.sqrt()) / 4.0;
    let found = try_sphere_cone(&b, &k, ((-100.0, 100.0), (0.0, 5.0)), &tol).expect("on the axis");
    assert_circles(&found, &b, &k, o, &[(v, v)]);
}

/// The ice-cream cone (apex 25 below a ball of radius 15, half-angle
/// atan(3/4)) touches its ball along a circle: nothing is emitted, near
/// the origin and far from it. A ball through the apex meets the cone
/// in its other circle only; a ball off the axis keeps the marcher; a
/// ball clear of the cone meets nothing.
#[test]
fn touches_apex_points_and_off_axis_balls() {
    let tol = Tolerance::default();
    let alpha = 0.75f64.atan();
    for o in [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1234.5, -678.25, 910.125),
    ] {
        let found = try_sphere_cone(
            &ball(o, 15.0),
            &cone(o - Vec3::new(0.0, 0.0, 25.0), alpha),
            WIDE,
            &tol,
        );
        assert_eq!(found.map(|v| v.len()), Some(0));
    }
    let o = Point3::new(0.0, 0.0, 0.0);
    let q = std::f64::consts::FRAC_PI_4;
    // Through the apex: 2v² − 20v = 0, v = 0 (a point) and v = 10.
    let b = ball(Point3::new(0.0, 0.0, 10.0), 10.0);
    let k = cone(o, q);
    let found = try_sphere_cone(&b, &k, WIDE, &tol).expect("on the axis");
    assert_circles(&found, &b, &k, o, &[(10.0, 10.0)]);
    assert!(try_sphere_cone(&ball(Point3::new(0.5, 0.0, 10.0), 8.0), &k, WIDE, &tol).is_none());
    let clear = try_sphere_cone(&ball(Point3::new(0.0, 0.0, 30.0), 5.0), &k, WIDE, &tol);
    assert_eq!(clear.map(|v| v.len()), Some(0));
}
