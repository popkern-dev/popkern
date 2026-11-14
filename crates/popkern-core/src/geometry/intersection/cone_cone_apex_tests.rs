use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;

fn cone(axis: Vec3) -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis,
        half_angle: std::f64::consts::FRAC_PI_4,
    }
}

/// Two cones of half-angle 45° about z and about z turned 60° about y,
/// apex at the origin, share the rulings along (1, ±√2, √3)/√6 (the
/// spherical sectors case): each is exact on both charts.
#[test]
fn cones_sharing_an_apex_share_their_rulings() {
    let (s, c) = (60f64.to_radians().psin(), 60f64.to_radians().pcos());
    let (k1, k2) = (cone(Vec3::new(0.0, 0.0, 1.0)), cone(Vec3::new(s, 0.0, c)));
    let w = ((0.0, std::f64::consts::TAU), (0.0, 10.0));
    let tol = Tolerance::default();
    let found = try_cone_cone_apex(&k1, w, &k2, w, &tol).expect("two cones, one apex");
    assert_eq!(found.len(), 2, "{found:?}");
    let want = [
        Vec3::new(1.0, 2f64.sqrt(), 3f64.sqrt()),
        Vec3::new(1.0, -(2f64.sqrt()), 3f64.sqrt()),
    ];
    for f in &found {
        let SurfaceSurfaceIntersection::Curve(cos) = f else {
            panic!("a curve");
        };
        let (t0, t1) = cos.curve_3d.parameter_range();
        let dir = (cos.curve_3d.point_at(t1) - cos.curve_3d.point_at(t0))
            .normalized()
            .expect("dir");
        assert!(
            want.iter()
                .any(|w| w.normalized().expect("w").cross(&dir).length() < 1e-12),
            "{dir:?}"
        );
        for k in 0..=8 {
            let t = t0 + (t1 - t0) * f64::from(k) / 8.0;
            let p = cos.curve_3d.point_at(t);
            let (q1, q2) = (cos.curve_uv1.point_at(t), cos.curve_uv2.point_at(t));
            assert!(
                k1.point_at(q1.x, q1.y).distance_to(&p) < 1e-12,
                "uv1 at {t}"
            );
            assert!(
                k2.point_at(q2.x, q2.y).distance_to(&p) < 1e-12,
                "uv2 at {t}"
            );
        }
    }
}

/// Cones with distinct apexes, or one axis, are not this arm's.
#[test]
fn other_cone_pairs_fall_through() {
    let tol = Tolerance::default();
    let w = ((0.0, std::f64::consts::TAU), (0.0, 10.0));
    let k1 = cone(Vec3::new(0.0, 0.0, 1.0));
    let moved = SurfaceGeometry::Cone {
        apex: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 1.0),
        half_angle: 0.5,
    };
    assert!(try_cone_cone_apex(&k1, w, &moved, w, &tol).is_none());
    assert!(try_cone_cone_apex(&k1, w, &cone(Vec3::new(0.0, 0.0, 1.0)), w, &tol).is_none());
}

/// PK-1138: a window padded past the apex reaches the other nappe,
/// where a cone face never is: the rulings run from the apex out, not
/// through it to the padding.
#[test]
fn a_window_padded_past_the_apex_keeps_its_own_nappe() {
    let (s, c) = (60f64.to_radians().psin(), 60f64.to_radians().pcos());
    let (k1, k2) = (cone(Vec3::new(0.0, 0.0, 1.0)), cone(Vec3::new(s, 0.0, c)));
    let padded = ((0.0, std::f64::consts::TAU), (-3.5, 10.6));
    let tol = Tolerance::default();
    let found = try_cone_cone_apex(&k1, padded, &k2, padded, &tol).expect("one apex");
    assert_eq!(found.len(), 2, "{found:?}");
    for f in &found {
        let SurfaceSurfaceIntersection::Curve(cos) = f else {
            panic!("a curve");
        };
        let (t0, t1) = cos.curve_3d.parameter_range();
        let (p, q) = (cos.curve_3d.point_at(t0), cos.curve_3d.point_at(t1));
        let origin = Point3::new(0.0, 0.0, 0.0);
        let near = p.distance_to(&origin).min(q.distance_to(&origin));
        assert!(near < 1e-12, "a ruling from the apex: {p:?} to {q:?}");
        assert!(
            p.z >= -1e-12 && q.z >= -1e-12,
            "on the upper nappe: {p:?} to {q:?}"
        );
    }
    // A window on the lower nappe keeps that one.
    assert_eq!(super::one_nappe((-10.6, 3.5)), (-10.6, 0.0));
    assert_eq!(super::one_nappe((10.6, -3.5)), (0.0, 10.6));
}
