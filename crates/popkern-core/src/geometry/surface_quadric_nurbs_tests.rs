use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{PI, TAU};

/// Every point of the conversion lies on the analytic surface, and at
/// the knots' angles the parameters are the chart's.
fn assert_exact(s: &SurfaceGeometry, u: (f64, f64), v: (f64, f64), on: impl Fn(Point3) -> f64) {
    let n = s.to_nurbs_in_range(u, v);
    let twin = SurfaceGeometry::Nurbs(n);
    for i in 0..=24 {
        for j in 0..=12 {
            let (a, b) = (
                u.0 + (u.1 - u.0) * f64::from(i) / 24.0,
                v.0 + (v.1 - v.0) * f64::from(j) / 12.0,
            );
            let p = twin.point_at(a, b);
            assert!(on(p) < 1e-12 * 10.0, "({a}, {b}): {:e}", on(p));
        }
    }
    // At the window's corners the parameters are the chart's.
    for (a, b) in [(u.0, v.0), (u.1, v.0), (u.0, v.1), (u.1, v.1)] {
        assert!(
            twin.point_at(a, b).distance_to(&s.point_at(a, b)) < 1e-9,
            "corner ({a}, {b})"
        );
    }
}

/// A full-turn cylinder of radius 5 (the sampled net read 4.394 from
/// its axis at a quarter turn), a partial cone and a sphere with a
/// tilted pole and with a +z pole convert exactly.
#[test]
fn cylinders_cones_and_spheres_convert_exactly() {
    let axis = Vec3::new(0.3, -0.2, 0.9).normalized().expect("axis");
    let o = Point3::new(5.0, -1.0, 2.0);
    let cyl = SurfaceGeometry::Cylinder {
        origin: o,
        axis,
        radius: 5.0,
    };
    assert_exact(&cyl, (0.0, TAU), (-12.0, 12.0), |p| {
        let d = p - o;
        ((d - axis * d.dot(&axis)).length() - 5.0).abs()
    });
    let alpha: f64 = 0.6;
    let cone = SurfaceGeometry::Cone {
        apex: o,
        axis,
        half_angle: alpha,
    };
    assert_exact(&cone, (0.4, 2.9), (1.0, 7.0), |p| {
        let d = p - o;
        let z = d.dot(&axis);
        ((d - axis * z).length() - z.abs() * alpha.tan()).abs()
    });
    for pole in [axis, Vec3::Z] {
        let ball = SurfaceGeometry::Sphere {
            center: o,
            radius: 4.0,
            axis: pole,
        };
        assert_exact(&ball, (0.0, TAU), (0.0, PI), |p| {
            (p.distance_to(&o) - 4.0).abs()
        });
        assert_exact(&ball, (1.0, 2.0), (0.3, 1.2), |p| {
            (p.distance_to(&o) - 4.0).abs()
        });
    }
}
