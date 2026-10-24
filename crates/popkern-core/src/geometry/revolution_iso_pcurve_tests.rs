use std::f64::consts::{FRAC_PI_2, PI, TAU};

use super::iso_circle_pcurve;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::pcurve::PcurveGeometry;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

fn arc(center: Point3, normal: Vec3, radius: f64, start: f64, end: f64) -> CurveGeometry {
    CurveGeometry::Arc {
        center,
        normal,
        radius,
        start_angle: start,
        end_angle: end,
    }
}

/// `surface(pcurve(t)) == arc(t)` at 33 stations, to machine precision.
fn assert_exact(surface: &SurfaceGeometry, curve: &CurveGeometry, range: (f64, f64)) {
    let pcurve = iso_circle_pcurve(surface, curve, range).expect("an iso arm applies");
    assert!(matches!(pcurve.geometry, PcurveGeometry::Linear { .. }));
    let mut worst = 0.0_f64;
    for i in 0..=32 {
        let t = range.0 + (range.1 - range.0) * i as f64 / 32.0;
        let uv = pcurve.point_at(t);
        let on_surface = surface.point_at(uv.x, uv.y);
        worst = worst.max(on_surface.distance_to(&curve.point_at(t)));
    }
    assert!(worst < 1.0e-12, "worst deviation {worst:.3e}");
}

#[test]
fn sphere_latitude_circles_are_u_lines() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.5, 0.5, 0.5),
        radius: 0.6,
        axis: crate::geometry::Vec3::Z,
    };
    let polar = 0.8_f64;
    let center = Point3::new(0.5, 0.5, 0.5 + 0.6 * polar.cos());
    for normal in [Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0)] {
        let curve = arc(center, normal, 0.6 * polar.sin(), -1.0, 2.5);
        assert_exact(&surface, &curve, (-1.0, 2.5));
    }
}

#[test]
fn sphere_meridian_arcs_are_v_lines_until_they_cross_a_pole() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
        axis: crate::geometry::Vec3::Z,
    };
    // A great circle in the xz-plane: normal along y. Where its own
    // parameter reaches the north pole depends on the arc's basis, so
    // locate it: z(t) = R·cos(t − φ).
    let curve = arc(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        2.0,
        0.0,
        TAU,
    );
    let phi = curve.point_at(FRAC_PI_2).z.atan2(curve.point_at(0.0).z);
    // Both endpoints away from the poles, no pole inside: exact.
    assert_exact(&surface, &curve, (phi + 0.2, phi + 1.4));
    // Ending exactly at a pole is still one-sided.
    assert_exact(&surface, &curve, (phi - 1.3, phi));
    // The same circle with normal along -y runs the other way.
    let reversed = arc(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        2.0,
        0.0,
        TAU,
    );
    let phi_r = reversed
        .point_at(FRAC_PI_2)
        .z
        .atan2(reversed.point_at(0.0).z);
    assert_exact(&surface, &reversed, (phi_r + 0.2, phi_r + 1.4));
    // Straddling the pole folds v: declined, the caller fits instead.
    assert!(iso_circle_pcurve(&surface, &curve, (phi - 0.5, phi + 1.0)).is_none());
}

#[test]
fn a_tilted_sphere_meridian_ending_on_its_pole_is_exact() {
    // An intersection solver's arc centre sits a few ulps off the
    // sphere's, so its end misses the pole by ~1e-16 in z / r; acos
    // there lost half its digits and the fold check declined the
    // meridian (PK-930). The polar angle is atan2 now.
    let axis = Vec3::new(0.3, -0.4, 0.8).normalized().expect("axis");
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 1.0, 1.0),
        radius: 0.5,
        axis,
    };
    let normal = axis
        .cross(&Vec3::new(1.0, 0.0, 0.0))
        .normalized()
        .expect("normal");
    let centre = Point3::new(1.0 + 2e-15, 1.0 - 1e-15, 1.0 + 1e-15);
    let curve = arc(centre, normal, 0.5, 0.0, TAU);
    let height = |t: f64| (curve.point_at(t) - centre).dot(&axis);
    let phi = height(FRAC_PI_2).atan2(height(0.0));
    for range in [(phi - 1.3, phi), (phi + PI, phi + PI + 1.1)] {
        assert_exact(&surface, &curve, range);
    }
}

#[test]
fn torus_latitude_and_meridian_circles_are_chart_lines() {
    let axis = Vec3::new(0.0, 1.0, 0.0);
    let surface = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis,
        major_radius: 1.0,
        minor_radius: 0.3,
    };
    // Latitude: the tube's top circle (v = π/2) and its outer equator.
    let top = arc(Point3::new(0.0, 0.3, 0.0), axis, 1.0, 0.3, 4.0);
    assert_exact(&surface, &top, (0.3, 4.0));
    let equator = arc(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        1.3,
        0.0,
        PI,
    );
    assert_exact(&surface, &equator, (0.0, PI));
    // Meridian: a tube circle at azimuth u = 1.1 (its plane contains
    // the axis), either sense.
    let spine = Point3::new(1.1_f64.cos(), 0.0, 1.1_f64.sin());
    let tangent = Vec3::new(-(1.1_f64.sin()), 0.0, 1.1_f64.cos());
    for n in [tangent, tangent * -1.0] {
        let tube = arc(spine, n, 0.3, -2.0, 1.0);
        assert_exact(&surface, &tube, (-2.0, 1.0));
    }
}

#[test]
fn cylinder_and_cone_stations_are_u_lines_and_everything_else_declines() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.5,
    };
    let rim = arc(
        Point3::new(0.0, 0.0, 0.7),
        Vec3::new(0.0, 0.0, 1.0),
        0.5,
        0.0,
        TAU,
    );
    assert_exact(&cylinder, &rim, (0.0, TAU));
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.4,
    };
    let station = arc(
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(0.0, 0.0, -1.0),
        2.0 * 0.4_f64.tan(),
        1.0,
        2.0,
    );
    assert_exact(&cone, &station, (1.0, 2.0));
    // Off-axis circle on the cylinder: not an iso-line.
    let tilted = arc(
        Point3::new(0.1, 0.0, 0.7),
        Vec3::new(0.0, 0.0, 1.0),
        0.5,
        0.0,
        TAU,
    );
    assert!(iso_circle_pcurve(&cylinder, &tilted, (0.0, TAU)).is_none());
    // A plane never has iso-circles; a line is never an arc.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(iso_circle_pcurve(&plane, &rim, (0.0, TAU)).is_none());
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(iso_circle_pcurve(&cylinder, &line, (0.0, 1.0)).is_none());
}
