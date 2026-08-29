//! Unit tests for exact strip end-section curves (PK-262).

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{Point3, Vec3};

use super::{verified_circle_section, verified_ellipse_section};

const AXIS: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};

fn origin() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

fn radial(p: &Point3) -> f64 {
    (p.x * p.x + p.y * p.y).sqrt()
}

#[test]
fn perpendicular_section_is_a_short_circle_arc() {
    let n_hat = Vec3::new(0.0, 0.0, 1.0);
    let plane_origin = Point3::new(0.0, 0.0, 2.0);
    let c0 = Point3::new(0.0, 0.0, 2.0);
    let p1 = Point3::new(1.0, 0.0, 2.0);
    let p2 = Point3::new(0.0, 1.0, 2.0);
    let (curve, s1, s2) =
        verified_circle_section(&origin(), &AXIS, 1.0, c0, p1, p2, &plane_origin, &n_hat)
            .expect("perpendicular section must author");
    assert!(matches!(curve.geometry, CurveGeometry::Arc { .. }));
    assert!(curve.t_end - curve.t_start <= std::f64::consts::PI + 1.0e-9);
    let ceval: &dyn CurveEval = &curve.geometry;
    assert!(ceval.point_at(curve.t_start).distance_to(&s1) < 1.0e-9);
    assert!(ceval.point_at(curve.t_end).distance_to(&s2) < 1.0e-9);
}

/// The section of the unit cylinder by an oblique plane, with feet
/// generated ON the true section, authors an exact ellipse arc: feet
/// round-trip, every sample lies on the cylinder and in the plane,
/// and the sweep takes the short side.
#[test]
fn oblique_section_is_an_exact_ellipse_arc() {
    let alpha: f64 = 0.6;
    let n_hat = Vec3::new(alpha.sin(), 0.0, alpha.cos());
    let c0 = Point3::new(0.0, 0.0, 1.5);
    let plane_origin = c0;
    let axis_cos = alpha.cos();
    // Feet on the true section: the ellipse with semi-major
    // r/cos(alpha) along the axis's in-plane projection.
    let major_dir = (AXIS - n_hat * axis_cos).normalized().unwrap();
    let minor_dir = n_hat.cross(&major_dir);
    let foot = |theta: f64| c0 + major_dir * (theta.cos() / axis_cos) + minor_dir * theta.sin();
    let p1 = foot(0.3);
    let p2 = foot(1.8);
    let (curve, s1, s2) = verified_ellipse_section(
        &origin(),
        &AXIS,
        1.0,
        c0,
        p1,
        p2,
        &plane_origin,
        &n_hat,
        axis_cos,
    )
    .expect("oblique section must author");
    assert!(matches!(curve.geometry, CurveGeometry::Ellipse { .. }));
    assert!(curve.t_end - curve.t_start <= std::f64::consts::PI + 1.0e-9);
    assert!(s1.distance_to(&p1) < 1.0e-9, "foot 1 already on section");
    assert!(s2.distance_to(&p2) < 1.0e-9, "foot 2 already on section");
    let ceval: &dyn CurveEval = &curve.geometry;
    assert!(ceval.point_at(curve.t_start).distance_to(&s1) < 1.0e-9);
    assert!(ceval.point_at(curve.t_end).distance_to(&s2) < 1.0e-9);
    for k in 0..=16 {
        let t = curve.t_start + (curve.t_end - curve.t_start) * k as f64 / 16.0;
        let q = ceval.point_at(t);
        assert!(
            (radial(&q) - 1.0).abs() < 1.0e-9,
            "sample {k} must lie exactly on the cylinder; radial={}",
            radial(&q)
        );
        assert!(
            (q - plane_origin).dot(&n_hat).abs() < 1.0e-9,
            "sample {k} must lie exactly in the end plane"
        );
    }
}

/// Feet subtending more than π in parameter still author — via the
/// flipped orientation, whose positive sweep is the short side. The
/// author never returns a long-way arc (the complementary-side class
/// measured on pyramid3_extreme's far runout).
#[test]
fn long_way_feet_author_the_short_side() {
    let alpha: f64 = 0.4;
    let n_hat = Vec3::new(alpha.sin(), 0.0, alpha.cos());
    let c0 = Point3::new(0.0, 0.0, -0.5);
    let axis_cos = alpha.cos();
    let major_dir = (AXIS - n_hat * axis_cos).normalized().unwrap();
    let minor_dir = n_hat.cross(&major_dir);
    let foot = |theta: f64| c0 + major_dir * (theta.cos() / axis_cos) + minor_dir * theta.sin();
    let p1 = foot(0.2);
    let p2 = foot(4.5);
    let (curve, _, _) = verified_ellipse_section(
        &origin(),
        &AXIS,
        1.0,
        c0,
        p1,
        p2,
        &c0.clone(),
        &n_hat,
        axis_cos,
    )
    .expect("short-side orientation must author");
    assert!(
        curve.t_end - curve.t_start <= std::f64::consts::PI + 1.0e-9,
        "sweep must be the short side; got {}",
        curve.t_end - curve.t_start
    );
}
