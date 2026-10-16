//! Projected correspondence (PK-499 exploration E2): closed seams
//! whose pcurves start at different phases, and a rational-arc
//! parameter against an angle, certify under foot-point
//! correspondence where fraction correspondence fails.

use super::*;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{FRAC_PI_2, TAU};

const R: f64 = 1.5;

fn sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: R,
        axis: crate::geometry::Vec3::Z,
    }
}

/// Sphere equator from longitude `phase` once round.
fn equator(phase: f64) -> Pcurve {
    Pcurve::linear((phase, FRAC_PI_2), (phase + TAU, FRAC_PI_2), 0.0, 1.0)
}

fn side(surface: &SurfaceGeometry, pcurve: &Pcurve) -> SeamSideSource {
    SeamSideSource::new(surface, pcurve, false).unwrap_or_else(|e| panic!("{e}"))
}

fn assert_projected_rail(cert: &SeamCertificate, correspondence: SeamCorrespondence) {
    assert!(cert.is_certified(), "{cert:?}");
    assert_eq!(cert.correspondence, correspondence);
    let options = SeamCertifyOptions::for_model(R, 1e-6);
    assert!(
        cert.g0_upper <= options.correspondence_tolerance,
        "{cert:?}"
    );
    assert!(
        cert.g1_upper <= 2.0 * options.targets.normal_sine,
        "{cert:?}"
    );
    // The cylinder is flat along its axis, the sphere curves 1/R: the
    // transverse jump is exactly 1/R.
    let slack = options.targets.curvature + options.targets.relative / R;
    assert!(
        cert.g2_upper >= 1.0 / R && cert.g2_upper <= 1.0 / R + 2.0 * slack,
        "{cert:?}"
    );
}

#[test]
fn a_phase_shifted_closed_rail_certifies_under_its_phase_shift() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: R,
    };
    let rim = Pcurve::linear((0.0, 0.0), (TAU, 0.0), 0.0, 1.0);
    let (s, c) = (sphere(), cylinder);
    let mut a = side(&s, &equator(1.1));
    let mut b = side(&c, &rim);
    let cert = certify_seam(&mut a, &mut b, &SeamCertifyOptions::for_model(R, 1e-6));
    assert_projected_rail(&cert, SeamCorrespondence::PhaseShifted);
    // E2 cost gate: a full-circle rail under the order-4 Taylor forms.
    assert!(cert.intervals <= 2048, "{} intervals", cert.intervals);
}

/// A cylinder as a rational-quadratic NURBS (nine-point circle in u,
/// linear in v): its u is a rational-arc parameter, not an angle.
fn rational_cylinder() -> NurbsSurface {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let ring = [
        (1.0, 0.0, 1.0),
        (1.0, 1.0, h),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, h),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, h),
        (0.0, -1.0, 1.0),
        (1.0, -1.0, h),
        (1.0, 0.0, 1.0),
    ];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for &(x, y, w) in &ring {
        for z in [0.0, 2.0] {
            control_points.push(Point3::new(R * x, R * y, z));
            weights.push(w);
        }
    }
    NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 9,
        num_v: 2,
        weights,
    }
}

#[test]
fn a_rational_arc_rail_against_an_angle_rail_certifies_under_projection() {
    let (s, c) = (sphere(), SurfaceGeometry::Nurbs(rational_cylinder()));
    let rim = Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0);
    let mut a = side(&s, &equator(0.0));
    let mut b = side(&c, &rim);
    let options = SeamCertifyOptions::for_model(R, 1e-6);
    let cert = certify_seam(&mut a, &mut b, &options);
    assert_projected_rail(&cert, SeamCorrespondence::Projected);
    // E2 cost gate: a full-circle rail under the order-4 Taylor forms.
    assert!(cert.intervals <= 2048, "{} intervals", cert.intervals);
}

#[test]
fn a_rail_that_misses_its_partner_certifies_the_miss() {
    // The cylinder is 1e-4 wider than the sphere: a real 1e-4 gap the
    // projected certificate must bound from both sides.
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: R + 1e-4,
    };
    let rim = Pcurve::linear((0.0, 0.0), (TAU, 0.0), 0.0, 1.0);
    let s = sphere();
    let mut a = side(&s, &equator(2.0));
    let mut b = side(&cylinder, &rim);
    let cert = certify_seam(&mut a, &mut b, &SeamCertifyOptions::for_model(R, 1e-6));
    assert!(
        cert.lower[0] >= 1e-4 * 0.999 && cert.g0_upper <= 1e-4 * 1.01,
        "{cert:?}"
    );
    assert!(
        !cert.is_certified(),
        "a 1e-4 miss certified as a seam: {cert:?}"
    );
}

/// PK-543: a seam running into a sphere's chart pole. The meridian
/// `u = 0` from the equator to the north pole (`v → 0`, where
/// `S_u × S_v` vanishes) against a cylinder of radius `R` about the y
/// axis, tangent to the sphere along that great circle. The sphere
/// side is evaluated chart-free (normal `(C − P)/R`, umbilic second
/// form), so the pole costs nothing: G0/G1 at roundoff, the transverse
/// curvature jump exactly `1/R`, nothing uncertified.
#[test]
fn a_seam_into_a_sphere_pole_certifies_chart_free() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: R,
    };
    let meridian = Pcurve::linear((0.0, FRAC_PI_2), (0.0, 0.0), 0.0, 1.0);
    // Cylinder frame: perp1 = x, perp2 = y × x = −z, so the great
    // circle point R(sin v, 0, cos v) sits at angle θ = v − π/2.
    let rail = Pcurve::linear((0.0, 0.0), (-FRAC_PI_2, 0.0), 0.0, 1.0);
    let (s, c) = (sphere(), cylinder);
    let mut a = side(&s, &meridian);
    let mut b = side(&c, &rail);
    let options = SeamCertifyOptions::for_model(R, 1e-6);
    let cert = certify_seam(&mut a, &mut b, &options);
    assert!(cert.is_certified(), "{cert:?}");
    assert_eq!(cert.uncertified_length, 0.0);
    assert!(
        cert.g0_upper <= options.targets.gap && cert.g1_upper < 1e-6,
        "{cert:?}"
    );
    let slack = options.targets.curvature + options.targets.relative / R;
    assert!(
        cert.g2_upper >= 1.0 / R && cert.g2_upper <= 1.0 / R + 2.0 * slack,
        "{cert:?}"
    );
}
