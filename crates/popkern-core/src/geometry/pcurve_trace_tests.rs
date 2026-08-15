//! Tests for refined pcurve-trace nearest-point measurement (PK-287).

use super::nearest_on_trace;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

fn unit_z_cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

/// Full-period rim circle on a cylinder as a Linear pcurve in UV.
fn rim_pcurve(v: f64) -> Pcurve {
    Pcurve {
        geometry: PcurveGeometry::Linear {
            start: (0.0, v),
            end: (std::f64::consts::TAU, v),
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// The old nearest-of-samples scan was floored at half the sample
/// spacing (~π·r/128 on a full rim). The refined measurement must
/// certify an on-trace point to far below any stale tolerance,
/// including between coarse samples.
#[test]
fn on_trace_points_measure_at_machine_precision() {
    let radius = 0.3;
    let surface = unit_z_cylinder(radius);
    let pcurve = rim_pcurve(1.0);
    for i in 0..37 {
        // Angles deliberately incommensurate with the 128 coarse samples.
        let angle = 0.171 + std::f64::consts::TAU * (i as f64 / 37.0);
        let point = surface.point_at(angle, 1.0);
        let nearest = nearest_on_trace(&pcurve, &surface, &point, 128);
        assert!(
            nearest.distance < 1.0e-9,
            "angle {angle}: refined distance {:.3e} should be at machine precision",
            nearest.distance
        );
    }
}

/// Off-trace points must measure their true analytic distance, not a
/// sample-floored overestimate.
#[test]
fn off_trace_distance_matches_analytic_value() {
    let radius = 0.3;
    let surface = unit_z_cylinder(radius);
    let pcurve = rim_pcurve(0.5);
    // Radially outward from the rim at an arbitrary angle: distance is
    // exactly the radial offset.
    let angle = 2.31_f64;
    let offset = 0.017;
    let dir = Vec3::new(angle.cos(), angle.sin(), 0.0);
    let point = Point3::new(0.0, 0.0, 0.5) + dir * (radius + offset);
    let nearest = nearest_on_trace(&pcurve, &surface, &point, 128);
    assert!(
        (nearest.distance - offset).abs() < 1.0e-9,
        "distance {:.6e} should equal the radial offset {offset}",
        nearest.distance
    );
}

/// A query point whose nearest trace point sits just inside the far
/// end of the seam: the coarse best lands on the first sample, and the
/// refinement must still find the minimum near the opposite end.
#[test]
fn seam_straddling_minimum_is_found() {
    let radius = 0.3;
    let surface = unit_z_cylinder(radius);
    let pcurve = rim_pcurve(0.0);
    // Just below the seam angle 0 (i.e. at TAU - small): nearest trace
    // parameter is near t = 1.
    let angle = std::f64::consts::TAU - 1.0e-3;
    let point = surface.point_at(angle, 0.0);
    let nearest = nearest_on_trace(&pcurve, &surface, &point, 128);
    assert!(
        nearest.distance < 1.0e-9,
        "seam-adjacent point measured {:.3e}",
        nearest.distance
    );
}

/// The returned parameter must locate the nearest point, not merely
/// its distance: consumers use the parameter sequence for direction
/// agreement.
#[test]
fn returned_parameter_locates_the_nearest_point() {
    let radius = 0.3;
    let surface = unit_z_cylinder(radius);
    let pcurve = rim_pcurve(0.0);
    let angle = 1.234_f64;
    let point = surface.point_at(angle, 0.0);
    let nearest = nearest_on_trace(&pcurve, &surface, &point, 128);
    let expected_t = angle / std::f64::consts::TAU;
    assert!(
        (nearest.t - expected_t).abs() < 1.0e-6,
        "nearest.t {:.8} should be the angle fraction {:.8}",
        nearest.t,
        expected_t
    );
}

/// Non-finite parameter ranges decline with an infinite distance
/// rather than sampling garbage.
#[test]
fn non_finite_range_declines() {
    let pcurve = Pcurve {
        geometry: PcurveGeometry::Linear {
            start: (0.0, 0.0),
            end: (1.0, 0.0),
        },
        t_start: 0.0,
        t_end: f64::INFINITY,
    };
    let surface = unit_z_cylinder(1.0);
    let nearest = nearest_on_trace(&pcurve, &surface, &Point3::new(0.0, 0.0, 0.0), 32);
    assert!(nearest.distance.is_infinite());
}

/// Descending parameter ranges (twin-side and STEP descending-domain
/// pcurves) must refine exactly like ascending ones: the coarse best
/// is always kept as a candidate, and bracketing runs on the
/// ascending hull.
#[test]
fn descending_range_refines_exactly() {
    let radius = 0.3;
    let surface = unit_z_cylinder(radius);
    let pcurve = Pcurve {
        geometry: PcurveGeometry::Linear {
            start: (0.0, 1.0),
            end: (std::f64::consts::TAU, 1.0),
        },
        t_start: 1.0,
        t_end: 0.0,
    };
    for i in 0..17 {
        let angle = 0.113 + std::f64::consts::TAU * (i as f64 / 17.0);
        let point = surface.point_at(angle, 1.0);
        let nearest = nearest_on_trace(&pcurve, &surface, &point, 128);
        assert!(
            nearest.distance < 1.0e-9,
            "angle {angle}: descending-range refined distance {:.3e}",
            nearest.distance
        );
    }
}
