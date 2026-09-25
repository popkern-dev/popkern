//! SameParameter repair: a foreign pcurve longer than its edge, or
//! parameterized differently, is rebuilt on the edge's parameterization;
//! a pcurve whose image leaves the edge is not.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, SQRT_2};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::PcurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;

use super::*;

const TOL: f64 = 1e-6;

fn xy_plane() -> SurfaceOnFace {
    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    }
}

/// The standard 9-point rational full circle of radius `r` about the
/// origin, on [0, 1] (as OpenCASCADE writes a plane circle's pcurve):
/// exact, but not parameterized by angle.
fn rational_full_circle(r: f64) -> Pcurve {
    let w = SQRT_2 / 2.0;
    let cps = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    Pcurve::from_nurbs(NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: cps
            .iter()
            .map(|&(x, y)| Point3::new(r * x, r * y, 0.0))
            .collect(),
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    })
}

/// An arc edge of radius 3 about the origin in the XY plane.
fn arc_edge(start_angle: f64, end_angle: f64) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    }
}

/// Same-fraction distance between the rebuilt pcurve's image and the
/// edge in half-edge direction, over a dense probe set.
fn same_fraction_gap(pcurve: &Pcurve, curve: &CurveOnEdge, forward: bool) -> f64 {
    let (e0, e1) = if forward {
        (curve.t_start, curve.t_end)
    } else {
        (curve.t_end, curve.t_start)
    };
    let plane = xy_plane();
    (0..=500)
        .map(|k| {
            let f = k as f64 / 500.0;
            let uv = pcurve.point_at(pcurve.t_start + f * (pcurve.t_end - pcurve.t_start));
            let on_surface = plane.geometry.point_at(uv.x, uv.y);
            on_surface.distance_to(&curve.geometry.point_at(e0 + f * (e1 - e0)))
        })
        .fold(0.0, f64::max)
}

#[test]
fn a_full_rational_circle_is_trimmed_and_reparameterized_onto_a_quarter_arc() {
    let edge = arc_edge(0.3, 0.3 + FRAC_PI_2);
    for forward in [true, false] {
        let repair =
            same_parameter_repair(&rational_full_circle(3.0), &edge, &xy_plane(), forward, TOL)
                .unwrap_or_else(|| panic!("forward={forward}: repair declined"));
        assert!(
            repair.max_gap <= TOL,
            "forward={forward}: {}",
            repair.max_gap
        );
        assert!(matches!(repair.pcurve.geometry, PcurveGeometry::Nurbs(_)));
        let dense = same_fraction_gap(&repair.pcurve, &edge, forward);
        assert!(dense <= 2.0 * TOL, "forward={forward}: dense {dense:e}");
    }
    // Negative control: the file's pcurve as written fails the
    // same-fraction check by the full arc length's mismatch.
    let raw = same_fraction_gap(&rational_full_circle(3.0), &edge, true);
    assert!(
        raw > 1.0,
        "raw pcurve unexpectedly co-parameterized ({raw:e})"
    );
}

#[test]
fn an_edge_across_the_pcurve_seam_is_followed_through_the_wrap() {
    // The rational circle starts at angle 0; this edge crosses it.
    let edge = arc_edge(-FRAC_PI_4, FRAC_PI_4);
    let repair = same_parameter_repair(&rational_full_circle(3.0), &edge, &xy_plane(), true, TOL)
        .unwrap_or_else(|| panic!("repair declined across the seam"));
    assert!(repair.max_gap <= TOL, "{}", repair.max_gap);
    assert!(same_fraction_gap(&repair.pcurve, &edge, true) <= 2.0 * TOL);
}

#[test]
fn a_pcurve_off_the_edge_is_not_repaired() {
    // Radius 3.1 against a radius-3 edge: 0.1 away everywhere.
    let edge = arc_edge(0.3, 1.2);
    assert!(
        same_parameter_repair(&rational_full_circle(3.1), &edge, &xy_plane(), true, TOL).is_none()
    );
}

#[test]
fn the_tally_reports_once_with_the_worst_residual() {
    let mut tally = SameParameterTally::default();
    let mut log = DiagnosticLog::new();
    tally.report(&mut log);
    assert!(log.is_empty(), "nothing rebuilt, nothing reported");
    tally.record(2e-8);
    tally.record(5e-7);
    tally.report(&mut log);
    let entries = log.entries();
    assert_eq!(entries.len(), 1);
    assert!(
        entries[0].message.contains("rebuilt 2 STEP pcurve(s)"),
        "{}",
        entries[0].message
    );
    assert!(
        entries[0].message.contains("5.000e-7"),
        "{}",
        entries[0].message
    );
}

#[test]
fn a_unit_segment_line_pcurve_is_extended_along_a_long_edge() {
    // A STEP LINE resolves to its unit segment; the edge runs 60 units
    // along it, starting 5 units past the segment.
    let pcurve = Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0);
    let edge = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(5.0, 0.0, 0.0),
            direction: Vec3::new(60.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    for forward in [true, false] {
        let repair = same_parameter_repair(&pcurve, &edge, &xy_plane(), forward, TOL)
            .unwrap_or_else(|| panic!("forward={forward}: repair declined"));
        assert!(repair.max_gap <= TOL, "{}", repair.max_gap);
        assert!(same_fraction_gap(&repair.pcurve, &edge, forward) <= 2.0 * TOL);
    }
}
