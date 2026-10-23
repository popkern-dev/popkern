//! Wrapped-window extraction (`clip_curve_wrapped`) and window
//! extraction (`clip_curve_to_crossings`) tests.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, PortableMath};
use std::f64::consts::TAU;

use super::{clip_curve_to_crossings, clip_curve_wrapped};

/// A periodic chart pcurve: u sweeps one full period while the 3D
/// curve closes, so the pcurve's END sits one period from its START.
fn azimuth_line(v: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, v, 0.0), Point3::new(TAU, v, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

#[test]
fn a_wrapped_window_of_a_periodic_pcurve_stays_continuous_through_the_seam() {
    // Window [0.8, 1.2] of the closed record runs through the seam:
    // the kept pcurve must sweep u ∈ [0.8·2π, 1.2·2π] monotonically,
    // not fold its tail back across the chart (it used to run
    // 2π → 0.2·2π backwards).
    let Some(piece) = clip_curve_wrapped(&azimuth_line(2.0), 0.8, 1.2) else {
        panic!("wrapped window must extract");
    };
    let (t0, t1) = piece.parameter_range();
    let mut prev = f64::NEG_INFINITY;
    for i in 0..=40 {
        let p = piece.point_at(t0 + (t1 - t0) * i as f64 / 40.0);
        assert!(
            p.x >= prev - 1e-12,
            "u must not fold back: {} after {prev}",
            p.x
        );
        assert!((p.y - 2.0).abs() < 1e-12);
        prev = p.x;
    }
    assert!((piece.point_at(t0).x - 0.8 * TAU).abs() < 1e-9);
    assert!((piece.point_at(t1).x - 1.2 * TAU).abs() < 1e-9);
}

#[test]
fn a_wrapped_window_of_a_closed_3d_curve_is_not_translated() {
    // A closed 3D curve's closure is round-off, never a period: the
    // spliced tail keeps its place on the curve.
    let circle = |t: f64| Point3::new(t.pcos(), t.psin(), 0.0);
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let ring = NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };
    let Some(piece) = clip_curve_wrapped(&ring, 0.75, 1.25) else {
        panic!("wrapped window must extract");
    };
    let (t0, t1) = piece.parameter_range();
    for i in 0..=20 {
        let p = piece.point_at(t0 + (t1 - t0) * i as f64 / 20.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 1e-12, "piece left the ring: r = {r}");
    }
    // Quarter-span knots: t = 1.25 wraps to t = 0.25, the +y point.
    assert!(piece.point_at(t1).distance_to(&circle(0.25 * TAU)) < 1e-9);
}

#[test]
fn a_window_of_a_polyline_keeps_its_vertices() {
    // A fitted sphere-chart image is a degree-1 POLYLINE. Its window
    // used to come back as the chord between the window's ends: a
    // meridian from the pole read as a diagonal across the chart, 0.35
    // off its curve at mid-span (PK-930). The window keeps the corner.
    let polyline = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0; 3],
    };
    let Some(piece) = clip_curve_to_crossings(&polyline, 0.25, 0.75) else {
        panic!("window must extract");
    };
    let (t0, t1) = piece.parameter_range();
    let ends = (piece.point_at(t0), piece.point_at(t1));
    assert!(ends.0.distance_to(&Point3::new(0.5, 0.0, 0.0)) < 1e-12);
    assert!(ends.1.distance_to(&Point3::new(1.0, 0.5, 0.0)) < 1e-12);
    let corner = piece.point_at(0.5 * (t0 + t1));
    assert!(
        corner.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-12,
        "window mid-span {corner:?} left the polyline's corner"
    );
}
