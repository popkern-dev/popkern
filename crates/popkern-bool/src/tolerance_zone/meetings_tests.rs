use std::f64::consts::FRAC_1_SQRT_2;

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::*;
use crate::tolerance_zone::tests::record;

const R: f64 = 2e-6;

/// The exact rational circle of radius `r` about `c` in the xy plane.
fn circle(c: (f64, f64), r: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    let p = |x: f64, y: f64| Point3::new(c.0 + r * x, c.1 + r * y, 0.0);
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            p(1.0, 0.0),
            p(1.0, 1.0),
            p(0.0, 1.0),
            p(-1.0, 1.0),
            p(-1.0, 0.0),
            p(-1.0, -1.0),
            p(0.0, -1.0),
            p(1.0, -1.0),
            p(1.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Villarceau's circles in their plane: radius 5 about (0, ±3), crossing
/// at (±4, 0). Two meetings, both crossings.
#[test]
fn two_crossing_circles_meet_twice() {
    let records = [
        record(circle((0.0, 3.0), 5.0)),
        record(circle((0.0, -3.0), 5.0)),
    ];
    let m = record_meetings(&records, R);
    assert_eq!(m.len(), 2, "{m:?}");
    for x in [4.0, -4.0] {
        let at = Point3::new(x, 0.0, 0.0);
        let found = m
            .iter()
            .find(|k| k.point.distance_to(&at) < 1e-7)
            .expect("a meeting");
        assert!(found.sin_angle > 0.5, "{found:?}");
        assert!(
            found.zone < 1e-5,
            "a crossing's zone is the radius over its sine"
        );
        assert_eq!(found.records, [0, 1]);
    }
}

/// Circles of radius 5 about (0, 0) and (10, 0) touch at (5, 0): one
/// meeting, tangential, whose zone is the touch window √(2r/κ) wide,
/// far wider than the radius.
#[test]
fn touching_circles_meet_once_with_a_wide_zone() {
    let records = [
        record(circle((0.0, 0.0), 5.0)),
        record(circle((10.0, 0.0), 5.0)),
    ];
    let m = record_meetings(&records, R);
    assert_eq!(m.len(), 1, "{m:?}");
    let k = &m[0];
    assert!(
        k.point.distance_to(&Point3::new(5.0, 0.0, 0.0)) < 1e-3,
        "{k:?}"
    );
    assert!(k.sin_angle < 1e-2, "{k:?}");
    // Relative curvature 2/5: the window is about √(2 · 2e-6 / 0.4) = 3e-3.
    assert!(k.zone > 1e-3 && k.zone < 1e-2, "{k:?}");
}

/// A figure eight (a polyline through (sin t, sin 2t)) crosses itself
/// once, at the origin: one meeting of the record with itself.
#[test]
fn a_figure_eight_meets_itself_once() {
    use popkern_core::geometry::PortableMath;
    let n = 400;
    let pts: Vec<Point3> = (0..=n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            Point3::new(5.0 * t.psin(), 5.0 * (2.0 * t).psin(), 0.0)
        })
        .collect();
    let mut knots = vec![0.0];
    knots.extend((0..=n).map(|k| k as f64 / n as f64));
    knots.push(1.0);
    let curve = NurbsCurve {
        degree: 1,
        knots,
        weights: vec![1.0; pts.len()],
        control_points: pts,
    };
    let m = record_meetings(&[record(curve)], R);
    assert_eq!(m.len(), 1, "{m:?}");
    assert_eq!(m[0].records, [0, 0]);
    assert!(
        m[0].point.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-6,
        "{m:?}"
    );
    assert!(m[0].params[0] < m[0].params[1]);
}

/// Circles apart do not meet; records on no common face are not
/// compared at all.
#[test]
fn apart_or_on_no_common_face_no_meeting() {
    let apart = [
        record(circle((0.0, 0.0), 5.0)),
        record(circle((20.0, 0.0), 5.0)),
    ];
    assert!(record_meetings(&apart, R).is_empty());
    let mut other = record(circle((0.0, -3.0), 5.0));
    let key =
        |n: u64| -> popkern_core::topology::FaceHandle { slotmap::KeyData::from_ffi(n).into() };
    other.face_a = key((1 << 32) | 7);
    other.face_b = key((1 << 32) | 9);
    let crossing = [record(circle((0.0, 3.0), 5.0)), other];
    assert!(record_meetings(&crossing, R).is_empty());
}

/// A closed record whose ends differ by its fit (a circle rotated to
/// start at a vertex and refitted) is closed for separation along it:
/// its two ends are no circumference apart. Negative control: with the
/// closing distance zero, the ends read a whole circle apart.
#[test]
fn a_refit_closed_record_closes_within_the_radius() {
    let c = circle((0.0, 0.0), 5.0);
    let (a, b) = (0.0, 1.0);
    // The end moved off the start by 1e-8, as a refit leaves it.
    let point_at = |s: f64| {
        let p = c.point_at(s);
        if s >= 1.0 {
            Point3::new(p.x, p.y + 1e-8, p.z)
        } else {
            p
        }
    };
    let close = super::super::resolve_record::separation(point_at, (a, b), (0.0, 1.0), 2e-6);
    assert!(close < 1e-6, "{close}");
    let open = super::super::resolve_record::separation(point_at, (a, b), (0.0, 1.0), 0.0);
    assert!(open > 31.0, "{open}");
}
