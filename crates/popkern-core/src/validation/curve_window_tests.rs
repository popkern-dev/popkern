use super::*;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn a_window_inside_its_span_passes_up_to_whole_turns_and_either_order() {
    let span = (-FRAC_PI_2, PI);
    assert!(window_within_span((0.0, FRAC_PI_2), span, 1e-9));
    assert!(window_within_span((FRAC_PI_2, 0.0), span, 1e-9));
    assert!(window_within_span((TAU, TAU + FRAC_PI_2), span, 1e-9));
    assert!(window_within_span(
        (-TAU - FRAC_PI_2, -TAU + PI),
        span,
        1e-9
    ));
    // A span declared backward.
    let backward = (0.0, -FRAC_PI_2);
    assert!(window_within_span((0.0, -FRAC_PI_2), backward, 1e-9));
    assert!(!window_within_span((0.0, FRAC_PI_2), backward, 1e-9));
}

#[test]
fn a_window_on_the_complement_is_flagged() {
    // The pocketed tray's fillet end arc: span [-pi/2, pi], window the
    // missing quarter [pi, 3pi/2].
    let span = (-FRAC_PI_2, PI);
    assert!(!window_within_span((PI, 1.5 * PI), span, 1e-9));
    // Half inside, half out.
    assert!(!window_within_span((FRAC_PI_2, 1.25 * PI), span, 1e-9));
    // Out by more than the slack only.
    assert!(window_within_span((-FRAC_PI_2 - 1e-10, 0.0), span, 1e-9));
    assert!(!window_within_span((-FRAC_PI_2 - 1e-6, 0.0), span, 1e-9));
    assert!(!window_within_span((f64::NAN, 0.0), span, 1e-9));
}

#[test]
fn spanning_declares_the_window_on_angular_curves_only() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.5,
        start_angle: -FRAC_PI_2,
        end_angle: PI,
    };
    let moved = arc.spanning((1.5 * PI, PI));
    assert_eq!(moved.angular_span(), Some((PI, 1.5 * PI)));
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(matches!(
        line.spanning((0.0, 1.0)),
        CurveGeometry::Line { origin, .. } if origin == Point3::new(0.0, 0.0, 0.0)
    ));
}
