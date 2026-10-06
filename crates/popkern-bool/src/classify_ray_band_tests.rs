//! The ray edge band is capped by the local knot span (PK-679).

use super::knot_span_width;

#[test]
fn span_width_is_the_distinct_span_holding_the_parameter() {
    // Clamped cubic over three spans.
    let knots = [0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 1.0, 1.0, 1.0, 1.0];
    assert_eq!(knot_span_width(&knots, 0.1), 0.25);
    assert_eq!(knot_span_width(&knots, 0.3), 0.25);
    assert_eq!(knot_span_width(&knots, 0.75), 0.5);
    // Out of range and on the last knot: the end spans.
    assert_eq!(knot_span_width(&knots, -0.001), 0.25);
    assert_eq!(knot_span_width(&knots, 1.0), 0.5);
    assert_eq!(knot_span_width(&knots, 1.2), 0.5);
}

#[test]
fn a_long_coil_gets_a_band_far_under_half_a_percent() {
    // A 112-turn wall fitted at 24 spans per turn: the tuned 0.5%
    // band was over half a turn; the span cap keeps it to an eighth
    // of a span (~1/100 of a turn).
    let spans = 112 * 24;
    let mut knots = vec![0.0; 4];
    knots.extend((1..spans).map(|i| i as f64 / spans as f64));
    knots.extend([1.0; 4]);
    let band = (0.005f64).min(0.25 * knot_span_width(&knots, 0.0));
    assert!(
        band * 112.0 < 0.011,
        "band {band} covers {} turns",
        band * 112.0
    );
    // A face of few spans keeps the tuned band.
    let few = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    assert_eq!((0.005f64).min(0.25 * knot_span_width(&few, 0.3)), 0.005);
}

#[test]
fn degenerate_knot_vectors_impose_no_cap() {
    assert!(knot_span_width(&[], 0.5).is_infinite());
    assert!(knot_span_width(&[0.0, 0.0], 0.5).is_infinite());
}
