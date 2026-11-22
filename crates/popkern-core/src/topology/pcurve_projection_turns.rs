//! Turns of a projected pcurve's chart track that local refinement
//! waits through (split from `pcurve_projection_fit.rs`, 500-line cap).

use crate::geometry::Point3;

/// Whether the chart step from `a` to `b` turns a periodic coordinate
/// (periods `periods`) past a quarter period.
pub(super) fn turns_past_quarter(
    periods: (Option<f64>, Option<f64>),
    a: (f64, f64),
    b: (f64, f64),
) -> bool {
    [(b.0 - a.0, periods.0), (b.1 - a.1, periods.1)]
        .iter()
        .any(|&(d, p)| p.is_some_and(|p| p > 0.0 && d.abs() > 0.25 * p))
}

/// Whether the track through three samples turns past 60° (not across a period).
pub(super) fn turns_sharply(periods: (Option<f64>, Option<f64>), w: &[Point3]) -> bool {
    let (a, b, c) = ((w[0].x, w[0].y), (w[1].x, w[1].y), (w[2].x, w[2].y));
    if turns_past_quarter(periods, a, b) || turns_past_quarter(periods, b, c) {
        return false;
    }
    let (s, t) = ((b.0 - a.0, b.1 - a.1), (c.0 - b.0, c.1 - b.1));
    let lens = (s.0 * s.0 + s.1 * s.1).sqrt() * (t.0 * t.0 + t.1 * t.1).sqrt();
    lens > 0.0 && (s.0 * t.0 + s.1 * t.1) < 0.5 * lens
}
