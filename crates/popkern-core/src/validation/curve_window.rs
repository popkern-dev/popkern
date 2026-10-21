//! Edge windows inside their curves' declared spans (PK-906).
//!
//! An Arc or Ellipse declares an angular span, and the edge's window
//! `[t_start, t_end]` must lie inside it, up to whole turns. Evaluation
//! is periodic, so a window on the complement still traces the edge and
//! every other check passes; but `closest_point` clamps into the
//! declared span and `parameter_range` reports it, so a projection onto
//! the edge lands on the wrong end (measured by the PK-505 chordal
//! certificate on the pocketed tray's fillet end arcs).

use std::f64::consts::TAU;

use crate::error::ValidationError;
use crate::geometry::CurveGeometry;
use crate::topology::BRep;

impl BRep {
    /// Every edge whose Arc or Ellipse window leaves the curve's
    /// declared span by more than the body's linear tolerance (as an
    /// angle at the curve's smallest radius).
    pub fn validate_curve_window_all(&self) -> Vec<ValidationError> {
        let linear = self.tolerance.linear();
        let mut errors = Vec::new();
        for (edge, e) in self.edges.iter() {
            let Some(curve) = e.curve.as_ref() else {
                continue;
            };
            let Some((span_start, span_end)) = curve.geometry.angular_span() else {
                continue;
            };
            let radius = match &curve.geometry {
                CurveGeometry::Arc { radius, .. } => *radius,
                CurveGeometry::Ellipse {
                    u_basis, v_basis, ..
                } => u_basis.length().min(v_basis.length()),
                _ => continue,
            };
            let slack = linear / radius.max(linear);
            let window = (curve.t_start, curve.t_end);
            if !window_within_span(window, (span_start, span_end), slack) {
                errors.push(ValidationError::CurveWindowOutsideSpan {
                    edge,
                    t_start: window.0,
                    t_end: window.1,
                    span_start,
                    span_end,
                });
            }
        }
        errors
    }
}

/// Whether `window` lies inside `span` (each in either order) shifted by some
/// whole number of turns, within `slack` radians at each end.
fn window_within_span(window: (f64, f64), span: (f64, f64), slack: f64) -> bool {
    // Before min/max, which drop a NaN.
    if ![window.0, window.1, span.0, span.1]
        .iter()
        .all(|x| x.is_finite())
    {
        return false;
    }
    let (lo, hi) = (window.0.min(window.1), window.0.max(window.1));
    // A span may be declared either way round (a producer sweeping
    // backward writes [0, -pi/2]).
    let span = (span.0.min(span.1), span.0.max(span.1));
    // The smallest turn lifting the window's start into the span; a
    // larger one only pushes its end further out.
    let turns = ((span.0 - slack - lo) / TAU).ceil();
    hi + turns * TAU <= span.1 + slack
}

#[cfg(test)]
#[path = "curve_window_tests.rs"]
mod tests;
