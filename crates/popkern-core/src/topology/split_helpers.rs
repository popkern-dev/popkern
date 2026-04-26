//! Local helpers for split.rs — kept in a sibling file to fit the
//! 500-line per-file structural cap. Used only by `split.rs` and
//! `split_boundary.rs` siblings; not part of the public API.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};

/// Split a pcurve at parameter `t` (in [0,1] of the edge).
///
/// Returns `(left, right)` where `left` covers the first portion and
/// `right` covers the second. If the pcurve is `None` or splitting
/// fails, returns `(None, None)`. Handles both `Linear` and `Nurbs`
/// pcurve variants.
pub(crate) fn split_pcurve_at(pcurve: &Option<Pcurve>, t: f64) -> (Option<Pcurve>, Option<Pcurve>) {
    let Some(pc) = pcurve else {
        return (None, None);
    };
    let (p_start, p_end) = (pc.t_start, pc.t_end);
    let split_param = p_start + t * (p_end - p_start);
    match &pc.geometry {
        PcurveGeometry::Linear { start, end } => {
            let span = p_end - p_start;
            let s = if span.abs() > 0.0 { t } else { 0.0 };
            let mid = (
                start.0 + s * (end.0 - start.0),
                start.1 + s * (end.1 - start.1),
            );
            (
                Some(Pcurve::linear(*start, mid, p_start, split_param)),
                Some(Pcurve::linear(mid, *end, split_param, p_end)),
            )
        }
        PcurveGeometry::Elliptical {
            center,
            u_radius,
            v_radius,
            rotation,
            start_angle,
            end_angle,
        } => {
            // Angle interpolates linearly with normalised parameter.
            let mid_angle = start_angle + t * (end_angle - start_angle);
            let left = Pcurve::elliptical(
                *center,
                *u_radius,
                *v_radius,
                *rotation,
                *start_angle,
                mid_angle,
                p_start,
                split_param,
            );
            let right = Pcurve::elliptical(
                *center,
                *u_radius,
                *v_radius,
                *rotation,
                mid_angle,
                *end_angle,
                split_param,
                p_end,
            );
            (Some(left), Some(right))
        }
        PcurveGeometry::Sinusoidal {
            u_offset,
            v_const,
            v_cos,
            v_sin,
        } => {
            // Sinusoidal form is closed-form on any sub-range —
            // restrict the parameter window without changing the
            // analytic formula.
            let left =
                Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, p_start, split_param);
            let right = Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, split_param, p_end);
            (Some(left), Some(right))
        }
        PcurveGeometry::Nurbs(curve) => match curve.split_at(split_param) {
            Ok((left, right)) => (
                Some(Pcurve::from_nurbs(left)),
                Some(Pcurve::from_nurbs(right)),
            ),
            Err(_) => (None, None),
        },
    }
}

/// Reverse a NURBS curve (swap direction of traversal).
///
/// The reversed curve evaluates the same geometry but in the opposite
/// direction: `reversed.point_at(t_start)` equals
/// `original.point_at(t_end)`.
pub(crate) fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    use crate::geometry::curve::CurveEval;
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();

    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();

    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}
