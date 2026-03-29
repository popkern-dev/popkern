//! Residual and Jacobian evaluation for spline endpoint constraints.
//!
//! Implements PointOnSpline, SplineEndTangent, SplineEndTangentToArc,
//! SplineFixedEndTangent, and SplineFixedEndCurvature.

use super::spline_deriv::{eval_curvature, eval_tangent};
use super::spline_eval::eval_point;
use super::{ar, ln, pt, EquationOutput};
use crate::constraint::ConstraintHandle;
use crate::elements::{ArcHandle, LineHandle, PointHandle, SplineHandle};
use crate::Sketch;

/// Compute the parameter value for a spline endpoint.
///
/// `end_idx` 0 = start (t = knots[degree]), `end_idx` 1 = end (t = knots[n_cp]).
fn endpoint_t(sketch: &Sketch, spline_handle: SplineHandle, end_idx: usize) -> f64 {
    let sp = &sketch.splines[spline_handle];
    if end_idx == 0 {
        sp.knots[sp.degree]
    } else {
        sp.knots[sp.control_points.len()]
    }
}

// ── 1. PointOnSpline ──────────────────────────────────────────────────────

/// `PointOnSpline(Point, Spline)` — 2 equations, 1 aux param (t).
///
/// ```text
/// F1 = S_x(t) - px = 0
/// F2 = S_y(t) - py = 0
/// ```
pub(crate) fn eval_point_on_spline(
    sketch: &Sketch,
    point: PointHandle,
    spline: SplineHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, point);
    let sp = &sketch.splines[spline];

    // Look up aux parameter (t).
    let aux_idx = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t = sketch.params[aux_idx];

    let eval = eval_point(sp, t, sketch);

    let mut jacobian = Vec::new();

    // dF1/dpx = -1, dF2/dpy = -1
    jacobian.push((row, ipx, -1.0));
    jacobian.push((row + 1, ipy, -1.0));

    // Control point contributions from eval_point Jacobian.
    for &(param_idx, dx_dp, dy_dp) in &eval.jacobian {
        if dx_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dx_dp));
        }
        if dy_dp.abs() > f64::EPSILON {
            jacobian.push((row + 1, param_idx, dy_dp));
        }
    }

    // Auxiliary parameter t: dF1/dt = S'_x(t), dF2/dt = S'_y(t)
    let tang = eval_tangent(sp, t, sketch);
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row, aux_idx, tang.x));
    }
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_idx, tang.y));
    }

    EquationOutput {
        residuals: vec![eval.x - px, eval.y - py],
        jacobian,
    }
}

// ── 2. SplineEndTangent ───────────────────────────────────────────────────

/// `SplineEndTangent(Spline, end_idx, Line)` — 1 equation.
///
/// The spline tangent at the endpoint is parallel to the line direction:
/// ```text
/// F = S'_x(t_end) * dy - S'_y(t_end) * dx = 0
/// ```
/// where `(dx, dy) = (x2 - x1, y2 - y1)` from the line.
pub(crate) fn eval_spline_end_tangent(
    sketch: &Sketch,
    spline: SplineHandle,
    end_idx: usize,
    line: LineHandle,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let t_end = endpoint_t(sketch, spline, end_idx);
    let tang = eval_tangent(sp, t_end, sketch);

    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, line);
    let dx = x2 - x1;
    let dy = y2 - y1;

    // F = tang.x * dy - tang.y * dx
    let residual = tang.x * dy - tang.y * dx;

    let mut jacobian = Vec::new();

    // dF/d(control_points): from tangent Jacobian
    for &(param_idx, dtx_dp, dty_dp) in &tang.jacobian {
        let val = dtx_dp * dy - dty_dp * dx;
        if val.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, val));
        }
    }

    // dF/d(line params):
    // dF/dx1 = tang.x * (-1) ... no wait:
    // dx = x2 - x1, dy = y2 - y1
    // dF/dx1 = tang.x * d(dy)/dx1 - tang.y * d(dx)/dx1 = 0 - tang.y * (-1) = tang.y
    // dF/dy1 = tang.x * (-1) - tang.y * 0 = -tang.x
    // dF/dx2 = 0 - tang.y * 1 = -tang.y
    // dF/dy2 = tang.x * 1 - 0 = tang.x
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row, ix1, tang.y));
    }
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row, iy1, -tang.x));
    }
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row, ix2, -tang.y));
    }
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row, iy2, tang.x));
    }

    EquationOutput {
        residuals: vec![residual],
        jacobian,
    }
}

// ── 3. SplineEndTangentToArc ──────────────────────────────────────────────

/// `SplineEndTangentToArc(Spline, end_idx, Arc)` — 1 equation.
///
/// The spline tangent at the endpoint is perpendicular to the radius
/// from the arc centre to the spline endpoint position:
/// ```text
/// F = S'_x(t_end) * (py - cy) - S'_y(t_end) * (px - cx) = 0
/// ```
/// where `(px, py)` is the spline endpoint position and `(cx, cy)` is the
/// arc centre.
pub(crate) fn eval_spline_end_tangent_to_arc(
    sketch: &Sketch,
    spline: SplineHandle,
    end_idx: usize,
    arc: ArcHandle,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let t_end = endpoint_t(sketch, spline, end_idx);
    let tang = eval_tangent(sp, t_end, sketch);
    let pos = eval_point(sp, t_end, sketch);

    let (cx, cy, _r, icx, icy, _ir) = ar(sketch, arc);

    let rx = pos.x - cx;
    let ry = pos.y - cy;

    // F = tang.x * ry - tang.y * rx
    let residual = tang.x * ry - tang.y * rx;

    let mut jacobian = Vec::new();

    // dF/d(control_points): contributions from both tangent and position.
    // From tangent: dF = d(tang.x) * ry - d(tang.y) * rx
    // From position (via rx, ry): dF = tang.x * d(pos.y) - tang.y * d(pos.x)
    for &(param_idx, dtx_dp, dty_dp) in &tang.jacobian {
        let val = dtx_dp * ry - dty_dp * rx;
        if val.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, val));
        }
    }
    for &(param_idx, dpx_dp, dpy_dp) in &pos.jacobian {
        let val = tang.x * dpy_dp - tang.y * dpx_dp;
        if val.abs() > f64::EPSILON {
            // Merge with existing entry or add new.
            merge_or_push(&mut jacobian, row, param_idx, val);
        }
    }

    // dF/d(arc centre):
    // dF/dcx = tang.x * d(ry)/dcx - tang.y * d(rx)/dcx
    //        = tang.x * 0 - tang.y * (-1) = tang.y
    // dF/dcy = tang.x * (-1) - tang.y * 0 = -tang.x
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row, icx, tang.y));
    }
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row, icy, -tang.x));
    }

    EquationOutput {
        residuals: vec![residual],
        jacobian,
    }
}

// ── 4. SplineFixedEndTangent ──────────────────────────────────────────────

/// `SplineFixedEndTangent(Spline, end_idx, tx, ty)` — 2 equations.
///
/// ```text
/// F1 = S'_x(t_end) - tx = 0
/// F2 = S'_y(t_end) - ty = 0
/// ```
pub(crate) fn eval_spline_fixed_end_tangent(
    sketch: &Sketch,
    spline: SplineHandle,
    end_idx: usize,
    fixed_values: &[f64],
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let t_end = endpoint_t(sketch, spline, end_idx);
    let tang = eval_tangent(sp, t_end, sketch);

    let tx = fixed_values[0];
    let ty = fixed_values[1];

    let mut jacobian = Vec::new();

    // dF/d(control_points) from tangent Jacobian.
    for &(param_idx, dtx_dp, dty_dp) in &tang.jacobian {
        if dtx_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dtx_dp));
        }
        if dty_dp.abs() > f64::EPSILON {
            jacobian.push((row + 1, param_idx, dty_dp));
        }
    }

    EquationOutput {
        residuals: vec![tang.x - tx, tang.y - ty],
        jacobian,
    }
}

// ── 5. SplineFixedEndCurvature ────────────────────────────────────────────

/// `SplineFixedEndCurvature(Spline, end_idx, value)` — 1 equation.
///
/// ```text
/// F = kappa(t_end) - value = 0
/// ```
pub(crate) fn eval_spline_fixed_end_curvature(
    sketch: &Sketch,
    spline: SplineHandle,
    end_idx: usize,
    value: f64,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let t_end = endpoint_t(sketch, spline, end_idx);
    let (kappa, kappa_jac) = eval_curvature(sp, t_end, sketch);

    let mut jacobian = Vec::new();
    for &(param_idx, dk_dp) in &kappa_jac {
        if dk_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dk_dp));
        }
    }

    EquationOutput {
        residuals: vec![kappa - value],
        jacobian,
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────

/// Merge a Jacobian entry into the vector, adding to an existing entry for
/// the same (row, col) pair or pushing a new one.
fn merge_or_push(jac: &mut Vec<(usize, usize, f64)>, row: usize, col: usize, val: f64) {
    if let Some(entry) = jac.iter_mut().find(|(r, c, _)| *r == row && *c == col) {
        entry.2 += val;
    } else {
        jac.push((row, col, val));
    }
}
