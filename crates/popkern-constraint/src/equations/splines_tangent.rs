//! Residual and Jacobian evaluation for spline-to-line and spline-to-arc
//! interior tangency constraints.
//!
//! Spline-to-conic tangency (ellipse, parabola, hyperbola) is in the
//! sibling [`splines_tangent_conic`] module, split for file-size reasons.

use super::spline_deriv::eval_tangent;
use super::spline_eval::eval_point;
use super::splines::merge_or_push;
use super::{ar, ln, EquationOutput};
use crate::constraint::ConstraintHandle;
use crate::elements::{ArcHandle, LineHandle, SplineHandle};
use crate::Sketch;

// ── 1. SplineTangentToLine ───────────────────────────────────────────────

/// `SplineTangentToLine(Spline, Line)` -- 2 equations, 1 aux param (t).
///
/// ```text
/// F1 = (S_x(t) - x1)*(y2-y1) - (S_y(t) - y1)*(x2-x1) = 0   (on line)
/// F2 = S'_x(t)*(y2-y1) - S'_y(t)*(x2-x1) = 0                (parallel)
/// ```
pub(crate) fn eval_spline_tangent_to_line(
    sketch: &Sketch,
    spline: SplineHandle,
    line: LineHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, line);
    let dx = x2 - x1;
    let dy = y2 - y1;

    // Aux param t.
    let aux_idx = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t = sketch.params[aux_idx];

    let pos = eval_point(sp, t, sketch);
    let tang = eval_tangent(sp, t, sketch);

    let sx = pos.x;
    let sy = pos.y;

    // F1 = (sx - x1)*dy - (sy - y1)*dx
    let f1 = (sx - x1) * dy - (sy - y1) * dx;
    // F2 = tang.x * dy - tang.y * dx
    let f2 = tang.x * dy - tang.y * dx;

    // -- F1 Jacobian: line parameters --
    let mut jacobian = vec![
        (row, ix1, sy - y2),
        (row, iy1, x2 - sx),
        (row, ix2, -(sy - y1)),
        (row, iy2, sx - x1),
    ];

    // dF1/d(spline CPs): via pos jacobian
    for &(param_idx, dsx_dp, dsy_dp) in &pos.jacobian {
        let val = dsx_dp * dy - dsy_dp * dx;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row, param_idx, val);
        }
    }

    // dF1/dt = tang.x * dy - tang.y * dx = f2
    if f2.abs() > f64::EPSILON {
        jacobian.push((row, aux_idx, f2));
    }

    // -- F2 Jacobian: line parameters --
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, ix1, tang.y));
    }
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row + 1, iy1, -tang.x));
    }
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, ix2, -tang.y));
    }
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row + 1, iy2, tang.x));
    }

    // dF2/d(spline CPs): via tangent jacobian
    for &(param_idx, dtx_dp, dty_dp) in &tang.jacobian {
        let val = dtx_dp * dy - dty_dp * dx;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 1, param_idx, val);
        }
    }

    // dF2/dt = S''_x(t) * dy - S''_y(t) * dx
    let accel = super::spline_deriv::eval_second_deriv(sp, t, sketch);
    let df2_dt = accel.x * dy - accel.y * dx;
    if df2_dt.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 1, aux_idx, df2_dt);
    }

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian,
    }
}

// ── 2. SplineTangentToArc ────────────────────────────────────────────────

/// `SplineTangentToArc(Spline, Arc)` -- 2 equations, 1 aux param (t).
///
/// ```text
/// F1 = (S_x(t) - cx)^2 + (S_y(t) - cy)^2 - r^2 = 0   (on circle)
/// F2 = (S_x(t) - cx)*S'_x(t) + (S_y(t) - cy)*S'_y(t) = 0  (tangent perp radius)
/// ```
pub(crate) fn eval_spline_tangent_to_arc(
    sketch: &Sketch,
    spline: SplineHandle,
    arc: ArcHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let (cx, cy, r, icx, icy, ir) = ar(sketch, arc);

    let aux_idx = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t = sketch.params[aux_idx];

    let pos = eval_point(sp, t, sketch);
    let tang = eval_tangent(sp, t, sketch);

    let rx = pos.x - cx;
    let ry = pos.y - cy;

    // F1 = rx^2 + ry^2 - r^2
    let f1 = rx * rx + ry * ry - r * r;
    // F2 = rx * tang.x + ry * tang.y
    let f2 = rx * tang.x + ry * tang.y;

    let mut jacobian = Vec::new();

    // -- F1 Jacobian --
    jacobian.push((row, icx, -2.0 * rx));
    jacobian.push((row, icy, -2.0 * ry));
    jacobian.push((row, ir, -2.0 * r));

    for &(param_idx, dsx_dp, dsy_dp) in &pos.jacobian {
        let val = 2.0 * rx * dsx_dp + 2.0 * ry * dsy_dp;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row, param_idx, val);
        }
    }

    // dF1/dt = 2*f2
    let df1_dt = 2.0 * f2;
    if df1_dt.abs() > f64::EPSILON {
        jacobian.push((row, aux_idx, df1_dt));
    }

    // -- F2 Jacobian --
    if tang.x.abs() > f64::EPSILON {
        jacobian.push((row + 1, icx, -tang.x));
    }
    if tang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, icy, -tang.y));
    }

    // dF2/d(spline CPs): from pos and tang
    for &(param_idx, dsx_dp, dsy_dp) in &pos.jacobian {
        let val = dsx_dp * tang.x + dsy_dp * tang.y;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 1, param_idx, val);
        }
    }
    for &(param_idx, dtx_dp, dty_dp) in &tang.jacobian {
        let val = rx * dtx_dp + ry * dty_dp;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 1, param_idx, val);
        }
    }

    // dF2/dt = |S'|^2 + rx*S''_x + ry*S''_y
    let accel = super::spline_deriv::eval_second_deriv(sp, t, sketch);
    let df2_dt = tang.x * tang.x + tang.y * tang.y + rx * accel.x + ry * accel.y;
    if df2_dt.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 1, aux_idx, df2_dt);
    }

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian,
    }
}
