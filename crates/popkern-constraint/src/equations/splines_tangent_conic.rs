//! Spline-to-conic tangency constraints: ellipse, parabola, hyperbola.
//!
//! Split from `splines_tangent.rs` to stay within the 500-line file limit.

use super::EquationOutput;
use super::conics::{el, hyp, par};
use super::spline_deriv::eval_tangent;
use super::spline_eval::eval_point;
use super::splines::merge_or_push;
use crate::Sketch;
use crate::constraint::ConstraintHandle;
use crate::elements::{EllipseHandle, HyperbolaHandle, ParabolaHandle, SplineHandle};
use popkern_math::PortableMath;

// ── 1. SplineTangentToEllipse ────────────────────────────────────────────

/// `SplineTangentToEllipse(Spline, Ellipse)` -- 3 eqs, 2 aux (t_s, t_e).
///
/// ```text
/// F1 = S_x(t_s) - E_x(t_e) = 0
/// F2 = S_y(t_s) - E_y(t_e) = 0
/// F3 = S'_x(t_s)*E'_y(t_e) - S'_y(t_s)*E'_x(t_e) = 0
/// ```
pub(crate) fn eval_spline_tangent_to_ellipse(
    sketch: &Sketch,
    spline: SplineHandle,
    ellipse: EllipseHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = el(sketch, ellipse);

    let aux_start = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t_s = sketch.params[aux_start];
    let t_e = sketch.params[aux_start + 1];

    let cos_th = theta.pcos();
    let sin_th = theta.psin();
    let cos_te = t_e.pcos();
    let sin_te = t_e.psin();

    // Spline point and tangent
    let spos = eval_point(sp, t_s, sketch);
    let stang = eval_tangent(sp, t_s, sketch);

    // Ellipse point and tangent
    let ex = cx + a * cos_te * cos_th - b * sin_te * sin_th;
    let ey = cy + a * cos_te * sin_th + b * sin_te * cos_th;
    let etx = -a * sin_te * cos_th - b * cos_te * sin_th;
    let ety = -a * sin_te * sin_th + b * cos_te * cos_th;

    let f1 = spos.x - ex;
    let f2 = spos.y - ey;
    let f3 = stang.x * ety - stang.y * etx;

    let mut jacobian = Vec::new();

    // ── F1, F2 Jacobian: spline CPs ──
    for &(param_idx, dsx_dp, dsy_dp) in &spos.jacobian {
        if dsx_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dsx_dp));
        }
        if dsy_dp.abs() > f64::EPSILON {
            jacobian.push((row + 1, param_idx, dsy_dp));
        }
    }

    // ── F1, F2 Jacobian: ellipse params ──
    jacobian.push((row, icx, -1.0));
    jacobian.push((row + 1, icy, -1.0));

    let dex_da = cos_te * cos_th;
    let dey_da = cos_te * sin_th;
    let dex_db = -sin_te * sin_th;
    let dey_db = sin_te * cos_th;
    let dex_dtheta = -(a * cos_te * sin_th + b * sin_te * cos_th);
    let dey_dtheta = a * cos_te * cos_th - b * sin_te * sin_th;

    jacobian.push((row, ia, -dex_da));
    jacobian.push((row + 1, ia, -dey_da));
    jacobian.push((row, ib, -dex_db));
    jacobian.push((row + 1, ib, -dey_db));
    jacobian.push((row, itheta, -dex_dtheta));
    jacobian.push((row + 1, itheta, -dey_dtheta));

    // dF1/dt_s = S'_x(t_s), dF2/dt_s = S'_y(t_s)
    if stang.x.abs() > f64::EPSILON {
        jacobian.push((row, aux_start, stang.x));
    }
    if stang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start, stang.y));
    }

    // dF1/dt_e = -etx, dF2/dt_e = -ety
    if etx.abs() > f64::EPSILON {
        jacobian.push((row, aux_start + 1, -etx));
    }
    if ety.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start + 1, -ety));
    }

    // ── F3 Jacobian ──
    // F3 = stang.x * ety - stang.y * etx
    for &(param_idx, dtx_dp, dty_dp) in &stang.jacobian {
        let val = dtx_dp * ety - dty_dp * etx;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 2, param_idx, val);
        }
    }

    let detx_da = -sin_te * cos_th;
    let dety_da = -sin_te * sin_th;
    let detx_db = -cos_te * sin_th;
    let dety_db = cos_te * cos_th;
    let detx_dtheta = a * sin_te * sin_th - b * cos_te * cos_th;
    let dety_dtheta = -a * sin_te * cos_th - b * cos_te * sin_th;

    let df3_da = stang.x * dety_da - stang.y * detx_da;
    let df3_db = stang.x * dety_db - stang.y * detx_db;
    let df3_dtheta = stang.x * dety_dtheta - stang.y * detx_dtheta;

    if df3_da.abs() > f64::EPSILON {
        jacobian.push((row + 2, ia, df3_da));
    }
    if df3_db.abs() > f64::EPSILON {
        jacobian.push((row + 2, ib, df3_db));
    }
    if df3_dtheta.abs() > f64::EPSILON {
        jacobian.push((row + 2, itheta, df3_dtheta));
    }

    // dF3/dt_s = S''_x(t_s)*ety - S''_y(t_s)*etx
    let accel = super::spline_deriv::eval_second_deriv(sp, t_s, sketch);
    let df3_dts = accel.x * ety - accel.y * etx;
    if df3_dts.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start, df3_dts);
    }

    // dF3/dt_e
    let detx_dte = -a * cos_te * cos_th + b * sin_te * sin_th;
    let dety_dte = -a * cos_te * sin_th - b * sin_te * cos_th;
    let df3_dte = stang.x * dety_dte - stang.y * detx_dte;
    if df3_dte.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start + 1, df3_dte);
    }

    EquationOutput {
        residuals: vec![f1, f2, f3],
        jacobian,
    }
}

// ── 2. SplineTangentToParabola ───────────────────────────────────────────

/// `SplineTangentToParabola(Spline, Parabola)` -- 3 eqs, 2 aux (t_s, t_p).
///
/// Same structure as ellipse: position match + tangent cross-product = 0.
pub(crate) fn eval_spline_tangent_to_parabola(
    sketch: &Sketch,
    spline: SplineHandle,
    parabola: ParabolaHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let (vx, vy, f, theta, ivx, ivy, if_, itheta) = par(sketch, parabola);

    let aux_start = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t_s = sketch.params[aux_start];
    let t_p = sketch.params[aux_start + 1];

    let cos_th = theta.pcos();
    let sin_th = theta.psin();

    // Parabola point and tangent at t_p
    let u = t_p * t_p / (4.0 * f);
    let du_dt = t_p / (2.0 * f);

    let px = vx + u * cos_th - t_p * sin_th;
    let py = vy + u * sin_th + t_p * cos_th;
    let ptx = du_dt * cos_th - sin_th;
    let pty = du_dt * sin_th + cos_th;

    let spos = eval_point(sp, t_s, sketch);
    let stang = eval_tangent(sp, t_s, sketch);

    let f1 = spos.x - px;
    let f2 = spos.y - py;
    let f3 = stang.x * pty - stang.y * ptx;

    let mut jacobian = Vec::new();

    // ── F1, F2: spline CPs ──
    for &(param_idx, dsx_dp, dsy_dp) in &spos.jacobian {
        if dsx_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dsx_dp));
        }
        if dsy_dp.abs() > f64::EPSILON {
            jacobian.push((row + 1, param_idx, dsy_dp));
        }
    }

    // ── F1, F2: parabola params ──
    jacobian.push((row, ivx, -1.0));
    jacobian.push((row + 1, ivy, -1.0));

    let du_df = -t_p * t_p / (4.0 * f * f);
    let dpx_df = du_df * cos_th;
    let dpy_df = du_df * sin_th;
    let dpx_dtheta = -u * sin_th - t_p * cos_th;
    let dpy_dtheta = u * cos_th - t_p * sin_th;

    jacobian.push((row, if_, -dpx_df));
    jacobian.push((row + 1, if_, -dpy_df));
    jacobian.push((row, itheta, -dpx_dtheta));
    jacobian.push((row + 1, itheta, -dpy_dtheta));

    if stang.x.abs() > f64::EPSILON {
        jacobian.push((row, aux_start, stang.x));
    }
    if stang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start, stang.y));
    }

    if ptx.abs() > f64::EPSILON {
        jacobian.push((row, aux_start + 1, -ptx));
    }
    if pty.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start + 1, -pty));
    }

    // ── F3: spline CPs via tangent Jacobian ──
    for &(param_idx, dtx_dp, dty_dp) in &stang.jacobian {
        let val = dtx_dp * pty - dty_dp * ptx;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 2, param_idx, val);
        }
    }

    // ── F3: parabola params ──
    let ddu_dt_df = -t_p / (2.0 * f * f);
    let dptx_df = ddu_dt_df * cos_th;
    let dpty_df = ddu_dt_df * sin_th;
    let dptx_dtheta = -du_dt * sin_th - cos_th;
    let dpty_dtheta = du_dt * cos_th - sin_th;

    let df3_df = stang.x * dpty_df - stang.y * dptx_df;
    let df3_dtheta = stang.x * dpty_dtheta - stang.y * dptx_dtheta;

    if df3_df.abs() > f64::EPSILON {
        jacobian.push((row + 2, if_, df3_df));
    }
    if df3_dtheta.abs() > f64::EPSILON {
        jacobian.push((row + 2, itheta, df3_dtheta));
    }

    let accel = super::spline_deriv::eval_second_deriv(sp, t_s, sketch);
    let df3_dts = accel.x * pty - accel.y * ptx;
    if df3_dts.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start, df3_dts);
    }

    let d2u_dt2 = 1.0 / (2.0 * f);
    let dptx_dtp = d2u_dt2 * cos_th;
    let dpty_dtp = d2u_dt2 * sin_th;
    let df3_dtp = stang.x * dpty_dtp - stang.y * dptx_dtp;
    if df3_dtp.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start + 1, df3_dtp);
    }

    EquationOutput {
        residuals: vec![f1, f2, f3],
        jacobian,
    }
}

// ── 3. SplineTangentToHyperbola ──────────────────────────────────────────

/// `SplineTangentToHyperbola(Spline, Hyperbola)` -- 3 eqs, 2 aux (t_s, t_h).
///
/// Same structure: position match + tangent cross-product = 0.
pub(crate) fn eval_spline_tangent_to_hyperbola(
    sketch: &Sketch,
    spline: SplineHandle,
    hyperbola: HyperbolaHandle,
    handle: Option<ConstraintHandle>,
    row: usize,
) -> EquationOutput {
    let sp = &sketch.splines[spline];
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = hyp(sketch, hyperbola);

    let aux_start = handle
        .and_then(|h| sketch.aux_param_map.get(&h).copied())
        .unwrap_or(0);
    let t_s = sketch.params[aux_start];
    let t_h = sketch.params[aux_start + 1];

    let cos_th = theta.pcos();
    let sin_th = theta.psin();
    let cosh_t = t_h.pcosh();
    let sinh_t = t_h.psinh();

    // Hyperbola point and tangent at t_h
    let hx = cx + a * cosh_t * cos_th - b * sinh_t * sin_th;
    let hy = cy + a * cosh_t * sin_th + b * sinh_t * cos_th;
    let htx = a * sinh_t * cos_th - b * cosh_t * sin_th;
    let hty = a * sinh_t * sin_th + b * cosh_t * cos_th;

    let spos = eval_point(sp, t_s, sketch);
    let stang = eval_tangent(sp, t_s, sketch);

    let f1 = spos.x - hx;
    let f2 = spos.y - hy;
    let f3 = stang.x * hty - stang.y * htx;

    let mut jacobian = Vec::new();

    // ── F1, F2: spline CPs ──
    for &(param_idx, dsx_dp, dsy_dp) in &spos.jacobian {
        if dsx_dp.abs() > f64::EPSILON {
            jacobian.push((row, param_idx, dsx_dp));
        }
        if dsy_dp.abs() > f64::EPSILON {
            jacobian.push((row + 1, param_idx, dsy_dp));
        }
    }

    // ── F1, F2: hyperbola params ──
    jacobian.push((row, icx, -1.0));
    jacobian.push((row + 1, icy, -1.0));

    let dhx_da = cosh_t * cos_th;
    let dhy_da = cosh_t * sin_th;
    let dhx_db = -sinh_t * sin_th;
    let dhy_db = sinh_t * cos_th;
    let dhx_dtheta = -(a * cosh_t * sin_th + b * sinh_t * cos_th);
    let dhy_dtheta = a * cosh_t * cos_th - b * sinh_t * sin_th;

    jacobian.push((row, ia, -dhx_da));
    jacobian.push((row + 1, ia, -dhy_da));
    jacobian.push((row, ib, -dhx_db));
    jacobian.push((row + 1, ib, -dhy_db));
    jacobian.push((row, itheta, -dhx_dtheta));
    jacobian.push((row + 1, itheta, -dhy_dtheta));

    if stang.x.abs() > f64::EPSILON {
        jacobian.push((row, aux_start, stang.x));
    }
    if stang.y.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start, stang.y));
    }

    if htx.abs() > f64::EPSILON {
        jacobian.push((row, aux_start + 1, -htx));
    }
    if hty.abs() > f64::EPSILON {
        jacobian.push((row + 1, aux_start + 1, -hty));
    }

    // ── F3: spline CPs via tangent Jacobian ──
    for &(param_idx, dtx_dp, dty_dp) in &stang.jacobian {
        let val = dtx_dp * hty - dty_dp * htx;
        if val.abs() > f64::EPSILON {
            merge_or_push(&mut jacobian, row + 2, param_idx, val);
        }
    }

    // ── F3: hyperbola params ──
    let dhtx_da = sinh_t * cos_th;
    let dhty_da = sinh_t * sin_th;
    let dhtx_db = -cosh_t * sin_th;
    let dhty_db = cosh_t * cos_th;
    let dhtx_dtheta = -(a * sinh_t * sin_th + b * cosh_t * cos_th);
    let dhty_dtheta = a * sinh_t * cos_th - b * cosh_t * sin_th;

    let df3_da = stang.x * dhty_da - stang.y * dhtx_da;
    let df3_db = stang.x * dhty_db - stang.y * dhtx_db;
    let df3_dtheta = stang.x * dhty_dtheta - stang.y * dhtx_dtheta;

    if df3_da.abs() > f64::EPSILON {
        jacobian.push((row + 2, ia, df3_da));
    }
    if df3_db.abs() > f64::EPSILON {
        jacobian.push((row + 2, ib, df3_db));
    }
    if df3_dtheta.abs() > f64::EPSILON {
        jacobian.push((row + 2, itheta, df3_dtheta));
    }

    let accel = super::spline_deriv::eval_second_deriv(sp, t_s, sketch);
    let df3_dts = accel.x * hty - accel.y * htx;
    if df3_dts.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start, df3_dts);
    }

    let dhtx_dth = a * cosh_t * cos_th - b * sinh_t * sin_th;
    let dhty_dth = a * cosh_t * sin_th + b * sinh_t * cos_th;
    let df3_dth = stang.x * dhty_dth - stang.y * dhtx_dth;
    if df3_dth.abs() > f64::EPSILON {
        merge_or_push(&mut jacobian, row + 2, aux_start + 1, df3_dth);
    }

    EquationOutput {
        residuals: vec![f1, f2, f3],
        jacobian,
    }
}
