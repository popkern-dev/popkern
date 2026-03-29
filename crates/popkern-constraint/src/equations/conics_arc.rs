//! Arc-to-conic tangency constraint evaluators.
//!
//! Each constraint uses one auxiliary parameter `t` (the contact parameter
//! on the conic) and produces two equations:
//!
//! 1. Distance from arc centre to conic point at `t` equals arc radius.
//! 2. Conic tangent at `t` is perpendicular to the centre-to-contact
//!    vector (dot product = 0).
//!
//! Split from `conics_tangent.rs` to stay within the 500-line limit.

use super::conics::{el, hyp, par};
use super::{ar, EquationOutput};
use crate::elements::{ArcHandle, EllipseHandle, HyperbolaHandle, ParabolaHandle};
use crate::Sketch;

/// TangentArcEllipse: arc tangent to ellipse via auxiliary parameter `t`.
///
/// F1 = (acx - Ex(t))^2 + (acy - Ey(t))^2 - r^2
/// F2 = (acx - Ex(t)) * Tx(t) + (acy - Ey(t)) * Ty(t)
///
/// Ellipse parameterisation:
///   Ex(t) = cx + a*cos(t)*cos(th) - b*sin(t)*sin(th)
///   Ey(t) = cy + a*cos(t)*sin(th) + b*sin(t)*cos(th)
///   Tx(t) = -a*sin(t)*cos(th) - b*cos(t)*sin(th)
///   Ty(t) = -a*sin(t)*sin(th) + b*cos(t)*cos(th)
pub(crate) fn eval_tangent_arc_ellipse(
    sketch: &Sketch,
    arc_h: ArcHandle,
    ellipse_h: EllipseHandle,
    aux_idx: usize,
    row: usize,
) -> EquationOutput {
    let (acx, acy, r, iacx, iacy, ir) = ar(sketch, arc_h);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = el(sketch, ellipse_h);
    let t = sketch.params[aux_idx];

    let cos_th = theta.cos();
    let sin_th = theta.sin();
    let cos_t = t.cos();
    let sin_t = t.sin();

    // Conic point
    let ex = cx + a * cos_t * cos_th - b * sin_t * sin_th;
    let ey = cy + a * cos_t * sin_th + b * sin_t * cos_th;

    // Conic tangent (dE/dt)
    let tx = -a * sin_t * cos_th - b * cos_t * sin_th;
    let ty = -a * sin_t * sin_th + b * cos_t * cos_th;

    let dx = acx - ex;
    let dy = acy - ey;

    // F1: distance squared minus r squared
    let f1 = dx * dx + dy * dy - r * r;
    // F2: tangent perpendicular to centre-to-contact (dot product)
    let f2 = dx * tx + dy * ty;

    // ── Jacobian ─────────────────────────────────────────────────────

    // d^2E/dt^2 (second derivative of conic point w.r.t. t):
    let dtx_dt = -a * cos_t * cos_th + b * sin_t * sin_th;
    let dty_dt = -a * cos_t * sin_th - b * sin_t * cos_th;

    // Partial derivatives of E w.r.t. conic params:
    let dex_da = cos_t * cos_th;
    let dey_da = cos_t * sin_th;
    let dex_db = -sin_t * sin_th;
    let dey_db = sin_t * cos_th;
    let dex_dtheta = -(a * cos_t * sin_th + b * sin_t * cos_th);
    let dey_dtheta = a * cos_t * cos_th - b * sin_t * sin_th;

    // Partial derivatives of T w.r.t. conic params:
    let dtx_da = -sin_t * cos_th;
    let dty_da = -sin_t * sin_th;
    let dtx_db = -cos_t * sin_th;
    let dty_db = cos_t * cos_th;
    let dtx_dtheta = a * sin_t * sin_th - b * cos_t * cos_th;
    let dty_dtheta = -a * sin_t * cos_th - b * cos_t * sin_th;

    // --- F1 derivatives ---
    let df1_acx = 2.0 * dx;
    let df1_acy = 2.0 * dy;
    let df1_r = -2.0 * r;
    let df1_cx = -2.0 * dx;
    let df1_cy = -2.0 * dy;
    let df1_a = -2.0 * (dx * dex_da + dy * dey_da);
    let df1_b = -2.0 * (dx * dex_db + dy * dey_db);
    let df1_theta = -2.0 * (dx * dex_dtheta + dy * dey_dtheta);
    let df1_dt = -2.0 * (dx * tx + dy * ty);

    // --- F2 derivatives ---
    // F2 = dx*tx + dy*ty (dot product)
    // dF2/d(acx) = tx,  dF2/d(acy) = ty
    // dF2/d(cx) = -tx,  dF2/d(cy) = -ty
    // dF2/d(r) = 0
    let df2_acx = tx;
    let df2_acy = ty;
    let df2_r = 0.0;
    let df2_cx = -tx;
    let df2_cy = -ty;
    // dF2/d(a) = -dex_da*tx + dx*dtx_da - dey_da*ty + dy*dty_da
    //          = ... rearranged for clarity:
    let df2_a = -dex_da * tx - dey_da * ty + dx * dtx_da + dy * dty_da;
    let df2_b = -dex_db * tx - dey_db * ty + dx * dtx_db + dy * dty_db;
    let df2_theta = -dex_dtheta * tx - dey_dtheta * ty + dx * dtx_dtheta + dy * dty_dtheta;
    // dF2/dt = d(dx)/dt*tx + dx*dtx_dt + d(dy)/dt*ty + dy*dty_dt
    //        = -tx*tx + dx*dtx_dt - ty*ty + dy*dty_dt
    let df2_dt = -(tx * tx + ty * ty) + dx * dtx_dt + dy * dty_dt;

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, iacx, df1_acx),
            (row, iacy, df1_acy),
            (row, ir, df1_r),
            (row, icx, df1_cx),
            (row, icy, df1_cy),
            (row, ia, df1_a),
            (row, ib, df1_b),
            (row, itheta, df1_theta),
            (row, aux_idx, df1_dt),
            (row + 1, iacx, df2_acx),
            (row + 1, iacy, df2_acy),
            (row + 1, ir, df2_r),
            (row + 1, icx, df2_cx),
            (row + 1, icy, df2_cy),
            (row + 1, ia, df2_a),
            (row + 1, ib, df2_b),
            (row + 1, itheta, df2_theta),
            (row + 1, aux_idx, df2_dt),
        ],
    }
}

/// TangentArcParabola: arc tangent to parabola via auxiliary parameter `t`.
///
/// Parabola parameterisation (t is local v-coordinate):
///   u(t) = t^2 / (4f)
///   Ex(t) = vx + u(t)*cos(th) - t*sin(th)
///   Ey(t) = vy + u(t)*sin(th) + t*cos(th)
///   Tx(t) = (t/(2f))*cos(th) - sin(th)
///   Ty(t) = (t/(2f))*sin(th) + cos(th)
pub(crate) fn eval_tangent_arc_parabola(
    sketch: &Sketch,
    arc_h: ArcHandle,
    parabola_h: ParabolaHandle,
    aux_idx: usize,
    row: usize,
) -> EquationOutput {
    let (acx, acy, r, iacx, iacy, ir) = ar(sketch, arc_h);
    let (vx, vy, f, theta, ivx, ivy, if_, itheta) = par(sketch, parabola_h);
    let t = sketch.params[aux_idx];

    let cos_th = theta.cos();
    let sin_th = theta.sin();

    let u = t * t / (4.0 * f);
    let du_dt = t / (2.0 * f);
    let du_df = -t * t / (4.0 * f * f);

    // Conic point
    let ex = vx + u * cos_th - t * sin_th;
    let ey = vy + u * sin_th + t * cos_th;

    // Conic tangent (dE/dt)
    let tx = du_dt * cos_th - sin_th;
    let ty = du_dt * sin_th + cos_th;

    let dx = acx - ex;
    let dy = acy - ey;

    let f1 = dx * dx + dy * dy - r * r;
    let f2 = dx * tx + dy * ty;

    // ── F1 Jacobian ──────────────────────────────────────────────────

    let dex_df = du_df * cos_th;
    let dey_df = du_df * sin_th;
    let dex_dtheta = -u * sin_th - t * cos_th;
    let dey_dtheta = u * cos_th - t * sin_th;

    let df1_acx = 2.0 * dx;
    let df1_acy = 2.0 * dy;
    let df1_r = -2.0 * r;
    let df1_vx = -2.0 * dx;
    let df1_vy = -2.0 * dy;
    let df1_f = -2.0 * (dx * dex_df + dy * dey_df);
    let df1_theta = -2.0 * (dx * dex_dtheta + dy * dey_dtheta);
    let df1_dt = -2.0 * (dx * tx + dy * ty);

    // ── F2 Jacobian ──────────────────────────────────────────────────

    let ddu_dt_df = -t / (2.0 * f * f);
    let dtx_df = ddu_dt_df * cos_th;
    let dty_df = ddu_dt_df * sin_th;
    let dtx_dtheta = -du_dt * sin_th - cos_th;
    let dty_dtheta = du_dt * cos_th - sin_th;
    let d2u_dt2 = 1.0 / (2.0 * f);
    let dtx_dt = d2u_dt2 * cos_th;
    let dty_dt = d2u_dt2 * sin_th;

    // F2 = dx*tx + dy*ty
    let df2_acx = tx;
    let df2_acy = ty;
    let df2_r = 0.0;
    let df2_vx = -tx;
    let df2_vy = -ty;
    let df2_f = -dex_df * tx - dey_df * ty + dx * dtx_df + dy * dty_df;
    let df2_theta = -dex_dtheta * tx - dey_dtheta * ty + dx * dtx_dtheta + dy * dty_dtheta;
    // dF2/dt = -tx*tx - ty*ty + dx*dtx_dt + dy*dty_dt
    let df2_dt = -(tx * tx + ty * ty) + dx * dtx_dt + dy * dty_dt;

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, iacx, df1_acx),
            (row, iacy, df1_acy),
            (row, ir, df1_r),
            (row, ivx, df1_vx),
            (row, ivy, df1_vy),
            (row, if_, df1_f),
            (row, itheta, df1_theta),
            (row, aux_idx, df1_dt),
            (row + 1, iacx, df2_acx),
            (row + 1, iacy, df2_acy),
            (row + 1, ir, df2_r),
            (row + 1, ivx, df2_vx),
            (row + 1, ivy, df2_vy),
            (row + 1, if_, df2_f),
            (row + 1, itheta, df2_theta),
            (row + 1, aux_idx, df2_dt),
        ],
    }
}

/// TangentArcHyperbola: arc tangent to hyperbola via auxiliary parameter `t`.
///
/// Hyperbola parameterisation (using sinh/cosh):
///   Ex(t) = cx + a*cosh(t)*cos(th) - b*sinh(t)*sin(th)
///   Ey(t) = cy + a*cosh(t)*sin(th) + b*sinh(t)*cos(th)
///   Tx(t) = a*sinh(t)*cos(th) - b*cosh(t)*sin(th)
///   Ty(t) = a*sinh(t)*sin(th) + b*cosh(t)*cos(th)
pub(crate) fn eval_tangent_arc_hyperbola(
    sketch: &Sketch,
    arc_h: ArcHandle,
    hyperbola_h: HyperbolaHandle,
    aux_idx: usize,
    row: usize,
) -> EquationOutput {
    let (acx, acy, r, iacx, iacy, ir) = ar(sketch, arc_h);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = hyp(sketch, hyperbola_h);
    let t = sketch.params[aux_idx];

    let cos_th = theta.cos();
    let sin_th = theta.sin();
    let cosh_t = t.cosh();
    let sinh_t = t.sinh();

    // Conic point
    let ex = cx + a * cosh_t * cos_th - b * sinh_t * sin_th;
    let ey = cy + a * cosh_t * sin_th + b * sinh_t * cos_th;

    // Conic tangent (dE/dt)
    let tx = a * sinh_t * cos_th - b * cosh_t * sin_th;
    let ty = a * sinh_t * sin_th + b * cosh_t * cos_th;

    let dx = acx - ex;
    let dy = acy - ey;

    let f1 = dx * dx + dy * dy - r * r;
    let f2 = dx * tx + dy * ty;

    // ── F1 Jacobian ──────────────────────────────────────────────────

    let dex_da = cosh_t * cos_th;
    let dey_da = cosh_t * sin_th;
    let dex_db = -sinh_t * sin_th;
    let dey_db = sinh_t * cos_th;
    let dex_dtheta = -(a * cosh_t * sin_th + b * sinh_t * cos_th);
    let dey_dtheta = a * cosh_t * cos_th - b * sinh_t * sin_th;

    let df1_acx = 2.0 * dx;
    let df1_acy = 2.0 * dy;
    let df1_r = -2.0 * r;
    let df1_cx = -2.0 * dx;
    let df1_cy = -2.0 * dy;
    let df1_a = -2.0 * (dx * dex_da + dy * dey_da);
    let df1_b = -2.0 * (dx * dex_db + dy * dey_db);
    let df1_theta = -2.0 * (dx * dex_dtheta + dy * dey_dtheta);
    let df1_dt = -2.0 * (dx * tx + dy * ty);

    // ── F2 Jacobian ──────────────────────────────────────────────────

    let dtx_dt = a * cosh_t * cos_th - b * sinh_t * sin_th;
    let dty_dt = a * cosh_t * sin_th + b * sinh_t * cos_th;

    let dtx_da = sinh_t * cos_th;
    let dty_da = sinh_t * sin_th;
    let dtx_db = -cosh_t * sin_th;
    let dty_db = cosh_t * cos_th;
    let dtx_dtheta = -(a * sinh_t * sin_th + b * cosh_t * cos_th);
    let dty_dtheta = a * sinh_t * cos_th - b * cosh_t * sin_th;

    // F2 = dx*tx + dy*ty
    let df2_acx = tx;
    let df2_acy = ty;
    let df2_r = 0.0;
    let df2_cx = -tx;
    let df2_cy = -ty;
    let df2_a = -dex_da * tx - dey_da * ty + dx * dtx_da + dy * dty_da;
    let df2_b = -dex_db * tx - dey_db * ty + dx * dtx_db + dy * dty_db;
    let df2_theta = -dex_dtheta * tx - dey_dtheta * ty + dx * dtx_dtheta + dy * dty_dtheta;
    let df2_dt = -(tx * tx + ty * ty) + dx * dtx_dt + dy * dty_dt;

    EquationOutput {
        residuals: vec![f1, f2],
        jacobian: vec![
            (row, iacx, df1_acx),
            (row, iacy, df1_acy),
            (row, ir, df1_r),
            (row, icx, df1_cx),
            (row, icy, df1_cy),
            (row, ia, df1_a),
            (row, ib, df1_b),
            (row, itheta, df1_theta),
            (row, aux_idx, df1_dt),
            (row + 1, iacx, df2_acx),
            (row + 1, iacy, df2_acy),
            (row + 1, ir, df2_r),
            (row + 1, icx, df2_cx),
            (row + 1, icy, df2_cy),
            (row + 1, ia, df2_a),
            (row + 1, ib, df2_b),
            (row + 1, itheta, df2_theta),
            (row + 1, aux_idx, df2_dt),
        ],
    }
}
