//! Line-to-conic tangency constraint evaluators.
//!
//! Uses the dual conic formulation to avoid normalisation of the line
//! normal. Split from `conics.rs` to stay within the 500-line limit.

use super::conics::{el, hyp, par};
use super::{EquationOutput, ln};
use crate::Sketch;
use crate::elements::{EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle};

/// TangentLineEllipse(L, E): dual conic tangency condition.
///
/// In the ellipse local frame, the line through two endpoints is
/// `l*u + m*v + n = 0`. The dual conic tangency for x^2/a^2 + y^2/b^2 = 1
/// is `a^2*l^2 + b^2*m^2 = n^2`, i.e., residual = a^2*l^2 + b^2*m^2 - n^2.
pub(crate) fn eval_tangent_line_ellipse(
    sketch: &Sketch,
    line: LineHandle,
    ellipse: EllipseHandle,
    row: usize,
) -> EquationOutput {
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, line);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = el(sketch, ellipse);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    // Line direction in global frame
    let gdx = x2 - x1;
    let gdy = y2 - y1;

    // Line direction in local frame
    let dx_l = gdx * cos_t + gdy * sin_t;
    let dy_l = -gdx * sin_t + gdy * cos_t;

    // First point in local frame
    let ux1 = (x1 - cx) * cos_t + (y1 - cy) * sin_t;
    let vx1 = -(x1 - cx) * sin_t + (y1 - cy) * cos_t;

    // Line equation in local: l*u + m*v + n = 0
    let l = dy_l;
    let m = -dx_l;
    let n = dx_l * vx1 - dy_l * ux1;

    let residual = a * a * l * l + b * b * m * m - n * n;

    let df_dl = 2.0 * a * a * l;
    let df_dm = 2.0 * b * b * m;
    let df_dn = -2.0 * n;

    let (jac_line, jac_center, jac_theta) =
        line_lmn_jacobian(cos_t, sin_t, dx_l, dy_l, ux1, vx1, df_dl, df_dm, df_dn);

    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ix1, jac_line[0]),
            (row, iy1, jac_line[1]),
            (row, ix2, jac_line[2]),
            (row, iy2, jac_line[3]),
            (row, icx, jac_center[0]),
            (row, icy, jac_center[1]),
            (row, ia, 2.0 * a * l * l),
            (row, ib, 2.0 * b * m * m),
            (row, itheta, jac_theta),
        ],
    }
}

/// TangentLineParabola(L, Par): dual conic tangency condition.
///
/// For a parabola v^2 = 4*f*u in local frame, a line l*u + m*v + n = 0
/// is tangent when f*m^2 - l*n = 0.
pub(crate) fn eval_tangent_line_parabola(
    sketch: &Sketch,
    line: LineHandle,
    parabola: ParabolaHandle,
    row: usize,
) -> EquationOutput {
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, line);
    let (vx, vy, f, theta, ivx, ivy, if_, itheta) = par(sketch, parabola);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let gdx = x2 - x1;
    let gdy = y2 - y1;

    let dx_l = gdx * cos_t + gdy * sin_t;
    let dy_l = -gdx * sin_t + gdy * cos_t;

    let ux1 = (x1 - vx) * cos_t + (y1 - vy) * sin_t;
    let vx1 = -(x1 - vx) * sin_t + (y1 - vy) * cos_t;

    let l = dy_l;
    let m = -dx_l;
    let n = dx_l * vx1 - dy_l * ux1;

    let residual = f * m * m - l * n;

    let df_dl = -n;
    let df_dm = 2.0 * f * m;
    let df_dn = -l;

    let (jac_line, jac_center, jac_theta) =
        line_lmn_jacobian(cos_t, sin_t, dx_l, dy_l, ux1, vx1, df_dl, df_dm, df_dn);

    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ix1, jac_line[0]),
            (row, iy1, jac_line[1]),
            (row, ix2, jac_line[2]),
            (row, iy2, jac_line[3]),
            (row, ivx, jac_center[0]),
            (row, ivy, jac_center[1]),
            (row, if_, m * m),
            (row, itheta, jac_theta),
        ],
    }
}

/// TangentLineHyperbola(L, H): dual conic tangency condition.
///
/// For a hyperbola u^2/a^2 - v^2/b^2 = 1, a line l*u + m*v + n = 0
/// is tangent when a^2*l^2 - b^2*m^2 = n^2.
pub(crate) fn eval_tangent_line_hyperbola(
    sketch: &Sketch,
    line: LineHandle,
    hyperbola: HyperbolaHandle,
    row: usize,
) -> EquationOutput {
    let (x1, y1, ix1, iy1, x2, y2, ix2, iy2) = ln(sketch, line);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = hyp(sketch, hyperbola);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let gdx = x2 - x1;
    let gdy = y2 - y1;

    let dx_l = gdx * cos_t + gdy * sin_t;
    let dy_l = -gdx * sin_t + gdy * cos_t;

    let ux1 = (x1 - cx) * cos_t + (y1 - cy) * sin_t;
    let vx1 = -(x1 - cx) * sin_t + (y1 - cy) * cos_t;

    let l = dy_l;
    let m = -dx_l;
    let n = dx_l * vx1 - dy_l * ux1;

    let residual = a * a * l * l - b * b * m * m - n * n;

    let df_dl = 2.0 * a * a * l;
    let df_dm = -2.0 * b * b * m;
    let df_dn = -2.0 * n;

    let (jac_line, jac_center, jac_theta) =
        line_lmn_jacobian(cos_t, sin_t, dx_l, dy_l, ux1, vx1, df_dl, df_dm, df_dn);

    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ix1, jac_line[0]),
            (row, iy1, jac_line[1]),
            (row, ix2, jac_line[2]),
            (row, iy2, jac_line[3]),
            (row, icx, jac_center[0]),
            (row, icy, jac_center[1]),
            (row, ia, 2.0 * a * l * l),
            (row, ib, -2.0 * b * m * m),
            (row, itheta, jac_theta),
        ],
    }
}

// ── Shared Jacobian helper ────────────────────────────────────────────────

/// Compute the chain-rule Jacobian of F(l, m, n) through the line-to-local
/// frame transformation, for line endpoint, conic centre, and rotation
/// parameters.
///
/// Returns `(line_jac, center_jac, theta_jac)` where:
/// - `line_jac = [dF/dx1, dF/dy1, dF/dx2, dF/dy2]`
/// - `center_jac = [dF/dcx, dF/dcy]`
/// - `theta_jac = dF/dtheta`
///
/// The caller supplies `dF/dl`, `dF/dm`, `dF/dn` and the local-frame
/// intermediates.
#[allow(clippy::too_many_arguments)]
fn line_lmn_jacobian(
    cos_t: f64,
    sin_t: f64,
    dx_l: f64,
    dy_l: f64,
    ux1: f64,
    vx1: f64,
    df_dl: f64,
    df_dm: f64,
    df_dn: f64,
) -> ([f64; 4], [f64; 2], f64) {
    let df = |dl: f64, dm: f64, dn: f64| -> f64 { df_dl * dl + df_dm * dm + df_dn * dn };

    // n = dx_l * vx1 - dy_l * ux1
    // dn/dvar = d(dx_l)/dvar * vx1 + dx_l * d(vx1)/dvar
    //         - d(dy_l)/dvar * ux1 - dy_l * d(ux1)/dvar

    // Variable: x1 (gdx changes by -1)
    // d(dx_l)/dx1 = -cos, d(dy_l)/dx1 = sin, d(ux1)/dx1 = cos, d(vx1)/dx1 = -sin
    let dn_dx1 = (-cos_t) * vx1 + dx_l * (-sin_t) - sin_t * ux1 - dy_l * cos_t;
    let df_dx1 = df(sin_t, cos_t, dn_dx1);

    // Variable: y1 (gdy changes by -1)
    // d(dx_l)/dy1 = -sin, d(dy_l)/dy1 = -cos, d(ux1)/dy1 = sin, d(vx1)/dy1 = cos
    let dn_dy1 = (-sin_t) * vx1 + dx_l * cos_t - (-cos_t) * ux1 - dy_l * sin_t;
    let df_dy1 = df(-cos_t, sin_t, dn_dy1);

    // Variable: x2 (gdx changes by +1, ux1/vx1 unaffected)
    // d(dx_l)/dx2 = cos, d(dy_l)/dx2 = -sin
    let dn_dx2 = cos_t * vx1 - (-sin_t) * ux1;
    let df_dx2 = df(-sin_t, -cos_t, dn_dx2);

    // Variable: y2 (gdy changes by +1)
    // d(dx_l)/dy2 = sin, d(dy_l)/dy2 = cos
    let dn_dy2 = sin_t * vx1 - cos_t * ux1;
    let df_dy2 = df(cos_t, -sin_t, dn_dy2);

    // Variable: cx (only ux1, vx1 change; l, m, dx_l, dy_l unaffected)
    // d(ux1)/dcx = -cos, d(vx1)/dcx = sin
    let dn_dcx = dx_l * sin_t - dy_l * (-cos_t);
    let df_dcx = df(0.0, 0.0, dn_dcx);

    // Variable: cy
    // d(ux1)/dcy = -sin, d(vx1)/dcy = -cos
    let dn_dcy = dx_l * (-cos_t) - dy_l * (-sin_t);
    let df_dcy = df(0.0, 0.0, dn_dcy);

    // Variable: theta
    // d(dx_l)/dtheta = dy_l, d(dy_l)/dtheta = -dx_l
    // d(ux1)/dtheta = vx1, d(vx1)/dtheta = -ux1
    let dn_dtheta = dy_l * vx1 + dx_l * (-ux1) - (-dx_l) * ux1 - dy_l * vx1;
    let df_dtheta = df(-dx_l, -dy_l, dn_dtheta);

    (
        [df_dx1, df_dy1, df_dx2, df_dy2],
        [df_dcx, df_dcy],
        df_dtheta,
    )
}
