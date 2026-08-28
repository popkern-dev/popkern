//! Residual and Jacobian evaluators for conic constraints (ellipse,
//! parabola, hyperbola): point-on-conic, fixed, concentric, dimensional,
//! and eccentricity.
//!
//! Line-to-conic tangency is in the sibling [`conics_tangent`] module.

use super::{EquationOutput, ar, pt};
use crate::Sketch;
use crate::elements::{ArcHandle, EllipseHandle, HyperbolaHandle, ParabolaHandle, PointHandle};

// ── Conic helper functions ────────────────────────────────────────────────

/// Read ellipse parameters and their indices from the sketch.
///
/// Returns `(cx, cy, a, b, theta, icx, icy, ia, ib, itheta)`.
pub(crate) fn el(
    sketch: &Sketch,
    h: EllipseHandle,
) -> (f64, f64, f64, f64, f64, usize, usize, usize, usize, usize) {
    let e = &sketch.ellipses[h];
    let [icx, icy, ia, ib, itheta] = e.param_indices();
    (
        sketch.params[icx],
        sketch.params[icy],
        sketch.params[ia],
        sketch.params[ib],
        sketch.params[itheta],
        icx,
        icy,
        ia,
        ib,
        itheta,
    )
}

/// Read parabola parameters and their indices from the sketch.
///
/// Returns `(vx, vy, f, theta, ivx, ivy, if_, itheta)`.
pub(crate) fn par(
    sketch: &Sketch,
    h: ParabolaHandle,
) -> (f64, f64, f64, f64, usize, usize, usize, usize) {
    let p = &sketch.parabolas[h];
    let [ivx, ivy, if_, itheta] = p.param_indices();
    (
        sketch.params[ivx],
        sketch.params[ivy],
        sketch.params[if_],
        sketch.params[itheta],
        ivx,
        ivy,
        if_,
        itheta,
    )
}

/// Read hyperbola parameters and their indices from the sketch.
///
/// Returns `(cx, cy, a, b, theta, icx, icy, ia, ib, itheta)`.
pub(crate) fn hyp(
    sketch: &Sketch,
    h: HyperbolaHandle,
) -> (f64, f64, f64, f64, f64, usize, usize, usize, usize, usize) {
    let hp = &sketch.hyperbolas[h];
    let [icx, icy, ia, ib, itheta] = hp.param_indices();
    (
        sketch.params[icx],
        sketch.params[icy],
        sketch.params[ia],
        sketch.params[ib],
        sketch.params[itheta],
        icx,
        icy,
        ia,
        ib,
        itheta,
    )
}

/// Transform a point into the local frame of a conic with centre
/// `(cx, cy)` and rotation `theta`.
///
/// Returns `(u, v)` where:
/// - `u = (px - cx) * cos(theta) + (py - cy) * sin(theta)`
/// - `v = -(px - cx) * sin(theta) + (py - cy) * cos(theta)`
fn to_local(px: f64, py: f64, cx: f64, cy: f64, theta: f64) -> (f64, f64) {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dx = px - cx;
    let dy = py - cy;
    (dx * cos_t + dy * sin_t, -dx * sin_t + dy * cos_t)
}

// ── Point-on-conic constraints ────────────────────────────────────────────

/// PointOnEllipse(P, E): F = (u/a)^2 + (v/b)^2 - 1
pub(crate) fn eval_point_on_ellipse(
    sketch: &Sketch,
    p: PointHandle,
    e: EllipseHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = el(sketch, e);
    let (u, v) = to_local(px, py, cx, cy, theta);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let ua2 = 2.0 * u / (a * a);
    let vb2 = 2.0 * v / (b * b);

    let residual = (u / a).powi(2) + (v / b).powi(2) - 1.0;

    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ipx, ua2 * cos_t - vb2 * sin_t),
            (row, ipy, ua2 * sin_t + vb2 * cos_t),
            (row, icx, -ua2 * cos_t + vb2 * sin_t),
            (row, icy, -ua2 * sin_t - vb2 * cos_t),
            (row, ia, -2.0 * u * u / (a * a * a)),
            (row, ib, -2.0 * v * v / (b * b * b)),
            (row, itheta, 2.0 * u * v * (1.0 / (a * a) - 1.0 / (b * b))),
        ],
    }
}

/// PointOnParabola(P, Par): F = v^2 - 4*f*u
pub(crate) fn eval_point_on_parabola(
    sketch: &Sketch,
    p: PointHandle,
    pa: ParabolaHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (vx, vy, f, theta, ivx, ivy, if_, itheta) = par(sketch, pa);
    let (u, v) = to_local(px, py, vx, vy, theta);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let residual = v * v - 4.0 * f * u;

    // dF/du = -4f, dF/dv = 2v
    // du/dpx = cos, dv/dpx = -sin
    // du/dpy = sin, dv/dpy = cos
    // du/dvx = -cos, dv/dvx = sin
    // du/dvy = -sin, dv/dvy = -cos
    // du/dtheta = v, dv/dtheta = -u

    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ipx, -4.0 * f * cos_t - 2.0 * v * sin_t),
            (row, ipy, -4.0 * f * sin_t + 2.0 * v * cos_t),
            (row, ivx, 4.0 * f * cos_t + 2.0 * v * sin_t),
            (row, ivy, 4.0 * f * sin_t - 2.0 * v * cos_t),
            (row, if_, -4.0 * u),
            (row, itheta, -2.0 * v * (2.0 * f + u)),
        ],
    }
}

/// PointOnHyperbola(P, H): F = (u/a)^2 - (v/b)^2 - 1
pub(crate) fn eval_point_on_hyperbola(
    sketch: &Sketch,
    p: PointHandle,
    h: HyperbolaHandle,
    row: usize,
) -> EquationOutput {
    let (px, py, ipx, ipy) = pt(sketch, p);
    let (cx, cy, a, b, theta, icx, icy, ia, ib, itheta) = hyp(sketch, h);
    let (u, v) = to_local(px, py, cx, cy, theta);
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let ua2 = 2.0 * u / (a * a);
    let vb2 = 2.0 * v / (b * b);

    let residual = (u / a).powi(2) - (v / b).powi(2) - 1.0;

    // dF/du = 2u/a^2, dF/dv = -2v/b^2
    EquationOutput {
        residuals: vec![residual],
        jacobian: vec![
            (row, ipx, ua2 * cos_t + vb2 * sin_t),
            (row, ipy, ua2 * sin_t - vb2 * cos_t),
            (row, icx, -ua2 * cos_t - vb2 * sin_t),
            (row, icy, -ua2 * sin_t + vb2 * cos_t),
            (row, ia, -2.0 * u * u / (a * a * a)),
            (row, ib, 2.0 * v * v / (b * b * b)),
            (row, itheta, 2.0 * u * v * (1.0 / (a * a) + 1.0 / (b * b))),
        ],
    }
}

// ── Fixed conic constraints ───────────────────────────────────────────────

/// FixedEllipse(E): 5 equations, each param minus its captured value.
pub(crate) fn eval_fixed_ellipse(
    sketch: &Sketch,
    e: EllipseHandle,
    fixed: &[f64],
    row: usize,
) -> EquationOutput {
    let el = &sketch.ellipses[e];
    let indices = el.param_indices();
    let mut residuals = Vec::with_capacity(5);
    let mut jacobian = Vec::with_capacity(5);
    for (i, &idx) in indices.iter().enumerate() {
        residuals.push(sketch.params[idx] - fixed[i]);
        jacobian.push((row + i, idx, 1.0));
    }
    EquationOutput {
        residuals,
        jacobian,
    }
}

/// FixedParabola(Par): 4 equations, each param minus its captured value.
pub(crate) fn eval_fixed_parabola(
    sketch: &Sketch,
    p: ParabolaHandle,
    fixed: &[f64],
    row: usize,
) -> EquationOutput {
    let pa = &sketch.parabolas[p];
    let indices = pa.param_indices();
    let mut residuals = Vec::with_capacity(4);
    let mut jacobian = Vec::with_capacity(4);
    for (i, &idx) in indices.iter().enumerate() {
        residuals.push(sketch.params[idx] - fixed[i]);
        jacobian.push((row + i, idx, 1.0));
    }
    EquationOutput {
        residuals,
        jacobian,
    }
}

/// FixedHyperbola(H): 5 equations, each param minus its captured value.
pub(crate) fn eval_fixed_hyperbola(
    sketch: &Sketch,
    h: HyperbolaHandle,
    fixed: &[f64],
    row: usize,
) -> EquationOutput {
    let hp = &sketch.hyperbolas[h];
    let indices = hp.param_indices();
    let mut residuals = Vec::with_capacity(5);
    let mut jacobian = Vec::with_capacity(5);
    for (i, &idx) in indices.iter().enumerate() {
        residuals.push(sketch.params[idx] - fixed[i]);
        jacobian.push((row + i, idx, 1.0));
    }
    EquationOutput {
        residuals,
        jacobian,
    }
}

// ── Concentric conic constraints ──────────────────────────────────────────

/// ConcentricEllipse(E1, E2): F0 = cx1-cx2, F1 = cy1-cy2
pub(crate) fn eval_concentric_ellipse(
    sketch: &Sketch,
    e1: EllipseHandle,
    e2: EllipseHandle,
    row: usize,
) -> EquationOutput {
    let (cx1, cy1, _, _, _, icx1, icy1, _, _, _) = el(sketch, e1);
    let (cx2, cy2, _, _, _, icx2, icy2, _, _, _) = el(sketch, e2);
    EquationOutput {
        residuals: vec![cx1 - cx2, cy1 - cy2],
        jacobian: vec![
            (row, icx1, 1.0),
            (row, icx2, -1.0),
            (row + 1, icy1, 1.0),
            (row + 1, icy2, -1.0),
        ],
    }
}

/// ConcentricArcEllipse(A, E): F0 = arc_cx - ellipse_cx, F1 = arc_cy - ellipse_cy
pub(crate) fn eval_concentric_arc_ellipse(
    sketch: &Sketch,
    a: ArcHandle,
    e: EllipseHandle,
    row: usize,
) -> EquationOutput {
    let (acx, acy, _, iacx, iacy, _) = ar(sketch, a);
    let (ecx, ecy, _, _, _, iecx, iecy, _, _, _) = el(sketch, e);
    EquationOutput {
        residuals: vec![acx - ecx, acy - ecy],
        jacobian: vec![
            (row, iacx, 1.0),
            (row, iecx, -1.0),
            (row + 1, iacy, 1.0),
            (row + 1, iecy, -1.0),
        ],
    }
}

// ── Dimensional conic constraints ─────────────────────────────────────────

/// EqualSemiAxis(E1, E2): F = a1 - a2
pub(crate) fn eval_equal_semi_axis(
    sketch: &Sketch,
    e1: EllipseHandle,
    e2: EllipseHandle,
    row: usize,
) -> EquationOutput {
    let (_, _, a1, _, _, _, _, ia1, _, _) = el(sketch, e1);
    let (_, _, a2, _, _, _, _, ia2, _, _) = el(sketch, e2);
    EquationOutput {
        residuals: vec![a1 - a2],
        jacobian: vec![(row, ia1, 1.0), (row, ia2, -1.0)],
    }
}

/// FocalDistance(E, d): F = a^2 - b^2 - d^2
pub(crate) fn eval_focal_distance(
    sketch: &Sketch,
    e: EllipseHandle,
    d: f64,
    row: usize,
) -> EquationOutput {
    let (_, _, a, b, _, _, _, ia, ib, _) = el(sketch, e);
    EquationOutput {
        residuals: vec![a * a - b * b - d * d],
        jacobian: vec![(row, ia, 2.0 * a), (row, ib, -2.0 * b)],
    }
}

/// EllipseEccentricity(E, ecc): F = a^2 * ecc^2 - a^2 + b^2
///
/// Eccentricity of an ellipse: e = sqrt(1 - b^2/a^2), so e^2 = 1 - b^2/a^2,
/// hence a^2 * e^2 = a^2 - b^2, rearranged: a^2 * e^2 - a^2 + b^2 = 0.
pub(crate) fn eval_ellipse_eccentricity(
    sketch: &Sketch,
    e: EllipseHandle,
    ecc: f64,
    row: usize,
) -> EquationOutput {
    let (_, _, a, b, _, _, _, ia, ib, _) = el(sketch, e);
    let ecc2 = ecc * ecc;
    EquationOutput {
        residuals: vec![a * a * ecc2 - a * a + b * b],
        jacobian: vec![(row, ia, 2.0 * a * (ecc2 - 1.0)), (row, ib, 2.0 * b)],
    }
}

/// HyperbolaEccentricity(H, ecc): F = a^2 * ecc^2 - a^2 - b^2
///
/// Eccentricity of a hyperbola: e = sqrt(1 + b^2/a^2), so e^2 = 1 + b^2/a^2,
/// hence a^2 * e^2 = a^2 + b^2, rearranged: a^2 * e^2 - a^2 - b^2 = 0.
pub(crate) fn eval_hyperbola_eccentricity(
    sketch: &Sketch,
    h: HyperbolaHandle,
    ecc: f64,
    row: usize,
) -> EquationOutput {
    let (_, _, a, b, _, _, _, ia, ib, _) = hyp(sketch, h);
    let ecc2 = ecc * ecc;
    EquationOutput {
        residuals: vec![a * a * ecc2 - a * a - b * b],
        jacobian: vec![(row, ia, 2.0 * a * (ecc2 - 1.0)), (row, ib, -2.0 * b)],
    }
}
