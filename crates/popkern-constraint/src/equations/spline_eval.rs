//! B-spline evaluation with analytical Jacobians — point evaluation.
//!
//! Provides the core data structures, basis-function computation, and
//! point evaluation for the [`Spline`] element.  Derivative-level
//! functions (tangent, second derivative, curvature, arc length) are in
//! the sibling [`spline_deriv`](super::spline_deriv) module.

use crate::Sketch;
use crate::elements::Spline;

/// Result of evaluating a spline quantity (point, tangent, etc.).
pub(crate) struct SplineEvalResult {
    /// x component of the evaluated quantity.
    pub x: f64,
    /// y component of the evaluated quantity.
    pub y: f64,
    /// Sparse Jacobian: `(param_index, dx/dparam, dy/dparam)`.
    /// Only nonzero for control points whose basis functions are active at `t`.
    pub jacobian: Vec<(usize, f64, f64)>,
}

// ── Knot-vector helpers ──────────────────────────────────────────────────

/// Build the extended knot vector for a periodic spline.
///
/// For a periodic spline with `n` control points and degree `d`, the
/// effective knot vector has `n + 2d + 1` entries obtained by wrapping
/// the interior knots with the period.
pub(crate) fn periodic_knot_vector(knots: &[f64], n: usize, degree: usize) -> Vec<f64> {
    let period = knots[n] - knots[0];
    let mut ext = Vec::with_capacity(n + 2 * degree + 1);
    // Prepend `degree` knots from the end of the previous period.
    // These are knots[n-degree..n] shifted back by one period.
    for &k in &knots[(n - degree)..n] {
        ext.push(k - period);
    }
    // Original knots (all n + 1 of them).
    ext.extend_from_slice(&knots[..=n]);
    // Append `degree` knots from the start of the next period.
    // These are knots[1..=degree] shifted forward by one period.
    for &k in &knots[1..=degree] {
        ext.push(k + period);
    }
    ext
}

/// Find the knot span index `i` such that `knots[i] <= t < knots[i+1]`.
///
/// The parameter is clamped to the active range
/// `[knots[degree], knots[n]]` where `n` is the number of control
/// points.  When `t` equals the upper boundary, the last valid span is
/// returned.
pub(crate) fn find_span(knots: &[f64], degree: usize, n_cp: usize, t: f64) -> (usize, f64) {
    let lo = knots[degree];
    let hi = knots[n_cp];
    let t = t.clamp(lo, hi);

    // Special case: t at the upper boundary.
    if (t - hi).abs() <= f64::EPSILON * hi.abs().max(1.0) {
        // Walk back to find the last span with a nonzero interval.
        let mut i = n_cp;
        while i > degree && (knots[i] - knots[i - 1]).abs() < f64::EPSILON {
            i -= 1;
        }
        return (i - 1, t);
    }

    // Binary search.
    let mut low = degree;
    let mut high = n_cp;
    while low + 1 < high {
        let mid = (low + high) / 2;
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
    }
    (low, t)
}

/// Compute the `degree + 1` nonzero basis function values at `t` for the
/// given knot span `span`.
///
/// Returns `N[span - degree], ..., N[span]` (indices relative to span).
pub(crate) fn basis_functions(knots: &[f64], span: usize, degree: usize, t: f64) -> Vec<f64> {
    let p = degree;
    let mut n = vec![0.0; p + 1];
    n[0] = 1.0;
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    for j in 1..=p {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            if denom.abs() < f64::EPSILON {
                // Degenerate span: 0/0 = 0.
            } else {
                let temp = n[r] / denom;
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
        }
        n[j] = saved;
    }
    n
}

/// Accumulate a 2D Jacobian entry, merging with an existing entry if the
/// parameter index already exists.
pub(super) fn accumulate_jac(jac: &mut Vec<(usize, f64, f64)>, idx: usize, dx: f64, dy: f64) {
    if let Some(entry) = jac.iter_mut().find(|(i, _, _)| *i == idx) {
        entry.1 += dx;
        entry.2 += dy;
    } else {
        jac.push((idx, dx, dy));
    }
}

/// Accumulate a scalar Jacobian entry, merging with an existing entry if
/// the parameter index already exists.
pub(super) fn accumulate_scalar_jac(jac: &mut Vec<(usize, f64)>, idx: usize, val: f64) {
    if let Some(entry) = jac.iter_mut().find(|(i, _)| *i == idx) {
        entry.1 += val;
    } else {
        jac.push((idx, val));
    }
}

// ── Public evaluation functions ──────────────────────────────────────────

/// Evaluate the B-spline at parameter `t`, returning position + Jacobian.
pub(crate) fn eval_point(spline: &Spline, t: f64, sketch: &Sketch) -> SplineEvalResult {
    let n = spline.control_points.len();
    let degree = spline.degree;

    if spline.periodic {
        return eval_point_periodic(spline, t, sketch);
    }

    let knots = &spline.knots;
    let (span, t) = find_span(knots, degree, n, t);
    let basis = basis_functions(knots, span, degree, t);

    let is_rational = spline
        .weights
        .iter()
        .any(|&w| (w - 1.0).abs() > f64::EPSILON);

    if is_rational {
        eval_point_rational(spline, sketch, span, &basis)
    } else {
        eval_point_non_rational(spline, sketch, span, &basis)
    }
}

/// Non-rational point evaluation.
fn eval_point_non_rational(
    spline: &Spline,
    sketch: &Sketch,
    span: usize,
    basis: &[f64],
) -> SplineEvalResult {
    let degree = spline.degree;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut jac = Vec::with_capacity(2 * (degree + 1));
    for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
        let cp_idx = span - degree + k;
        let pt = &sketch.points[spline.control_points[cp_idx]];
        let [ix, iy] = pt.param_indices();
        let px = sketch.params[ix];
        let py = sketch.params[iy];
        x += nk * px;
        y += nk * py;
        if nk.abs() > f64::EPSILON {
            jac.push((ix, nk, 0.0));
            jac.push((iy, 0.0, nk));
        }
    }
    SplineEvalResult {
        x,
        y,
        jacobian: jac,
    }
}

/// Rational (NURBS) point evaluation.
fn eval_point_rational(
    spline: &Spline,
    sketch: &Sketch,
    span: usize,
    basis: &[f64],
) -> SplineEvalResult {
    let degree = spline.degree;
    // Compute weighted sum and total weight.
    let mut wx = 0.0;
    let mut wy = 0.0;
    let mut w_total = 0.0;
    for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
        let cp_idx = span - degree + k;
        let pt = &sketch.points[spline.control_points[cp_idx]];
        let [ix, _iy] = pt.param_indices();
        let px = sketch.params[ix];
        let py = sketch.params[ix + 1];
        let wk = spline.weights[cp_idx];
        wx += nk * wk * px;
        wy += nk * wk * py;
        w_total += nk * wk;
    }
    let inv_w = if w_total.abs() > f64::EPSILON {
        1.0 / w_total
    } else {
        0.0
    };
    let x = wx * inv_w;
    let y = wy * inv_w;

    // Jacobian: dS_x/dP_jx = R_j(t) = N_j * w_j / W
    let mut jac = Vec::with_capacity(2 * (degree + 1));
    for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
        let cp_idx = span - degree + k;
        let pt = &sketch.points[spline.control_points[cp_idx]];
        let [ix, iy] = pt.param_indices();
        let wk = spline.weights[cp_idx];
        let rk = nk * wk * inv_w;
        if rk.abs() > f64::EPSILON {
            jac.push((ix, rk, 0.0));
            jac.push((iy, 0.0, rk));
        }
    }
    SplineEvalResult {
        x,
        y,
        jacobian: jac,
    }
}

/// Periodic spline point evaluation.
fn eval_point_periodic(spline: &Spline, t: f64, sketch: &Sketch) -> SplineEvalResult {
    let n = spline.control_points.len();
    let degree = spline.degree;
    let ext_knots = periodic_knot_vector(&spline.knots, n, degree);
    let n_ext = n + degree; // effective CP count for evaluation
    let (span, t) = find_span(&ext_knots, degree, n_ext, t);
    let basis = basis_functions(&ext_knots, span, degree, t);

    let mut x = 0.0;
    let mut y = 0.0;
    let mut jac: Vec<(usize, f64, f64)> = Vec::with_capacity(2 * (degree + 1));
    for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
        let raw_idx = span - degree + k;
        let cp_idx = raw_idx % n;
        let pt = &sketch.points[spline.control_points[cp_idx]];
        let [ix, iy] = pt.param_indices();
        let px = sketch.params[ix];
        let py = sketch.params[iy];
        x += nk * px;
        y += nk * py;
        if nk.abs() > f64::EPSILON {
            accumulate_jac(&mut jac, ix, nk, 0.0);
            accumulate_jac(&mut jac, iy, 0.0, nk);
        }
    }
    SplineEvalResult {
        x,
        y,
        jacobian: jac,
    }
}
