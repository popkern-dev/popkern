//! B-spline basis function evaluation: values plus first and second
//! derivatives at a knot span (PK-281 / PK-285).
//!
//! Values come from the standard left/right Cox–de Boor triangle;
//! derivatives from the lower-degree rows via the exact recurrences
//! `N′_i^p = p·(N_i^{p−1}/(k[i+p]−k[i]) − N_{i+1}^{p−1}/(k[i+p+1]−k[i+1]))`
//! (applied once for first, twice for second derivatives), with
//! zero-multiplicity denominators contributing zero. Local index `r`
//! corresponds to global basis index `span - degree + r` throughout.

use super::MAX_DE_BOOR_DEGREE;

/// The Cox–de Boor triangle up to `degree`, returning the final row
/// plus the two lower rows needed for derivative recurrences.
///
/// `rows.0` = degree-p values, `rows.1` = degree-(p−1) row,
/// `rows.2` = degree-(p−2) row (all zero-filled beyond their length;
/// lower rows are the trivial `[1.0]` row when the degree is small).
fn triangle_rows(
    degree: usize,
    knots: &[f64],
    span: usize,
    t: f64,
) -> (
    [f64; MAX_DE_BOOR_DEGREE],
    [f64; MAX_DE_BOOR_DEGREE],
    [f64; MAX_DE_BOOR_DEGREE],
) {
    debug_assert!(
        degree < MAX_DE_BOOR_DEGREE,
        "basis triangle: degree {degree} exceeds MAX_DE_BOOR_DEGREE"
    );
    let mut values = [0.0; MAX_DE_BOOR_DEGREE];
    let mut lower = [0.0; MAX_DE_BOOR_DEGREE];
    let mut lower2 = [0.0; MAX_DE_BOOR_DEGREE];
    let mut left = [0.0; MAX_DE_BOOR_DEGREE];
    let mut right = [0.0; MAX_DE_BOOR_DEGREE];

    values[0] = 1.0;
    lower[0] = 1.0;
    lower2[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() < 1e-30 {
                0.0
            } else {
                values[r] / denom
            };
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
        if degree >= 2 && j == degree - 2 {
            lower2[..=j].copy_from_slice(&values[..=j]);
        }
        if j == degree - 1 {
            lower[..=j].copy_from_slice(&values[..=j]);
        }
    }
    (values, lower, lower2)
}

/// One derivative-recurrence step: given the degree-`(d−1)` row
/// `low` (local index base `span - (d-1)`), produce the derivative
/// contribution array for degree `d` (local index base `span - d`).
fn derivative_step(
    d: usize,
    knots: &[f64],
    span: usize,
    low: &[f64; MAX_DE_BOOR_DEGREE],
) -> [f64; MAX_DE_BOOR_DEGREE] {
    let mut out = [0.0; MAX_DE_BOOR_DEGREE];
    if d == 0 {
        return out;
    }
    let df = d as f64;
    for (r, slot) in out.iter_mut().enumerate().take(d + 1) {
        let i = span - d + r;
        let mut val = 0.0;
        if r >= 1 {
            let denom = knots[i + d] - knots[i];
            if denom.abs() > 1e-30 {
                val += low[r - 1] / denom;
            }
        }
        if r < d {
            let denom = knots[i + d + 1] - knots[i + 1];
            if denom.abs() > 1e-30 {
                val -= low[r] / denom;
            }
        }
        *slot = df * val;
    }
    out
}

/// Basis function values and first derivatives at `t` on `span`.
///
/// Returns `(values, first_derivs)` for local indices `0..=degree`.
pub(crate) fn basis_funs_and_first_derivs(
    degree: usize,
    knots: &[f64],
    span: usize,
    t: f64,
) -> ([f64; MAX_DE_BOOR_DEGREE], [f64; MAX_DE_BOOR_DEGREE]) {
    let (values, lower, _) = triangle_rows(degree, knots, span, t);
    let derivs = derivative_step(degree, knots, span, &lower);
    (values, derivs)
}

/// Basis function values, first, and second derivatives at `t`.
///
/// Second derivatives apply the derivative recurrence twice: the
/// degree-(p−2) row lifts to first derivatives of the degree-(p−1)
/// functions, which lift to second derivatives of the degree-p
/// functions. Degrees 0 and 1 have identically zero second
/// derivatives.
pub(crate) fn basis_funs_and_two_derivs(
    degree: usize,
    knots: &[f64],
    span: usize,
    t: f64,
) -> (
    [f64; MAX_DE_BOOR_DEGREE],
    [f64; MAX_DE_BOOR_DEGREE],
    [f64; MAX_DE_BOOR_DEGREE],
) {
    let (values, lower, lower2) = triangle_rows(degree, knots, span, t);
    let derivs = derivative_step(degree, knots, span, &lower);
    let mut derivs2 = [0.0; MAX_DE_BOOR_DEGREE];
    if degree >= 2 {
        // First derivatives of the degree-(p−1) functions from the
        // degree-(p−2) row, then the same (linear, knots-constant)
        // recurrence lifts them to second derivatives of the degree-p
        // functions — identical index bookkeeping both times.
        let lower_derivs = derivative_step(degree - 1, knots, span, &lower2);
        derivs2 = derivative_step(degree, knots, span, &lower_derivs);
    }
    (values, derivs, derivs2)
}
