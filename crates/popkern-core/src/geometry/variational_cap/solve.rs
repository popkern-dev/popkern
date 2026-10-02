//! The level solves: the stacked Householder QR in the null space's
//! structure, the dense stacked QR behind it, and a ridged Cholesky of
//! the normal equations as the last resort — see [`solve`].

use faer::Mat;

use super::VariationalCapError;
use super::factor::EnergyFactor;
use super::linalg::{cholesky, lstsq, solve_lower};

/// Minimise `λ ‖energy rows‖² + ‖A y − b‖²` by the stacked QR in the
/// energy factor's block structure ([`banded`](super::banded)), with
/// the dense stacked QR behind it for a row that breaks the structure.
///
/// There is no normal-equation fast path (PK-620). Its guard read the
/// diagonal spread of `λG + AᵀA`, which depends on the null-space basis:
/// the dense basis mixed every column and always read past 1e16 on the
/// fits that mattered, but the structured basis's unit interior columns
/// pass it, and the Cholesky then squared a condition number of ~1e7
/// (measured: the octagon h0.4 variable apex's corner-cell G1 step
/// landed 8% over the minimum objective, and the fit bent to 3.57 /
/// r_min instead of 1.59). With the structured QR at ~0.1 s per step
/// the shortcut buys nothing.
pub(super) fn solve(
    energy: &EnergyFactor,
    lambda: f64,
    rows: &Mat<f64>,
    rhs: &[f64],
) -> Result<[Vec<f64>; 3], VariationalCapError> {
    if let Some(y) = super::banded::solve(energy, lambda, rows, rhs) {
        return Ok(y);
    }
    if let Ok(y) = solve_dense_qr(energy, lambda, rows, rhs) {
        return Ok(y);
    }
    // Extreme target weights (rows scaled 1e6..1e7 against λG ~ 1e-2)
    // can defeat both factorisations; a relative ridge, grown until
    // the Cholesky succeeds, still yields a step (it shortens the
    // directions the rows do not see) for the merit to judge.
    solve_ridged(energy, lambda, rows, rhs).ok_or(VariationalCapError::Singular { stage: "fit" })
}

/// The normal equations `(λG + AᵀA) y = Aᵀb − λh` by a Cholesky with a
/// relative ridge (see [`ridged_cholesky`]).
fn solve_ridged(
    energy: &EnergyFactor,
    lambda: f64,
    rows: &Mat<f64>,
    rhs: &[f64],
) -> Option<[Vec<f64>; 3]> {
    use faer::linalg::matmul::matmul;
    let k = energy.k;
    let dim = 3 * k;
    let mut h = energy.dense_gram();
    matmul(
        h.as_mut(),
        rows.transpose(),
        rows.as_ref(),
        Some(lambda),
        1.0,
        faer::Parallelism::None,
    );
    let mut g = vec![0.0; dim];
    for (j, gj) in g.iter_mut().enumerate() {
        // Aᵀb − λ L r, with L r = h.
        *gj = rows
            .col_as_slice(j)
            .iter()
            .zip(rhs)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            - lambda * energy.h[j];
    }
    let l = cholesky(h.clone()).or_else(|| ridged_cholesky(&h))?;
    let z = solve_lower(&l, &g);
    // Back substitution with Lᵀ.
    let mut y = vec![0.0; dim];
    for i in (0..dim).rev() {
        let mut s = z[i];
        for (j, &yj) in y.iter().enumerate().skip(i + 1) {
            s -= l.read(j, i) * yj;
        }
        y[i] = s / l.read(i, i);
    }
    if y.iter().any(|v| !v.is_finite()) {
        return None;
    }
    Some(std::array::from_fn(|c| y[c * k..(c + 1) * k].to_vec()))
}

/// Minimise `λ ‖Lᵀ y + r‖² + ‖A y − b‖²` by the dense stacked QR.
fn solve_dense_qr(
    energy: &EnergyFactor,
    lambda: f64,
    rows: &Mat<f64>,
    rhs: &[f64],
) -> Result<[Vec<f64>; 3], VariationalCapError> {
    let k = energy.k;
    let dim = 3 * k;
    let sl = lambda.sqrt();
    let (factor_rows, r) = energy.dense_rows();
    let top = factor_rows.nrows();
    let mut a = Mat::<f64>::zeros(top + rows.nrows(), dim);
    let mut b = Mat::<f64>::zeros(top + rows.nrows(), 1);
    for (i, &ri) in r.iter().enumerate() {
        for j in 0..dim {
            a.write(i, j, sl * factor_rows.read(i, j));
        }
        b.write(i, 0, -sl * ri);
    }
    for (i, &value) in rhs.iter().enumerate().take(rows.nrows()) {
        for j in 0..dim {
            a.write(top + i, j, rows.read(i, j));
        }
        b.write(top + i, 0, value);
    }
    let x = lstsq(a, b, 1e-12).ok_or(VariationalCapError::Singular { stage: "fit" })?;
    Ok(std::array::from_fn(|c| {
        (0..k).map(|i| x.read(c * k + i, 0)).collect()
    }))
}

/// Cholesky of `h + ρI` for the smallest `ρ = 10^{2k} · 1e-14 · max diag`
/// (up to 1e-6 · max diag) that factorises.
fn ridged_cholesky(h: &Mat<f64>) -> Option<Mat<f64>> {
    let dim = h.nrows();
    let diag_max = (0..dim).map(|i| h.read(i, i)).fold(0.0, f64::max);
    let mut ridge = 1e-14 * diag_max;
    while ridge <= 1e-6 * diag_max {
        let mut hr = h.clone();
        for i in 0..dim {
            hr.write(i, i, hr.read(i, i) + ridge);
        }
        if let Some(l) = cholesky(hr) {
            return Some(l);
        }
        ridge *= 100.0;
    }
    None
}
