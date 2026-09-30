//! The level solves: the normal equations by Cholesky (fast path),
//! the stacked Householder QR (accurate fallback), and a ridged
//! Cholesky (last resort) — see [`solve`].

use faer::Mat;

use super::VariationalCapError;
use super::level::EnergyFactor;
use super::linalg::{cholesky, lstsq, solve_lower};

/// Minimise `λ ‖Lᵀ y + r‖² + ‖A y − b‖² + Σ μ_j (y_j − y_k,j)²`.
///
/// With rim rows the normal matrix `λ G + AᵀA (+ diag(μ))` is
/// symmetric positive definite (G is) and a Cholesky solve costs ~20×
/// less than the stacked QR (measured need: 2000 unknowns, 5–10 s per QR
/// step). The stacked QR remains the path for the energy-only solve and
/// the fallback when the Cholesky factorisation fails.
pub(super) fn solve(
    energy: &EnergyFactor,
    lambda: f64,
    rim: Option<(&Mat<f64>, &[f64])>,
    damping: Option<(&[f64], &[Vec<f64>; 3])>,
) -> Result<[Vec<f64>; 3], VariationalCapError> {
    if let Some((rows, rhs)) = rim {
        let dim = energy.l.nrows();
        let zeros_mu = vec![0.0; dim];
        let zeros_y: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; dim / 3]);
        let (mu, yk) = damping.unwrap_or((&zeros_mu, &zeros_y));
        if let Some(y) = solve_normal(energy, lambda, rows, rhs, mu, yk, false) {
            return Ok(y);
        }
        if let Ok(y) = solve_qr(energy, lambda, rim, damping) {
            return Ok(y);
        }
        // Extreme target weights (rows scaled 1e6..1e7 against λG ~ 1e-2)
        // can defeat both factorisations; a relative ridge, grown until
        // the Cholesky succeeds, still yields a step (it shortens the
        // directions the rows do not see) for the merit to judge.
        if let Some(y) = solve_normal(energy, lambda, rows, rhs, mu, yk, true) {
            return Ok(y);
        }
    }
    solve_qr(energy, lambda, rim, damping)
}

/// Diagonal spread past which the normal equations cannot factorise
/// (see [`solve_normal`]).
const NORMAL_EQUATIONS_MAX_DIAGONAL_RATIO: f64 = 1e16;

fn solve_normal(
    energy: &EnergyFactor,
    lambda: f64,
    rows: &Mat<f64>,
    rhs: &[f64],
    mu: &[f64],
    yk: &[Vec<f64>; 3],
    ridge_allowed: bool,
) -> Option<[Vec<f64>; 3]> {
    use faer::linalg::matmul::matmul;
    let dim = energy.l.nrows();
    let k = dim / 3;
    // κ(H) ≥ max Hᵢᵢ / min Hᵢᵢ for symmetric positive definite H. Past
    // 1/ε the factor's smallest pivot keeps no correct digit and the
    // factorisation fails (measured on every n = 5..8 pyramid fit:
    // ratios 7e16..1e20, all declined): read the diagonal from the
    // row norms and skip straight to the QR without forming AᵀA
    // (PK-598: ~0.05 s per step at 666 unknowns).
    if !ridge_allowed {
        let diag: Vec<f64> = (0..dim)
            .map(|j| {
                (0..rows.nrows())
                    .map(|i| rows.read(i, j).powi(2))
                    .sum::<f64>()
                    + lambda * energy.gram.read(j, j)
                    + mu[j]
            })
            .collect();
        let (lo, hi) = diag.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &d| {
            (lo.min(d), hi.max(d))
        });
        if !(lo > 0.0 && hi / lo <= NORMAL_EQUATIONS_MAX_DIAGONAL_RATIO) {
            return None;
        }
    }
    let mut h = Mat::<f64>::zeros(dim, dim);
    matmul(
        h.as_mut(),
        rows.transpose(),
        rows.as_ref(),
        None,
        1.0,
        faer::Parallelism::None,
    );
    let mut g = vec![0.0; dim];
    for (j, gj) in g.iter_mut().enumerate() {
        *gj = (0..rows.nrows())
            .map(|i| rows.read(i, j) * rhs[i])
            .sum::<f64>();
    }
    // −λ L r
    for (i, gi) in g.iter_mut().enumerate() {
        let lr: f64 = (0..=i).map(|j| energy.l.read(i, j) * energy.r[j]).sum();
        *gi -= lambda * lr;
    }
    for i in 0..dim {
        for j in 0..dim {
            h.write(i, j, h.read(i, j) + lambda * energy.gram.read(i, j));
        }
        h.write(i, i, h.read(i, i) + mu[i]);
        g[i] += mu[i] * yk[i / k][i % k];
    }
    let l = cholesky(h.clone()).or_else(|| {
        if ridge_allowed {
            ridged_cholesky(&h)
        } else {
            None
        }
    })?;
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

/// Minimise `λ ‖Lᵀ y + r‖² + ‖A y − b‖² (+ damping)` by stacked QR.
fn solve_qr(
    energy: &EnergyFactor,
    lambda: f64,
    rim: Option<(&Mat<f64>, &[f64])>,
    damping: Option<(&[f64], &[Vec<f64>; 3])>,
) -> Result<[Vec<f64>; 3], VariationalCapError> {
    let dim = energy.l.nrows();
    let k = dim / 3;
    let rim_rows = rim.map_or(0, |(a, _)| a.nrows());
    let extra = rim_rows + damping.map_or(0, |_| dim);
    let sl = lambda.sqrt();
    let mut a = Mat::<f64>::zeros(dim + extra, dim);
    let mut b = Mat::<f64>::zeros(dim + extra, 1);
    for i in 0..dim {
        for j in 0..=i {
            // Row j of Lᵀ holds L[i][j] at column i >= j.
            a.write(j, i, sl * energy.l.read(i, j));
        }
        b.write(i, 0, -sl * energy.r[i]);
    }
    if let Some((rows, rhs)) = rim {
        for (i, &value) in rhs.iter().enumerate().take(rows.nrows()) {
            for j in 0..dim {
                a.write(dim + i, j, rows.read(i, j));
            }
            b.write(dim + i, 0, value);
        }
    }
    // Levenberg–Marquardt rows: √(μ_j) (y_j − y_k,j).
    if let Some((mu, yk)) = damping {
        for (j, &m) in mu.iter().enumerate().take(dim) {
            let s = m.sqrt();
            a.write(dim + rim_rows + j, j, s);
            b.write(dim + rim_rows + j, 0, s * yk[j / k][j % k]);
        }
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
