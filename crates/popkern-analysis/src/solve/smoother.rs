//! Chebyshev smoothing on D⁻¹A (Arbenz §5.2).
//!
//! λ_max of D⁻¹A comes from 10 Lanczos steps on the symmetric
//! D^(−1/2) A D^(−1/2); the polynomial targets [λ_max/30, 1.1 λ_max].
//! The same polynomial is applied before and after the coarse
//! correction, so the V-cycle is symmetric and serves as a CG
//! preconditioner. If Lanczos gives no usable estimate the smoother is
//! damped Jacobi with ω = 0.6, and the hierarchy counts it.

use crate::linalg::sparse::{Csr, SparseError};

/// Chebyshev polynomial degree.
pub const DEGREE: usize = 4;
/// Lanczos steps for the λ_max estimate.
pub const LANCZOS_STEPS: usize = 10;
/// Lower end of the smoothing interval as a fraction of λ_max.
pub const LOWER_FRACTION: f64 = 1.0 / 30.0;
/// Safety factor on the λ_max estimate.
pub const UPPER_SAFETY: f64 = 1.1;
/// Damping of the Jacobi fallback.
pub const JACOBI_OMEGA: f64 = 0.6;

/// A Chebyshev smoother (or its Jacobi fallback) for one operator.
#[derive(Debug, Clone, PartialEq)]
pub struct Chebyshev {
    inv_diag: Vec<f64>,
    lambda_max: f64,
    jacobi: bool,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// λ_max of D^(−1/2) A D^(−1/2) by `LANCZOS_STEPS` Lanczos steps,
/// then power iteration on the tridiagonal matrix. `None` on breakdown
/// before one full step or a non-finite result.
fn lanczos_max(a: &Csr, inv_sqrt_d: &[f64]) -> Option<f64> {
    let n = a.rows();
    if n == 0 {
        return None;
    }
    let apply = |v: &[f64]| -> Option<Vec<f64>> {
        let scaled: Vec<f64> = v.iter().zip(inv_sqrt_d).map(|(x, s)| x * s).collect();
        let mut y = vec![0.0; n];
        a.mul_vec(&scaled, &mut y).ok()?;
        Some(y.iter().zip(inv_sqrt_d).map(|(x, s)| x * s).collect())
    };
    // Deterministic start: alternating-weight vector, normalised.
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 + 0.5 * ((i % 7) as f64) / 7.0).collect();
    let norm = dot(&v, &v).sqrt();
    v.iter_mut().for_each(|x| *x /= norm);
    let mut v_prev = vec![0.0; n];
    let (mut alpha, mut beta) = (Vec::new(), Vec::new());
    let mut b_prev = 0.0;
    for _ in 0..LANCZOS_STEPS.min(n) {
        let mut w = apply(&v)?;
        let a_k = dot(&w, &v);
        for i in 0..n {
            w[i] -= a_k * v[i] + b_prev * v_prev[i];
        }
        alpha.push(a_k);
        let b_k = dot(&w, &w).sqrt();
        if b_k.is_nan() || b_k <= 1e-14 * a_k.abs() {
            break;
        }
        beta.push(b_k);
        v_prev = v;
        v = w.into_iter().map(|x| x / b_k).collect();
        b_prev = b_k;
    }
    let m = alpha.len();
    if m == 0 {
        return None;
    }
    // Largest eigenvalue of the m×m tridiagonal: power iteration.
    let mut x = vec![1.0; m];
    let mut lambda = 0.0;
    for _ in 0..200 {
        let mut y = vec![0.0; m];
        for i in 0..m {
            y[i] = alpha[i] * x[i];
            if i > 0 {
                y[i] += beta[i - 1] * x[i - 1];
            }
            if i + 1 < m {
                y[i] += beta[i] * x[i + 1];
            }
        }
        let ny = dot(&y, &y).sqrt();
        if ny.is_nan() || ny <= 0.0 {
            return None;
        }
        lambda = dot(&y, &x) / dot(&x, &x);
        x = y.into_iter().map(|v| v / ny).collect();
    }
    (lambda.is_finite() && lambda > 0.0).then_some(lambda)
}

impl Chebyshev {
    /// The smoother for `a` (SPD with a positive diagonal).
    pub fn new(a: &Csr) -> Chebyshev {
        let d = a.diagonal();
        let positive = d.iter().all(|&x| x > 0.0 && x.is_finite());
        let inv_diag: Vec<f64> = d
            .iter()
            .map(|&x| if x > 0.0 { 1.0 / x } else { 0.0 })
            .collect();
        let estimate = if positive {
            let inv_sqrt: Vec<f64> = d.iter().map(|x| 1.0 / x.sqrt()).collect();
            lanczos_max(a, &inv_sqrt)
        } else {
            None
        };
        match estimate {
            Some(l) => Chebyshev {
                inv_diag,
                lambda_max: UPPER_SAFETY * l,
                jacobi: false,
            },
            None => Chebyshev {
                inv_diag,
                lambda_max: 0.0,
                jacobi: true,
            },
        }
    }

    /// True when this level fell back to damped Jacobi.
    pub fn is_jacobi(&self) -> bool {
        self.jacobi
    }

    /// The (safety-scaled) upper end of the smoothing interval.
    pub fn lambda_max(&self) -> f64 {
        self.lambda_max
    }

    /// Improves `x` towards A⁻¹ b with one application of the
    /// polynomial (or `DEGREE` damped Jacobi sweeps).
    pub fn smooth(&self, a: &Csr, b: &[f64], x: &mut [f64]) -> Result<(), SparseError> {
        let n = b.len();
        let mut ax = vec![0.0; n];
        a.mul_vec(x, &mut ax)?;
        let r: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
        self.smooth_residual(a, b, x, r)
    }

    /// [`Chebyshev::smooth`] from x = 0, whose residual is b: one
    /// product fewer (b − 0 is b, so the result is the same).
    pub fn smooth_from_zero(&self, a: &Csr, b: &[f64], x: &mut [f64]) -> Result<(), SparseError> {
        x.fill(0.0);
        self.smooth_residual(a, b, x, b.to_vec())
    }

    /// The smoothing step from `x` with residual `r` = b − A x.
    fn smooth_residual(
        &self,
        a: &Csr,
        b: &[f64],
        x: &mut [f64],
        mut r: Vec<f64>,
    ) -> Result<(), SparseError> {
        let n = b.len();
        let mut ax = vec![0.0; n];
        if self.jacobi {
            for _ in 0..DEGREE {
                for i in 0..n {
                    x[i] += JACOBI_OMEGA * self.inv_diag[i] * r[i];
                }
                a.mul_vec(x, &mut ax)?;
                for i in 0..n {
                    r[i] = b[i] - ax[i];
                }
            }
            return Ok(());
        }
        let upper = self.lambda_max;
        let lower = LOWER_FRACTION * upper;
        let theta = 0.5 * (upper + lower);
        let delta = 0.5 * (upper - lower);
        let sigma = theta / delta;
        let mut rho = 1.0 / sigma;
        let mut d: Vec<f64> = r
            .iter()
            .zip(&self.inv_diag)
            .map(|(ri, di)| ri * di / theta)
            .collect();
        let mut ad = vec![0.0; n];
        for _ in 1..DEGREE {
            for (xi, di) in x.iter_mut().zip(&d) {
                *xi += di;
            }
            a.mul_vec(&d, &mut ad)?;
            for (ri, adi) in r.iter_mut().zip(&ad) {
                *ri -= adi;
            }
            let rho_next = 1.0 / (2.0 * sigma - rho);
            for i in 0..n {
                d[i] = rho_next * rho * d[i] + 2.0 * rho_next / delta * self.inv_diag[i] * r[i];
            }
            rho = rho_next;
        }
        for (xi, di) in x.iter_mut().zip(&d) {
            *xi += di;
        }
        Ok(())
    }
}
