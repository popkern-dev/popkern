//! Dense one-sided Jacobi SVD for the mate solver's small Jacobians
//! (rows: mate residuals; columns: six per free instance). Accurate to
//! working precision on rank-deficient matrices, which is exactly the
//! case the diagnostics need (redundant mates, free motions).

/// `A = U Σ Vᵀ` of an `m × n` matrix given row-major.
pub struct Svd {
    /// Rows of the input.
    pub m: usize,
    /// Columns of the input.
    pub n: usize,
    /// `U Σ` columns (`m × n`, column-major): the orthogonalised columns.
    us: Vec<Vec<f64>>,
    /// Singular values, one per column.
    pub sigma: Vec<f64>,
    /// `V` (`n × n`, column-major).
    pub v: Vec<Vec<f64>>,
}

impl Svd {
    /// Decompose `a` (`m` rows of `n` entries).
    pub fn new(a: &[Vec<f64>], n: usize) -> Svd {
        let m = a.len();
        let mut cols: Vec<Vec<f64>> = (0..n)
            .map(|j| a.iter().map(|row| row[j]).collect())
            .collect();
        let mut v: Vec<Vec<f64>> = (0..n)
            .map(|j| (0..n).map(|i| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        for _sweep in 0..60 {
            let mut rotated = false;
            for p in 0..n {
                for q in p + 1..n {
                    let alpha: f64 = cols[p].iter().map(|x| x * x).sum();
                    let beta: f64 = cols[q].iter().map(|x| x * x).sum();
                    let gamma: f64 = cols[p].iter().zip(&cols[q]).map(|(x, y)| x * y).sum();
                    if gamma.abs() <= 1e-15 * (alpha * beta).sqrt() || gamma == 0.0 {
                        continue;
                    }
                    rotated = true;
                    let zeta = (beta - alpha) / (2.0 * gamma);
                    let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                    let t = if zeta == 0.0 { 1.0 } else { t };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = c * t;
                    rotate(&mut cols, p, q, c, s);
                    rotate(&mut v, p, q, c, s);
                }
            }
            if !rotated {
                break;
            }
        }
        let sigma = cols
            .iter()
            .map(|c| c.iter().map(|x| x * x).sum::<f64>().sqrt())
            .collect();
        Svd {
            m,
            n,
            us: cols,
            sigma,
            v,
        }
    }

    /// The largest singular value.
    pub fn max_sigma(&self) -> f64 {
        self.sigma.iter().copied().fold(0.0, f64::max)
    }

    /// Numerical rank at relative threshold `rel`.
    pub fn rank(&self, rel: f64) -> usize {
        let cut = rel * self.max_sigma();
        self.sigma.iter().filter(|&&s| s > cut).count()
    }

    /// The Levenberg–Marquardt step for `A x = b` with damping `lambda`:
    /// `x = V (Σ² + λ²)⁻¹ Σ Uᵀ b` — between the Gauss–Newton step
    /// (λ → 0) and a short gradient step (λ large).
    pub fn damped_solve(&self, b: &[f64], lambda: f64) -> Vec<f64> {
        let mut x = vec![0.0; self.n];
        for j in 0..self.n {
            let s = self.sigma[j];
            if s == 0.0 {
                continue;
            }
            // u_j = us_j / s; coefficient = s (u_j · b) / (s² + λ²).
            let coefficient: f64 = self.us[j].iter().zip(b).map(|(u, bi)| u * bi).sum::<f64>()
                / (s * s + lambda * lambda);
            for (k, xk) in x.iter_mut().enumerate() {
                *xk += coefficient * self.v[j][k];
            }
        }
        x
    }

    /// The truncated pseudo-inverse applied to `b` (`x = V Σ⁺ Uᵀ b`):
    /// the minimum-norm least-squares solution of `A x = b`.
    pub fn solve(&self, b: &[f64], rel: f64) -> Vec<f64> {
        let cut = rel * self.max_sigma();
        let mut x = vec![0.0; self.n];
        for j in 0..self.n {
            let s = self.sigma[j];
            if s <= cut || s == 0.0 {
                continue;
            }
            // u_j = us_j / s; coefficient = (u_j · b) / s.
            let coefficient: f64 =
                self.us[j].iter().zip(b).map(|(u, bi)| u * bi).sum::<f64>() / (s * s);
            for (k, xk) in x.iter_mut().enumerate() {
                *xk += coefficient * self.v[j][k];
            }
        }
        x
    }

    /// Basis of the (numerical) null space: the `V` columns whose
    /// singular value is at or below the threshold.
    pub fn null_space(&self, rel: f64) -> Vec<&[f64]> {
        let cut = rel * self.max_sigma();
        (0..self.n)
            .filter(|&j| self.sigma[j] <= cut)
            .map(|j| self.v[j].as_slice())
            .collect()
    }
}

/// Apply the plane rotation `(c, s)` to columns `p < q` of `cols`.
fn rotate(cols: &mut [Vec<f64>], p: usize, q: usize, c: f64, s: f64) {
    let (head, tail) = cols.split_at_mut(q);
    for (x, y) in head[p].iter_mut().zip(tail[0].iter_mut()) {
        let (a, b) = (*x, *y);
        *x = c * a - s * b;
        *y = s * a + c * b;
    }
}

/// Numerical rank of the row vectors `rows` (each `n` long).
pub fn rank_of(rows: &[Vec<f64>], n: usize, rel: f64) -> usize {
    if rows.is_empty() || n == 0 {
        return 0;
    }
    Svd::new(rows, n).rank(rel)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_and_least_squares_on_a_deficient_matrix() {
        // Third row = first + second: rank 2.
        let a = vec![
            vec![1.0, 2.0, 0.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 3.0, 1.0],
        ];
        let svd = Svd::new(&a, 3);
        assert_eq!(svd.rank(1e-10), 2);
        assert_eq!(svd.null_space(1e-10).len(), 1);
        // Consistent right-hand side: exact residual zero.
        let b = [1.0, 2.0, 3.0];
        let x = svd.solve(&b, 1e-10);
        for (row, bi) in a.iter().zip(b) {
            let r: f64 = row.iter().zip(&x).map(|(p, q)| p * q).sum::<f64>() - bi;
            assert!(r.abs() < 1e-12, "{r}");
        }
    }
}
