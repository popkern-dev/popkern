//! Small dense Cholesky factorisation A = L·Lᵀ for symmetric positive
//! definite matrices: the coarsest multigrid level and the EAS
//! condensation.
//!
//! No pivoting, and a fixed loop order (row by row, inner sums in
//! increasing column order), so the factor's bits depend only on the
//! input. A pivot that is not safely positive fails loudly with its
//! index, value and original diagonal; it is never clamped. `sqrt` is
//! correctly rounded by IEEE 754, so it needs no portable wrapper.

/// Why a matrix could not be factored or solved.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CholeskyError {
    /// The value slice is not n×n.
    #[error("{len} values do not form a {n}x{n} matrix")]
    Shape {
        /// Number of values given.
        len: usize,
        /// Requested order.
        n: usize,
    },
    /// An entry of the lower triangle is NaN or infinite.
    #[error("entry ({row}, {col}) is not finite")]
    NonFinite {
        /// Row.
        row: usize,
        /// Column.
        col: usize,
    },
    /// The pivot at `index` is not above `min_pivot_ratio` times the
    /// original diagonal: the matrix is not (numerically) positive
    /// definite.
    #[error("pivot {index} is {value} (diagonal {diagonal}); matrix is not positive definite")]
    NonPositivePivot {
        /// Row and column of the failing pivot.
        index: usize,
        /// The pivot value before the square root.
        value: f64,
        /// The original diagonal entry a_ii.
        diagonal: f64,
    },
    /// A right-hand side has the wrong length.
    #[error("right-hand side has {len} entries, expected {n}")]
    RhsLength {
        /// Entries given.
        len: usize,
        /// Matrix order.
        n: usize,
    },
}

/// The Cholesky factor L of an SPD matrix, stored row-major (the
/// strict upper triangle is zero).
#[derive(Debug, Clone, PartialEq)]
pub struct Cholesky {
    n: usize,
    l: Vec<f64>,
}

impl Cholesky {
    /// Factors the n×n row-major matrix `a`, reading only its lower
    /// triangle (the caller guarantees symmetry).
    ///
    /// Fails when a pivot d_j = a_jj − Σ_k l_jk² is not finite or is at
    /// most `min_pivot_ratio · a_jj`. A ratio of 0 rejects exactly the
    /// non-positive pivots; a small positive ratio (e.g. n·ε) also
    /// rejects numerically singular matrices whose rounded pivot is
    /// tiny but positive.
    pub fn factor(a: &[f64], n: usize, min_pivot_ratio: f64) -> Result<Cholesky, CholeskyError> {
        if n.checked_mul(n) != Some(a.len()) {
            return Err(CholeskyError::Shape { len: a.len(), n });
        }
        for row in 0..n {
            if let Some(col) = (0..=row).find(|&col| !a[row * n + col].is_finite()) {
                return Err(CholeskyError::NonFinite { row, col });
            }
        }
        let mut l = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..i {
                let mut s = a[i * n + j];
                for k in 0..j {
                    s -= l[i * n + k] * l[j * n + k];
                }
                l[i * n + j] = s / l[j * n + j];
            }
            let diagonal = a[i * n + i];
            let mut d = diagonal;
            for k in 0..i {
                d -= l[i * n + k] * l[i * n + k];
            }
            if !(d.is_finite() && d > min_pivot_ratio * diagonal && d > 0.0) {
                return Err(CholeskyError::NonPositivePivot {
                    index: i,
                    value: d,
                    diagonal,
                });
            }
            l[i * n + i] = d.sqrt();
        }
        Ok(Cholesky { n, l })
    }

    /// Order of the factored matrix.
    pub fn order(&self) -> usize {
        self.n
    }

    /// The factor L, row-major n×n, upper triangle zero.
    pub fn factor_l(&self) -> &[f64] {
        &self.l
    }

    /// Solves A·x = b by forward then back substitution, overwriting
    /// `b` with x.
    pub fn solve_in_place(&self, b: &mut [f64]) -> Result<(), CholeskyError> {
        let n = self.n;
        if b.len() != n {
            return Err(CholeskyError::RhsLength { len: b.len(), n });
        }
        let l = &self.l;
        for i in 0..n {
            let mut s = b[i];
            for k in 0..i {
                s -= l[i * n + k] * b[k];
            }
            b[i] = s / l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut s = b[i];
            for k in i + 1..n {
                s -= l[k * n + i] * b[k];
            }
            b[i] = s / l[i * n + i];
        }
        Ok(())
    }

    /// Solves A·x = b, returning x.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, CholeskyError> {
        let mut x = b.to_vec();
        self.solve_in_place(&mut x)?;
        Ok(x)
    }
}
