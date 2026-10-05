//! Dense linear algebra with the same bits on every target (PK-676).
//!
//! faer and gemm choose their SIMD kernels by runtime CPU detection
//! (FMA, AVX2, AVX-512 natively; simd128 in wasm), so a solve's last
//! bits depended on the machine: the gallery's G2 corner caps differed
//! natively and in the browser by 1e-12..7e-10 although everything
//! upstream was bit-identical. These kernels are plain `f64` loops in a
//! fixed order, with no fused multiply-add and no lane-width-dependent
//! reductions, so a result depends only on its inputs.
//!
//! Matrices are column-major ([`DenseMat`]). The kernels are the
//! textbook algorithms at the sizes the kernel solves (hundreds of
//! unknowns): partial-pivot LU, Cholesky, triangular solves,
//! Householder QR, one-sided Jacobi SVD, and products with one fixed
//! summation order.

mod eigen;
mod factor;
mod sparse_ldlt;
mod svd;

pub use eigen::symmetric_eigenvalues;
pub use factor::{
    ColPivQr, HouseholderQr, cholesky, lu_solve, solve_lower_in_place,
    solve_lower_transposed_in_place, solve_upper_in_place,
};
pub use sparse_ldlt::{LdltFactor, LdltSymbolic};
pub use svd::{Svd, jacobi_svd};

/// A dense column-major matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMat {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl DenseMat {
    /// `rows × cols` zeros.
    pub fn zeros(rows: usize, cols: usize) -> DenseMat {
        DenseMat {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Entry `(i, j)` from `f`.
    pub fn from_fn(rows: usize, cols: usize, f: impl Fn(usize, usize) -> f64) -> DenseMat {
        let mut m = DenseMat::zeros(rows, cols);
        for j in 0..cols {
            for i in 0..rows {
                m.data[j * rows + i] = f(i, j);
            }
        }
        m
    }

    /// The identity of size `n`.
    pub fn identity(n: usize) -> DenseMat {
        DenseMat::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 })
    }

    /// Row count.
    pub fn nrows(&self) -> usize {
        self.rows
    }

    /// Column count.
    pub fn ncols(&self) -> usize {
        self.cols
    }

    /// Entry `(i, j)`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.rows + i]
    }

    /// Set entry `(i, j)`.
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.data[j * self.rows + i] = value;
    }

    /// Column `j`.
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// Column `j`, mutably.
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// The transpose (a copy).
    pub fn transpose(&self) -> DenseMat {
        DenseMat::from_fn(self.cols, self.rows, |i, j| self.get(j, i))
    }

    /// The column-major entries.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
}

/// `dst ← keep·dst + beta·(lhs · rhs)` (`keep = None` overwrites),
/// summing over the inner dimension in ascending order for every entry.
///
/// # Panics
///
/// Never for conforming shapes; the shapes are the caller's contract
/// (checked with `debug_assert`).
pub fn matmul(dst: &mut DenseMat, keep: Option<f64>, lhs: &DenseMat, rhs: &DenseMat, beta: f64) {
    debug_assert_eq!(lhs.cols, rhs.rows);
    debug_assert_eq!((dst.rows, dst.cols), (lhs.rows, rhs.cols));
    let mut sum = vec![0.0; lhs.rows];
    for j in 0..rhs.cols {
        sum.iter_mut().for_each(|s| *s = 0.0);
        for k in 0..lhs.cols {
            let b = rhs.get(k, j);
            if b == 0.0 {
                continue;
            }
            for (s, &a) in sum.iter_mut().zip(lhs.col(k)) {
                *s += a * b;
            }
        }
        let column = dst.col_mut(j);
        match keep {
            None => column
                .iter_mut()
                .zip(&sum)
                .for_each(|(d, &s)| *d = beta * s),
            Some(k) => column
                .iter_mut()
                .zip(&sum)
                .for_each(|(d, &s)| *d = k * *d + beta * s),
        }
    }
}

/// `lhsᵀ · rhs` (no transpose copy): entry `(i, j)` is the dot product
/// of column `i` of `lhs` and column `j` of `rhs`, in ascending order.
pub fn matmul_tn(lhs: &DenseMat, rhs: &DenseMat) -> DenseMat {
    debug_assert_eq!(lhs.rows, rhs.rows);
    DenseMat::from_fn(lhs.cols, rhs.cols, |i, j| dot(lhs.col(i), rhs.col(j)))
}

/// The upper triangle of `q qᵀ` (the strict lower triangle zero), each
/// entry summed over the columns of `q` in ascending order.
pub fn syrk_upper(q: &DenseMat) -> DenseMat {
    let n = q.rows;
    let mut g = DenseMat::zeros(n, n);
    for k in 0..q.cols {
        let column = q.col(k);
        for (j, &qj) in column.iter().enumerate() {
            if qj == 0.0 {
                continue;
            }
            let target = &mut g.data[j * n..j * n + j + 1];
            for (t, &qi) in target.iter_mut().zip(column) {
                *t += qi * qj;
            }
        }
    }
    g
}

/// `Σ a_i b_i` in ascending order.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0;
    for (x, y) in a.iter().zip(b) {
        s += x * y;
    }
    s
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
