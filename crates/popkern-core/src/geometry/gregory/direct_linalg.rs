//! Small dense linear algebra helpers shared across the
//! direct-construction modules.
//!
//! The direct-construction pipeline (corner-C Shirman-Séquin,
//! blend-coefficient 2×2 solves, row-1 KKT min-norm) all need
//! the same small-system primitives. Collecting them here
//! keeps the math-layer modules focused on the construction
//! logic and makes the numerical primitives auditable in one
//! place.
//!
//! All helpers are deterministic, allocation-disciplined (no
//! hidden heap allocations beyond the input sizes), and use
//! partial-pivoting Gaussian elimination suitable for systems
//! up to ~30×30 — appropriate for this stage of the pipeline
//! where per-spoke constraint systems are at most 27×36.

use crate::geometry::Vec3;
use faer::Mat;

/// Singular-value condition report for a dense matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DenseConditionReport {
    /// Largest singular value.
    pub sigma_max: f64,
    /// Smallest singular value.
    pub sigma_min: f64,
    /// `sigma_max / sigma_min`, or infinity when the matrix is singular.
    pub condition_number: f64,
}

/// Numerical rank of a dense matrix via SVD: count of singular
/// values strictly greater than `tol`. Returns 0 for empty input
/// or rows of inconsistent length.
#[cfg(test)]
pub(crate) fn dense_rank(m: &[Vec<f64>], tol: f64) -> usize {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    if rows == 0 || cols == 0 || m.iter().any(|row| row.len() != cols) {
        return 0;
    }
    let mat = Mat::from_fn(rows, cols, |r, c| m[r][c]);
    let singular_values = mat.singular_values();
    singular_values.iter().filter(|&&s| s > tol).count()
}

/// Measure a dense matrix condition number via SVD.
pub(crate) fn condition_number_dense(m: &[Vec<f64>]) -> Option<DenseConditionReport> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    if rows == 0 || cols == 0 || m.iter().any(|row| row.len() != cols) {
        return None;
    }
    let mat = Mat::from_fn(rows, cols, |r, c| m[r][c]);
    let singular_values = mat.singular_values();
    let sigma_max = singular_values.first().copied().unwrap_or(0.0);
    let sigma_min = singular_values.last().copied().unwrap_or(0.0);
    let condition_number = if sigma_min > 0.0 {
        sigma_max / sigma_min
    } else {
        f64::INFINITY
    };
    Some(DenseConditionReport {
        sigma_max,
        sigma_min,
        condition_number,
    })
}

/// In-place Gaussian elimination with partial pivoting. Solves
/// `M · x = b` for a dense square system. Returns `None` if the
/// matrix is singular to machine precision.
pub(crate) fn solve_dense(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for (r, row) in m.iter().enumerate().skip(i + 1) {
            let v = row[i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != i {
            m.swap(i, pivot);
            b.swap(i, pivot);
        }
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (c, xc) in x.iter().enumerate().skip(i + 1) {
            s -= m[i][c] * xc;
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

/// Deterministic choice of an orthonormal basis `(t1, t2)` in
/// the plane perpendicular to `n`. Picks the seed axis least
/// aligned with `n` so the result is stable across runs
/// regardless of `n`'s orientation.
pub(crate) fn orthonormal_tangent_basis(n: Vec3) -> (Vec3, Vec3) {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let dot = seed.dot(&n);
    let raw_t1 = Vec3::new(seed.x - dot * n.x, seed.y - dot * n.y, seed.z - dot * n.z);
    let t1 = raw_t1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let t2 = n.cross(&t1);
    (t1, t2)
}

/// Binomial coefficient `C(n, k)` as `f64`. Multiplicative form
/// (no factorial overflow). Exact through `n ≈ 53` under IEEE
/// 754 doubles — well beyond our bidegree-7 / elevation-to-8
/// needs.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut out = 1.0_f64;
    for i in 0..k {
        out *= (n - i) as f64 / (i + 1) as f64;
    }
    out
}

#[cfg(test)]
#[path = "direct_linalg_tests.rs"]
mod tests;
