//! Dense linear algebra primitives shared across the
//! direct-construction and KP20-recovery modules.
//!
//! The direct-construction pipeline (corner-C Shirman-Séquin,
//! blend-coefficient 2×2 solves, row-1 KKT min-norm) and the KP20
//! §4.2/§4.3 scalar τ solves all route through the same primitives.
//! Collecting them here keeps the math-layer modules focused on the
//! construction logic and makes the numerical primitives auditable
//! in one place.
//!
//! All helpers are deterministic partial-pivoting Gaussian
//! elimination. The module originally served per-spoke constraint
//! systems (≤ ~30×36); the KP20 scalar KKT systems now push order
//! ~500–800 matrices through [`solve_dense`], so its elimination
//! loop is allocation-free (PK-289) and its arithmetic order is
//! pinned bit-for-bit by the clone-based reference oracle in the
//! test module — any future solver rework (blocked elimination,
//! factorization reuse) must either preserve that pin or replace it
//! deliberately alongside a battery-wide ranking/verdict re-check.

use crate::geometry::Vec3;
use faer::Mat;

/// Singular values of a dense matrix, in nonincreasing order,
/// computed strictly on the calling thread.
///
/// `Mat::singular_values()` runs under faer's global parallelism,
/// which hands the divide-and-conquer bidiagonal SVD to rayon's
/// global worker pool. Those workers keep the default ~2 MiB stack,
/// and the debug-profile recursion overflows it (PK-120) even when
/// the caller runs on a deliberately roomy thread. Pinning
/// `Parallelism::None` keeps every frame on the caller's stack; the
/// matrices here are small enough that parallel SVD buys nothing.
///
/// Returns `None` only if the workspace size computation overflows
/// (not reachable for our systems — the largest are the ~800×800
/// KKT matrices measured by the `*_diagnosed` solve entry points).
fn sequential_singular_values(mat: &Mat<f64>) -> Option<Vec<f64>> {
    use faer::dyn_stack::{GlobalPodBuffer, PodStack};
    use faer::linalg::svd::{ComputeVectors, SvdParams, compute_svd, compute_svd_req};
    let dim = usize::min(mat.nrows(), mat.ncols());
    let params = SvdParams::default();
    let req = compute_svd_req::<f64>(
        mat.nrows(),
        mat.ncols(),
        ComputeVectors::No,
        ComputeVectors::No,
        faer::Parallelism::None,
        params,
    )
    .ok()?;
    let mut buffer = GlobalPodBuffer::new(req);
    let mut s = faer::Col::<f64>::zeros(dim);
    compute_svd(
        mat.as_ref(),
        s.as_mut(),
        None,
        None,
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
        params,
    );
    Some((0..dim).map(|i| s.read(i)).collect())
}

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
    let Some(singular_values) = sequential_singular_values(&mat) else {
        return 0;
    };
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
    let singular_values = sequential_singular_values(&mat)?;
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
///
/// The elimination's floating-point operation order is a contract:
/// results are bit-identical to the pre-PK-289 clone-based loop
/// (see the reference-oracle test), so KP20 battery rankings and
/// every downstream verdict are unchanged by the allocation-free
/// rewrite.
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
        let (pivot_rows, elimination_rows) = m.split_at_mut(i + 1);
        let pivot_row = &pivot_rows[i];
        for (offset, row) in elimination_rows.iter_mut().enumerate() {
            let factor = row[i] / pivot_row[i];
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                row[c] -= factor * pc;
            }
            b[i + 1 + offset] -= factor * b[i];
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

/// Finite-guarded dense solve: faer partial-pivot LU first, falling
/// back to a normalized Gauss–Jordan elimination when the LU
/// solution carries non-finite entries. Returns `None` on any
/// non-finite input or intermediate, or when a pivot falls below
/// `1e-18`.
///
/// This is the regularized-least-squares route of the KP20 τ
/// assembly (normal equations of campaign-diagnostic fits), moved
/// here verbatim (PK-289) so the dense primitives live in one
/// auditable place. It intentionally differs from [`solve_dense`]:
/// the faer LU fast path and the explicit finiteness guards suit
/// normal-equation systems that may carry non-finite contamination,
/// while `solve_dense`'s fixed elimination order is the
/// bit-identity contract for the KKT solves.
pub(crate) fn solve_dense_finite_guarded(
    mut matrix: Vec<Vec<f64>>,
    mut rhs: Vec<f64>,
) -> Option<Vec<f64>> {
    use faer::prelude::SpSolver;
    let n = rhs.len();
    if matrix.iter().flatten().any(|value| !value.is_finite())
        || rhs.iter().any(|value| !value.is_finite())
    {
        return None;
    }
    let mut faer_matrix = Mat::<f64>::zeros(n, n);
    let mut faer_rhs = Mat::<f64>::zeros(n, 1);
    for row in 0..n {
        faer_rhs.write(row, 0, rhs[row]);
        for (column, &value) in matrix[row][..n].iter().enumerate() {
            faer_matrix.write(row, column, value);
        }
    }
    let lu = faer_matrix.partial_piv_lu();
    let solution = lu.solve(&faer_rhs);
    let out = (0..n).map(|row| solution.read(row, 0)).collect::<Vec<_>>();
    if out.iter().all(|value| value.is_finite()) {
        return Some(out);
    }

    for pivot in 0..n {
        let mut best = pivot;
        let mut best_abs = matrix[pivot][pivot].abs();
        if !best_abs.is_finite() {
            return None;
        }
        for (candidate, row) in matrix.iter().enumerate().skip(pivot + 1) {
            let value = row[pivot].abs();
            if !value.is_finite() {
                return None;
            }
            if value > best_abs {
                best = candidate;
                best_abs = value;
            }
        }
        if best_abs < 1.0e-18 {
            return None;
        }
        if best != pivot {
            matrix.swap(best, pivot);
            rhs.swap(best, pivot);
        }
        let pivot_value = matrix[pivot][pivot];
        if !pivot_value.is_finite() {
            return None;
        }
        for value in matrix[pivot].iter_mut().skip(pivot) {
            *value /= pivot_value;
            if !value.is_finite() {
                return None;
            }
        }
        rhs[pivot] /= pivot_value;
        if !rhs[pivot].is_finite() {
            return None;
        }
        let pivot_row = matrix[pivot].clone();
        for row_index in 0..n {
            if row_index == pivot {
                continue;
            }
            let factor = matrix[row_index][pivot];
            if !factor.is_finite() {
                return None;
            }
            if factor == 0.0 {
                continue;
            }
            for column in pivot..n {
                matrix[row_index][column] -= factor * pivot_row[column];
                if !matrix[row_index][column].is_finite() {
                    return None;
                }
            }
            rhs[row_index] -= factor * rhs[pivot];
            if !rhs[row_index].is_finite() {
                return None;
            }
        }
    }
    if rhs.iter().any(|value| !value.is_finite()) {
        return None;
    }
    Some(rhs)
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
