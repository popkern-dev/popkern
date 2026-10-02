//! Dense linear algebra for the fit, strictly single-threaded
//! (`Parallelism::None`): the result must be bit-identical run to run,
//! and faer's global parallelism hands work to a rayon pool whose
//! partitioning is not part of the contract (see `gregory::direct_linalg`
//! for the stack-depth reason as well).

use faer::Mat;
use faer::dyn_stack::{GlobalPodBuffer, PodStack};

/// Singular value decomposition `A = U Σ Vᵀ`.
pub(crate) struct Svd {
    /// Singular values, nonincreasing (length `min(rows, cols)`).
    pub s: Vec<f64>,
    /// Left vectors, `rows × min(rows, cols)` (columns of zero singular
    /// values are zero).
    pub u: Mat<f64>,
    /// Right vectors, `cols × cols`.
    pub v: Mat<f64>,
}

/// SVD of `a` by one-sided Jacobi (see `jacobi_svd`): `U` with
/// `min(rows, cols)` columns, `V` square. `None` when the sweeps do not
/// converge or the factors fail to reproduce `a`
/// ([`SVD_RECONSTRUCTION_TOL`]).
pub(crate) fn svd(a: &Mat<f64>) -> Option<Svd> {
    factor_checked(a)
}

fn factor_checked(a: &Mat<f64>) -> Option<Svd> {
    let (rows, cols) = (a.nrows(), a.ncols());
    let dim = rows.min(cols);
    let (sigma, w, v) = super::jacobi_svd::jacobi_svd(a)?;
    let d = Svd {
        s: sigma[..dim].to_vec(),
        u: w.submatrix(0, 0, rows, dim).to_owned(),
        v,
    };
    // A − U Σ Vᵀ, entrywise against A's largest entry.
    let us = Mat::from_fn(rows, dim, |i, t| d.u.read(i, t) * d.s[t]);
    let mut residual = a.clone();
    faer::linalg::matmul::matmul(
        residual.as_mut(),
        us.as_ref(),
        d.v.submatrix(0, 0, cols, dim).transpose(),
        Some(1.0),
        -1.0,
        faer::Parallelism::None,
    );
    let largest = |m: &Mat<f64>| {
        (0..m.ncols())
            .flat_map(|j| m.col_as_slice(j).iter().copied())
            .fold(0.0_f64, |acc, v| acc.max(v.abs()))
    };
    let scale = largest(a).max(f64::MIN_POSITIVE);
    (largest(&residual) <= SVD_RECONSTRUCTION_TOL * scale).then_some(d)
}

/// Largest entry of `A − U Σ Vᵀ` an SVD may leave, relative to the
/// largest entry of `A` (a correct factorisation leaves ~1e-15).
const SVD_RECONSTRUCTION_TOL: f64 = 1e-10;

/// Least-squares solution of `A X = B` by Householder QR (no normal
/// equations). `A` must be tall (rows >= cols). Returns `None` when
/// `R` is numerically rank deficient: a diagonal entry at or below
/// `rank_tol` times its column's norm in `A`.
pub(crate) fn lstsq(a: Mat<f64>, b: Mat<f64>, rank_tol: f64) -> Option<Mat<f64>> {
    // Rank is judged per column: the rows carry very different weights
    // (target-scaled rim rows against a small fairness weight), so a
    // global diagonal ratio would read a well-posed weighted system as
    // singular.
    let col_norms: Vec<f64> = (0..a.ncols())
        .map(|j| {
            (0..a.nrows())
                .map(|i| a.read(i, j) * a.read(i, j))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    lstsq_with_norms(a, b, rank_tol, &col_norms)
}

/// [`lstsq`] with the rank test against given column norms (a reduced
/// system's columns, measured in the system it was reduced from).
pub(crate) fn lstsq_with_norms(
    mut a: Mat<f64>,
    mut b: Mat<f64>,
    rank_tol: f64,
    col_norms: &[f64],
) -> Option<Mat<f64>> {
    use faer::linalg::qr::no_pivoting::compute::{
        QrComputeParams, qr_in_place, qr_in_place_req, recommended_blocksize,
    };
    use faer::linalg::qr::no_pivoting::solve::{solve_in_place, solve_in_place_req};
    let (rows, cols) = (a.nrows(), a.ncols());
    if rows < cols || b.nrows() != rows || cols == 0 || col_norms.len() != cols {
        return None;
    }
    let blocksize = recommended_blocksize::<f64>(rows, cols);
    let mut householder = Mat::<f64>::zeros(blocksize, cols);
    let params = QrComputeParams::default();
    let req = qr_in_place_req::<f64>(rows, cols, blocksize, faer::Parallelism::None, params)
        .ok()?
        .or(solve_in_place_req::<f64>(cols, blocksize, b.ncols()).ok()?);
    let mut buffer = GlobalPodBuffer::new(req);
    qr_in_place(
        a.as_mut(),
        householder.as_mut(),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
        params,
    );
    if (0..cols).any(|i| {
        let d = a.read(i, i).abs();
        d.is_nan() || d <= rank_tol * col_norms[i]
    }) {
        return None;
    }
    solve_in_place(
        a.as_ref(),
        householder.as_ref(),
        faer::Conj::No,
        b.as_mut(),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
    );
    Some(Mat::from_fn(cols, b.ncols(), |i, j| b.read(i, j)))
}

/// Householder QR of `m`'s first `left` columns, applied to the rest:
/// on return `m` is `Qᵀ m`: its first `left` columns are `R` over zeros,
/// every other column is `Qᵀ` times its input. With `rank = Some((norms,
/// tol))`, `None` when `R` is numerically rank deficient against those
/// column norms (as [`lstsq_with_norms`]'s test); `None` too when `m`
/// has fewer rows than `left`.
pub(crate) fn qr_eliminate(
    m: &mut Mat<f64>,
    left: usize,
    rank: Option<(&[f64], f64)>,
) -> Option<()> {
    use faer::linalg::householder::{
        apply_block_householder_sequence_transpose_on_the_left_in_place_req,
        apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj,
    };
    use faer::linalg::qr::no_pivoting::compute::{
        QrComputeParams, qr_in_place, qr_in_place_req, recommended_blocksize,
    };
    let rows = m.nrows();
    let right = m.ncols() - left;
    if rows < left {
        return None;
    }
    if left == 0 {
        return Some(());
    }
    let blocksize = recommended_blocksize::<f64>(rows, left);
    let mut householder = Mat::<f64>::zeros(blocksize, left);
    let params = QrComputeParams::default();
    let req = qr_in_place_req::<f64>(rows, left, blocksize, faer::Parallelism::None, params)
        .ok()?
        .or(
            apply_block_householder_sequence_transpose_on_the_left_in_place_req::<f64>(
                rows, blocksize, right,
            )
            .ok()?,
        );
    let mut buffer = GlobalPodBuffer::new(req);
    let (mut a, mut rest) = m.as_mut().split_at_col_mut(left);
    qr_in_place(
        a.as_mut(),
        householder.as_mut(),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
        params,
    );
    if let Some((col_norms, rank_tol)) = rank
        && (col_norms.len() != left
            || (0..left).any(|i| {
                let d = a.read(i, i).abs();
                d.is_nan() || d <= rank_tol * col_norms[i]
            }))
    {
        return None;
    }
    apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj(
        a.as_ref(),
        householder.as_ref(),
        faer::Conj::No,
        rest.as_mut(),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
    );
    // The QR leaves its Householder vectors below the diagonal; the
    // caller reads the transformed rows, where they are zero.
    for j in 0..left {
        for i in j + 1..rows {
            a.write(i, j, 0.0);
        }
    }
    Some(())
}

/// `(0..n).map(f)` in parallel when the feature allows, collected in
/// order: each item is independent, so the result is the serial one.
pub(super) fn par_map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        (0..n).into_par_iter().map(f).collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        (0..n).map(f).collect()
    }
}

/// Lower Cholesky factor `L` of the symmetric positive definite `a`
/// (`a = L Lᵀ`, lower triangle only); `None` when a pivot is not
/// positive.
pub(super) fn cholesky(mut a: Mat<f64>) -> Option<Mat<f64>> {
    use faer::linalg::cholesky::llt::compute::{
        LltParams, LltRegularization, cholesky_in_place, cholesky_in_place_req,
    };
    let n = a.nrows();
    let params = LltParams::default();
    let req = cholesky_in_place_req::<f64>(n, faer::Parallelism::None, params).ok()?;
    let mut buffer = GlobalPodBuffer::new(req);
    cholesky_in_place(
        a.as_mut(),
        LltRegularization::default(),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
        params,
    )
    .ok()?;
    Some(Mat::from_fn(
        n,
        n,
        |i, j| if j <= i { a.read(i, j) } else { 0.0 },
    ))
}

/// Solve `L x = b` for lower-triangular `L`.
pub(super) fn solve_lower(l: &Mat<f64>, b: &[f64]) -> Vec<f64> {
    let n = l.nrows();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for (k, &xk) in x.iter().enumerate().take(i) {
            s -= l.read(i, k) * xk;
        }
        x[i] = s / l.read(i, i);
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lstsq_recovers_an_overdetermined_consistent_solution() {
        let a = Mat::from_fn(7, 3, |i, j| {
            ((i * 3 + j) as f64 * 0.7).sin() + (i == j) as u8 as f64
        });
        let x = [0.5, -1.25, 2.0];
        let b = Mat::from_fn(7, 1, |i, _| (0..3).map(|j| a.read(i, j) * x[j]).sum());
        let sol = lstsq(a, b, 1e-12).expect("full rank");
        for (j, xj) in x.iter().enumerate() {
            assert!((sol.read(j, 0) - xj).abs() < 1e-12);
        }
    }

    #[test]
    fn svd_reproduces_rank_deficient_matrices_with_zero_rows() {
        // A rank-21 product with a zero first row and column, as a spoke
        // block's: faer's SVD did not reproduce a zero first column.
        for &(r, c) in &[(5usize, 7usize), (7, 5), (34, 36), (36, 34), (40, 40)] {
            let rank = 21.min(r).min(c);
            let b = Mat::from_fn(r, rank, |i, j| {
                ((i * 7 + j * 13) as f64 * 0.37).sin() + if i == j { 0.5 } else { 0.0 }
            });
            let f = Mat::from_fn(rank, c, |i, j| ((i * 3 + j * 5) as f64 * 0.71).cos());
            let mut a = &b * &f;
            for j in 0..c {
                a.write(0, j, 0.0);
            }
            for i in 0..r {
                a.write(i, 0, 0.0);
            }
            let d = svd(&a).expect("svd");
            let mut worst = 0.0_f64;
            for i in 0..r {
                for j in 0..c {
                    let v: f64 = (0..d.s.len())
                        .map(|t| d.u.read(i, t) * d.s[t] * d.v.read(j, t))
                        .sum();
                    worst = worst.max((v - a.read(i, j)).abs());
                }
            }
            assert!(worst < 1e-12, "{r} x {c}: reconstruction error {worst:.2e}");
            assert_eq!((d.u.nrows(), d.v.nrows(), d.v.ncols()), (r, c, c));
        }
    }

    #[test]
    fn svd_exposes_the_null_space_of_a_wide_matrix() {
        // Rows are e0 - e1 and e1 - e2: null space spanned by (1,1,1).
        let a = Mat::from_fn(2, 3, |i, j| match (i, j) {
            (0, 0) | (1, 1) => 1.0,
            (0, 1) | (1, 2) => -1.0,
            _ => 0.0,
        });
        let d = svd(&a).expect("svd");
        let rank = d.s.iter().filter(|&&s| s > 1e-12).count();
        assert_eq!(rank, 2);
        let z: Vec<f64> = (0..3).map(|i| d.v.read(i, 2)).collect();
        assert!((z[0] - z[1]).abs() < 1e-14 && (z[1] - z[2]).abs() < 1e-14);
    }
}
