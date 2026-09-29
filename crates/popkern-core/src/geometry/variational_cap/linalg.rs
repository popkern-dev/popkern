//! Dense linear algebra for the fit, strictly single-threaded
//! (`Parallelism::None`): the result must be bit-identical run to run,
//! and faer's global parallelism hands work to a rayon pool whose
//! partitioning is not part of the contract (see `gregory::direct_linalg`
//! for the stack-depth reason as well).

use faer::Mat;
use faer::dyn_stack::{GlobalPodBuffer, PodStack};

/// Thin-U, full-V singular value decomposition `A = U Σ Vᵀ`.
pub(crate) struct Svd {
    /// Singular values, nonincreasing (length `min(rows, cols)`).
    pub s: Vec<f64>,
    /// Left vectors, `rows × min(rows, cols)`.
    pub u: Mat<f64>,
    /// Right vectors, `cols × cols`.
    pub v: Mat<f64>,
}

/// SVD of `a`; `None` only if the workspace size overflows.
pub(crate) fn svd(a: &Mat<f64>) -> Option<Svd> {
    use faer::linalg::svd::{ComputeVectors, SvdParams, compute_svd, compute_svd_req};
    let (rows, cols) = (a.nrows(), a.ncols());
    let dim = rows.min(cols);
    let params = SvdParams::default();
    let req = compute_svd_req::<f64>(
        rows,
        cols,
        ComputeVectors::Thin,
        ComputeVectors::Full,
        faer::Parallelism::None,
        params,
    )
    .ok()?;
    let mut buffer = GlobalPodBuffer::new(req);
    let mut s = faer::Col::<f64>::zeros(dim);
    let mut u = Mat::<f64>::zeros(rows, dim);
    let mut v = Mat::<f64>::zeros(cols, cols);
    compute_svd(
        a.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        Some(v.as_mut()),
        faer::Parallelism::None,
        PodStack::new(&mut buffer),
        params,
    );
    Some(Svd {
        s: (0..dim).map(|i| s.read(i)).collect(),
        u,
        v,
    })
}

/// Least-squares solution of `A X = B` by Householder QR (no normal
/// equations). `A` must be tall (rows >= cols). Returns `None` when
/// `R` is numerically rank deficient: a diagonal entry at or below
/// `rank_tol` times its column's norm in `A`.
pub(crate) fn lstsq(mut a: Mat<f64>, mut b: Mat<f64>, rank_tol: f64) -> Option<Mat<f64>> {
    use faer::linalg::qr::no_pivoting::compute::{
        QrComputeParams, qr_in_place, qr_in_place_req, recommended_blocksize,
    };
    use faer::linalg::qr::no_pivoting::solve::{solve_in_place, solve_in_place_req};
    let (rows, cols) = (a.nrows(), a.ncols());
    if rows < cols || b.nrows() != rows || cols == 0 {
        return None;
    }
    // Rank is judged per column: the rows carry very different weights
    // (target-scaled rim rows against a small fairness weight), so a
    // global diagonal ratio would read a well-posed weighted system as
    // singular.
    let col_norms: Vec<f64> = (0..cols)
        .map(|j| {
            (0..rows)
                .map(|i| a.read(i, j) * a.read(i, j))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
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
