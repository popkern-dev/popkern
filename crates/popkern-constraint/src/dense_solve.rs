//! The sketch solver's dense algebra on the deterministic kernels of
//! `popkern_math::linalg` (PK-676): faer's products, LU, pivoted QR and
//! SVD chose SIMD kernels by runtime CPU detection, so a sketch could
//! solve to different last bits on two machines. faer `Mat` remains the
//! container.

use faer::Mat;
use popkern_math::linalg::{self, ColPivQr, DenseMat};

fn dense(m: &Mat<f64>) -> DenseMat {
    DenseMat::from_fn(m.nrows(), m.ncols(), |i, j| m.read(i, j))
}

fn mat(d: &DenseMat) -> Mat<f64> {
    Mat::from_fn(d.nrows(), d.ncols(), |i, j| d.get(i, j))
}

/// `a · b`.
pub(crate) fn mul(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    let mut out = DenseMat::zeros(a.nrows(), b.ncols());
    linalg::matmul(&mut out, None, &dense(a), &dense(b), 1.0);
    mat(&out)
}

/// The least-squares solution of `a x = b` by column-pivoted QR (the
/// full `R` solved, as faer's `solve_lstsq` did).
pub(crate) fn lstsq(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    mat(&ColPivQr::new(dense(a)).solve_lstsq(dense(b)))
}

/// `a x = b` by partial-pivot LU. A singular `a` gives NaN entries, as
/// faer's LU did, so a step built from it fails its line search and the
/// solver falls back instead of reading a zero step as converged.
pub(crate) fn lu(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    match linalg::lu_solve(dense(a), dense(b)) {
        Some(x) => mat(&x),
        None => Mat::from_fn(b.nrows(), b.ncols(), |_, _| f64::NAN),
    }
}

/// The thin SVD's singular values (`min(rows, cols)`, nonincreasing)
/// and right vectors (`cols × min(rows, cols)`); `None` when the Jacobi
/// sweeps do not converge.
pub(crate) fn thin_svd(a: &Mat<f64>) -> Option<(Vec<f64>, Mat<f64>)> {
    let k = a.nrows().min(a.ncols());
    let svd = linalg::jacobi_svd(&dense(a))?;
    let v = Mat::from_fn(a.ncols(), k, |i, j| svd.v.get(i, j));
    Some((svd.sigma[..k].to_vec(), v))
}
