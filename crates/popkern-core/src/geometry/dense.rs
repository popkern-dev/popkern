//! `faer::Mat` storage ↔ [`DenseMat`] for the deterministic kernels
//! (PK-676). faer stays the container in the fit code (its storage and
//! element access are plain memory); every product, factorisation and
//! solve goes through `popkern_math::linalg`, whose results do not
//! depend on the CPU.

use faer::Mat;
use popkern_math::linalg::DenseMat;

/// `m` as a [`DenseMat`].
pub(crate) fn dense(m: &Mat<f64>) -> DenseMat {
    DenseMat::from_fn(m.nrows(), m.ncols(), |i, j| m.read(i, j))
}

/// `d` as a faer `Mat`.
pub(crate) fn mat(d: &DenseMat) -> Mat<f64> {
    Mat::from_fn(d.nrows(), d.ncols(), |i, j| d.get(i, j))
}
