//! Small, deterministic linear algebra: a dense Cholesky factorisation,
//! small symmetric eigenproblems and a 3×3-block sparse matrix.

pub mod csr;
pub mod dense_chol;
pub mod dense_eigen;
pub mod sparse;

pub use csr::{Block, BlockCsr, BlockTriplet, CsrError};
pub use dense_chol::{Cholesky, CholeskyError};

#[cfg(test)]
#[path = "csr_tests.rs"]
mod csr_tests;
#[cfg(test)]
#[path = "dense_chol_tests.rs"]
mod dense_chol_tests;
#[cfg(test)]
#[path = "dense_eigen_tests.rs"]
mod dense_eigen_tests;
#[cfg(test)]
#[path = "sparse_tests.rs"]
mod sparse_tests;
#[cfg(test)]
#[path = "spd_test_helpers.rs"]
pub(crate) mod spd_test_helpers;
