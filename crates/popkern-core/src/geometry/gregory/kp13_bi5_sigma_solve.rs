//! σ F_3 constrained minimisation (KP13 §4.3).
//!
//! Builds the F_3 quadratic form Hessian as a 36×36 matrix on σ's
//! symmetric-reduced scalar parameters, pins the four corners, and
//! solves the resulting KKT system to get the F_3-optimal σ.
//!
//! # Parameter layout (36 scalars)
//!
//! After enforcing the KP13 §4.3 symmetry `σ(u, v) = swap(σ(v, u))`,
//! σ has 36 unique scalar DoFs:
//!
//! | range       | what                                                 |
//! |-------------|------------------------------------------------------|
//! | `[0..6]`    | 6 diagonal cells `cps[i][i]`, each 1 scalar with `x = y` |
//! | `[6..36]`   | 15 off-diagonal pairs `cps[i][j]` for `i < j`, 2 scalars (x, y) each; `cps[j][i] = swap(cps[i][j])` derived |
//!
//! The off-diagonal block has 15 × 2 = 30 scalars in lex order
//! `(0,1), (0,2), ..., (0,5), (1,2), ..., (4,5)`.
//!
//! # Corner constraints
//!
//! Four flat-vector positions are pinned to known values (σ is a
//! "sector-to-sector" reparameterisation whose unit-square corners
//! are fixed by construction):
//!
//! | flat | scalar                | value | reason                  |
//! |------|-----------------------|-------|-------------------------|
//! | 0    | `d[0] = cps[0][0].x`  | 0     | σ(0, 0) = (0, 0), star   |
//! | 5    | `d[5] = cps[5][5].x`  | 1     | σ(1, 1) = (1, 1), vₖ     |
//! | 14   | `o[0][5].x`           | 0     | σ(0, 1) = (0, 1), m_{k+1} |
//! | 15   | `o[0][5].y`           | 1     |                          |
//!
//! By the symmetry identity, `cps[5][0] = swap(cps[0][5]) = (1, 0)`
//! without a separate constraint.
//!
//! That's **32 free scalars** for the F_3 solve. The paper's
//! additional "boundary on sector partition lines + G² via ρ"
//! constraints would reduce this to 15 DoFs per KP13 §4.3; those
//! are still deferred (they involve ρ / β-reparameterised χ
//! boundary data which this solver doesn't consume).

#![allow(dead_code)]

use super::kp13_bi5_sigma::{SigmaPatch, f3_functional};
use crate::geometry::gregory::direct_linalg::solve_dense;

/// Total scalar DoFs after symmetry (no corner pinning yet).
pub const SIGMA_FLAT_LEN: usize = 36;

/// Indices in the 36-flat vector pinned by the corner constraints.
pub const FIXED_FLAT_INDICES: [usize; 4] = [0, 5, 14, 15];

/// Values corresponding to [`FIXED_FLAT_INDICES`].
pub const FIXED_FLAT_VALUES: [f64; 4] = [0.0, 1.0, 0.0, 1.0];

/// Encode the symmetric-reduced 36 scalar DoFs into a full
/// `SigmaPatch` (6×6 2-D CP grid).
pub fn decode_flat(x: &[f64; SIGMA_FLAT_LEN]) -> SigmaPatch {
    let mut cps = [[(0.0, 0.0); 6]; 6];
    // Diagonal cells: x = y.
    #[allow(clippy::needless_range_loop)]
    for i in 0..6 {
        cps[i][i] = (x[i], x[i]);
    }
    // Off-diagonal pairs (i < j).
    let mut idx = 6;
    #[allow(clippy::needless_range_loop)]
    for i in 0..6 {
        for j in (i + 1)..6 {
            let px = x[idx];
            let py = x[idx + 1];
            cps[i][j] = (px, py);
            cps[j][i] = (py, px); // swap
            idx += 2;
        }
    }
    SigmaPatch { cps }
}

/// Flat-index of the off-diagonal pair `(i, j)` with `i < j`.
/// Returns `(x_index, y_index)`.
pub fn off_diag_flat_index(i: usize, j: usize) -> (usize, usize) {
    assert!(i < j, "off_diag_flat_index needs i < j, got ({i}, {j})");
    assert!(j < 6);
    // Lex order: iterate i = 0..5, for each i iterate j = i+1..6.
    // pair_count before (i, j) = sum_{a < i} (5 - a) + (j - i - 1)
    //                         = i*(11 - i)/2 + (j - i - 1)
    let pair_count = i * (11 - i) / 2 + (j - i - 1);
    (6 + 2 * pair_count, 6 + 2 * pair_count + 1)
}

/// Build the F_3 Hessian `H` as a 36×36 matrix via central-difference
/// finite differences on the `f3_functional`. `H_{ij} = ∂²F_3/∂x_i∂x_j`.
///
/// Because F_3 is a quadratic form (σ enters F_3 linearly through the
/// Bernstein basis, so F_3 is bilinear in σ pairs), the Hessian is
/// **exact**: `F_3(x) = ½ x^T H x` for any `x`.
pub fn build_f3_hessian() -> [[f64; SIGMA_FLAT_LEN]; SIGMA_FLAT_LEN] {
    // Diagonal entries: H_ii = 2·F_3(e_i).
    // Off-diagonal:     H_ij = F_3(e_i + e_j) - F_3(e_i) - F_3(e_j).
    let mut h = [[0.0_f64; SIGMA_FLAT_LEN]; SIGMA_FLAT_LEN];
    let mut f_ei = [0.0_f64; SIGMA_FLAT_LEN];

    // Precompute F_3(e_i) for each basis vector.
    for i in 0..SIGMA_FLAT_LEN {
        let mut x = [0.0_f64; SIGMA_FLAT_LEN];
        x[i] = 1.0;
        let patch = decode_flat(&x);
        f_ei[i] = f3_functional(&patch);
    }

    for i in 0..SIGMA_FLAT_LEN {
        h[i][i] = 2.0 * f_ei[i];
        for j in (i + 1)..SIGMA_FLAT_LEN {
            let mut x = [0.0_f64; SIGMA_FLAT_LEN];
            x[i] = 1.0;
            x[j] = 1.0;
            let patch = decode_flat(&x);
            let f_sum = f3_functional(&patch);
            let off = f_sum - f_ei[i] - f_ei[j];
            h[i][j] = off;
            h[j][i] = off;
        }
    }
    h
}

/// Solve the corner-pinned σ F_3 minimisation.
///
/// Returns the `SigmaPatch` minimising `F_3` subject to symmetry +
/// the four corner constraints listed in [`FIXED_FLAT_INDICES`] /
/// [`FIXED_FLAT_VALUES`].
///
/// # Algorithm
///
/// Let `x_free` be the 32 free scalars and `x_fixed` be the 4 pinned
/// ones. Block-partitioning `H`:
///
/// ```text
///   F_3(x) = ½ x^T H x
///          = ½ x_free^T H_ff x_free + x_free^T H_fc x_fixed
///          + ½ x_fixed^T H_cc x_fixed
/// ```
///
/// Minimising over `x_free`:
///
/// ```text
///   H_ff · x_free = − H_fc · x_fixed
/// ```
///
/// which we solve via the existing `direct_linalg::solve_dense`.
///
/// Falls back to the identity reparameterisation if the reduced
/// Hessian `H_ff` is singular. That branch is not expected for the
/// well-posed symmetric F_3 system, but identity preserves a valid
/// sector map better than partially-populated coefficients.
pub fn solve_sigma_f3_minimal() -> SigmaPatch {
    let h = build_f3_hessian();

    // Split indices.
    let mut is_fixed = [false; SIGMA_FLAT_LEN];
    for &i in &FIXED_FLAT_INDICES {
        is_fixed[i] = true;
    }
    let free_indices: Vec<usize> = (0..SIGMA_FLAT_LEN).filter(|i| !is_fixed[*i]).collect();
    let n_free = free_indices.len();
    let n_fixed = FIXED_FLAT_INDICES.len();
    assert_eq!(n_free + n_fixed, SIGMA_FLAT_LEN);

    // Build H_ff (n_free × n_free) and H_fc · x_fixed (n_free vector).
    let mut h_ff = vec![vec![0.0_f64; n_free]; n_free];
    let mut rhs = vec![0.0_f64; n_free];
    for (r, &i_free) in free_indices.iter().enumerate() {
        for (c, &j_free) in free_indices.iter().enumerate() {
            h_ff[r][c] = h[i_free][j_free];
        }
        // - H_fc · x_fixed.
        let mut s = 0.0;
        for (c_fixed_pos, &j_fixed) in FIXED_FLAT_INDICES.iter().enumerate() {
            s += h[i_free][j_fixed] * FIXED_FLAT_VALUES[c_fixed_pos];
        }
        rhs[r] = -s;
    }

    let Some(x_free) = solve_dense(h_ff, rhs) else {
        return SigmaPatch::identity();
    };

    // Reassemble full flat vector.
    let mut x_full = [0.0_f64; SIGMA_FLAT_LEN];
    for (r, &i_free) in free_indices.iter().enumerate() {
        x_full[i_free] = x_free[r];
    }
    for (c_fixed_pos, &i_fixed) in FIXED_FLAT_INDICES.iter().enumerate() {
        x_full[i_fixed] = FIXED_FLAT_VALUES[c_fixed_pos];
    }

    decode_flat(&x_full)
}

#[cfg(test)]
#[path = "kp13_bi5_sigma_solve_tests.rs"]
mod tests;
