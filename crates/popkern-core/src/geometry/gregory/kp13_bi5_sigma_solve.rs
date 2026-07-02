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

/// One scalar component of one σ BB coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SigmaScalarCoeff {
    /// Tensor-product BB `u` index.
    pub i: usize,
    /// Tensor-product BB `v` index.
    pub j: usize,
    /// Coordinate axis: `0` for x, `1` for y.
    pub axis: usize,
}

impl SigmaScalarCoeff {
    /// Construct a scalar σ coefficient handle.
    pub const fn new(i: usize, j: usize, axis: usize) -> Self {
        Self { i, j, axis }
    }
}

/// One scalar linear row over the symmetry-reduced σ coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct SigmaScalarConstraint {
    /// Sparse scalar terms.
    pub terms: Vec<(SigmaScalarCoeff, f64)>,
    /// Target scalar value.
    pub target: f64,
}

impl SigmaScalarConstraint {
    /// Pin one scalar σ coefficient.
    pub fn pin(coeff: SigmaScalarCoeff, target: f64) -> Self {
        Self {
            terms: vec![(coeff, 1.0)],
            target,
        }
    }

    /// Construct a homogeneous scalar row.
    pub fn homogeneous(terms: Vec<(SigmaScalarCoeff, f64)>) -> Self {
        Self { terms, target: 0.0 }
    }
}

/// Result of a constrained σ solve.
#[derive(Debug, Clone)]
pub struct SigmaSolveOutput {
    /// Solved σ patch.
    pub sigma: SigmaPatch,
    /// Maximum absolute residual over the supplied scalar rows.
    pub max_constraint_residual: f64,
    /// Number of independent scalar rows used in the KKT solve.
    pub independent_constraint_count: usize,
}

/// Error from the constrained σ solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SigmaSolveError {
    /// A coefficient row referenced an invalid BB index or axis.
    InvalidCoeff(SigmaScalarCoeff),
    /// The KKT system was singular.
    SingularSystem,
}

impl std::fmt::Display for SigmaSolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCoeff(coeff) => write!(f, "invalid sigma coefficient {coeff:?}"),
            Self::SingularSystem => write!(f, "sigma constrained F3 KKT system singular"),
        }
    }
}

impl std::error::Error for SigmaSolveError {}

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

/// Scalar rows pinning the four σ corners in the symmetry-reduced basis.
pub fn sigma_corner_constraints() -> Vec<SigmaScalarConstraint> {
    FIXED_FLAT_INDICES
        .iter()
        .zip(FIXED_FLAT_VALUES)
        .map(|(&flat_index, value)| SigmaScalarConstraint {
            terms: vec![(
                scalar_coeff_for_flat_index(flat_index),
                scalar_coeff_weight_for_flat_index(flat_index),
            )],
            target: value,
        })
        .collect()
}

/// KP20/KP13 §4.3 boundary-line rows for the two sector partition axes.
///
/// In the canonical sector chart, the bottom edge lies on the first sector
/// axis (`sigma_y = 0`) and, by symmetry, the left edge lies on the second
/// sector axis (`sigma_x = 0`). These rows intentionally do not pin the
/// opposite outer tensor-border edges; those must come from the characteristic
/// tensor-border/β construction.
pub fn sigma_sector_axis_boundary_constraints() -> Vec<SigmaScalarConstraint> {
    let mut rows = Vec::new();
    for i in 0..6 {
        rows.push(SigmaScalarConstraint::pin(
            SigmaScalarCoeff::new(i, 0, 1),
            0.0,
        ));
        rows.push(SigmaScalarConstraint::pin(
            SigmaScalarCoeff::new(0, i, 0),
            0.0,
        ));
    }
    rows
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
    solve_sigma_f3_with_constraints(&sigma_corner_constraints())
        .map(|output| output.sigma)
        .unwrap_or_else(|_| SigmaPatch::identity())
}

/// Solve σ F3 minimization subject to scalar rows over the
/// symmetry-reduced coefficient basis.
pub fn solve_sigma_f3_with_constraints(
    constraints: &[SigmaScalarConstraint],
) -> Result<SigmaSolveOutput, SigmaSolveError> {
    let h = build_f3_hessian();
    let dense_rows = dense_sigma_constraints(constraints)?;
    let targets = constraints.iter().map(|row| row.target).collect::<Vec<_>>();
    let active_indices = independent_constraint_indices(&dense_rows);
    let active_rows = active_indices
        .iter()
        .map(|&index| dense_rows[index].clone())
        .collect::<Vec<_>>();
    let active_targets = active_indices
        .iter()
        .map(|&index| targets[index])
        .collect::<Vec<_>>();

    let variable_count = SIGMA_FLAT_LEN;
    let row_count = active_rows.len();
    let size = variable_count + row_count;
    let mut matrix = vec![vec![0.0_f64; size]; size];
    let mut rhs = vec![0.0_f64; size];
    for r in 0..variable_count {
        for c in 0..variable_count {
            matrix[r][c] = h[r][c];
        }
    }
    for (row_index, row) in active_rows.iter().enumerate() {
        let kkt_row = variable_count + row_index;
        rhs[kkt_row] = active_targets[row_index];
        for (var, &weight) in row.iter().enumerate() {
            if weight != 0.0 {
                matrix[var][kkt_row] += weight;
                matrix[kkt_row][var] += weight;
            }
        }
    }
    let Some(solution) = solve_dense(matrix, rhs) else {
        return Err(SigmaSolveError::SingularSystem);
    };

    let mut x_full = [0.0_f64; SIGMA_FLAT_LEN];
    x_full.copy_from_slice(&solution[..SIGMA_FLAT_LEN]);
    let sigma = decode_flat(&x_full);
    Ok(SigmaSolveOutput {
        sigma,
        max_constraint_residual: max_constraint_residual(&dense_rows, &targets, &x_full),
        independent_constraint_count: active_rows.len(),
    })
}

fn dense_sigma_constraints(
    constraints: &[SigmaScalarConstraint],
) -> Result<Vec<Vec<f64>>, SigmaSolveError> {
    constraints
        .iter()
        .map(|constraint| {
            let mut row = vec![0.0_f64; SIGMA_FLAT_LEN];
            for &(coeff, weight) in &constraint.terms {
                let (flat_index, flat_weight) = flat_term_for_coeff(coeff)?;
                row[flat_index] += weight * flat_weight;
            }
            Ok(row)
        })
        .collect()
}

fn flat_term_for_coeff(coeff: SigmaScalarCoeff) -> Result<(usize, f64), SigmaSolveError> {
    if coeff.i >= 6 || coeff.j >= 6 || coeff.axis >= 2 {
        return Err(SigmaSolveError::InvalidCoeff(coeff));
    }
    if coeff.i == coeff.j {
        return Ok((coeff.i, 1.0));
    }
    let (lower_i, lower_j, swapped) = if coeff.i < coeff.j {
        (coeff.i, coeff.j, false)
    } else {
        (coeff.j, coeff.i, true)
    };
    let (x_index, y_index) = off_diag_flat_index(lower_i, lower_j);
    let effective_axis = if swapped { 1 - coeff.axis } else { coeff.axis };
    Ok((
        if effective_axis == 0 {
            x_index
        } else {
            y_index
        },
        1.0,
    ))
}

fn scalar_coeff_for_flat_index(flat_index: usize) -> SigmaScalarCoeff {
    if flat_index < 6 {
        return SigmaScalarCoeff::new(flat_index, flat_index, 0);
    }
    let mut next = 6;
    for i in 0..6 {
        for j in (i + 1)..6 {
            if flat_index == next {
                return SigmaScalarCoeff::new(i, j, 0);
            }
            if flat_index == next + 1 {
                return SigmaScalarCoeff::new(i, j, 1);
            }
            next += 2;
        }
    }
    panic!("flat sigma index {flat_index} outside 0..{SIGMA_FLAT_LEN}");
}

fn scalar_coeff_weight_for_flat_index(_flat_index: usize) -> f64 {
    1.0
}

fn independent_constraint_indices(rows: &[Vec<f64>]) -> Vec<usize> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut indices = Vec::new();
    for (row_idx, row) in rows.iter().enumerate() {
        let mut residual = row.clone();
        for basis_row in &basis {
            let projection = dot(&residual, basis_row);
            for (value, basis_value) in residual.iter_mut().zip(basis_row) {
                *value -= projection * basis_value;
            }
        }
        let norm = dot(&residual, &residual).sqrt();
        if norm > 1.0e-10 {
            for value in &mut residual {
                *value /= norm;
            }
            basis.push(residual);
            indices.push(row_idx);
        }
    }
    indices
}

fn max_constraint_residual(rows: &[Vec<f64>], targets: &[f64], solved: &[f64]) -> f64 {
    rows.iter()
        .zip(targets)
        .map(|(row, target)| (dot(row, solved) - target).abs())
        .fold(0.0_f64, f64::max)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
#[path = "kp13_bi5_sigma_solve_tests.rs"]
mod tests;
