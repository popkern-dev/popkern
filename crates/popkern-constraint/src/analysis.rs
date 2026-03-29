//! DOF analysis and constraint status detection.
//!
//! After the Newton–Raphson solver converges, these routines compute the
//! null space of the Jacobian (via SVD) to identify remaining free degrees
//! of freedom, and check each constraint's residual to determine per-
//! constraint satisfaction status.

use crate::constraint::{ConstraintHandle, ConstraintKind, ConstraintStatus};
use crate::equations::{equation_count, evaluate_constraint};
use crate::error::FreedomInfo;
use crate::jacobian::assemble;
use crate::Sketch;
use faer::Mat;

/// Singular-value threshold: values below this are treated as zero (null
/// space).
const SV_THRESHOLD: f64 = 1e-8;

/// Analyse the degrees of freedom remaining in a solved sketch.
///
/// Computes the Jacobian's null space via thin SVD. Each right-singular
/// vector whose singular value is below [`SV_THRESHOLD`] represents one
/// free DOF. The non-zero entries of that vector are mapped back to sketch
/// elements to produce [`FreedomInfo`] records.
pub(crate) fn analyze_dof(sketch: &Sketch) -> Vec<FreedomInfo> {
    let sys = assemble(sketch);
    let n = sys.num_params;
    let m = sys.num_equations;

    if n == 0 {
        return Vec::new();
    }
    if m == 0 {
        // No equations at all — every parameter is free.
        return params_to_freedoms(sketch, n);
    }

    // Build dense Jacobian (m × n).
    let mut j_mat = Mat::<f64>::zeros(m, n);
    for &(row, col, val) in &sys.jacobian_triplets {
        let existing = j_mat.read(row, col);
        j_mat.write(row, col, existing + val);
    }

    // Thin SVD: U (m×k) S (k) V (n×k) where k = min(m, n).
    let svd = j_mat.thin_svd();
    let s = svd.s_diagonal();
    let v = svd.v(); // n × k

    let k = s.nrows();

    // Collect null-space column indices (singular value ≈ 0).
    let mut null_cols: Vec<usize> = Vec::new();
    for i in 0..k {
        if s.read(i).abs() < SV_THRESHOLD {
            null_cols.push(i);
        }
    }

    // Also count any columns of V beyond k (if n > m, there are n-k
    // implicit null-space directions that thin SVD omits). In that case
    // we report them as free parameters without specific directions.
    let implicit_null = n.saturating_sub(k);

    if null_cols.is_empty() && implicit_null == 0 {
        return Vec::new();
    }

    // Map null-space vectors back to elements.
    let mut freedoms: Vec<FreedomInfo> = Vec::new();
    let entry_threshold = 0.1; // component magnitude threshold

    for &col_idx in &null_cols {
        add_freedom_from_vector(sketch, v.as_ref(), col_idx, entry_threshold, &mut freedoms);
    }

    // For implicit null-space directions (n > k), identify which parameters
    // have no Jacobian rows at all by checking the V matrix coverage.
    if implicit_null > 0 {
        add_implicit_freedoms(sketch, &j_mat, n, &mut freedoms);
    }

    // Deduplicate: merge directions for the same element.
    merge_freedoms(&mut freedoms);

    freedoms
}

/// Analyse constraint satisfaction status after solving.
///
/// Sets each constraint's [`ConstraintStatus`] and detects over-constrained
/// systems by comparing the Jacobian rank to the equation count.
///
/// Returns the list of conflicting constraint handles when the system is
/// over-constrained (more equations than Jacobian rank with unsatisfied
/// residuals).
pub(crate) fn analyze_constraints(sketch: &mut Sketch, tolerance: f64) -> Vec<ConstraintHandle> {
    let sys = assemble(sketch);
    let m = sys.num_equations;
    let n = sys.num_params;

    // Phase 1: compute per-constraint status into a temporary vec (read-only).
    let mut status_updates: Vec<(ConstraintHandle, ConstraintStatus)> = Vec::new();
    let mut row_offset: usize = 0;
    for (handle, entry) in sketch.constraints.iter() {
        let neq = equation_count(&entry.constraint);
        if entry.kind != ConstraintKind::Driving {
            let output = evaluate_constraint(sketch, entry, 0);
            let max_residual = output
                .residuals
                .iter()
                .map(|r| r.abs())
                .fold(0.0_f64, f64::max);
            let status = if max_residual < tolerance {
                ConstraintStatus::Satisfied
            } else {
                ConstraintStatus::Violated
            };
            status_updates.push((handle, status));
            continue;
        }
        let max_residual = sys.residuals[row_offset..row_offset + neq]
            .iter()
            .map(|r| r.abs())
            .fold(0.0_f64, f64::max);
        let status = if max_residual < tolerance {
            ConstraintStatus::Satisfied
        } else {
            ConstraintStatus::Violated
        };
        status_updates.push((handle, status));
        row_offset += neq;
    }

    // Phase 2: apply statuses (mutable).
    for (handle, status) in &status_updates {
        if let Some(entry) = sketch.constraints.get_mut(*handle) {
            entry.status = *status;
        }
    }

    // Detect over-constrained: compute Jacobian rank via SVD.
    if m == 0 || n == 0 {
        return Vec::new();
    }

    let mut j_mat = Mat::<f64>::zeros(m, n);
    for &(row, col, val) in &sys.jacobian_triplets {
        let existing = j_mat.read(row, col);
        j_mat.write(row, col, existing + val);
    }

    let svd = j_mat.thin_svd();
    let s = svd.s_diagonal();
    let rank = (0..s.nrows())
        .filter(|&i| s.read(i).abs() >= SV_THRESHOLD)
        .count();

    if m <= rank {
        return Vec::new();
    }

    // Over-constrained: m equations but only `rank` independent ones.
    find_conflicting_constraints(sketch, &sys.residuals, &j_mat, rank, tolerance)
}

/// Identify constraints whose rows are linearly dependent on other rows.
///
/// Strategy: incrementally add constraint row-blocks to a matrix and check
/// if each addition increases the rank. If it does not, that constraint is
/// conflicting.
fn find_conflicting_constraints(
    sketch: &Sketch,
    _residuals: &[f64],
    full_j: &Mat<f64>,
    _full_rank: usize,
    _tolerance: f64,
) -> Vec<ConstraintHandle> {
    let n = full_j.ncols();
    let mut conflicting = Vec::new();
    let mut accumulated_rows: Vec<Vec<f64>> = Vec::new();

    for (handle, entry) in sketch.constraints.iter() {
        if entry.kind != ConstraintKind::Driving {
            continue;
        }
        let neq = equation_count(&entry.constraint);
        let output = evaluate_constraint(sketch, entry, 0);

        let prev_rank = matrix_rank(&accumulated_rows, n);

        // Add this constraint's Jacobian rows.
        for eq_idx in 0..neq {
            let mut row = vec![0.0; n];
            for &(r, col, val) in &output.jacobian {
                if r == eq_idx {
                    row[col] += val;
                }
            }
            accumulated_rows.push(row);
        }

        let new_rank = matrix_rank(&accumulated_rows, n);
        let rank_increase = new_rank - prev_rank;

        if rank_increase < neq {
            conflicting.push(handle);
        }
    }

    conflicting
}

/// Compute the rank of a row-major matrix given as a vector of row vectors.
fn matrix_rank(rows: &[Vec<f64>], ncols: usize) -> usize {
    let m = rows.len();
    if m == 0 || ncols == 0 {
        return 0;
    }
    let mut mat = Mat::<f64>::zeros(m, ncols);
    for (i, row) in rows.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            mat.write(i, j, val);
        }
    }
    let svd = mat.thin_svd();
    let s = svd.s_diagonal();
    (0..s.nrows())
        .filter(|&i| s.read(i).abs() >= SV_THRESHOLD)
        .count()
}

#[path = "analysis_freedoms.rs"]
mod freedom_helpers;
use freedom_helpers::{
    add_freedom_from_vector, add_implicit_freedoms, merge_freedoms, params_to_freedoms,
};

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
