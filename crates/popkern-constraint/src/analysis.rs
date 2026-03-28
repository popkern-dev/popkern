//! DOF analysis and constraint status detection.
//!
//! After the Newton–Raphson solver converges, these routines compute the
//! null space of the Jacobian (via SVD) to identify remaining free degrees
//! of freedom, and check each constraint's residual to determine per-
//! constraint satisfaction status.

use crate::constraint::{ConstraintHandle, ConstraintKind, ConstraintStatus};
use crate::elements::ElementHandle;
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

/// Map all parameters to freedom entries (used when there are no equations).
fn params_to_freedoms(sketch: &Sketch, _n: usize) -> Vec<FreedomInfo> {
    let mut freedoms = Vec::new();
    for (handle, point) in sketch.points.iter() {
        let _ = point;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Point(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, arc) in sketch.arcs.iter() {
        let _ = arc;
        // Arcs have 5 DOF but we report as element-level freedom.
        freedoms.push(FreedomInfo {
            element: ElementHandle::Arc(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    freedoms
}

/// Extract freedom info from a single null-space vector column of V.
fn add_freedom_from_vector(
    sketch: &Sketch,
    v: faer::MatRef<'_, f64>,
    col_idx: usize,
    threshold: f64,
    freedoms: &mut Vec<FreedomInfo>,
) {
    // Check which elements have significant components in this null-space
    // vector.
    for (handle, point) in sketch.points.iter() {
        let [ix, iy] = point.param_indices();
        let vx = if ix < v.nrows() {
            v.read(ix, col_idx)
        } else {
            0.0
        };
        let vy = if iy < v.nrows() {
            v.read(iy, col_idx)
        } else {
            0.0
        };
        let mag = (vx * vx + vy * vy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Point(handle),
                directions: vec![(vx * inv, vy * inv)],
            });
        }
    }
    for (handle, arc) in sketch.arcs.iter() {
        let [icx, icy, _ir, _is, _ie] = arc.param_indices();
        let vcx = if icx < v.nrows() {
            v.read(icx, col_idx)
        } else {
            0.0
        };
        let vcy = if icy < v.nrows() {
            v.read(icy, col_idx)
        } else {
            0.0
        };
        let mag = (vcx * vcx + vcy * vcy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Arc(handle),
                directions: vec![(vcx * inv, vcy * inv)],
            });
        }
    }
}

/// Add freedoms for parameters that have zero Jacobian coverage (implicit
/// null-space directions when n > min(m, n)).
fn add_implicit_freedoms(
    sketch: &Sketch,
    j_mat: &Mat<f64>,
    n: usize,
    freedoms: &mut Vec<FreedomInfo>,
) {
    // Find parameters with all-zero Jacobian columns.
    let mut col_has_entry = vec![false; n];
    for (col, has_entry) in col_has_entry.iter_mut().enumerate() {
        for row in 0..j_mat.nrows() {
            if j_mat.read(row, col).abs() > 1e-15 {
                *has_entry = true;
                break;
            }
        }
    }

    for (handle, point) in sketch.points.iter() {
        let [ix, iy] = point.param_indices();
        let free_x = ix < n && !col_has_entry[ix];
        let free_y = iy < n && !col_has_entry[iy];
        if free_x || free_y {
            let mut dirs = Vec::new();
            if free_x {
                dirs.push((1.0, 0.0));
            }
            if free_y {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Point(handle),
                directions: dirs,
            });
        }
    }
    for (handle, arc) in sketch.arcs.iter() {
        let [icx, icy, _ir, _is, _ie] = arc.param_indices();
        let free_cx = icx < n && !col_has_entry[icx];
        let free_cy = icy < n && !col_has_entry[icy];
        if free_cx || free_cy {
            let mut dirs = Vec::new();
            if free_cx {
                dirs.push((1.0, 0.0));
            }
            if free_cy {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Arc(handle),
                directions: dirs,
            });
        }
    }
}

/// Merge [`FreedomInfo`] entries that reference the same element.
fn merge_freedoms(freedoms: &mut Vec<FreedomInfo>) {
    if freedoms.len() <= 1 {
        return;
    }
    let mut merged: Vec<FreedomInfo> = Vec::new();
    for f in freedoms.drain(..) {
        if let Some(existing) = merged.iter_mut().find(|e| e.element == f.element) {
            for dir in f.directions {
                // Only add if not a near-duplicate direction.
                let is_dup = existing.directions.iter().any(|&(dx, dy)| {
                    let dot = (dx * dir.0 + dy * dir.1).abs();
                    dot > 0.99
                });
                if !is_dup {
                    existing.directions.push(dir);
                }
            }
        } else {
            merged.push(f);
        }
    }
    *freedoms = merged;
}

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
