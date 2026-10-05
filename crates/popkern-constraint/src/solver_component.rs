//! Newton-Raphson solver loop for a single connected component.
//!
//! Extracted from `solver.rs` to stay within the 500-line file limit.

use super::SolverOptions;
use crate::Sketch;
use crate::decompose;
use crate::error::{SolveResult, SolveStatus};
use crate::jacobian::assemble_component;
use crate::rigid_set;
use faer::Mat;

/// Solve a single connected component of the constraint graph.
///
/// Assembles only the component's constraints with remapped indices,
/// runs Newton-Raphson on the local parameter subset, and writes
/// results back to the global parameter vector.
pub(crate) fn solve_component(
    sketch: &mut Sketch,
    component: &decompose::SubProblem,
    options: &SolverOptions,
) -> SolveResult {
    let (iterations, converged) = newton_loop_component(sketch, component, options);
    let norm = if converged {
        0.0
    } else {
        let sys = assemble_component(sketch, component);
        residual_norm(&sys.residuals)
    };
    SolveResult {
        status: if converged {
            SolveStatus::FullyConstrained
        } else {
            SolveStatus::Failed(format!(
                "component did not converge after {} iterations",
                iterations
            ))
        },
        iterations,
        residual_norm: norm,
    }
}

/// Newton-Raphson iteration loop for a single component.
///
/// Returns `(iterations, converged)`.
fn newton_loop_component(
    sketch: &mut Sketch,
    component: &decompose::SubProblem,
    options: &SolverOptions,
) -> (u32, bool) {
    let n_local = component.global_indices.len();

    for iteration in 0..options.max_iterations {
        let sys = assemble_component(sketch, component);

        if sys.num_equations == 0 {
            return (iteration as u32, true);
        }

        let norm = residual_norm(&sys.residuals);
        if norm < options.tolerance {
            return (iteration as u32, true);
        }

        let m = sys.num_equations;
        let n = sys.num_params;

        let mut j_mat = Mat::<f64>::zeros(m, n);
        for &(row, col, val) in &sys.jacobian_triplets {
            let existing = j_mat.read(row, col);
            j_mat.write(row, col, existing + val);
        }

        let mut rhs = Mat::<f64>::zeros(m, 1);
        for (i, &r) in sys.residuals.iter().enumerate() {
            rhs.write(i, 0, -r);
        }

        let dx: Vec<f64> = if m >= n {
            let dx_mat = crate::dense_solve::lstsq(&j_mat, &rhs);
            (0..n).map(|i| dx_mat.read(i, 0)).collect()
        } else {
            super::solve_underdetermined(&j_mat, &rhs, n)
        };

        let step_norm = dx.iter().map(|d| d * d).sum::<f64>().sqrt();
        if step_norm < options.tolerance * 1e-2 {
            return (iteration as u32, true);
        }

        // Save local params for backtracking.
        let saved_local: Vec<f64> = component
            .global_indices
            .iter()
            .map(|&gi| sketch.params[gi])
            .collect();

        let mut step_accepted = false;
        let mut alpha = 1.0;
        for _ in 0..options.line_search_max {
            for (local_i, &global_i) in component.global_indices.iter().enumerate() {
                sketch.params[global_i] = saved_local[local_i] + alpha * dx[local_i];
            }
            rigid_set::sync_rigid_sets(sketch);

            let new_sys = assemble_component(sketch, component);
            let new_norm = residual_norm(&new_sys.residuals);

            if new_norm < norm {
                step_accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if step_accepted {
            continue;
        }

        // Restore and try LM.
        for (local_i, &global_i) in component.global_indices.iter().enumerate() {
            sketch.params[global_i] = saved_local[local_i];
        }

        let lm_accepted = try_levenberg_marquardt_component(
            sketch,
            component,
            &j_mat,
            &sys.residuals,
            norm,
            n_local,
            options,
        );

        if !lm_accepted {
            for (local_i, &global_i) in component.global_indices.iter().enumerate() {
                sketch.params[global_i] = saved_local[local_i];
            }
            return ((iteration + 1) as u32, false);
        }
    }

    let final_norm = residual_norm(&assemble_component(sketch, component).residuals);
    (
        options.max_iterations as u32,
        final_norm < options.tolerance,
    )
}

/// Levenberg-Marquardt fallback for a component.
fn try_levenberg_marquardt_component(
    sketch: &mut Sketch,
    component: &decompose::SubProblem,
    j_mat: &Mat<f64>,
    residuals: &[f64],
    current_norm: f64,
    n: usize,
    options: &SolverOptions,
) -> bool {
    let m = j_mat.nrows();
    let saved_local: Vec<f64> = component
        .global_indices
        .iter()
        .map(|&gi| sketch.params[gi])
        .collect();

    let jt = j_mat.transpose().to_owned();
    let jtj = crate::dense_solve::mul(&jt, j_mat);

    let mut neg_f = Mat::<f64>::zeros(m, 1);
    for (i, &r) in residuals.iter().enumerate() {
        neg_f.write(i, 0, -r);
    }
    let jtf = crate::dense_solve::mul(&jt, &neg_f);

    let mut lambda = options.lm_lambda_init;
    while lambda <= options.lm_lambda_max {
        let mut damped = jtj.to_owned();
        for i in 0..n {
            let existing = damped.read(i, i);
            damped.write(i, i, existing + lambda);
        }

        let dx_mat = crate::dense_solve::lu(&damped, &jtf);
        let dx: Vec<f64> = (0..n).map(|i| dx_mat.read(i, 0)).collect();

        for (local_i, &global_i) in component.global_indices.iter().enumerate() {
            sketch.params[global_i] = saved_local[local_i] + dx[local_i];
        }
        rigid_set::sync_rigid_sets(sketch);

        let new_norm = residual_norm(&assemble_component(sketch, component).residuals);
        if new_norm < current_norm {
            return true;
        }

        lambda *= 10.0;
    }

    for (local_i, &global_i) in component.global_indices.iter().enumerate() {
        sketch.params[global_i] = saved_local[local_i];
    }
    false
}

/// L2 norm of a residual vector.
fn residual_norm(residuals: &[f64]) -> f64 {
    residuals.iter().map(|r| r * r).sum::<f64>().sqrt()
}
