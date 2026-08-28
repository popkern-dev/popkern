//! Newton-Raphson solver with backtracking line search and
//! Levenberg-Marquardt fallback.
//!
//! The solver assembles the equation system from all driving constraints,
//! builds a dense Jacobian, and iteratively updates the sketch parameter
//! vector until the residual norm drops below a tolerance.
//!
//! When the constraint graph decomposes into multiple independent connected
//! components, each sub-problem is solved separately to reduce system size.

#[path = "solver_component.rs"]
mod component;

use crate::Sketch;
use crate::analysis;
use crate::decompose;
use crate::error::{ConstraintError, SolveResult, SolveStatus};
use crate::jacobian::assemble;
use crate::rigid_set;
use faer::Mat;
use faer::prelude::{SpSolver, SpSolverLstsq};

/// Options controlling the Newton-Raphson solver.
#[derive(Debug, Clone)]
pub struct SolverOptions {
    /// Convergence tolerance on the L2 residual norm.
    pub tolerance: f64,
    /// Maximum number of Newton iterations.
    pub max_iterations: usize,
    /// Maximum number of backtracking line-search halvings.
    pub line_search_max: usize,
    /// Initial Levenberg-Marquardt damping parameter.
    pub lm_lambda_init: f64,
    /// Maximum Levenberg-Marquardt damping parameter before giving up.
    pub lm_lambda_max: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-10,
            max_iterations: 50,
            line_search_max: 10,
            lm_lambda_init: 1e-3,
            lm_lambda_max: 1e6,
        }
    }
}

/// Solve the constraint system, updating `sketch.params` in place.
///
/// When the constraint graph decomposes into multiple independent
/// components, each is solved separately. A single-component sketch
/// takes the fast path with zero decomposition overhead.
///
/// Returns the solve result with status, iteration count, and final
/// residual norm.
pub(crate) fn solve(
    sketch: &mut Sketch,
    options: &SolverOptions,
) -> Result<SolveResult, ConstraintError> {
    // Evaluate dimension expressions before the Newton loop so the solver
    // sees up-to-date constant values.
    crate::expression::evaluate_all(sketch).map_err(ConstraintError::ExpressionError)?;

    let components = decompose::find_components(sketch);

    // Fast path: 0 or 1 component means no decomposition overhead.
    if components.len() <= 1 {
        return solve_single(sketch, options);
    }

    // Multi-component path: solve each sub-problem independently.
    let total_iterations = {
        let mut max_iter = 0u32;
        for comp in &components {
            let result = component::solve_component(sketch, comp, options);
            if result.iterations > max_iter {
                max_iter = result.iterations;
            }
        }
        max_iter
    };

    // Run post-solve analysis on the full sketch after all components.
    let final_sys = assemble(sketch);
    let final_norm = residual_norm(&final_sys.residuals);
    let converged = final_norm < options.tolerance || final_sys.num_equations == 0;

    finalize(sketch, options, total_iterations, converged)
}

/// Solve the full (monolithic) constraint system.
///
/// This is the original solver path: assemble the complete system and
/// run the Newton-Raphson loop on all parameters at once.
fn solve_single(
    sketch: &mut Sketch,
    options: &SolverOptions,
) -> Result<SolveResult, ConstraintError> {
    let (iterations, converged) = newton_loop(sketch, options);
    finalize(sketch, options, iterations, converged)
}

/// Run the Newton-Raphson iteration loop (monolithic system).
///
/// Returns `(iterations, converged)`.
fn newton_loop(sketch: &mut Sketch, options: &SolverOptions) -> (u32, bool) {
    for iteration in 0..options.max_iterations {
        let sys = assemble(sketch);

        if sys.num_equations == 0 {
            return (iteration as u32, true);
        }

        let norm = residual_norm(&sys.residuals);
        if norm < options.tolerance {
            return (iteration as u32, true);
        }

        let m = sys.num_equations;
        let n = sys.num_params;

        // Build dense Jacobian matrix and RHS from triplets.
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
            let qr = j_mat.col_piv_qr();
            let dx_mat = qr.solve_lstsq(&rhs);
            (0..n).map(|i| dx_mat.read(i, 0)).collect()
        } else {
            solve_underdetermined(&j_mat, &rhs, n)
        };

        // If the step is tiny, we are at a (possibly non-zero) minimum.
        let step_norm = dx.iter().map(|d| d * d).sum::<f64>().sqrt();
        if step_norm < options.tolerance * 1e-2 {
            return (iteration as u32, true);
        }

        // Backtracking line search.
        let saved_params = sketch.params.clone();
        let mut step_accepted = false;

        let mut alpha = 1.0;
        for _ in 0..options.line_search_max {
            for (i, p) in sketch.params.iter_mut().enumerate() {
                *p = saved_params[i] + alpha * dx[i];
            }
            rigid_set::sync_rigid_sets(sketch);

            let new_sys = assemble(sketch);
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

        // Line search failed -- try Levenberg-Marquardt.
        sketch.params.copy_from_slice(&saved_params);

        let lm_accepted = try_levenberg_marquardt(sketch, &j_mat, &sys.residuals, norm, options);

        if !lm_accepted {
            sketch.params.copy_from_slice(&saved_params);
            return ((iteration + 1) as u32, false);
        }
    }

    // Exhausted iterations -- check if we converged on the last step.
    let final_norm = residual_norm(&assemble(sketch).residuals);
    (
        options.max_iterations as u32,
        final_norm < options.tolerance,
    )
}

/// Post-solve analysis: DOF analysis + constraint status.
fn finalize(
    sketch: &mut Sketch,
    options: &SolverOptions,
    iterations: u32,
    converged: bool,
) -> Result<SolveResult, ConstraintError> {
    let final_norm = residual_norm(&assemble(sketch).residuals);

    // Run constraint status detection (sets per-constraint Satisfied/Violated).
    let conflicting = analysis::analyze_constraints(sketch, options.tolerance);

    if !converged && final_norm >= options.tolerance {
        return Ok(SolveResult {
            status: SolveStatus::Failed(format!(
                "did not converge after {} iterations (residual {:.2e})",
                iterations, final_norm
            )),
            iterations,
            residual_norm: final_norm,
        });
    }

    // Check for over-constrained system.
    if !conflicting.is_empty() {
        return Ok(SolveResult {
            status: SolveStatus::OverConstrained(conflicting),
            iterations,
            residual_norm: final_norm,
        });
    }

    // DOF analysis.
    let freedoms = analysis::analyze_dof(sketch);

    let status = if freedoms.is_empty() {
        SolveStatus::FullyConstrained
    } else {
        SolveStatus::UnderConstrained(freedoms)
    };

    Ok(SolveResult {
        status,
        iterations,
        residual_norm: final_norm,
    })
}

/// Minimum-norm solution for underdetermined systems (m < n):
///   dx = J^T (J J^T)^{-1} (-F)
///
/// This finds the smallest step that satisfies the linearised equations.
pub(super) fn solve_underdetermined(j_mat: &Mat<f64>, neg_f: &Mat<f64>, n: usize) -> Vec<f64> {
    let jt = j_mat.transpose().to_owned();
    let jjt = j_mat * &jt; // m x m
    let lu = jjt.partial_piv_lu();
    let y = lu.solve(neg_f); // m x 1
    let dx_mat = jt * &y; // n x 1
    (0..n).map(|i| dx_mat.read(i, 0)).collect()
}

/// Attempt a Levenberg-Marquardt step: solve (J^T J + lambda I) dx = -J^T F
/// with increasing lambda until the residual decreases.
///
/// Returns `true` if a step was accepted (and sketch.params is updated).
fn try_levenberg_marquardt(
    sketch: &mut Sketch,
    j_mat: &Mat<f64>,
    residuals: &[f64],
    current_norm: f64,
    options: &SolverOptions,
) -> bool {
    let m = j_mat.nrows();
    let n = j_mat.ncols();
    let saved_params = sketch.params.clone();

    let jt = j_mat.transpose().to_owned();
    let jtj = &jt * j_mat;

    let mut neg_f = Mat::<f64>::zeros(m, 1);
    for (i, &r) in residuals.iter().enumerate() {
        neg_f.write(i, 0, -r);
    }
    let jtf = jt * &neg_f;

    let mut lambda = options.lm_lambda_init;
    while lambda <= options.lm_lambda_max {
        let mut damped = jtj.to_owned();
        for i in 0..n {
            let existing = damped.read(i, i);
            damped.write(i, i, existing + lambda);
        }

        let lu = damped.partial_piv_lu();
        let dx_mat = lu.solve(&jtf);
        let dx: Vec<f64> = (0..n).map(|i| dx_mat.read(i, 0)).collect();

        for (i, p) in sketch.params.iter_mut().enumerate() {
            *p = saved_params[i] + dx[i];
        }
        rigid_set::sync_rigid_sets(sketch);

        let new_norm = residual_norm(&assemble(sketch).residuals);
        if new_norm < current_norm {
            return true;
        }

        lambda *= 10.0;
    }

    sketch.params.copy_from_slice(&saved_params);
    false
}

/// L2 norm of a residual vector.
fn residual_norm(residuals: &[f64]) -> f64 {
    residuals.iter().map(|r| r * r).sum::<f64>().sqrt()
}
