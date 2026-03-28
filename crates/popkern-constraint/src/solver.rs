//! Newton–Raphson solver with backtracking line search and
//! Levenberg–Marquardt fallback.
//!
//! The solver assembles the equation system from all driving constraints,
//! builds a dense Jacobian, and iteratively updates the sketch parameter
//! vector until the residual norm drops below a tolerance.

use crate::error::{ConstraintError, SolveResult, SolveStatus};
use crate::jacobian::assemble;
use crate::Sketch;
use faer::prelude::{SpSolver, SpSolverLstsq};
use faer::Mat;

/// Options controlling the Newton–Raphson solver.
#[derive(Debug, Clone)]
pub struct SolverOptions {
    /// Convergence tolerance on the L2 residual norm.
    pub tolerance: f64,
    /// Maximum number of Newton iterations.
    pub max_iterations: usize,
    /// Maximum number of backtracking line-search halvings.
    pub line_search_max: usize,
    /// Initial Levenberg–Marquardt damping parameter.
    pub lm_lambda_init: f64,
    /// Maximum Levenberg–Marquardt damping parameter before giving up.
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
/// Returns the solve result with status, iteration count, and final
/// residual norm.
pub(crate) fn solve(
    sketch: &mut Sketch,
    options: &SolverOptions,
) -> Result<SolveResult, ConstraintError> {
    for iteration in 0..options.max_iterations {
        // (a) Assemble equation system.
        let sys = assemble(sketch);

        if sys.num_equations == 0 {
            return Ok(SolveResult {
                status: SolveStatus::UnderConstrained(Vec::new()),
                iterations: iteration as u32,
                residual_norm: 0.0,
            });
        }

        // (b) Check convergence.
        let norm = residual_norm(&sys.residuals);
        if norm < options.tolerance {
            let status = if sys.num_equations == sys.num_params {
                SolveStatus::FullyConstrained
            } else if sys.num_equations < sys.num_params {
                SolveStatus::UnderConstrained(Vec::new())
            } else {
                // More equations than params but still satisfied — redundant.
                SolveStatus::FullyConstrained
            };
            return Ok(SolveResult {
                status,
                iterations: iteration as u32,
                residual_norm: norm,
            });
        }

        let m = sys.num_equations;
        let n = sys.num_params;

        // (c) Build dense Jacobian matrix and RHS from triplets.
        let mut j_mat = Mat::<f64>::zeros(m, n);
        for &(row, col, val) in &sys.jacobian_triplets {
            // Accumulate in case multiple constraints write to the same entry.
            let existing = j_mat.read(row, col);
            j_mat.write(row, col, existing + val);
        }

        let mut rhs = Mat::<f64>::zeros(m, 1);
        for (i, &r) in sys.residuals.iter().enumerate() {
            rhs.write(i, 0, -r);
        }

        // (d) Solve J * dx = -F.
        //     Overdetermined (m >= n): QR least-squares.
        //     Underdetermined (m < n): minimum-norm via J^T (J J^T)^{-1} (-F).
        let dx: Vec<f64> = if m >= n {
            let qr = j_mat.col_piv_qr();
            let dx_mat = qr.solve_lstsq(&rhs);
            (0..n).map(|i| dx_mat.read(i, 0)).collect()
        } else {
            solve_underdetermined(&j_mat, &rhs, n)
        };

        // (e) Backtracking line search.
        let saved_params = sketch.params.clone();
        let mut step_accepted = false;

        let mut alpha = 1.0;
        for _ in 0..options.line_search_max {
            // Apply step: params = saved + alpha * dx
            for (i, p) in sketch.params.iter_mut().enumerate() {
                *p = saved_params[i] + alpha * dx[i];
            }

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

        // (f) Line search failed — try Levenberg–Marquardt.
        // Restore params before LM attempts.
        sketch.params.copy_from_slice(&saved_params);

        let lm_accepted = try_levenberg_marquardt(sketch, &j_mat, &sys.residuals, norm, options);

        if !lm_accepted {
            // Neither line search nor LM improved the residual.
            sketch.params.copy_from_slice(&saved_params);
            let final_norm = residual_norm(&assemble(sketch).residuals);
            return Ok(SolveResult {
                status: SolveStatus::Failed(
                    "line search and Levenberg-Marquardt both failed to reduce residual"
                        .to_string(),
                ),
                iterations: (iteration + 1) as u32,
                residual_norm: final_norm,
            });
        }
    }

    // Exhausted iterations.
    let final_sys = assemble(sketch);
    let final_norm = residual_norm(&final_sys.residuals);
    if final_norm < options.tolerance {
        Ok(SolveResult {
            status: SolveStatus::FullyConstrained,
            iterations: options.max_iterations as u32,
            residual_norm: final_norm,
        })
    } else {
        Ok(SolveResult {
            status: SolveStatus::Failed(format!(
                "did not converge after {} iterations (residual {:.2e})",
                options.max_iterations, final_norm
            )),
            iterations: options.max_iterations as u32,
            residual_norm: final_norm,
        })
    }
}

/// Minimum-norm solution for underdetermined systems (m < n):
///   dx = J^T (J J^T)^{-1} (-F)
///
/// This finds the smallest step that satisfies the linearised equations.
fn solve_underdetermined(j_mat: &Mat<f64>, neg_f: &Mat<f64>, n: usize) -> Vec<f64> {
    let jt = j_mat.transpose().to_owned();
    let jjt = j_mat * &jt; // m x m
    let lu = jjt.partial_piv_lu();
    let y = lu.solve(neg_f); // m x 1
    let dx_mat = jt * &y; // n x 1
    (0..n).map(|i| dx_mat.read(i, 0)).collect()
}

/// Attempt a Levenberg–Marquardt step: solve (J^T J + lambda I) dx = -J^T F
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

    // J^T * J
    let jt = j_mat.transpose().to_owned();
    let jtj = &jt * j_mat;

    // J^T * (-F)
    let mut neg_f = Mat::<f64>::zeros(m, 1);
    for (i, &r) in residuals.iter().enumerate() {
        neg_f.write(i, 0, -r);
    }
    let jtf = jt * &neg_f;

    let mut lambda = options.lm_lambda_init;
    while lambda <= options.lm_lambda_max {
        // (J^T J + lambda I)
        let mut damped = jtj.to_owned();
        for i in 0..n {
            let existing = damped.read(i, i);
            damped.write(i, i, existing + lambda);
        }

        // Solve via LU.
        let lu = damped.partial_piv_lu();
        let dx_mat = lu.solve(&jtf);
        let dx: Vec<f64> = (0..n).map(|i| dx_mat.read(i, 0)).collect();

        // Apply step.
        for (i, p) in sketch.params.iter_mut().enumerate() {
            *p = saved_params[i] + dx[i];
        }

        let new_norm = residual_norm(&assemble(sketch).residuals);
        if new_norm < current_norm {
            return true;
        }

        lambda *= 10.0;
    }

    // Restore params — no improvement found.
    sketch.params.copy_from_slice(&saved_params);
    false
}

/// L2 norm of a residual vector.
fn residual_norm(residuals: &[f64]) -> f64 {
    residuals.iter().map(|r| r * r).sum::<f64>().sqrt()
}
