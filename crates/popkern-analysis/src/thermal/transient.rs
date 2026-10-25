//! Transient conduction by backward Euler (map Stage 2, item 5; Hughes
//! §8.1).
//!
//! With the steady system A x = b (conduction, convection, Nitsche,
//! sources) and the capacity matrix C = ∫ ρc N_a N_b dV on the same
//! moments and AgFEM extension, each step of size Δt solves
//! (C/Δt + A) x' = b + C x/Δt, as the increment d = x' − x from
//! (C/Δt + A) d = b − A x so that the solver's relative tolerance
//! applies to the step's change rather than to the absolute
//! temperature (PH-100: a 5e-4 K rise over 20 °C otherwise carried
//! 4e4 times the tolerance). The step is fixed, so one multigrid
//! hierarchy serves every step. The part starts at a uniform
//! temperature (the extension reproduces constants exactly).
//!
//! Radiating faces are linearised each step about the step's starting
//! temperature (semi-implicit, first order in Δt like backward Euler
//! itself), which reassembles the operator and its hierarchy per step.
//!
//! Energy: the constant field 1 tested against the step gives
//! 1ᵀC(x' − x)/Δt = heat in − heat out, so the stored energy and the
//! face heat balance each step to the solver's tolerance; the per-step
//! imbalance is reported.

use crate::error::AnalysisError;
use crate::linalg::sparse::Csr;
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::solve::pcg::pcg;

use super::assemble::{SolidInput, assemble_capacity};
use super::system::{assemble, field, finish, hierarchy, prepare, state};
use super::{ConductionProblem, FaceCondition, ThermalResult, internal};

/// A transient conduction problem.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientProblem {
    /// Conductivity, conditions and sources, as for a steady solve.
    pub conduction: ConductionProblem,
    /// Density, kg/m³.
    pub density: f64,
    /// Specific heat capacity, J/(kg·K).
    pub specific_heat: f64,
    /// Uniform initial temperature, °C.
    pub initial_c: f64,
    /// Time step, s.
    pub step_s: f64,
    /// Number of steps.
    pub steps: usize,
}

/// The state after one step.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct TransientStep {
    /// Time, s.
    pub time_s: f64,
    /// Highest vertex temperature, °C.
    pub max_temperature_c: f64,
    /// Energy stored since the start, J.
    pub stored_j: f64,
    /// Net heat out through the faces during the step, W.
    pub heat_out_w: f64,
    /// |stored rate + heat out − heat in| over the largest of the three.
    pub imbalance: f64,
}

/// A solved transient.
#[derive(Debug, Clone, PartialEq)]
pub struct TransientResult {
    /// Every step.
    pub steps: Vec<TransientStep>,
    /// The final state, as a steady result (field, vertex values, heat
    /// balance of the final step's faces).
    pub last: ThermalResult,
}

/// Solves `problem` on `model` with `options`.
pub fn solve_transient(
    model: &AnalysisModel,
    problem: &TransientProblem,
    options: &AnalysisOptions,
) -> Result<TransientResult, AnalysisError> {
    let conduction = &problem.conduction;
    if conduction
        .faces
        .iter()
        .any(|f| matches!(f, FaceCondition::Prescribed))
    {
        return Err(AnalysisError::Unsupported {
            what: "prescribed faces in a transient solve",
        });
    }
    let radiating = conduction
        .faces
        .iter()
        .any(|f| matches!(f, FaceCondition::Radiating { .. }));
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !(positive(problem.density)
        && positive(problem.specific_heat)
        && positive(problem.step_s)
        && problem.initial_c.is_finite())
        || problem.steps == 0
    {
        return Err(AnalysisError::Unsupported {
            what: "a transient needs positive density, specific heat, step and step count",
        });
    }
    let mut prepared = prepare(model, conduction, options, None, None)?;
    // ρc in J/(mm³·K): kg/m³ · J/(kg·K) · 1e-9 m³/mm³.
    let rho_c = 1e-9 * problem.density * problem.specific_heat;
    let mut blocks = Vec::new();
    for (clipped, moments, aggregation) in &prepared.solids {
        let input = SolidInput {
            grid: &prepared.grid,
            clipped,
            moments,
            aggregation,
            conditions: &prepared.conditions,
            prescribed: None,
            previous: None,
            cache: None,
        };
        blocks.push(
            assemble_capacity(&prepared.tables, &input, rho_c)
                .map_err(|e| internal("capacity", e))?,
        );
    }
    let c = Csr::block_diagonal(blocks).map_err(|e| internal("capacity", e))?;
    let n = prepared.b.len();
    let dt = problem.step_s;
    // C/Δt + A, and its hierarchy.
    let operator = |prepared: &mut super::system::Prepared| -> Result<_, AnalysisError> {
        let mut triplets = Vec::new();
        for r in 0..n {
            triplets.extend(c.row(r).map(|(col, v)| (r, col, v / dt)));
            triplets.extend(prepared.a.row(r).map(|(col, v)| (r, col, v)));
        }
        let a = Csr::from_triplets(n, n, &triplets).map_err(|e| internal("transient", e))?;
        let positions = prepared.positions.clone();
        hierarchy(a, positions, options, &mut prepared.report)
    };
    let mut hierarchy_now = operator(&mut prepared)?;
    let ones = vec![1.0; n];
    let mut c_ones = vec![0.0; n];
    c.mul_vec(&ones, &mut c_ones)
        .map_err(|e| internal("transient", e))?;
    let mut x = vec![problem.initial_c; n];
    let mut steps = Vec::with_capacity(problem.steps);
    let (mut stored, mut iterations) = (0.0, 0);
    let mut ax = vec![0.0; n];
    let mut linearised = None;
    for k in 1..=problem.steps {
        // Radiation linearised about this step's start (semi-implicit:
        // first order in Δt, as backward Euler is).
        if radiating {
            let f = field(&prepared, &x);
            assemble(&mut prepared, None, Some(&f))?;
            hierarchy_now = operator(&mut prepared)?;
            linearised = Some(f);
        }
        // The increment d = x' − x: (C/Δt + A) d = b − A x, the same
        // step with the tolerance on the change, not on the absolute
        // temperature (PH-100).
        prepared
            .a
            .mul_vec(&x, &mut ax)
            .map_err(|e| internal("transient", e))?;
        let rhs: Vec<f64> = prepared.b.iter().zip(&ax).map(|(b, v)| b - v).collect();
        let solution = pcg(
            &hierarchy_now,
            &rhs,
            options.cg_rel_tol,
            options.cg_max_iters,
        )
        .map_err(|e| internal("pcg", e))?;
        if !solution.converged {
            return Err(AnalysisError::SolverNotConverged {
                history: solution.history,
            });
        }
        iterations += solution.iterations;
        // 1ᵀC d: C symmetric.
        let gained: f64 = solution.x.iter().zip(&c_ones).map(|(d, w)| w * d).sum();
        stored += gained;
        x.iter_mut().zip(&solution.x).for_each(|(t, d)| *t += d);
        let (out, heat_in, max_c) = state(model, conduction, &prepared, &x, linearised.as_ref())?;
        let rate = gained / dt;
        let scale = rate
            .abs()
            .max(out.abs())
            .max(heat_in.abs())
            .max(f64::MIN_POSITIVE);
        steps.push(TransientStep {
            time_s: k as f64 * dt,
            max_temperature_c: max_c,
            stored_j: stored,
            heat_out_w: out,
            imbalance: (rate + out - heat_in).abs() / scale,
        });
    }
    prepared.report.cg_iterations = iterations;
    let mut last = finish(model, conduction, prepared, &x, (None, linearised.as_ref()))?;
    last.transient = steps.clone();
    Ok(TransientResult { steps, last })
}

#[cfg(test)]
#[path = "transient_tests.rs"]
mod tests;
