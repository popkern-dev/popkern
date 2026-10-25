//! Radiation by fixed point (map Stage 2, item 2).
//!
//! Grey radiation to surroundings at the ambient, q = εσ(T⁴ − T∞⁴), is
//! linearised per piece about the previous pass's temperature T₀ (the
//! ambient on the first) by Newton's tangent: q ≈ q(T₀) + 4εσT₀³(T − T₀),
//! a convection with coefficient 4εσT₀³ toward an effective ambient.
//! Passes repeat until the heat through every face changes by less than
//! [`TOLERANCE`] of the throughput; at the fixed point the solution
//! satisfies the nonlinear balance, to the discretisation. q is convex
//! and increasing in T, so after the first pass (which overshoots) the
//! iterates fall to it quadratically.

use crate::error::AnalysisError;
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::thermal_setup::{ABSOLUTE_ZERO_C, STEFAN_BOLTZMANN};

use super::{
    ConductionProblem, FaceCondition, TemperatureField, ThermalResult, solve_once, solve_prepared,
    system,
};

/// Face heat changes below this fraction of the throughput end the
/// iteration.
pub const TOLERANCE: f64 = 1.0e-6;

/// Most passes before the solve declines.
pub const MAX_PASSES: usize = 50;

/// Radiation linearised about `t0_c` (°C) to surroundings at
/// `ambient_c`: (coefficient W/(m²·K), effective ambient °C) with
/// q ≈ coefficient · (T − effective ambient) the tangent at T₀.
pub fn tangent(t0_c: f64, ambient_c: f64, emittance: f64) -> (f64, f64) {
    let (t, a) = (t0_c - ABSOLUTE_ZERO_C, ambient_c - ABSOLUTE_ZERO_C);
    let es = emittance * STEFAN_BOLTZMANN;
    let slope = 4.0 * es * t * t * t;
    let q = es * (t * t * t * t - a * a * a * a);
    (slope, t0_c - q / slope)
}

/// Solves `problem`, iterating radiating faces to their fixed point; a
/// problem without them takes one pass.
pub(super) fn solve_radiating(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    options: &AnalysisOptions,
    prescribed: Option<&TemperatureField>,
) -> Result<ThermalResult, AnalysisError> {
    let radiating = problem
        .faces
        .iter()
        .any(|f| matches!(f, FaceCondition::Radiating { .. }));
    if !radiating {
        return solve_once(model, problem, options, prescribed, None);
    }
    // The geometry is prepared once; each pass reassembles on it about
    // the previous pass's field (nothing else depends on it) and tests
    // the same heat balance the result reports.
    let mut prepared = system::prepare(model, problem, options, prescribed, None)?;
    let mut x = solve_prepared(&mut prepared, options)?;
    let mut heat = system::prepared_balance(model, problem, &prepared, &x, (prescribed, None));
    let mut history = Vec::new();
    for pass in 2..=MAX_PASSES {
        let field = system::field(&prepared, &x);
        system::assemble(&mut prepared, prescribed, Some(&field))?;
        x = solve_prepared(&mut prepared, options)?;
        let next =
            system::prepared_balance(model, problem, &prepared, &x, (prescribed, Some(&field)));
        let change = next
            .face_out
            .iter()
            .zip(&heat.face_out)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let scale = next.throughput().max(f64::MIN_POSITIVE);
        history.push(change / scale);
        heat = next;
        if change <= TOLERANCE * scale {
            let mut result =
                system::finish(model, problem, prepared, &x, (prescribed, Some(&field)))?;
            result.report.radiation_passes = pass;
            return Ok(result);
        }
    }
    Err(AnalysisError::SolverNotConverged { history })
}

#[cfg(test)]
#[path = "radiation_tests.rs"]
mod tests;
