//! The optimisation loop: evaluate, update, continue the projection,
//! until the design settles.

use serde::Serialize;

use crate::error::AnalysisError;
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::report::SolveReport;

use super::super::{ElasticProblem, internal};
use super::{DesignEntry, Evaluator, GenerativeOptions, UpdateMethod, mma, projection, update};

/// One update's record.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Iteration {
    /// Compliance of the design before the update, N·mm (the eroded
    /// design's under the robust formulation).
    pub compliance: f64,
    /// Its material volume over the part's.
    pub volume_fraction: f64,
    /// The largest change of a physical density in the update.
    pub change: f64,
    /// PCG iterations of its solve.
    pub cg_iterations: usize,
    /// The projection's sharpness at this update (0: off).
    pub beta: f64,
}

/// A finished design.
#[derive(Debug, Clone, PartialEq)]
pub struct Design {
    /// The entries, in (solid, entry) order.
    pub entries: Vec<DesignEntry>,
    /// Physical densities x̃ of the final design.
    pub densities: Vec<f64>,
    /// Each entry's density floor (the kept material).
    pub floors: Vec<f64>,
    /// Grid cell size, mm.
    pub cell_size: f64,
    /// Compliance of the final design, N·mm.
    pub compliance: f64,
    /// Its material volume over the part's.
    pub volume_fraction: f64,
    /// Every update.
    pub history: Vec<Iteration>,
    /// Whether the change fell below the tolerance.
    pub converged: bool,
    /// With a build direction: the fraction of the design's solid
    /// volume (density over ½, kept material aside) that would print
    /// unsupported by the exact layer rule.
    pub unsupported_fraction: Option<f64>,
    /// Its non-discreteness (Sigmund 2007, as Wang 2011 report it):
    /// Σ 4x̄(1 − x̄) v / Σ v, 0 for a solid-and-void design, 1 for one
    /// all at ½.
    pub grey_fraction: f64,
    /// Its thinnest members, mm: the 10% quantile by volume of the
    /// solid's (density ≥ ½) member thickness, capped at the filter's
    /// 2r (`thickness::thickness_quantile`). Measured on the crisp
    /// 24 × 8 × 4 cantilever (r = 1.5 mm, 0.5 mm cells): 0.5 mm, one
    /// cell, for the plain projection, 1.0 for the robust formulation
    /// (η = 0.3); 0 for a design without solid.
    pub member_thickness: f64,
    /// Under the robust formulation, the eroded design's compliance,
    /// N·mm (the intermediate design's is `compliance`).
    pub eroded_compliance: Option<f64>,
    /// The solve's report (the final solve's PCG fields).
    pub report: SolveReport,
}

/// Optimises `problem` on `model`, holding solid the material near the
/// faces in `keep`.
pub fn optimise(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    design: &GenerativeOptions,
    keep: &[crate::model::FaceId],
    options: &AnalysisOptions,
) -> Result<Design, AnalysisError> {
    optimise_observed(model, problem, design, keep, options, &mut |_| {
        std::ops::ControlFlow::Continue(())
    })
}

/// One update as an observer sees it: its record and the physical
/// densities it moved the design to (one per entry, the order of
/// [`Design::entries`]).
#[derive(Debug, Clone, Copy)]
pub struct Progress<'a> {
    /// The update's record.
    pub iteration: &'a Iteration,
    /// Updates so far, this one included.
    pub count: usize,
    /// The design after the update.
    pub densities: &'a [f64],
}

/// [`optimise`] reporting each update to `observer`, which may stop the
/// run: `Break` declines it as [`AnalysisError::Cancelled`]. The
/// observer reads only; the design is the same bits as without it.
pub fn optimise_observed(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    design: &GenerativeOptions,
    keep: &[crate::model::FaceId],
    options: &AnalysisOptions,
    observer: &mut dyn FnMut(Progress<'_>) -> std::ops::ControlFlow<()>,
) -> Result<Design, AnalysisError> {
    let mut evaluator = Evaluator::new(model, problem, *design, keep, options)?;
    let volumes: Vec<f64> = evaluator.entries().iter().map(|d| d.volume).collect();
    let total: f64 = volumes.iter().sum();
    let target = design.volume_fraction * total;
    let floors = evaluator.floors().to_vec();
    let kept: f64 = volumes.iter().zip(&floors).map(|(v, f)| v * f).sum();
    if kept > target {
        return Err(AnalysisError::KeptVolumeExceedsTarget { kept, target });
    }
    let mut x: Vec<f64> = floors
        .iter()
        .map(|f| f.max(design.volume_fraction))
        .collect();
    let mut history = Vec::new();
    let mut converged = false;
    let mut beta = if design.projection_beta > 0.0 {
        1.0f64.min(design.projection_beta)
    } else {
        0.0
    };
    let mut since_raise = 0;
    let mut mma = mma::Mma::new(floors.clone(), vec![1.0; floors.len()], design.move_limit);
    // MMA's objective is the compliance over the first design's.
    let mut scale = None;
    let project = |beta: f64| {
        (beta > 0.0).then_some(projection::Projection {
            beta,
            kind: design.projection_kind,
        })
    };
    evaluator.set_projection(project(beta));
    evaluator.set_robust(design.robust);
    evaluator.set_tolerance(design.update_tolerance);
    let mut eval = evaluator.evaluate(&x)?;
    for _ in 0..design.max_iterations {
        let volume: f64 = eval.physical.iter().zip(&volumes).map(|(p, v)| p * v).sum();
        let next = match design.method {
            UpdateMethod::OptimalityCriteria => {
                update::optimality_criteria(
                    &x,
                    &update::UpdateInput {
                        dc: &eval.gradient,
                        dv: &eval.volume_gradient,
                        floors: &floors,
                        volumes: &volumes,
                        filter: evaluator.filter(),
                        projection: project(beta),
                        print: evaluator.print(),
                        target,
                        move_limit: design.move_limit,
                    },
                )
                .0
            }
            UpdateMethod::Mma => {
                let c0 = *scale.get_or_insert(eval.compliance.abs().max(f64::MIN_POSITIVE));
                let df0: Vec<f64> = eval.gradient.iter().map(|g| g / c0).collect();
                let df1: Vec<f64> = eval.volume_gradient.iter().map(|g| g / target).collect();
                // The volume is explicit: met exactly.
                let filter = evaluator.filter();
                let print = evaluator.print();
                let projection = project(beta);
                let exact = |z: &[f64]| {
                    update::volume(z, &floors, &volumes, filter, projection, print) / target - 1.0
                };
                mma.update(&x, &df0, (volume / target - 1.0, &df1), Some(&exact))
            }
        };
        // A non-finite design would print as its floors (f64::max drops
        // NaN) and read as a collapse.
        if next.iter().any(|v| !v.is_finite()) {
            return Err(internal("generative update", "non-finite design variable"));
        }
        let trial = evaluator.evaluate(&next)?;
        // The change of the design itself, the physical densities: under
        // a saturated projection the variables can swing where the
        // design does not move (measured: ±0.2 at constant compliance).
        let change = trial
            .physical
            .iter()
            .zip(&eval.physical)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        history.push(Iteration {
            compliance: eval.compliance,
            volume_fraction: volume / total,
            change,
            cg_iterations: eval.cg_iterations,
            beta,
        });
        if let Some(iteration) = history.last() {
            let progress = Progress {
                iteration,
                count: history.len(),
                densities: &trial.physical,
            };
            if observer(progress).is_break() {
                return Err(AnalysisError::Cancelled {
                    iterations: history.len(),
                });
            }
        }
        x = next;
        eval = trial;
        since_raise += 1;
        let settled = change < design.tolerance;
        if beta > 0.0
            && beta < design.projection_beta
            && (settled || since_raise >= design.projection_every)
        {
            beta = (2.0 * beta).min(design.projection_beta);
            since_raise = 0;
            evaluator.set_projection(project(beta));
            eval = evaluator.evaluate(&x)?;
            continue;
        }
        if settled {
            converged = true;
            break;
        }
    }
    // Under the robust formulation the result is the intermediate
    // design, solved on its own; with a looser update tolerance, the
    // final design is solved again at the analysis's.
    let eroded_compliance = design.robust.map(|_| eval.compliance);
    if design.robust.is_some() || design.update_tolerance > options.cg_rel_tol {
        evaluator.set_robust(None);
        evaluator.set_tolerance(options.cg_rel_tol);
        eval = evaluator.evaluate(&x)?;
    }
    let last = eval;
    let volume: f64 = last.physical.iter().zip(&volumes).map(|(p, v)| p * v).sum();
    let grey_fraction = if total > 0.0 {
        last.physical
            .iter()
            .zip(&volumes)
            .map(|(p, v)| 4.0 * p * (1.0 - p) * v)
            .sum::<f64>()
            / total
    } else {
        0.0
    };
    let cell_size = evaluator.report().cell_size.unwrap_or(0.0);
    let member_thickness =
        super::cells::CellField::new(evaluator.entries(), &last.physical, cell_size)
            .and_then(|c| super::thickness::thickness_quantile(&c, 2.0 * design.filter_radius, 0.1))
            .unwrap_or(0.0);
    let unsupported_fraction = evaluator
        .print()
        .map(|f| f.unsupported(&last.physical, &floors, &volumes, 0.5));
    Ok(Design {
        entries: evaluator.entries().to_vec(),
        densities: last.physical,
        floors,
        cell_size,
        compliance: last.compliance,
        volume_fraction: volume / total,
        history,
        converged,
        unsupported_fraction,
        grey_fraction,
        eroded_compliance,
        member_thickness,
        report: evaluator.report().clone(),
    })
}

#[cfg(test)]
#[path = "optimise_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "robust_tests.rs"]
mod robust_tests;
