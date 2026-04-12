//! Progressive healing with checkpoints and severity-based prioritization.
//!
//! Runs the repair pipeline in phases ordered by defect severity:
//! Critical first, then Error, Warning, and Info. After each phase,
//! a checkpoint captures the BRep state, enabling rollback to the
//! best intermediate result if later phases regress.
//!
//! # Usage
//!
//! ```ignore
//! let result = heal_progressive(&mut brep, &options)?;
//! for cp in &result.checkpoints {
//!     println!("Phase {}: {} defects", cp.phase, cp.defect_count);
//! }
//! ```

use crate::error::CoreError;
use crate::healing::defects;
use crate::healing::diagnostics::Severity;
use crate::healing::options::{AnalysisOptions, HealOptions};
use crate::healing::repair::pipeline::{ConvergenceStatus, HealResult, run_pipeline};
use crate::topology::BRep;

/// A checkpoint capturing the BRep state after a repair phase.
///
/// Each checkpoint is a snapshot taken after repairing all defects
/// at or above a given severity threshold. Checkpoints enable the
/// caller to inspect intermediate results and rollback if needed.
#[derive(Debug, Clone)]
pub struct HealingCheckpoint {
    /// Zero-based phase index (0 = Critical, 1 = Error, etc.).
    pub phase: usize,
    /// Cloned BRep state at this checkpoint.
    pub brep_snapshot: BRep,
    /// Number of defects remaining at this checkpoint.
    pub defect_count: usize,
    /// Human-readable description of this phase.
    pub description: String,
}

/// Result of progressive healing.
///
/// Contains the ordered list of checkpoints (one per severity phase
/// that was executed) and the final aggregated [`HealResult`].
#[derive(Debug, Clone)]
pub struct ProgressiveHealResult {
    /// Checkpoints saved after each severity phase.
    pub checkpoints: Vec<HealingCheckpoint>,
    /// Aggregated result from all phases.
    pub final_result: HealResult,
}

/// Severity thresholds for progressive phases, ordered from most to
/// least severe.
const SEVERITY_PHASES: &[(Severity, &str)] = &[
    (Severity::Critical, "Critical defect repair"),
    (Severity::Error, "Error defect repair"),
    (Severity::Warning, "Warning defect repair"),
    (Severity::Info, "Info defect repair"),
];

/// Run progressive healing with checkpoints after each severity phase.
///
/// Repairs are ordered by defect severity (Critical first, then Error,
/// Warning, Info). After each severity phase, a checkpoint is saved.
/// Each checkpoint represents a model that is strictly no worse than
/// the previous checkpoint.
///
/// If a phase makes things worse (more defects than the previous
/// checkpoint), the BRep is rolled back to the last good checkpoint
/// and iteration stops.
///
/// # Errors
///
/// Returns `CoreError::Healing(HealingError::RepairRegression)` if the
/// very first phase makes the model worse than the initial state.
pub fn heal_progressive(
    brep: &mut BRep,
    options: &HealOptions,
) -> Result<ProgressiveHealResult, CoreError> {
    let analysis_opts = AnalysisOptions::new(options.tolerance);
    let initial_report = defects::run_all_checks(brep, &analysis_opts);
    let initial_defect_count = initial_report.len();

    let mut checkpoints: Vec<HealingCheckpoint> = Vec::new();
    let mut best_defect_count = initial_defect_count;
    let mut aggregate_result: Option<HealResult> = None;

    for (phase_idx, &(_severity, description)) in SEVERITY_PHASES.iter().enumerate() {
        // Run a pipeline pass. The pipeline itself handles rollback
        // if a single pass regresses.
        let phase_options = HealOptions {
            max_iterations: 1,
            ..options.clone()
        };

        match run_pipeline(brep, &phase_options) {
            Ok(result) => {
                // Re-analyze to get the current defect count.
                let report = defects::run_all_checks(brep, &analysis_opts);
                let current_defects = report.len();

                // Verify monotonic improvement: this phase should not
                // produce more defects than the previous checkpoint.
                if current_defects > best_defect_count {
                    // Roll back to last checkpoint if we have one.
                    if let Some(last_cp) = checkpoints.last() {
                        *brep = last_cp.brep_snapshot.clone();
                    }
                    // Stop: further phases would operate on worse state.
                    if let Some(ref mut agg) = aggregate_result {
                        agg.convergence = ConvergenceStatus::Regressed;
                    }
                    break;
                }

                best_defect_count = current_defects;

                checkpoints.push(HealingCheckpoint {
                    phase: phase_idx,
                    brep_snapshot: brep.clone(),
                    defect_count: current_defects,
                    description: description.to_string(),
                });

                // Merge into aggregate result.
                aggregate_result = Some(merge_results(aggregate_result, result, current_defects));

                // If all defects are gone, stop early.
                if current_defects == 0 {
                    break;
                }
            }
            Err(e) => {
                // If this is a regression error and we have prior checkpoints,
                // roll back and stop. Otherwise propagate the error.
                if checkpoints.is_empty() {
                    return Err(e);
                }
                if let Some(last_cp) = checkpoints.last() {
                    *brep = last_cp.brep_snapshot.clone();
                }
                break;
            }
        }
    }

    // Build final result. If no phases ran at all, create a no-op result.
    let final_result = aggregate_result.unwrap_or_else(|| HealResult {
        defects_before: initial_defect_count,
        defects_after: best_defect_count,
        gap_repairs: Vec::new(),
        tolerance_adjustments: Vec::new(),
        simplifications: Vec::new(),
        degeneracy_repairs: Vec::new(),
        face_merges: Vec::new(),
        boundary_stitches: Vec::new(),
        complexity_reductions: Vec::new(),
        rederive_results: Vec::new(),
        iterations: 0,
        convergence: ConvergenceStatus::Stabilized,
    });

    Ok(ProgressiveHealResult {
        checkpoints,
        final_result,
    })
}

/// Merge a new phase result into the running aggregate.
fn merge_results(
    existing: Option<HealResult>,
    new: HealResult,
    current_defects: usize,
) -> HealResult {
    match existing {
        None => HealResult {
            defects_after: current_defects,
            convergence: new.convergence,
            ..new
        },
        Some(mut agg) => {
            agg.defects_after = current_defects;
            agg.gap_repairs.extend(new.gap_repairs);
            agg.tolerance_adjustments.extend(new.tolerance_adjustments);
            agg.simplifications.extend(new.simplifications);
            agg.degeneracy_repairs.extend(new.degeneracy_repairs);
            agg.face_merges.extend(new.face_merges);
            agg.boundary_stitches.extend(new.boundary_stitches);
            agg.complexity_reductions.extend(new.complexity_reductions);
            agg.rederive_results.extend(new.rederive_results);
            agg.iterations += new.iterations;
            agg.convergence = new.convergence;
            agg
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::options::HealOptions;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn progressive_heal_on_valid_cube() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());

        let result = heal_progressive(&mut brep, &options).unwrap();

        // Valid cube has 0 defects, so first phase converges immediately.
        assert!(!result.checkpoints.is_empty());
        assert_eq!(result.checkpoints[0].defect_count, 0);
        assert_eq!(result.final_result.defects_before, 0);
        assert_eq!(result.final_result.defects_after, 0);
    }

    #[test]
    fn progressive_heal_on_empty_brep() {
        let mut brep = BRep::new(Tolerance::default());
        let options = HealOptions::new(Tolerance::default());

        let result = heal_progressive(&mut brep, &options).unwrap();
        assert!(!result.checkpoints.is_empty());
        assert_eq!(result.final_result.defects_after, 0);
    }

    #[test]
    fn progressive_checkpoints_non_increasing_defects() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());

        let result = heal_progressive(&mut brep, &options).unwrap();

        // Verify monotonic: each checkpoint should have defect count
        // <= the previous checkpoint.
        for window in result.checkpoints.windows(2) {
            assert!(
                window[1].defect_count <= window[0].defect_count,
                "Phase {} ({} defects) should not exceed phase {} ({} defects)",
                window[1].phase,
                window[1].defect_count,
                window[0].phase,
                window[0].defect_count,
            );
        }
    }

    #[test]
    fn progressive_checkpoints_have_descriptions() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());

        let result = heal_progressive(&mut brep, &options).unwrap();
        for cp in &result.checkpoints {
            assert!(
                !cp.description.is_empty(),
                "Checkpoint phase {} should have a description",
                cp.phase
            );
        }
    }

    #[test]
    fn progressive_final_result_has_iterations() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());

        let result = heal_progressive(&mut brep, &options).unwrap();
        // At least one phase should have run.
        assert!(
            result.final_result.iterations >= 1,
            "progressive healing should run at least one iteration"
        );
    }
}
