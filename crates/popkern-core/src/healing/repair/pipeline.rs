//! Multi-step repair pipeline with transactional rollback.
//!
//! Orchestrates gap repair, tolerance enforcement, and surface
//! simplification across multiple iterations.  If a repair pass makes
//! the model worse (more defects than before), the BRep is restored
//! from a backup and a `RepairRegression` error is returned.

use crate::error::{CoreError, HealingError};
use crate::healing::defects;
use crate::healing::options::{AnalysisOptions, HealOptions, HealStrategy};
use crate::healing::repair::complexity::{self, ComplexityReductionResult};
use crate::healing::repair::degeneracy::{self, DegeneracyRepairResult};
use crate::healing::repair::gap_repair::{self, GapRepairResult, GapRepairStrategy};
use crate::healing::repair::merge_faces::{self, FaceMergeResult};
use crate::healing::repair::rederive::{self, RederiveResult};
use crate::healing::repair::simplification::{self, SimplificationResult};
use crate::healing::repair::stitch::{self, StitchResult};
use crate::healing::repair::tolerance_mgmt::{self, ToleranceAdjustment};
use crate::topology::BRep;

/// A single step in the repair pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RepairStep {
    /// Close edge gaps using the given strategy.
    GapRepair(GapRepairStrategy),
    /// Enforce the vertex >= edge >= face tolerance hierarchy.
    ToleranceEnforcement,
    /// Replace over-parameterized NURBS with analytical surfaces.
    SurfaceSimplification,
    /// Remove zero-length edges and collapse small faces.
    DegeneracyRepair,
    /// Detect and merge adjacent coplanar faces.
    CoplanarFaceMerge,
    /// Identify and merge matching free boundary edges.
    BoundaryStitch,
    /// Remove unnecessary NURBS knots from edge curves.
    ComplexityReduction,
    /// Re-derive edge curves from adjacent surface intersections.
    Rederive,
}

/// Aggregate result of running the repair pipeline.
#[derive(Debug, Clone)]
pub struct HealResult {
    /// Number of defects before the first repair pass.
    pub defects_before: usize,
    /// Number of defects after the final repair pass.
    pub defects_after: usize,
    /// All gap repairs performed across all iterations.
    pub gap_repairs: Vec<GapRepairResult>,
    /// All tolerance adjustments made across all iterations.
    pub tolerance_adjustments: Vec<ToleranceAdjustment>,
    /// All surface simplifications performed across all iterations.
    pub simplifications: Vec<SimplificationResult>,
    /// All degeneracy repairs performed across all iterations.
    pub degeneracy_repairs: Vec<DegeneracyRepairResult>,
    /// All coplanar face merge results across all iterations.
    pub face_merges: Vec<FaceMergeResult>,
    /// All boundary stitch results across all iterations.
    pub boundary_stitches: Vec<StitchResult>,
    /// All complexity reduction results across all iterations.
    pub complexity_reductions: Vec<ComplexityReductionResult>,
    /// All edge re-derivation results across all iterations.
    pub rederive_results: Vec<RederiveResult>,
    /// Number of iterations the pipeline ran.
    pub iterations: usize,
}

/// Determine the ordered list of repair steps for a given strategy.
fn steps_for_strategy(strategy: HealStrategy) -> Vec<RepairStep> {
    match strategy {
        HealStrategy::Conservative => vec![
            RepairStep::ToleranceEnforcement,
            RepairStep::GapRepair(GapRepairStrategy::ToleranceAbsorption),
        ],
        HealStrategy::Standard => vec![
            RepairStep::ToleranceEnforcement,
            RepairStep::DegeneracyRepair,
            RepairStep::BoundaryStitch,
            RepairStep::GapRepair(GapRepairStrategy::GeometricAdjustment),
            RepairStep::SurfaceSimplification,
            RepairStep::CoplanarFaceMerge,
        ],
        HealStrategy::Aggressive => vec![
            RepairStep::ToleranceEnforcement,
            RepairStep::DegeneracyRepair,
            RepairStep::BoundaryStitch,
            RepairStep::GapRepair(GapRepairStrategy::GeometricAdjustment),
            RepairStep::Rederive,
            RepairStep::SurfaceSimplification,
            RepairStep::ComplexityReduction,
            RepairStep::CoplanarFaceMerge,
        ],
    }
}

/// Run the full repair pipeline on a BRep.
///
/// Takes a backup of the model before starting.  For each iteration
/// (up to `options.max_iterations`), runs the steps determined by the
/// healing strategy.  After each iteration, re-analyzes the model; if
/// the defect count did not improve, iteration stops early.
///
/// If the final defect count is worse than the initial count, the BRep
/// is restored from backup and a [`HealingError::RepairRegression`]
/// error is returned.
///
/// # Errors
///
/// - `CoreError::Healing(HealingError::RepairRegression)` if the
///   repair made the model worse.
pub fn run_pipeline(brep: &mut BRep, options: &HealOptions) -> Result<HealResult, CoreError> {
    let backup = brep.clone();
    let steps = steps_for_strategy(options.strategy);

    let analysis_opts = AnalysisOptions::new(options.tolerance);
    let initial_report = defects::run_all_checks(brep, &analysis_opts);
    let defects_before = initial_report.len();

    let mut all_gap_repairs: Vec<GapRepairResult> = Vec::new();
    let mut all_tol_adjustments: Vec<ToleranceAdjustment> = Vec::new();
    let mut all_simplifications: Vec<SimplificationResult> = Vec::new();
    let mut all_degeneracy_repairs: Vec<DegeneracyRepairResult> = Vec::new();
    let mut all_face_merges: Vec<FaceMergeResult> = Vec::new();
    let mut all_boundary_stitches: Vec<StitchResult> = Vec::new();
    let mut all_complexity_reductions: Vec<ComplexityReductionResult> = Vec::new();
    let mut all_rederive_results: Vec<RederiveResult> = Vec::new();
    let mut current_defects = defects_before;
    let mut iterations = 0;

    for _ in 0..options.max_iterations {
        iterations += 1;
        let defects_before_iter = current_defects;

        for step in &steps {
            match step {
                RepairStep::ToleranceEnforcement => {
                    let adjs =
                        tolerance_mgmt::auto_adjust_tolerances(brep, &options.tolerance_bounds);
                    all_tol_adjustments.extend(adjs);
                }
                RepairStep::GapRepair(strategy) => {
                    let repairs = gap_repair::repair_all_gaps(brep, *strategy, &options.tolerance);
                    all_gap_repairs.extend(repairs);
                }
                RepairStep::SurfaceSimplification => {
                    let simps =
                        simplification::simplify_all_faces(brep, options.tolerance.linear());
                    all_simplifications.extend(simps);
                }
                RepairStep::DegeneracyRepair => {
                    let result = degeneracy::repair_degeneracies(brep, &options.tolerance);
                    all_degeneracy_repairs.push(result);
                }
                RepairStep::CoplanarFaceMerge => {
                    let result = merge_faces::merge_coplanar_faces(brep, &options.tolerance);
                    all_face_merges.push(result);
                }
                RepairStep::BoundaryStitch => {
                    let result = stitch::stitch_boundaries(
                        brep,
                        options.max_displacement * 2.0,
                        options.max_displacement,
                        &options.tolerance,
                    );
                    all_boundary_stitches.push(result);
                }
                RepairStep::ComplexityReduction => {
                    let result = complexity::reduce_all_curves(brep, options.tolerance.linear());
                    all_complexity_reductions.push(result);
                }
                RepairStep::Rederive => {
                    let results = rederive::rederive_all_edges(brep, &options.tolerance);
                    all_rederive_results.extend(results);
                }
            }
        }

        // Re-analyze after this iteration.
        let report = defects::run_all_checks(brep, &analysis_opts);
        current_defects = report.len();

        // If no improvement, stop iterating.
        if current_defects >= defects_before_iter {
            break;
        }
    }

    // If we made the model worse, roll back.
    if current_defects > defects_before {
        *brep = backup;
        return Err(HealingError::RepairRegression {
            before: defects_before,
            after: current_defects,
        }
        .into());
    }

    Ok(HealResult {
        defects_before,
        defects_after: current_defects,
        gap_repairs: all_gap_repairs,
        tolerance_adjustments: all_tol_adjustments,
        simplifications: all_simplifications,
        degeneracy_repairs: all_degeneracy_repairs,
        face_merges: all_face_merges,
        boundary_stitches: all_boundary_stitches,
        complexity_reductions: all_complexity_reductions,
        rederive_results: all_rederive_results,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::options::{HealOptions, HealStrategy};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn pipeline_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let result = run_pipeline(&mut brep, &options).unwrap();
        assert_eq!(result.defects_before, 0);
        assert_eq!(result.defects_after, 0);
        assert!(result.gap_repairs.is_empty());
        assert!(result.simplifications.is_empty());
    }

    #[test]
    fn conservative_strategy_uses_tolerance_absorption() {
        let steps = steps_for_strategy(HealStrategy::Conservative);
        assert!(steps.contains(&RepairStep::ToleranceEnforcement));
        assert!(steps.contains(&RepairStep::GapRepair(
            GapRepairStrategy::ToleranceAbsorption
        )));
        assert!(!steps.contains(&RepairStep::SurfaceSimplification));
    }

    #[test]
    fn standard_strategy_includes_simplification() {
        let steps = steps_for_strategy(HealStrategy::Standard);
        assert!(steps.contains(&RepairStep::ToleranceEnforcement));
        assert!(steps.contains(&RepairStep::GapRepair(
            GapRepairStrategy::GeometricAdjustment
        )));
        assert!(steps.contains(&RepairStep::SurfaceSimplification));
    }

    #[test]
    fn aggressive_strategy_includes_simplification() {
        let steps = steps_for_strategy(HealStrategy::Aggressive);
        assert!(steps.contains(&RepairStep::SurfaceSimplification));
        assert!(steps.contains(&RepairStep::Rederive));
    }

    #[test]
    fn pipeline_on_empty_brep_succeeds() {
        let mut brep = BRep::new(Tolerance::default());
        let options = HealOptions::new(Tolerance::default());
        let result = run_pipeline(&mut brep, &options).unwrap();
        assert_eq!(result.defects_before, 0);
        assert_eq!(result.defects_after, 0);
        assert_eq!(result.iterations, 1);
    }

    #[test]
    fn pipeline_conservative_on_valid_cube() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions {
            strategy: HealStrategy::Conservative,
            ..HealOptions::new(Tolerance::default())
        };
        let result = run_pipeline(&mut brep, &options).unwrap();
        assert_eq!(result.defects_before, 0);
        assert_eq!(result.defects_after, 0);
    }

    #[test]
    fn heal_result_iterations_at_least_one() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let result = run_pipeline(&mut brep, &options).unwrap();
        assert!(
            result.iterations >= 1,
            "pipeline should run at least one iteration"
        );
    }

    #[test]
    fn pipeline_max_iterations_respected() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions {
            max_iterations: 1,
            ..HealOptions::new(Tolerance::default())
        };
        let result = run_pipeline(&mut brep, &options).unwrap();
        assert!(result.iterations <= 1, "should not exceed max_iterations");
    }
}
