//! Import repair pipeline — systematic repair of models imported from external formats.
//!
//! Provides a configurable sequence of repair operations tailored to the
//! characteristic defect patterns of imported models (Section 7.2).

use crate::healing::defects;
use crate::healing::options::AnalysisOptions;
use crate::healing::repair::gap_repair::{self, GapRepairStrategy};
use crate::healing::repair::simplification;
use crate::healing::repair::tolerance_mgmt;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A single import repair operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportRepairOp {
    /// Fix tolerance values that don't match the target system.
    ToleranceNormalization,
    /// Restore missing periodic surface closure topology.
    PeriodicTopologyRestore,
    /// Fix inconsistent face/edge orientations.
    OrientationFix,
    /// Close boundary gaps between adjacent faces.
    GapClosure,
    /// Simplify over-parameterized surfaces to exact analytical forms.
    SurfaceSimplification,
    /// Enforce tolerance hierarchy after all other repairs.
    ToleranceEnforcement,
}

/// Configuration for the import repair pipeline.
#[derive(Debug, Clone)]
pub struct ImportRepairConfig {
    /// Ordered list of repair operations to perform.
    pub operations: Vec<ImportRepairOp>,
    /// Maximum displacement allowed per repair.
    pub max_displacement: f64,
    /// Tolerance for geometric comparisons.
    pub tolerance: Tolerance,
}

impl ImportRepairConfig {
    /// Default import repair configuration (all operations in recommended order).
    ///
    /// The order is chosen to resolve dependency chains: orientation and
    /// topology first, then tolerance normalization, gap closure, surface
    /// simplification, and finally tolerance enforcement.
    pub fn default_config(tolerance: Tolerance) -> Self {
        Self {
            operations: vec![
                ImportRepairOp::OrientationFix,
                ImportRepairOp::PeriodicTopologyRestore,
                ImportRepairOp::ToleranceNormalization,
                ImportRepairOp::GapClosure,
                ImportRepairOp::SurfaceSimplification,
                ImportRepairOp::ToleranceEnforcement,
            ],
            max_displacement: tolerance.linear() * 10.0,
            tolerance,
        }
    }
}

/// Report from an import repair run.
#[derive(Debug, Clone)]
pub struct ImportRepairReport {
    /// Per-operation results: (operation, defects_before, defects_after).
    pub steps: Vec<(ImportRepairOp, usize, usize)>,
    /// Total defects before repair.
    pub total_defects_before: usize,
    /// Total defects after repair.
    pub total_defects_after: usize,
}

/// Count defects on the BRep using the standard analysis options.
fn count_defects(brep: &BRep, tolerance: &Tolerance) -> usize {
    let analysis_opts = AnalysisOptions::new(*tolerance);
    let report = defects::run_all_checks(brep, &analysis_opts);
    report.len()
}

/// Execute a single import repair operation on the BRep.
fn execute_op(brep: &mut BRep, op: ImportRepairOp, config: &ImportRepairConfig) {
    match op {
        ImportRepairOp::ToleranceNormalization => {
            let bounds = crate::healing::options::ToleranceBounds::default();
            tolerance_mgmt::auto_adjust_tolerances(brep, &bounds);
        }
        ImportRepairOp::GapClosure => {
            gap_repair::repair_all_gaps(
                brep,
                GapRepairStrategy::GeometricAdjustment,
                &config.tolerance,
            );
        }
        ImportRepairOp::SurfaceSimplification => {
            simplification::simplify_all_faces(brep, config.tolerance.linear());
        }
        ImportRepairOp::ToleranceEnforcement => {
            tolerance_mgmt::enforce_tolerance_hierarchy(brep);
        }
        // Stubs: orientation repair requires complex topology surgery,
        // and periodic topology restore needs Euler operators to add seam edges.
        ImportRepairOp::OrientationFix | ImportRepairOp::PeriodicTopologyRestore => {}
    }
}

/// Run the import repair pipeline on a BRep.
///
/// Executes each configured operation in order, re-analyzing between
/// steps to track progress. Returns a detailed report.
///
/// Each operation modifies the BRep in place. The report records the
/// defect count before and after each step so callers can evaluate
/// which operations had the most impact.
pub fn run_import_repair(brep: &mut BRep, config: &ImportRepairConfig) -> ImportRepairReport {
    let total_defects_before = count_defects(brep, &config.tolerance);

    let mut steps = Vec::with_capacity(config.operations.len());

    for &op in &config.operations {
        let before = count_defects(brep, &config.tolerance);
        execute_op(brep, op, config);
        let after = count_defects(brep, &config.tolerance);
        steps.push((op, before, after));
    }

    let total_defects_after = count_defects(brep, &config.tolerance);

    ImportRepairReport {
        steps,
        total_defects_before,
        total_defects_after,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn default_config_on_valid_cube_does_not_increase_defects() {
        let tol = Tolerance::default();
        let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
        let config = ImportRepairConfig::default_config(tol);
        let report = run_import_repair(&mut brep, &config);

        // The pipeline must not make things worse.
        assert!(
            report.total_defects_after <= report.total_defects_before,
            "import repair should not increase defects: before={}, after={}",
            report.total_defects_before,
            report.total_defects_after,
        );
        // Every step should also not increase defects.
        for (op, before, after) in &report.steps {
            assert!(
                *after <= *before,
                "{op:?} increased defects: {before} -> {after}"
            );
        }
    }

    #[test]
    fn empty_ops_returns_immediately() {
        let tol = Tolerance::default();
        let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
        let config = ImportRepairConfig {
            operations: vec![],
            max_displacement: tol.linear() * 10.0,
            tolerance: tol,
        };
        let report = run_import_repair(&mut brep, &config);

        assert!(
            report.steps.is_empty(),
            "empty operations should produce no steps"
        );
        assert_eq!(report.total_defects_before, report.total_defects_after);
    }

    #[test]
    fn default_config_has_correct_operation_order() {
        let tol = Tolerance::default();
        let config = ImportRepairConfig::default_config(tol);

        assert_eq!(config.operations.len(), 6);
        assert_eq!(config.operations[0], ImportRepairOp::OrientationFix);
        assert_eq!(
            config.operations[1],
            ImportRepairOp::PeriodicTopologyRestore
        );
        assert_eq!(config.operations[2], ImportRepairOp::ToleranceNormalization);
        assert_eq!(config.operations[3], ImportRepairOp::GapClosure);
        assert_eq!(config.operations[4], ImportRepairOp::SurfaceSimplification);
        assert_eq!(config.operations[5], ImportRepairOp::ToleranceEnforcement);
    }

    #[test]
    fn max_displacement_scales_with_tolerance() {
        let tol = Tolerance::default();
        let config = ImportRepairConfig::default_config(tol);
        let expected = tol.linear() * 10.0;
        assert!(
            (config.max_displacement - expected).abs() < f64::EPSILON,
            "max_displacement should be 10x linear tolerance"
        );
    }

    #[test]
    fn single_op_pipeline_runs() {
        let tol = Tolerance::default();
        let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
        let config = ImportRepairConfig {
            operations: vec![ImportRepairOp::ToleranceEnforcement],
            max_displacement: tol.linear() * 10.0,
            tolerance: tol,
        };
        let report = run_import_repair(&mut brep, &config);
        assert_eq!(report.steps.len(), 1);
        assert_eq!(report.steps[0].0, ImportRepairOp::ToleranceEnforcement);
    }

    #[test]
    fn empty_brep_returns_no_defects() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        let config = ImportRepairConfig::default_config(tol);
        let report = run_import_repair(&mut brep, &config);

        assert_eq!(report.total_defects_before, 0);
        assert_eq!(report.total_defects_after, 0);
    }
}
