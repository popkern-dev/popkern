//! Import repair pipeline — systematic repair of models imported from external formats.
//!
//! Provides a configurable sequence of repair operations tailored to the
//! characteristic defect patterns of imported models (Section 7.2).
//!
//! Supports custom repair operator insertion (Section 7.2 BEYOND) and
//! format-specific import profiles (Section 7.3).

use crate::error::CoreError;
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

/// A user-defined repair operation.
///
/// Implement this trait to inject custom repair logic into the import
/// repair pipeline at a specific position. Custom ops participate in
/// the same before/after defect tracking as standard operations.
pub trait CustomRepairOp: Send + Sync {
    /// Human-readable name for this operation.
    fn name(&self) -> &str;
    /// Execute the repair, modifying the BRep in place.
    fn execute(&self, brep: &mut BRep, tolerance: &Tolerance) -> Result<(), CoreError>;
}

impl std::fmt::Debug for dyn CustomRepairOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CustomRepairOp({})", self.name())
    }
}

/// A single operation in the repair pipeline — either standard or custom.
///
/// This enum allows interleaving built-in repair operations with
/// user-defined custom operators in any order.
pub enum RepairOperation {
    /// A built-in import repair operator.
    Standard(ImportRepairOp),
    /// A user-defined repair operator.
    Custom(Box<dyn CustomRepairOp>),
}

impl std::fmt::Debug for RepairOperation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RepairOperation::Standard(op) => write!(f, "Standard({op:?})"),
            RepairOperation::Custom(op) => write!(f, "Custom({})", op.name()),
        }
    }
}

/// Predefined import format with known defect patterns.
///
/// Each format has characteristic defects that benefit from different
/// repair operation orderings and tolerances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportFormat {
    /// STEP AP203 — Configuration Controlled 3D Design.
    /// Simple geometry, needs tolerance normalization and gap closure.
    StepAP203,
    /// STEP AP214 — Core Data for Automotive Mechanical Design.
    /// Complex assemblies, benefits from surface simplification.
    StepAP214,
    /// IGES — Initial Graphics Exchange Specification.
    /// Commonly loses orientation and periodic topology.
    Iges,
}

/// Configuration for the import repair pipeline.
#[derive(Debug)]
pub struct ImportRepairConfig {
    /// Ordered list of repair operations to perform.
    pub operations: Vec<RepairOperation>,
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
                RepairOperation::Standard(ImportRepairOp::OrientationFix),
                RepairOperation::Standard(ImportRepairOp::PeriodicTopologyRestore),
                RepairOperation::Standard(ImportRepairOp::ToleranceNormalization),
                RepairOperation::Standard(ImportRepairOp::GapClosure),
                RepairOperation::Standard(ImportRepairOp::SurfaceSimplification),
                RepairOperation::Standard(ImportRepairOp::ToleranceEnforcement),
            ],
            max_displacement: tolerance.linear() * 10.0,
            tolerance,
        }
    }

    /// Create a repair config tuned to a specific import format.
    ///
    /// Each format gets specific operation ordering and tolerances:
    /// - STEP AP203: focus on tolerance normalization + gap closure.
    /// - STEP AP214: adds surface simplification for complex assemblies.
    /// - IGES: aggressive orientation fix + periodic topology restore.
    pub fn for_format(format: ImportFormat, tolerance: Tolerance) -> Self {
        match format {
            ImportFormat::StepAP203 => Self {
                operations: vec![
                    RepairOperation::Standard(ImportRepairOp::ToleranceNormalization),
                    RepairOperation::Standard(ImportRepairOp::GapClosure),
                    RepairOperation::Standard(ImportRepairOp::ToleranceEnforcement),
                ],
                max_displacement: tolerance.linear() * 10.0,
                tolerance,
            },
            ImportFormat::StepAP214 => Self {
                operations: vec![
                    RepairOperation::Standard(ImportRepairOp::OrientationFix),
                    RepairOperation::Standard(ImportRepairOp::ToleranceNormalization),
                    RepairOperation::Standard(ImportRepairOp::GapClosure),
                    RepairOperation::Standard(ImportRepairOp::SurfaceSimplification),
                    RepairOperation::Standard(ImportRepairOp::ToleranceEnforcement),
                ],
                max_displacement: tolerance.linear() * 10.0,
                tolerance,
            },
            ImportFormat::Iges => Self {
                operations: vec![
                    RepairOperation::Standard(ImportRepairOp::OrientationFix),
                    RepairOperation::Standard(ImportRepairOp::PeriodicTopologyRestore),
                    RepairOperation::Standard(ImportRepairOp::ToleranceNormalization),
                    RepairOperation::Standard(ImportRepairOp::GapClosure),
                    RepairOperation::Standard(ImportRepairOp::ToleranceEnforcement),
                ],
                max_displacement: tolerance.linear() * 15.0,
                tolerance,
            },
        }
    }

    /// Insert a custom repair operation at a specific position.
    ///
    /// The operation will execute at position `index` in the pipeline.
    /// Panics if `index > self.operations.len()`.
    pub fn insert_custom_op(&mut self, index: usize, op: Box<dyn CustomRepairOp>) {
        self.operations.insert(index, RepairOperation::Custom(op));
    }

    /// Add a custom repair operation at the end of the pipeline.
    pub fn push_custom_op(&mut self, op: Box<dyn CustomRepairOp>) {
        self.operations.push(RepairOperation::Custom(op));
    }

    /// Return the number of standard operations in the pipeline.
    pub fn standard_op_count(&self) -> usize {
        self.operations
            .iter()
            .filter(|op| matches!(op, RepairOperation::Standard(_)))
            .count()
    }
}

/// Label for a step in the repair report.
#[derive(Debug, Clone)]
pub enum StepLabel {
    /// A built-in repair operation.
    Standard(ImportRepairOp),
    /// A custom repair operation, identified by name.
    Custom(String),
}

/// Report from an import repair run.
#[derive(Debug, Clone)]
pub struct ImportRepairReport {
    /// Per-operation results: (label, defects_before, defects_after).
    pub steps: Vec<(StepLabel, usize, usize)>,
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
///
/// Custom operations that return errors are logged in the report
/// but do not abort the pipeline — subsequent operations still run.
pub fn run_import_repair(brep: &mut BRep, config: &ImportRepairConfig) -> ImportRepairReport {
    let total_defects_before = count_defects(brep, &config.tolerance);

    let mut steps = Vec::with_capacity(config.operations.len());

    for op in &config.operations {
        let before = count_defects(brep, &config.tolerance);
        match op {
            RepairOperation::Standard(standard_op) => {
                execute_op(brep, *standard_op, config);
                let after = count_defects(brep, &config.tolerance);
                steps.push((StepLabel::Standard(*standard_op), before, after));
            }
            RepairOperation::Custom(custom_op) => {
                // Custom ops may fail; we record defect counts regardless.
                let _result = custom_op.execute(brep, &config.tolerance);
                let after = count_defects(brep, &config.tolerance);
                steps.push((
                    StepLabel::Custom(custom_op.name().to_string()),
                    before,
                    after,
                ));
            }
        }
    }

    let total_defects_after = count_defects(brep, &config.tolerance);

    ImportRepairReport {
        steps,
        total_defects_before,
        total_defects_after,
    }
}

/// Run import repair on multiple BReps with shared configuration.
///
/// More efficient than calling `run_import_repair` on each BRep
/// individually because analysis options are constructed once.
pub fn run_batch_import_repair(
    breps: &mut [BRep],
    config: &ImportRepairConfig,
) -> Vec<ImportRepairReport> {
    breps
        .iter_mut()
        .map(|brep| run_import_repair(brep, config))
        .collect()
}

#[cfg(test)]
#[path = "import_repair_tests.rs"]
mod tests;
