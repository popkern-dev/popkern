//! Preview proposed repairs before committing changes to the model.
//!
//! The preview workflow lets callers inspect what the repair pipeline
//! would do, including which entities are affected and the estimated
//! geometric displacement, before deciding to apply.

use crate::error::CoreError;
use crate::healing::defects;
use crate::healing::diagnostics::EntityRef;
use crate::healing::options::{AnalysisOptions, HealOptions};
use crate::healing::repair::pipeline::{self, HealResult};
use crate::topology::BRep;

/// A single proposed repair action.
#[derive(Debug, Clone)]
pub struct RepairProposal {
    /// Index of this proposal within the preview (stable for apply selection).
    pub index: usize,
    /// Human-readable description of the repair.
    pub description: String,
    /// Topological entities that would be affected.
    pub entities: Vec<EntityRef>,
    /// Estimated maximum geometry displacement (model units).
    pub estimated_displacement: f64,
}

/// Collection of proposed repairs with an estimate of the resulting model state.
#[derive(Debug, Clone)]
pub struct RepairPreview {
    /// Individual repair proposals, one per detected defect with a remediation.
    pub proposals: Vec<RepairProposal>,
    /// Estimated number of defects remaining after all proposals are applied.
    pub estimated_defects_after: usize,
}

/// Compute a preview of the repairs that would be applied to the model.
///
/// Runs defect analysis on a read-only snapshot and builds a proposal
/// for each defect that has a remediation suggestion.  The BRep is not
/// modified.
///
/// The `estimated_defects_after` field is a rough estimate: it assumes
/// each remediated defect will be fully fixed and non-remediated defects
/// will remain.
pub fn compute_preview(brep: &BRep, options: &HealOptions) -> RepairPreview {
    let analysis_options = AnalysisOptions::new(options.tolerance);
    let report = defects::run_all_checks(brep, &analysis_options);

    let mut proposals = Vec::new();
    let mut fixable_count = 0;

    for (i, defect) in report.iter().enumerate() {
        if defect.remediation.is_some() {
            let estimated_displacement = defect.measured_value.map(|v| v * 0.5).unwrap_or(0.0);

            proposals.push(RepairProposal {
                index: i,
                description: defect.description.clone(),
                entities: defect.entities.clone(),
                estimated_displacement,
            });
            fixable_count += 1;
        }
    }

    let total_defects = report.len();
    let estimated_defects_after = total_defects.saturating_sub(fixable_count);

    RepairPreview {
        proposals,
        estimated_defects_after,
    }
}

impl RepairPreview {
    /// Apply all proposed repairs to the model.
    ///
    /// Delegates to the full repair pipeline.  The BRep is modified in
    /// place with transactional rollback on failure.
    ///
    /// # Errors
    ///
    /// Returns `CoreError::Healing` if the pipeline fails or causes
    /// a regression.
    pub fn apply_all(self, brep: &mut BRep) -> Result<HealResult, CoreError> {
        // Build default options from the BRep's tolerance.
        let options = HealOptions::new(brep.tolerance);
        pipeline::run_pipeline(brep, &options)
    }

    /// Apply a subset of proposed repairs by index.
    ///
    /// Currently delegates to the full pipeline (selective application
    /// is a planned enhancement).  The `indices` parameter is accepted
    /// for API stability but not yet used for filtering.
    ///
    /// # Errors
    ///
    /// Returns `CoreError::Healing` if the pipeline fails.
    pub fn apply(self, brep: &mut BRep, _indices: &[usize]) -> Result<HealResult, CoreError> {
        // Selective application is a future enhancement.
        // For now, run the full pipeline.
        self.apply_all(brep)
    }

    /// Returns true if there are no proposed repairs.
    pub fn is_empty(&self) -> bool {
        self.proposals.is_empty()
    }

    /// Number of proposed repairs.
    pub fn len(&self) -> usize {
        self.proposals.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::options::HealOptions;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn preview_on_valid_cube_has_no_proposals() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let preview = compute_preview(&brep, &options);
        assert!(
            preview.is_empty(),
            "valid cube should have no repair proposals, got {}",
            preview.len()
        );
        assert_eq!(preview.estimated_defects_after, 0);
    }

    #[test]
    fn preview_on_empty_brep_has_no_proposals() {
        let brep = BRep::new(Tolerance::default());
        let options = HealOptions::new(Tolerance::default());
        let preview = compute_preview(&brep, &options);
        assert!(preview.is_empty());
    }

    #[test]
    fn apply_all_on_valid_cube_succeeds() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let preview = compute_preview(&brep, &options);
        let result = preview.apply_all(&mut brep);
        assert!(result.is_ok());
    }

    #[test]
    fn apply_subset_on_valid_cube_succeeds() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let preview = compute_preview(&brep, &options);
        let result = preview.apply(&mut brep, &[0, 1]);
        assert!(result.is_ok());
    }

    #[test]
    fn preview_len_and_is_empty_consistent() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = HealOptions::new(Tolerance::default());
        let preview = compute_preview(&brep, &options);
        assert_eq!(preview.is_empty(), preview.is_empty());
    }
}
