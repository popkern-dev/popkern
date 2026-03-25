//! Healing framework for repairing imported BRep models.
//!
//! After import, models may have small gaps, overlapping edges, or other
//! defects. The healing pipeline runs a configurable set of repair passes
//! to improve model quality.

pub mod sew_gaps;

use popkern_core::topology::BRep;

use crate::error::IoError;
use crate::step::diagnostics::{DiagnosticEntry, DiagnosticLog};

/// A single healing pass that can inspect and repair BRep defects.
pub trait HealingPass: Send + Sync {
    /// Human-readable name for logging.
    fn name(&self) -> &str;

    /// Returns true if this pass can attempt to repair the given diagnostic.
    fn applies_to(&self, entry: &DiagnosticEntry) -> bool;

    /// Attempt to repair the defect described by `entry`.
    ///
    /// Returns `Ok(true)` if repaired, `Ok(false)` if no change, `Err` on failure.
    fn attempt_repair(
        &self,
        brep: &mut BRep,
        entry: &DiagnosticEntry,
        log: &mut DiagnosticLog,
    ) -> Result<bool, IoError>;
}

/// How aggressively to heal imported models.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingLevel {
    /// Skip all healing.
    None,
    /// Apply only safe, conservative repairs.
    Basic,
    /// Apply aggressive repairs (may alter geometry).
    Aggressive,
}

/// Run healing passes on a set of imported solids.
///
/// Applies registered healing passes to each solid. At `Basic` level, only
/// the sew-gaps pass (vertex merging + edge curve reconciliation) runs.
/// At `Aggressive` level, the same passes run (more may be added later).
pub fn run_healing(
    solids: &mut [BRep],
    diagnostics: &DiagnosticLog,
    level: HealingLevel,
    log: &mut DiagnosticLog,
) -> Result<(), IoError> {
    if matches!(level, HealingLevel::None) {
        return Ok(());
    }

    let passes: Vec<Box<dyn HealingPass>> = vec![Box::new(sew_gaps::SewGapsPass)];

    for solid in solids.iter_mut() {
        for entry in diagnostics.entries() {
            for pass in &passes {
                if pass.applies_to(entry) {
                    pass.attempt_repair(solid, entry, log)?;
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn healing_none_is_noop() {
        let mut log = DiagnosticLog::new();
        let result = run_healing(&mut [], &DiagnosticLog::new(), HealingLevel::None, &mut log);
        assert!(result.is_ok());
        assert!(log.is_empty());
    }

    #[test]
    fn healing_basic_on_empty_is_noop() {
        let mut log = DiagnosticLog::new();
        let result = run_healing(
            &mut [],
            &DiagnosticLog::new(),
            HealingLevel::Basic,
            &mut log,
        );
        assert!(result.is_ok());
        // No solids and no diagnostics means no passes fire.
        assert!(log.is_empty());
    }

    #[test]
    fn healing_aggressive_on_empty_is_noop() {
        let mut log = DiagnosticLog::new();
        let result = run_healing(
            &mut [],
            &DiagnosticLog::new(),
            HealingLevel::Aggressive,
            &mut log,
        );
        assert!(result.is_ok());
        assert!(log.is_empty());
    }

    #[test]
    fn healing_basic_on_valid_brep_makes_no_changes() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut solids = vec![BRep::unit_cube(Tolerance::default()).unwrap()];
        let face_count_before = solids[0].faces.len();

        // A warning-level diagnostic with Validate phase triggers SewGapsPass.
        let mut diagnostics = DiagnosticLog::new();
        diagnostics.warn(crate::step::diagnostics::Phase::Validate, None, "test gap");

        let mut log = DiagnosticLog::new();
        let result = run_healing(&mut solids, &diagnostics, HealingLevel::Basic, &mut log);
        assert!(result.is_ok());
        // The cube is already valid — no structural changes expected.
        assert_eq!(solids[0].faces.len(), face_count_before);
    }

    #[test]
    fn healing_aggressive_on_valid_brep_same_as_basic() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut solids = vec![BRep::unit_cube(Tolerance::default()).unwrap()];
        let face_count_before = solids[0].faces.len();

        let mut diagnostics = DiagnosticLog::new();
        diagnostics.warn(crate::step::diagnostics::Phase::Validate, None, "test gap");

        let mut log = DiagnosticLog::new();
        let result = run_healing(
            &mut solids,
            &diagnostics,
            HealingLevel::Aggressive,
            &mut log,
        );
        assert!(result.is_ok());
        assert_eq!(solids[0].faces.len(), face_count_before);
    }

    #[test]
    fn healing_none_skips_all_levels() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut solids = vec![BRep::unit_cube(Tolerance::default()).unwrap()];
        let mut diagnostics = DiagnosticLog::new();
        diagnostics.warn(crate::step::diagnostics::Phase::Validate, None, "gap");

        let mut log = DiagnosticLog::new();
        // None level must return immediately without touching solids.
        let result = run_healing(&mut solids, &diagnostics, HealingLevel::None, &mut log);
        assert!(result.is_ok());
        // Log must remain empty — no passes ran.
        assert!(log.is_empty());
    }
}
