//! Defect detection for BRep models.
//!
//! Runs topological, geometric, degenerate, and import-specific checks,
//! producing a [`HealingReport`] with structured [`DefectRecord`]s. Each
//! check category lives in its own submodule and can be enabled or disabled
//! via [`AnalysisOptions`].

pub mod containment;
pub mod degenerate;
pub mod geometric;
pub mod import;
pub mod parameterization;
pub mod topological;

use crate::healing::diagnostics::{DefectCategory, DefectRecord, HealingReport, Severity};
use crate::healing::options::AnalysisOptions;
use crate::topology::BRep;

/// Run all enabled defect checks on a BRep.
///
/// Dispatches to the individual check modules based on the flags in
/// `options.checks`. Returns an aggregate report containing every
/// defect found across all enabled categories.
pub fn run_all_checks(brep: &BRep, options: &AnalysisOptions) -> HealingReport {
    let mut report = HealingReport::new();

    if options.checks.topological {
        topological::check_all(brep, &options.tolerance, &mut report);
    }
    if options.checks.geometric {
        geometric::check_all(
            brep,
            &options.tolerance,
            &options.sample_density,
            &mut report,
        );
    }
    if options.checks.degenerate {
        degenerate::check_all(brep, &options.tolerance, &mut report);
    }
    if options.checks.import {
        import::check_all(brep, &options.tolerance, &mut report);
    }

    if options.checks.quality_profile
        && let Some(tier) = options.quality_tier
    {
        let profile = crate::healing::quality::QualityProfile::from_tier(tier);
        let qr =
            crate::healing::quality::validate_against_profile(brep, &profile, &options.tolerance);
        for dim in &qr.dimensions {
            if !dim.passed {
                report.push(DefectRecord {
                    category: DefectCategory::Geometric,
                    severity: Severity::Warning,
                    description: format!("Quality check failed: {dim}"),
                    entities: Vec::new(),
                    location: None,
                    measured_value: Some(dim.measured),
                    threshold: Some(dim.threshold),
                    remediation: None,
                });
            }
        }
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::options::AnalysisOptions;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn run_all_checks_on_valid_cube_returns_empty_report() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = AnalysisOptions::new(Tolerance::default());
        let report = run_all_checks(&brep, &options);
        assert!(
            report.is_empty(),
            "Valid unit cube should produce no defects, got {} defects",
            report.len()
        );
    }

    #[test]
    fn run_all_checks_respects_disabled_categories() {
        use crate::healing::options::CheckSet;

        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let options = AnalysisOptions {
            checks: CheckSet {
                topological: false,
                geometric: false,
                degenerate: false,
                import: false,
                continuity: false,
                quality_profile: false,
            },
            ..AnalysisOptions::new(Tolerance::default())
        };
        let report = run_all_checks(&brep, &options);
        assert!(
            report.is_empty(),
            "All checks disabled should produce no defects"
        );
    }

    #[test]
    fn run_all_checks_on_empty_brep() {
        let brep = BRep::new(Tolerance::default());
        let options = AnalysisOptions::new(Tolerance::default());
        let report = run_all_checks(&brep, &options);
        // An empty BRep has no entities, so no defects.
        assert!(report.is_empty());
    }
}
