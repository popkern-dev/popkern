//! Surface quality analysis, defect detection, and geometry healing.
//!
//! This module provides a unified system for analyzing B-Rep model quality
//! and repairing geometric defects. It follows an analysis-then-repair
//! workflow: analyze to discover defects, review results, repair selectively,
//! then re-analyze to verify.
//!
//! # Progressive Disclosure
//!
//! Simple path — analyze or heal with sensible defaults:
//!
//! ```ignore
//! let report = healing::analyze(&brep, &tolerance);
//! let result = healing::heal(&mut brep, &tolerance)?;
//! ```
//!
//! Options path — full control via options structs:
//!
//! ```ignore
//! let report = healing::analyze_with(&brep, &AnalysisOptions {
//!     checks: CheckSet::FAST,
//!     quality_tier: Some(QualityTier::Tier2),
//!     ..AnalysisOptions::new(tolerance)
//! });
//! ```
//!
//! Preview path — review proposed repairs before committing:
//!
//! ```ignore
//! let preview = healing::preview_repairs(&brep, &HealOptions::new(tolerance));
//! preview.apply_all(&mut brep)?;
//! ```

pub mod continuity;
pub mod curvature;
pub mod defects;
pub mod deviation;
pub mod diagnostics;
pub mod fairing;
pub mod options;
pub mod parallel;
pub mod quality;
pub mod regional;
pub mod repair;
pub mod singularity;
pub mod telemetry;

// Re-export key types for ergonomic access.
pub use continuity::classification::{EdgeClassification, SmoothChain};
pub use continuity::enforcement::EnforcementResult;
pub use continuity::measurement::ContinuityMeasurement;
pub use curvature::curve_curvature::{curve_curvature_at, curve_torsion_at};
pub use curvature::fundamental_forms::{
    FirstFundamentalForm, SecondFundamentalForm, ShapeOperator,
};
pub use curvature::surface_curvature::{
    SurfaceCurvature, gaussian_curvature, mean_curvature, min_radius_of_curvature,
    principal_curvatures, surface_curvature_at,
};
pub use deviation::{
    DeviationResult, curve_curve_deviation, edge_surface_deviation, surface_surface_deviation,
};
pub use diagnostics::{
    DefectCategory, DefectLocation, DefectRecord, EntityRef, HealingReport, HealingSummary,
    Remediation, RemediationAction, Severity,
};
pub use fairing::FairingResult;
pub use options::{
    AnalysisOptions, CheckSet, HealOptions, HealStrategy, SampleDensity, ToleranceBounds,
};
pub use quality::{
    ContinuityLevel, QualityDimension, QualityProfile, QualityResult, QualityTier,
    validate_against_profile,
};
pub use regional::{CurvatureBounds, RegionalAnalysis};
pub use repair::complexity::{
    ComplexityReductionResult, CurveReductionResult, SurfaceReductionResult,
};
pub use repair::coons_patch::{CoonsBoundary, CoonsPatchResult};
pub use repair::degeneracy::DegeneracyRepairResult;
pub use repair::gap_repair::{GapRepairResult, GapRepairStrategy};
pub use repair::hole_fill::{DetectedHole, HoleFillResult, HoleFillStrategy};
pub use repair::import_repair::{ImportRepairConfig, ImportRepairOp, ImportRepairReport};
pub use repair::merge_faces::FaceMergeResult;
pub use repair::pipeline::HealResult;
pub use repair::preview::RepairPreview;
pub use repair::self_intersection::SelfIntersectionRepairResult;
pub use repair::simplification::{AnalyticalKind, SimplificationResult};
pub use repair::stitch::StitchResult;
pub use repair::tolerance_mgmt::ToleranceAdjustment;
pub use singularity::{SingularityKind, SingularityRecord};
pub use telemetry::{HealingTelemetry, OperationTiming};

use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Analyze a BRep model for defects using default settings.
///
/// Runs all defect checks at standard sampling density.
/// Returns a report with all detected defects.
///
/// For more control, use [`analyze_with`].
pub fn analyze(brep: &BRep, tolerance: &Tolerance) -> HealingReport {
    analyze_with(brep, &AnalysisOptions::new(*tolerance))
}

/// Analyze a BRep model for defects with custom options.
///
/// The options control which checks run, sampling density, and
/// optional quality tier validation.
pub fn analyze_with(brep: &BRep, options: &AnalysisOptions) -> HealingReport {
    defects::run_all_checks(brep, options)
}

/// Heal a BRep model using default settings.
///
/// Runs the Standard repair strategy: gap repair, tolerance enforcement,
/// and surface simplification. The BRep is modified in place.
///
/// Returns a summary of what was repaired. On failure, the BRep is
/// restored to its original state (transactional semantics).
///
/// For more control, use [`heal_with`].
pub fn heal(brep: &mut BRep, tolerance: &Tolerance) -> Result<HealResult, CoreError> {
    heal_with(brep, &HealOptions::new(*tolerance))
}

/// Heal a BRep model with custom options.
///
/// The options control repair strategy, max displacement, iteration
/// count, and tolerance bounds. The BRep is modified in place.
///
/// On failure, the BRep is restored to its original state.
pub fn heal_with(brep: &mut BRep, options: &HealOptions) -> Result<HealResult, CoreError> {
    repair::pipeline::run_pipeline(brep, options)
}

/// Preview proposed repairs without modifying the model.
///
/// Returns a [`RepairPreview`] describing what repairs would be
/// performed. The caller can review and selectively apply them:
///
/// ```ignore
/// let preview = healing::preview_repairs(&brep, &options);
/// for proposal in &preview.proposals {
///     println!("{}: displacement {:.6}", proposal.description, proposal.estimated_displacement);
/// }
/// preview.apply_all(&mut brep)?;
/// ```
pub fn preview_repairs(brep: &BRep, options: &HealOptions) -> RepairPreview {
    repair::preview::compute_preview(brep, options)
}

/// Analyze a BRep model for defects with telemetry tracking.
///
/// Wraps [`analyze_with`] and records per-category timing data.
/// Returns both the healing report and a telemetry collector with
/// one entry per enabled check category.
///
/// ```ignore
/// let (report, telemetry) = healing::analyze_with_telemetry(&brep, &options);
/// println!("{telemetry}");
/// println!("Slowest: {:?}", telemetry.slowest());
/// ```
pub fn analyze_with_telemetry(
    brep: &BRep,
    options: &AnalysisOptions,
) -> (HealingReport, HealingTelemetry) {
    let mut telemetry = HealingTelemetry::new();
    let mut report = HealingReport::new();

    let entity_count = brep.vertices.len() + brep.edges.len() + brep.faces.len();

    if options.checks.topological {
        let start = HealingTelemetry::start_timer();
        defects::topological::check_all(brep, &options.tolerance, &mut report);
        telemetry.record("topological_checks", start.elapsed(), entity_count);
    }

    if options.checks.geometric {
        let start = HealingTelemetry::start_timer();
        defects::geometric::check_all(
            brep,
            &options.tolerance,
            &options.sample_density,
            &mut report,
        );
        telemetry.record("geometric_checks", start.elapsed(), entity_count);
    }

    if options.checks.degenerate {
        let start = HealingTelemetry::start_timer();
        defects::degenerate::check_all(brep, &options.tolerance, &mut report);
        telemetry.record("degenerate_checks", start.elapsed(), entity_count);
    }

    if options.checks.import {
        let start = HealingTelemetry::start_timer();
        defects::import::check_all(brep, &options.tolerance, &mut report);
        telemetry.record("import_checks", start.elapsed(), entity_count);
    }

    if options.checks.quality_profile
        && let Some(tier) = options.quality_tier
    {
        let start = HealingTelemetry::start_timer();
        let profile = quality::QualityProfile::from_tier(tier);
        let qr = quality::validate_against_profile(brep, &profile, &options.tolerance);
        for dim in &qr.dimensions {
            if !dim.passed {
                report.push(diagnostics::DefectRecord {
                    category: diagnostics::DefectCategory::Geometric,
                    severity: diagnostics::Severity::Warning,
                    description: format!("Quality check failed: {dim}"),
                    entities: Vec::new(),
                    location: None,
                    measured_value: Some(dim.measured),
                    threshold: Some(dim.threshold),
                    remediation: None,
                });
            }
        }
        telemetry.record("quality_profile", start.elapsed(), entity_count);
    }

    (report, telemetry)
}

/// Heal a BRep model with telemetry tracking.
///
/// Wraps [`heal_with`] and records the total pipeline duration.
/// Returns the heal result and a telemetry collector.
///
/// ```ignore
/// let (result, telemetry) = healing::heal_with_telemetry(&mut brep, &options)?;
/// println!("{telemetry}");
/// ```
pub fn heal_with_telemetry(
    brep: &mut BRep,
    options: &HealOptions,
) -> Result<(HealResult, HealingTelemetry), CoreError> {
    let mut telemetry = HealingTelemetry::new();
    let entity_count = brep.vertices.len() + brep.edges.len() + brep.faces.len();

    let start = HealingTelemetry::start_timer();
    let result = repair::pipeline::run_pipeline(brep, options)?;
    telemetry.record("repair_pipeline", start.elapsed(), entity_count);

    Ok((result, telemetry))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn analyze_with_telemetry_returns_report_and_timings() {
        let brep = BRep::unit_cube(Tolerance::default())
            .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
        let options = AnalysisOptions::new(Tolerance::default());
        let (report, telemetry) = analyze_with_telemetry(&brep, &options);

        // With ALL checks enabled on a valid cube, report should be clean.
        assert!(report.is_clean());

        // Telemetry should have recorded operations for the enabled categories.
        assert!(
            telemetry.operation_count() > 0,
            "expected at least one timed operation"
        );
        assert!(
            telemetry.total_duration() >= std::time::Duration::ZERO,
            "total duration should be non-negative"
        );
    }

    #[test]
    fn analyze_with_telemetry_disabled_checks_records_nothing() {
        let brep = BRep::new(Tolerance::default());
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
        let (_report, telemetry) = analyze_with_telemetry(&brep, &options);
        assert_eq!(telemetry.operation_count(), 0);
    }

    #[test]
    fn analyze_with_telemetry_matches_analyze_with() {
        let brep = BRep::unit_cube(Tolerance::default())
            .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
        let options = AnalysisOptions::new(Tolerance::default());

        let plain_report = analyze_with(&brep, &options);
        let (tele_report, _telemetry) = analyze_with_telemetry(&brep, &options);

        // Both should find the same number of defects.
        assert_eq!(plain_report.len(), tele_report.len());
    }

    #[test]
    fn heal_with_telemetry_returns_result_and_timings() {
        let mut brep = BRep::unit_cube(Tolerance::default())
            .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
        let options = HealOptions::new(Tolerance::default());
        let result = heal_with_telemetry(&mut brep, &options);
        assert!(result.is_ok());
        let (_heal_result, telemetry) =
            result.unwrap_or_else(|_| panic!("heal_with_telemetry should succeed on a valid cube"));
        assert!(
            telemetry.operation_count() > 0,
            "expected repair_pipeline timing"
        );
    }
}
