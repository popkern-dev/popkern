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
pub mod options;
pub mod quality;
pub mod regional;
pub mod repair;
pub mod singularity;

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
pub use options::{
    AnalysisOptions, CheckSet, HealOptions, HealStrategy, SampleDensity, ToleranceBounds,
};
pub use quality::{
    ContinuityLevel, QualityDimension, QualityProfile, QualityResult, QualityTier,
    validate_against_profile,
};
pub use regional::{CurvatureBounds, RegionalAnalysis};
pub use repair::gap_repair::{GapRepairResult, GapRepairStrategy};
pub use repair::pipeline::HealResult;
pub use repair::preview::RepairPreview;
pub use repair::simplification::{AnalyticalKind, SimplificationResult};
pub use repair::tolerance_mgmt::ToleranceAdjustment;
pub use singularity::{SingularityKind, SingularityRecord};

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
