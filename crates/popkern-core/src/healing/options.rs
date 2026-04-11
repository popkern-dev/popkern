//! Configuration types for healing analysis and repair.
//!
//! Follows the Options struct pattern used by `FilletOptions`, `SweepOptions`,
//! and `ShellOptions` — named structs with `Default` implementations for
//! progressive disclosure.

use serde::{Deserialize, Serialize};

use crate::healing::quality::QualityTier;
use crate::tolerance::Tolerance;

/// Which categories of checks to run during analysis.
///
/// Each field enables or disables a category of defect checks.
/// Use the associated constants for common configurations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckSet {
    /// Topological checks: closure, manifold, orientation (Section 5.1).
    pub topological: bool,
    /// Geometric checks: edge-surface consistency, containment (Section 5.2).
    pub geometric: bool,
    /// Degenerate entity checks: small faces/edges, missing geometry (Section 5.3).
    pub degenerate: bool,
    /// Import-specific checks: tolerance mismatch, lost topology (Section 5.4).
    pub import: bool,
    /// Continuity analysis: G0/G1/G2 measurement across edges.
    pub continuity: bool,
    /// Quality profile validation against a tier.
    pub quality_profile: bool,
}

impl CheckSet {
    /// All checks enabled.
    pub const ALL: Self = Self {
        topological: true,
        geometric: true,
        degenerate: true,
        import: true,
        continuity: true,
        quality_profile: true,
    };

    /// Fast checks only — skips expensive continuity and quality analysis.
    pub const FAST: Self = Self {
        topological: true,
        geometric: true,
        degenerate: true,
        import: false,
        continuity: false,
        quality_profile: false,
    };

    /// Geometric checks only — useful for targeted analysis.
    pub const GEOMETRIC_ONLY: Self = Self {
        topological: false,
        geometric: true,
        degenerate: false,
        import: false,
        continuity: false,
        quality_profile: false,
    };
}

impl Default for CheckSet {
    fn default() -> Self {
        Self::ALL
    }
}

/// Sampling density for adaptive analysis operations.
///
/// Controls the number of sample points per edge or per surface patch
/// during deviation measurement and continuity analysis.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum SampleDensity {
    /// 8 samples per edge — fast but may miss narrow defects.
    Coarse,
    /// 16 samples per edge — good balance of speed and accuracy.
    #[default]
    Standard,
    /// 32 samples per edge — thorough, for final validation.
    Fine,
}

impl SampleDensity {
    /// Number of samples per edge for this density.
    pub fn samples_per_edge(self) -> usize {
        match self {
            SampleDensity::Coarse => 8,
            SampleDensity::Standard => 16,
            SampleDensity::Fine => 32,
        }
    }
}

/// Configuration for healing analysis.
///
/// Use `AnalysisOptions::new(tolerance)` for sensible defaults, or
/// customize individual fields via struct update syntax:
///
/// ```ignore
/// let options = AnalysisOptions {
///     checks: CheckSet::FAST,
///     quality_tier: Some(QualityTier::Tier2),
///     ..AnalysisOptions::new(tolerance)
/// };
/// ```
#[derive(Debug, Clone)]
pub struct AnalysisOptions {
    /// Which categories of checks to run.
    pub checks: CheckSet,
    /// Quality tier to validate against (if `checks.quality_profile` is true).
    pub quality_tier: Option<QualityTier>,
    /// Sampling density for deviation and continuity measurement.
    pub sample_density: SampleDensity,
    /// Tolerance for geometric comparisons.
    pub tolerance: Tolerance,
}

impl AnalysisOptions {
    /// Create options with sensible defaults and the given tolerance.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            checks: CheckSet::ALL,
            quality_tier: None,
            sample_density: SampleDensity::Standard,
            tolerance,
        }
    }
}

/// Strategy for automated healing — controls which repairs are attempted.
///
/// For a UI slider: Conservative < Standard < Aggressive.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum HealStrategy {
    /// Only tolerance absorption — no geometry changes.
    ///
    /// Widens entity tolerances to absorb small gaps. Fast and
    /// non-destructive, but does not improve actual geometry.
    Conservative,
    /// Gap repair + tolerance enforcement + surface simplification.
    ///
    /// Moves vertices and adjusts edge curves to close gaps.
    /// Enforces tolerance hierarchy. Detects and replaces
    /// over-parameterized NURBS with exact analytical forms.
    #[default]
    Standard,
    /// All repairs including aggressive geometry adjustment.
    ///
    /// Everything in Standard, plus larger displacements allowed
    /// and more aggressive simplification thresholds.
    Aggressive,
}

/// Bounds on entity tolerances during repair.
///
/// Slider-friendly: a UI can present `min_linear` and `max_linear`
/// as adjustable bounds with the defaults as starting values.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ToleranceBounds {
    /// Minimum allowed linear tolerance for any entity.
    pub min_linear: f64,
    /// Maximum allowed linear tolerance for any entity.
    pub max_linear: f64,
}

impl Default for ToleranceBounds {
    fn default() -> Self {
        Self {
            min_linear: 1e-10,
            max_linear: 1e-3,
        }
    }
}

/// Configuration for healing repair.
///
/// Use `HealOptions::new(tolerance)` for sensible defaults, or
/// customize individual fields via struct update syntax:
///
/// ```ignore
/// let options = HealOptions {
///     strategy: HealStrategy::Aggressive,
///     max_displacement: 0.1,
///     ..HealOptions::new(tolerance)
/// };
/// ```
#[derive(Debug, Clone)]
pub struct HealOptions {
    /// Which repair strategy to use.
    pub strategy: HealStrategy,
    /// Maximum geometry displacement allowed per repair (in model units).
    ///
    /// Repairs that would move geometry further than this are rejected.
    /// Default: 10x the linear tolerance.
    pub max_displacement: f64,
    /// Maximum repair iterations. Each iteration runs the pipeline
    /// and re-checks for remaining defects.
    pub max_iterations: usize,
    /// Bounds on entity tolerances during repair.
    pub tolerance_bounds: ToleranceBounds,
    /// Tolerance for geometric comparisons.
    pub tolerance: Tolerance,
}

impl HealOptions {
    /// Create options with sensible defaults and the given tolerance.
    ///
    /// Defaults: Standard strategy, max displacement = 10x linear tolerance,
    /// 3 iterations, default tolerance bounds.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            strategy: HealStrategy::Standard,
            max_displacement: tolerance.linear * 10.0,
            max_iterations: 3,
            tolerance_bounds: ToleranceBounds::default(),
            tolerance,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_set_all_enables_everything() {
        let cs = CheckSet::ALL;
        assert!(cs.topological);
        assert!(cs.geometric);
        assert!(cs.degenerate);
        assert!(cs.import);
        assert!(cs.continuity);
        assert!(cs.quality_profile);
    }

    #[test]
    fn check_set_fast_skips_expensive() {
        let cs = CheckSet::FAST;
        assert!(cs.topological);
        assert!(cs.geometric);
        assert!(!cs.continuity);
        assert!(!cs.quality_profile);
    }

    #[test]
    fn sample_density_values() {
        assert_eq!(SampleDensity::Coarse.samples_per_edge(), 8);
        assert_eq!(SampleDensity::Standard.samples_per_edge(), 16);
        assert_eq!(SampleDensity::Fine.samples_per_edge(), 32);
    }

    #[test]
    fn analysis_options_defaults() {
        let tol = Tolerance::default();
        let opts = AnalysisOptions::new(tol);
        assert!(opts.checks.topological);
        assert!(opts.quality_tier.is_none());
        assert_eq!(opts.sample_density, SampleDensity::Standard);
    }

    #[test]
    fn heal_options_defaults() {
        let tol = Tolerance::default();
        let opts = HealOptions::new(tol);
        assert_eq!(opts.strategy, HealStrategy::Standard);
        assert_eq!(opts.max_iterations, 3);
        assert!((opts.max_displacement - tol.linear * 10.0).abs() < 1e-15);
    }

    #[test]
    fn heal_strategy_ordering() {
        assert!(HealStrategy::Conservative < HealStrategy::Standard);
        assert!(HealStrategy::Standard < HealStrategy::Aggressive);
    }

    #[test]
    fn tolerance_bounds_defaults() {
        let bounds = ToleranceBounds::default();
        assert_eq!(bounds.min_linear, 1e-10);
        assert_eq!(bounds.max_linear, 1e-3);
    }

    #[test]
    fn struct_update_syntax_works() {
        let tol = Tolerance::default();
        let opts = HealOptions {
            strategy: HealStrategy::Aggressive,
            max_displacement: 0.1,
            ..HealOptions::new(tol)
        };
        assert_eq!(opts.strategy, HealStrategy::Aggressive);
        assert_eq!(opts.max_displacement, 0.1);
        assert_eq!(opts.max_iterations, 3); // inherited from default
    }
}
