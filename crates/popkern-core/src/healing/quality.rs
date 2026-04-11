//! Quality tiers and profiles for surface quality assessment.
//!
//! Defines three predefined quality tiers (Automotive Class-A, Industrial,
//! General) and a custom profile mechanism. Models are validated against
//! a profile to produce per-dimension pass/fail results.

use serde::{Deserialize, Serialize};

/// Geometric continuity level.
///
/// Each level subsumes the previous: G1 implies G0, G2 implies G1, etc.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum ContinuityLevel {
    /// Positional continuity — shared boundary within position tolerance.
    G0,
    /// Tangent-plane continuity — surface normals parallel along boundary.
    G1,
    /// Curvature continuity — normal curvature matches across boundary.
    G2,
    /// Curvature-rate continuity — curvature derivative matches.
    G3,
}

impl std::fmt::Display for ContinuityLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ContinuityLevel::G0 => write!(f, "G0"),
            ContinuityLevel::G1 => write!(f, "G1"),
            ContinuityLevel::G2 => write!(f, "G2"),
            ContinuityLevel::G3 => write!(f, "G3"),
        }
    }
}

/// Predefined quality tiers (Section 4.1).
///
/// Each tier maps to a `QualityProfile` with appropriate thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QualityTier {
    /// Automotive Class-A: G2+ continuity, tight curvature constraints,
    /// highlight-line and reflection-line continuity.
    Tier1,
    /// Industrial/functional: G1 continuity at smooth boundaries,
    /// functional curvature constraints, no self-intersections.
    Tier2,
    /// General/watertight: G0 continuity, manifold closure,
    /// no gaps, no self-intersections.
    Tier3,
}

impl std::fmt::Display for QualityTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QualityTier::Tier1 => write!(f, "Tier 1 (Automotive Class-A)"),
            QualityTier::Tier2 => write!(f, "Tier 2 (Industrial)"),
            QualityTier::Tier3 => write!(f, "Tier 3 (General)"),
        }
    }
}

/// A single check threshold within a custom quality profile.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CheckThreshold {
    /// Identifier for the check (e.g., "max_position_gap").
    pub check_id: String,
    /// Threshold value for the check.
    pub threshold: f64,
    /// Weight for scoring (higher = more important).
    pub weight: f64,
}

/// A quality profile: named collection of thresholds defining acceptable quality.
///
/// Use `QualityProfile::from_tier()` for predefined profiles, or construct
/// a custom profile by setting individual fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityProfile {
    /// Human-readable name for this profile.
    pub name: String,
    /// Minimum required continuity at smooth boundaries.
    pub min_continuity: ContinuityLevel,
    /// Maximum allowed position gap at G0 boundaries (model units).
    pub max_position_gap: f64,
    /// Maximum allowed inter-normal angle at smooth boundaries (degrees).
    pub max_normal_angle_deg: f64,
    /// Maximum allowed curvature difference at G2 boundaries.
    pub max_curvature_diff: f64,
    /// Minimum allowed edge length (model units).
    pub min_edge_length: f64,
    /// Minimum allowed face area (model units squared).
    pub min_face_area: f64,
    /// Additional custom check thresholds.
    pub custom_checks: Vec<CheckThreshold>,
}

impl QualityProfile {
    /// Create a quality profile from a predefined tier.
    pub fn from_tier(tier: QualityTier) -> Self {
        match tier {
            QualityTier::Tier1 => Self {
                name: "Automotive Class-A".into(),
                min_continuity: ContinuityLevel::G2,
                max_position_gap: 1e-6,
                max_normal_angle_deg: 0.1,
                max_curvature_diff: 1e-3,
                min_edge_length: 1e-4,
                min_face_area: 1e-8,
                custom_checks: Vec::new(),
            },
            QualityTier::Tier2 => Self {
                name: "Industrial".into(),
                min_continuity: ContinuityLevel::G1,
                max_position_gap: 1e-4,
                max_normal_angle_deg: 1.0,
                max_curvature_diff: 1.0,
                min_edge_length: 1e-4,
                min_face_area: 1e-6,
                custom_checks: Vec::new(),
            },
            QualityTier::Tier3 => Self {
                name: "General".into(),
                min_continuity: ContinuityLevel::G0,
                max_position_gap: 1e-3,
                max_normal_angle_deg: 10.0,
                max_curvature_diff: 100.0,
                min_edge_length: 1e-5,
                min_face_area: 1e-8,
                custom_checks: Vec::new(),
            },
        }
    }
}

/// Result of a single quality dimension check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityDimension {
    /// Name of the checked dimension (e.g., "Max position gap").
    pub name: String,
    /// Measured value.
    pub measured: f64,
    /// Threshold from the quality profile.
    pub threshold: f64,
    /// Whether this dimension passed.
    pub passed: bool,
}

impl std::fmt::Display for QualityDimension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = if self.passed { "PASS" } else { "FAIL" };
        write!(
            f,
            "[{status}] {}: {:.6} (threshold: {:.6})",
            self.name, self.measured, self.threshold
        )
    }
}

/// Aggregate quality assessment result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityResult {
    /// The tier that was checked against.
    pub tier: QualityTier,
    /// Whether all dimensions passed.
    pub passed: bool,
    /// Per-dimension results.
    pub dimensions: Vec<QualityDimension>,
}

impl QualityResult {
    /// Create from a tier and dimension results.
    pub fn new(tier: QualityTier, dimensions: Vec<QualityDimension>) -> Self {
        let passed = dimensions.iter().all(|d| d.passed);
        Self {
            tier,
            passed,
            dimensions,
        }
    }

    /// Get only the failed dimensions.
    pub fn failures(&self) -> Vec<&QualityDimension> {
        self.dimensions.iter().filter(|d| !d.passed).collect()
    }
}

impl std::fmt::Display for QualityResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = if self.passed { "PASSED" } else { "FAILED" };
        write!(f, "Quality {}: {status}", self.tier)?;
        for dim in &self.dimensions {
            write!(f, "\n  {dim}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn continuity_level_ordering() {
        assert!(ContinuityLevel::G0 < ContinuityLevel::G1);
        assert!(ContinuityLevel::G1 < ContinuityLevel::G2);
        assert!(ContinuityLevel::G2 < ContinuityLevel::G3);
    }

    #[test]
    fn tier1_profile_has_g2_continuity() {
        let profile = QualityProfile::from_tier(QualityTier::Tier1);
        assert_eq!(profile.min_continuity, ContinuityLevel::G2);
        assert!(profile.max_position_gap < 1e-5);
    }

    #[test]
    fn tier3_profile_has_g0_continuity() {
        let profile = QualityProfile::from_tier(QualityTier::Tier3);
        assert_eq!(profile.min_continuity, ContinuityLevel::G0);
        assert!(profile.max_position_gap > 1e-4);
    }

    #[test]
    fn quality_result_passes_when_all_dimensions_pass() {
        let result = QualityResult::new(
            QualityTier::Tier3,
            vec![QualityDimension {
                name: "gap".into(),
                measured: 1e-6,
                threshold: 1e-3,
                passed: true,
            }],
        );
        assert!(result.passed);
        assert!(result.failures().is_empty());
    }

    #[test]
    fn quality_result_fails_when_any_dimension_fails() {
        let result = QualityResult::new(
            QualityTier::Tier1,
            vec![
                QualityDimension {
                    name: "gap".into(),
                    measured: 1e-6,
                    threshold: 1e-3,
                    passed: true,
                },
                QualityDimension {
                    name: "angle".into(),
                    measured: 5.0,
                    threshold: 0.1,
                    passed: false,
                },
            ],
        );
        assert!(!result.passed);
        assert_eq!(result.failures().len(), 1);
    }

    #[test]
    fn tier_display() {
        assert!(QualityTier::Tier1.to_string().contains("Class-A"));
        assert!(QualityTier::Tier2.to_string().contains("Industrial"));
    }
}
