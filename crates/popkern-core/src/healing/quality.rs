//! Quality tiers and profiles for surface quality assessment.
//!
//! Defines three predefined quality tiers (Automotive Class-A, Industrial,
//! General) and a custom profile mechanism. Models are validated against
//! a profile to produce per-dimension pass/fail results.

use serde::{Deserialize, Serialize};

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

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
    /// The quality tier this profile belongs to.
    pub tier: QualityTier,
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
                tier: QualityTier::Tier1,
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
                tier: QualityTier::Tier2,
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
                tier: QualityTier::Tier3,
                min_continuity: ContinuityLevel::G0,
                max_position_gap: 1e-3,
                max_normal_angle_deg: 180.0,
                max_curvature_diff: 100.0,
                min_edge_length: 1e-5,
                min_face_area: 1e-8,
                custom_checks: Vec::new(),
            },
        }
    }
}

/// Validate a BRep model against a quality profile.
///
/// Measures each quality dimension and returns pass/fail results.
/// This runs the continuity and geometric checks relevant to the profile's
/// tier, comparing measured values against the profile thresholds.
pub fn validate_against_profile(
    brep: &BRep,
    profile: &QualityProfile,
    _tolerance: &Tolerance,
) -> QualityResult {
    let max_gap = measure_max_position_gap(brep);
    let max_angle = measure_max_normal_angle(brep);
    let min_len = measure_min_edge_length(brep);

    let dimensions = vec![
        QualityDimension {
            name: "Max position gap".into(),
            measured: max_gap,
            threshold: profile.max_position_gap,
            passed: max_gap <= profile.max_position_gap,
        },
        QualityDimension {
            name: "Max normal angle (deg)".into(),
            measured: max_angle,
            threshold: profile.max_normal_angle_deg,
            passed: max_angle <= profile.max_normal_angle_deg,
        },
        QualityDimension {
            name: "Min edge length".into(),
            measured: min_len,
            threshold: profile.min_edge_length,
            // For min edge length, pass when measured >= threshold.
            passed: min_len >= profile.min_edge_length,
        },
    ];

    QualityResult::new(profile.tier, dimensions)
}

/// Measure worst-case position gap across all edges.
///
/// For each edge with a curve, evaluates the curve at the start and end
/// parameters and compares against the corresponding vertex positions.
/// Returns the maximum gap found, or `0.0` if the model has no curved edges.
fn measure_max_position_gap(brep: &BRep) -> f64 {
    let mut worst = 0.0_f64;

    for (_eh, edge) in &brep.edges {
        let Some(ref curve) = edge.curve else {
            continue;
        };

        // Get the two endpoint vertices via the edge's half-edges.
        let Some(he1) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(he2) = brep.half_edges.get(he1.twin) else {
            continue;
        };

        // he1.origin is at t_start, he2.origin (twin's origin) is at t_end.
        let Some(v_start) = brep.vertices.get(he1.origin) else {
            continue;
        };
        let Some(v_end) = brep.vertices.get(he2.origin) else {
            continue;
        };

        let curve_start = curve.geometry.point_at(curve.t_start);
        let curve_end = curve.geometry.point_at(curve.t_end);

        let gap_start = v_start.position.distance_to(&curve_start);
        let gap_end = v_end.position.distance_to(&curve_end);
        worst = worst.max(gap_start).max(gap_end);
    }

    worst
}

/// Measure worst-case inter-normal angle across edges (degrees).
///
/// For each edge that has two distinct adjacent faces with surface geometry,
/// evaluates both surface normals at the edge midpoint and computes the
/// angle between them. Returns the maximum angle in degrees, or `0.0` if
/// no qualifying edges exist.
fn measure_max_normal_angle(brep: &BRep) -> f64 {
    let mut worst = 0.0_f64;

    for (_eh, edge) in &brep.edges {
        let Some(ref curve) = edge.curve else {
            continue;
        };

        // Get the two half-edges.
        let Some(he1) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(he2) = brep.half_edges.get(he1.twin) else {
            continue;
        };

        // Get the two adjacent faces via loop -> face.
        let Some(loop1) = brep.loops.get(he1.loop_) else {
            continue;
        };
        let Some(loop2) = brep.loops.get(he2.loop_) else {
            continue;
        };
        if loop1.face == loop2.face {
            continue; // seam edge — skip
        }

        let Some(face1) = brep.faces.get(loop1.face) else {
            continue;
        };
        let Some(face2) = brep.faces.get(loop2.face) else {
            continue;
        };
        let Some(ref surf1) = face1.surface else {
            continue;
        };
        let Some(ref surf2) = face2.surface else {
            continue;
        };

        // Get UV coordinates at the edge midpoint via pcurves or fallback.
        let (u1, v1) = uv_from_pcurve_or_center(brep, edge.half_edge, surf1);
        let (u2, v2) = uv_from_pcurve_or_center(brep, he1.twin, surf2);
        let _ = curve; // suppress unused warning

        let n1 = surf1.geometry.normal_at(u1, v1);
        let n2 = surf2.geometry.normal_at(u2, v2);

        // Compute angle between normals.
        let len1 = n1.length();
        let len2 = n2.length();
        if len1 < 1e-15 || len2 < 1e-15 {
            continue;
        }
        let cos_angle = (n1.dot(&n2) / (len1 * len2)).clamp(-1.0, 1.0);
        let angle_deg = cos_angle.acos().to_degrees();
        worst = worst.max(angle_deg);
    }

    worst
}

/// Get UV coordinates from a half-edge's pcurve midpoint, or fall back to
/// the surface parameter range center.
fn uv_from_pcurve_or_center(
    brep: &BRep,
    he_handle: crate::topology::HalfEdgeHandle,
    surface: &crate::geometry::SurfaceOnFace,
) -> (f64, f64) {
    if let Some(he) = brep.half_edges.get(he_handle)
        && let Some(ref pcurve) = he.pcurve
    {
        let (t0, t1) = pcurve.parameter_range();
        let t_mid = (t0 + t1) * 0.5;
        let uv = pcurve.point_at(t_mid);
        return (uv.x, uv.y);
    }
    // Fall back to surface center.
    let ((u_min, u_max), (v_min, v_max)) = surface.geometry.parameter_range();
    let u = if u_min.is_finite() && u_max.is_finite() {
        (u_min + u_max) * 0.5
    } else {
        0.0
    };
    let v = if v_min.is_finite() && v_max.is_finite() {
        (v_min + v_max) * 0.5
    } else {
        0.0
    };
    (u, v)
}

/// Measure the shortest edge in the model.
///
/// For edges with a curve, uses `length_between`. For edges without a curve,
/// falls back to the straight-line distance between endpoint vertices.
/// Returns `f64::INFINITY` if the model has no edges.
fn measure_min_edge_length(brep: &BRep) -> f64 {
    let mut shortest = f64::INFINITY;

    for (_eh, edge) in &brep.edges {
        let len = if let Some(ref curve) = edge.curve {
            curve.geometry.length_between(curve.t_start, curve.t_end)
        } else {
            // Fall back to vertex distance.
            let Some(he1) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(he2) = brep.half_edges.get(he1.twin) else {
                continue;
            };
            let Some(v1) = brep.vertices.get(he1.origin) else {
                continue;
            };
            let Some(v2) = brep.vertices.get(he2.origin) else {
                continue;
            };
            v1.position.distance_to(&v2.position)
        };

        if len < shortest {
            shortest = len;
        }
    }

    shortest
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
#[path = "quality_tests.rs"]
mod tests;
