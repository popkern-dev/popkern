//! Diagnostic types for healing analysis and repair.
//!
//! Provides structured defect records with severity, entity references,
//! geometric location, and remediation suggestions. Mirrors the
//! `BoolDiagnostics` pattern from popkern-bool for consistency.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::topology::{
    EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, ShellHandle, VertexHandle,
};

/// Severity level for diagnostic entries.
///
/// Ordered from least to most severe. Matches popkern-bool's `Severity`
/// for cross-crate consistency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Severity {
    /// Notable condition, no action needed.
    Info,
    /// Degraded quality — should be repaired for production use.
    Warning,
    /// Model is invalid without repair.
    Error,
    /// System-level failure or data corruption.
    Critical,
}

impl fmt::Display for Severity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Severity::Info => write!(f, "info"),
            Severity::Warning => write!(f, "warning"),
            Severity::Error => write!(f, "error"),
            Severity::Critical => write!(f, "critical"),
        }
    }
}

/// Category of a detected defect — maps to spec Sections 5.1–5.4.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DefectCategory {
    /// Broken topology: closure, manifold, orientation (Section 5.1).
    Topological,
    /// Geometric inconsistency: edge-surface deviation, containment (Section 5.2).
    Geometric,
    /// Degenerate entities: small faces/edges, missing geometry (Section 5.3).
    Degenerate,
    /// Tolerance hierarchy violations (Section 5.2.3).
    Tolerance,
    /// Import-specific defects: lost topology, representation errors (Section 5.4).
    Import,
}

impl fmt::Display for DefectCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DefectCategory::Topological => write!(f, "topological"),
            DefectCategory::Geometric => write!(f, "geometric"),
            DefectCategory::Degenerate => write!(f, "degenerate"),
            DefectCategory::Tolerance => write!(f, "tolerance"),
            DefectCategory::Import => write!(f, "import"),
        }
    }
}

/// Reference to one or more affected topological entities.
///
/// Used by UI for selection highlighting and by repair for targeting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EntityRef {
    /// A single vertex.
    Vertex(VertexHandle),
    /// A single edge.
    Edge(EdgeHandle),
    /// A single half-edge.
    HalfEdge(HalfEdgeHandle),
    /// A single face.
    Face(FaceHandle),
    /// A single loop.
    Loop(LoopHandle),
    /// A single shell.
    Shell(ShellHandle),
    /// A pair of edges (e.g., a gap between two edges).
    EdgePair(EdgeHandle, EdgeHandle),
    /// A pair of faces (e.g., self-intersection between two faces).
    FacePair(FaceHandle, FaceHandle),
}

/// Geometric location of a defect — for UI camera positioning and highlighting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DefectLocation {
    /// A specific 3D point.
    Point(Point3),
    /// A point in surface parameter space.
    ParameterPoint {
        /// U parameter.
        u: f64,
        /// V parameter.
        v: f64,
    },
    /// A region in surface parameter space.
    ParameterRegion {
        /// U parameter range.
        u_range: (f64, f64),
        /// V parameter range.
        v_range: (f64, f64),
    },
    /// A 3D bounding box enclosing the defect region.
    BoundingBox(Aabb),
}

/// Suggested repair action for a defect.
///
/// Mirrors popkern-bool's `Remediation` pattern. The UI can present
/// these as one-click repair options.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Remediation {
    /// What action to take.
    pub action: RemediationAction,
    /// Which parameter to adjust, if applicable.
    pub parameter: Option<String>,
    /// Suggested new value, if applicable.
    pub value: Option<f64>,
    /// Human-readable explanation of why this might help.
    pub rationale: String,
}

/// Possible remediation actions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RemediationAction {
    /// Widen tolerance to absorb the defect.
    AbsorbTolerance,
    /// Move geometry to close a gap or fix inconsistency.
    AdjustGeometry,
    /// Remove a degenerate entity.
    RemoveEntity,
    /// Split an entity at a defect location.
    SplitEntity,
    /// Replace a NURBS surface with an exact analytical form.
    SimplifySurface,
    /// Fix vertex >= edge >= face tolerance ordering.
    EnforceToleranceHierarchy,
    /// Reconstruct missing periodic-surface closure topology.
    RestorePeriodicTopology,
    /// Reconstruct missing singularity topology at poles.
    RestoreSingularityTopology,
    /// Re-derive edge curve from adjacent surfaces.
    RecomputeEdgeCurve,
    /// Reparameterize a surface to reduce UV distortion.
    ReparameterizeSurface,
    /// Recompute a half-edge's parameter curve from 3D geometry.
    RecomputePCurve,
}

/// A single detected defect with full diagnostic context.
///
/// Designed for both programmatic consumption and UI display.
/// The `description` field is a human-readable, UI-safe string.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DefectRecord {
    /// Defect category (topological, geometric, etc.).
    pub category: DefectCategory,
    /// How severe the defect is.
    pub severity: Severity,
    /// Human-readable description suitable for UI display.
    pub description: String,
    /// Affected entities — for selection highlighting.
    pub entities: Vec<EntityRef>,
    /// Geometric location — for camera positioning.
    pub location: Option<DefectLocation>,
    /// The measured value that triggered this defect.
    pub measured_value: Option<f64>,
    /// The threshold that was exceeded.
    pub threshold: Option<f64>,
    /// Suggested repair action.
    pub remediation: Option<Remediation>,
}

impl fmt::Display for DefectRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}] {}: {}",
            self.severity, self.category, self.description
        )?;
        if let (Some(measured), Some(threshold)) = (self.measured_value, self.threshold) {
            write!(f, " (measured: {measured:.6}, threshold: {threshold:.6})")?;
        }
        Ok(())
    }
}

/// Summary statistics for a healing report.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HealingSummary {
    /// Total number of defects.
    pub total: usize,
    /// Number of critical defects.
    pub critical: usize,
    /// Number of error defects.
    pub error: usize,
    /// Number of warning defects.
    pub warning: usize,
    /// Number of informational defects.
    pub info: usize,
}

/// Aggregate analysis report with filtering and summary methods.
///
/// Collects all defects found during analysis. Provides filtering
/// by severity, category, and entity for UI presentation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HealingReport {
    defects: Vec<DefectRecord>,
}

impl HealingReport {
    /// Create an empty report.
    pub fn new() -> Self {
        Self {
            defects: Vec::new(),
        }
    }

    /// Add a defect to the report.
    pub fn push(&mut self, defect: DefectRecord) {
        self.defects.push(defect);
    }

    /// Number of defects in the report.
    pub fn len(&self) -> usize {
        self.defects.len()
    }

    /// Returns true if no defects were found.
    pub fn is_empty(&self) -> bool {
        self.defects.is_empty()
    }

    /// Returns true if no defects of Warning severity or above were found.
    pub fn is_clean(&self) -> bool {
        !self.defects.iter().any(|d| d.severity >= Severity::Warning)
    }

    /// Returns the highest severity found, or `None` if the report is empty.
    pub fn worst_severity(&self) -> Option<Severity> {
        self.defects.iter().map(|d| d.severity).max()
    }

    /// Iterate over all defects.
    pub fn iter(&self) -> impl Iterator<Item = &DefectRecord> {
        self.defects.iter()
    }

    /// Filter defects by severity level (returns those at or above the given level).
    pub fn by_severity(&self, min_severity: Severity) -> Vec<&DefectRecord> {
        self.defects
            .iter()
            .filter(|d| d.severity >= min_severity)
            .collect()
    }

    /// Filter defects by category.
    pub fn by_category(&self, category: DefectCategory) -> Vec<&DefectRecord> {
        self.defects
            .iter()
            .filter(|d| d.category == category)
            .collect()
    }

    /// Get defects affecting a specific vertex.
    pub fn for_vertex(&self, handle: VertexHandle) -> Vec<&DefectRecord> {
        self.defects
            .iter()
            .filter(|d| {
                d.entities
                    .iter()
                    .any(|e| matches!(e, EntityRef::Vertex(v) if *v == handle))
            })
            .collect()
    }

    /// Get defects affecting a specific edge.
    pub fn for_edge(&self, handle: EdgeHandle) -> Vec<&DefectRecord> {
        self.defects
            .iter()
            .filter(|d| {
                d.entities
                    .iter()
                    .any(|e| matches!(e, EntityRef::Edge(eh) if *eh == handle))
            })
            .collect()
    }

    /// Get defects affecting a specific face.
    pub fn for_face(&self, handle: FaceHandle) -> Vec<&DefectRecord> {
        self.defects
            .iter()
            .filter(|d| {
                d.entities
                    .iter()
                    .any(|e| matches!(e, EntityRef::Face(fh) if *fh == handle))
            })
            .collect()
    }

    /// Compute summary statistics.
    pub fn summary(&self) -> HealingSummary {
        let mut s = HealingSummary {
            total: self.defects.len(),
            ..HealingSummary::default()
        };
        for d in &self.defects {
            match d.severity {
                Severity::Critical => s.critical += 1,
                Severity::Error => s.error += 1,
                Severity::Warning => s.warning += 1,
                Severity::Info => s.info += 1,
            }
        }
        s
    }

    /// Consume the report and return the defect vector.
    pub fn into_defects(self) -> Vec<DefectRecord> {
        self.defects
    }

    /// Access the defect slice.
    pub fn defects(&self) -> &[DefectRecord] {
        &self.defects
    }
}

impl fmt::Display for HealingReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = self.summary();
        write!(
            f,
            "Healing report: {} defects ({} critical, {} error, {} warning, {} info)",
            s.total, s.critical, s.error, s.warning, s.info
        )
    }
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod diagnostics_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_defect(severity: Severity, category: DefectCategory) -> DefectRecord {
        DefectRecord {
            category,
            severity,
            description: format!("Test {category} defect"),
            entities: Vec::new(),
            location: None,
            measured_value: None,
            threshold: None,
            remediation: None,
        }
    }

    #[test]
    fn empty_report_is_clean() {
        let report = HealingReport::new();
        assert!(report.is_empty());
        assert!(report.is_clean());
        assert!(report.worst_severity().is_none());
    }

    #[test]
    fn report_with_info_only_is_clean() {
        let mut report = HealingReport::new();
        report.push(sample_defect(Severity::Info, DefectCategory::Geometric));
        assert!(!report.is_empty());
        assert!(report.is_clean());
    }

    #[test]
    fn report_with_warning_is_not_clean() {
        let mut report = HealingReport::new();
        report.push(sample_defect(
            Severity::Warning,
            DefectCategory::Topological,
        ));
        assert!(!report.is_clean());
        assert_eq!(report.worst_severity(), Some(Severity::Warning));
    }

    #[test]
    fn filter_by_severity() {
        let mut report = HealingReport::new();
        report.push(sample_defect(Severity::Info, DefectCategory::Geometric));
        report.push(sample_defect(
            Severity::Warning,
            DefectCategory::Topological,
        ));
        report.push(sample_defect(Severity::Error, DefectCategory::Degenerate));
        assert_eq!(report.by_severity(Severity::Warning).len(), 2);
        assert_eq!(report.by_severity(Severity::Error).len(), 1);
    }

    #[test]
    fn filter_by_category() {
        let mut report = HealingReport::new();
        report.push(sample_defect(Severity::Warning, DefectCategory::Geometric));
        report.push(sample_defect(Severity::Warning, DefectCategory::Geometric));
        report.push(sample_defect(
            Severity::Warning,
            DefectCategory::Topological,
        ));
        assert_eq!(report.by_category(DefectCategory::Geometric).len(), 2);
        assert_eq!(report.by_category(DefectCategory::Topological).len(), 1);
    }

    #[test]
    fn summary_counts() {
        let mut report = HealingReport::new();
        report.push(sample_defect(Severity::Info, DefectCategory::Geometric));
        report.push(sample_defect(
            Severity::Warning,
            DefectCategory::Topological,
        ));
        report.push(sample_defect(Severity::Error, DefectCategory::Degenerate));
        report.push(sample_defect(Severity::Critical, DefectCategory::Tolerance));
        let s = report.summary();
        assert_eq!(s.total, 4);
        assert_eq!(s.info, 1);
        assert_eq!(s.warning, 1);
        assert_eq!(s.error, 1);
        assert_eq!(s.critical, 1);
    }

    #[test]
    fn severity_ordering() {
        assert!(Severity::Info < Severity::Warning);
        assert!(Severity::Warning < Severity::Error);
        assert!(Severity::Error < Severity::Critical);
    }

    #[test]
    fn defect_record_display() {
        let defect = DefectRecord {
            category: DefectCategory::Geometric,
            severity: Severity::Warning,
            description: "Gap between edges".into(),
            entities: Vec::new(),
            location: None,
            measured_value: Some(0.05),
            threshold: Some(0.001),
            remediation: None,
        };
        let s = defect.to_string();
        assert!(s.contains("warning"));
        assert!(s.contains("geometric"));
        assert!(s.contains("Gap between edges"));
        assert!(s.contains("0.050000"));
    }

    #[test]
    fn report_display() {
        let mut report = HealingReport::new();
        report.push(sample_defect(Severity::Error, DefectCategory::Geometric));
        let s = report.to_string();
        assert!(s.contains("1 defects"));
        assert!(s.contains("1 error"));
    }
}
