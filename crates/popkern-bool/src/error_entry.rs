//! Diagnostic entries: severity, the entities involved, and the
//! measurement and remediation a failure reports.

use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle;

use super::{Remediation, Stage};

/// Severity level for diagnostic entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// Succeeded, notable condition (e.g. "12 redundant edges removed").
    Info,
    /// Succeeded with degraded quality (e.g. "3 sliver faces in result").
    Warning,
    /// Failed, model unchanged.
    Error,
    /// System-level failure, possible resource leak.
    Critical,
}

/// A reference to a topological entity involved in a diagnostic.
#[derive(Debug, Clone)]
pub enum EntityRef {
    /// A face handle.
    Face(FaceHandle),
    /// A pair of faces (e.g. a coincident face pair).
    FacePair(FaceHandle, FaceHandle),
}

/// A structured diagnostic entry with severity, entity references, and
/// optional remediation suggestion.
#[derive(Debug, Clone)]
pub struct DiagnosticEntry {
    /// How severe this diagnostic is.
    pub severity: Severity,
    /// Which pipeline phase produced this diagnostic.
    pub phase: Stage,
    /// Human-readable description.
    pub message: String,
    /// Topological entities involved.
    pub entities: Vec<EntityRef>,
    /// Measured geometric value associated with the diagnostic (§13.1).
    /// For example a "gap of 0.003mm" diagnostic would set this to
    /// `Some(0.003)`. The unit is the model's linear unit; angular
    /// values are stored in radians. `None` when the diagnostic is
    /// not associated with a single scalar measurement.
    pub measured_value: Option<f64>,
    /// 3D location where the problem occurs (§13.1). Populated when the
    /// failure site has a known coordinate (e.g. an intersection endpoint,
    /// a face centroid, a measured gap point); `None` for diagnostics
    /// that are intrinsically topological.
    pub location: Option<Point3>,
    /// Machine-readable remediation suggestion.
    pub remediation: Option<Remediation>,
}
