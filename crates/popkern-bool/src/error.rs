//! Error and diagnostic types for popkern-bool.

use popkern_core::{
    error::CoreError, geometry::curve_on_surface::CurveOnSurface, topology::FaceHandle,
};
use thiserror::Error;

/// The type of boolean set operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOperation {
    /// Compute the union of two solids (A ∪ B).
    Union,
    /// Compute the intersection of two solids (A ∩ B).
    Intersection,
    /// Subtract the second solid from the first (A − B).
    Subtraction,
}

/// Stages of the boolean pipeline, used to identify where a failure occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// Face-face intersection curve computation.
    Intersection,
    /// Trimming faces along intersection curves.
    Split,
    /// Point-in-solid classification of face fragments.
    Classification,
    /// Selecting which fragments to keep.
    Selection,
    /// Stitching selected fragments into the result solid.
    Merge,
    /// Post-merge topological and geometric validation.
    Validation,
}

/// Classification of a face fragment relative to the opposite solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classification {
    /// The fragment lies entirely outside the opposite solid.
    Outside,
    /// The fragment lies entirely inside the opposite solid.
    Inside,
    /// The fragment lies on the boundary and the normals are co-oriented.
    OnBoundaryCooriented,
    /// The fragment lies on the boundary and the normals are anti-parallel (crossing).
    OnBoundaryCrossing,
}

/// Classification of a 3D point relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointClassification {
    /// The point is strictly inside the solid.
    Inside,
    /// The point is strictly outside the solid.
    Outside,
    /// The point lies on the boundary of the solid.
    OnBoundary,
}

/// Top-level error type for boolean operations.
#[derive(Debug, Error)]
pub enum BoolError {
    /// The two solids do not intersect; no boolean result is possible for this operation.
    #[error("solids do not intersect")]
    NoIntersection,

    /// Computing the face-face intersection curves failed.
    #[error("intersection curve computation failed: {0}")]
    IntersectionFailed(String),

    /// Splitting a face along an intersection curve failed.
    #[error("face split failed on face {face:?}: {reason}")]
    SplitFailed {
        /// The face that could not be split.
        face: FaceHandle,
        /// Description of the failure.
        reason: String,
    },

    /// Point-in-solid classification failed for a face fragment.
    #[error("classification failed on face {face:?}: {reason}")]
    ClassificationFailed {
        /// The face fragment that could not be classified.
        face: FaceHandle,
        /// Description of the failure.
        reason: String,
    },

    /// Stitching fragments into the result solid failed.
    #[error("merge failed at stage {stage:?}: {reason}")]
    MergeFailed {
        /// The merge sub-stage that failed.
        stage: String,
        /// Description of the failure.
        reason: String,
    },

    /// Post-merge validation of the result solid failed.
    #[error("result validation failed: {0}")]
    ValidationFailed(String),

    /// One or both inputs are geometrically degenerate (e.g. zero-volume solid).
    #[error("degenerate input: {0}")]
    DegenerateInput(String),

    /// Two faces of the input solids are coincident, creating an ambiguous boundary.
    #[error("coincident faces: {face_a:?} and {face_b:?}")]
    CoincidentFaces {
        /// Face from the first solid.
        face_a: FaceHandle,
        /// Face from the second solid.
        face_b: FaceHandle,
    },

    /// An error propagated from popkern-core.
    #[error(transparent)]
    CoreError(#[from] CoreError),
}

/// Diagnostic snapshot collected during a boolean operation.
///
/// Populated on both success and failure; useful for debugging and testing.
#[derive(Debug, Default)]
pub struct BoolDiagnostics {
    /// The error that caused the operation to fail, if any.
    pub error: Option<BoolError>,
    /// Intersection curves computed in Stage 1.
    pub intersection_curves: Vec<CurveOnSurface>,
    /// Per-face classification results from Stage 3.
    pub face_classifications: Vec<(FaceHandle, Classification)>,
    /// The pipeline stage at which the operation failed, if any.
    pub failed_stage: Option<Stage>,
    /// Structured diagnostic entries with severity and entity references.
    pub entries: Vec<DiagnosticEntry>,
}

impl BoolDiagnostics {
    /// Add a diagnostic entry.
    pub fn add(&mut self, entry: DiagnosticEntry) {
        self.entries.push(entry);
    }

    /// Add an informational diagnostic.
    pub fn info(&mut self, phase: Stage, message: impl Into<String>) {
        self.entries.push(DiagnosticEntry {
            severity: Severity::Info,
            phase,
            message: message.into(),
            entities: Vec::new(),
            remediation: None,
        });
    }

    /// Add a warning diagnostic.
    pub fn warn(&mut self, phase: Stage, message: impl Into<String>) {
        self.entries.push(DiagnosticEntry {
            severity: Severity::Warning,
            phase,
            message: message.into(),
            entities: Vec::new(),
            remediation: None,
        });
    }

    /// Returns true if any entry has `Error` or `Critical` severity.
    pub fn has_errors(&self) -> bool {
        self.entries
            .iter()
            .any(|e| matches!(e.severity, Severity::Error | Severity::Critical))
    }
}

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
    /// Machine-readable remediation suggestion.
    pub remediation: Option<Remediation>,
}

/// A suggested fix for a boolean operation failure.
#[derive(Debug, Clone)]
pub struct Remediation {
    /// What action to take.
    pub action: RemediationAction,
    /// Which tolerance parameter to adjust, if applicable.
    pub parameter: Option<String>,
    /// Suggested new value, if applicable.
    pub value: Option<f64>,
    /// Why this might help.
    pub rationale: String,
}

/// Possible remediation actions suggested by the engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemediationAction {
    /// Increase a tolerance parameter.
    IncreaseTolerance,
    /// Decrease a tolerance parameter.
    DecreaseTolerance,
    /// Heal the input geometry before retrying.
    HealInput,
    /// Use coincidence detection mode for the problematic face pair.
    UseCoincidenceMode,
    /// Simplify the input geometry.
    SimplifyGeometry,
    /// Retry with geometric perturbation.
    RetryWithPerturbation,
}
