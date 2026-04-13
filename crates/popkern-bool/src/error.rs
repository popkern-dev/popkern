//! Error and diagnostic types for popkern-bool.

use popkern_core::{
    error::CoreError, geometry::Point3, geometry::curve_on_surface::CurveOnSurface,
    topology::FaceHandle,
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

    /// The caller cancelled the operation via a cancellation token (§15).
    /// The input model is unchanged.
    #[error("boolean operation cancelled by caller")]
    Cancelled,

    /// At least one input is an open shell / sheet body and the
    /// requested operation does not yet have a sheet-aware code
    /// path (§2.3). The input model is unchanged. Callers can
    /// detect this case via [`crate::sheet::is_sheet_body`] before
    /// calling the boolean.
    #[error(
        "boolean operation on sheet body not yet supported: operand_a_is_sheet={a}, operand_b_is_sheet={b}"
    )]
    OpenShellNotSupported {
        /// True if operand A is an open shell.
        a: bool,
        /// True if operand B is an open shell.
        b: bool,
    },

    /// Interference detection failed (§13.2). The check could not
    /// determine whether the bodies overlap, typically because the
    /// underlying surface intersection or BVH lookup errored.
    #[error("interference detection failed: {0}")]
    InterferenceFailed(String),

    /// System-level resource failure (§13.2): out of memory,
    /// per-stage iteration limit exceeded, etc.
    #[error("system error: {0}")]
    SystemError(String),
}

/// Per-stage timing and counter sample (§16 telemetry).
///
/// One `PhaseMetric` is recorded per pipeline stage that runs to completion.
/// `entity_count` is stage-specific: pair count for intersection, fragment
/// count for split/classify, etc.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseMetric {
    /// Which pipeline stage this measurement covers.
    pub stage: Stage,
    /// Wall-clock duration of the stage.
    pub duration: std::time::Duration,
    /// Stage-specific entity count (pairs, fragments, etc.).
    pub entity_count: usize,
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
    /// Per-phase timing and counters (§16 telemetry).
    pub phase_metrics: Vec<PhaseMetric>,
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
            location: None,
            measured_value: None,
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
            location: None,
            measured_value: None,
            remediation: None,
        });
    }

    /// Add a warning diagnostic with a 3D location (§13.1).
    pub fn warn_at(&mut self, phase: Stage, message: impl Into<String>, location: Point3) {
        self.entries.push(DiagnosticEntry {
            severity: Severity::Warning,
            phase,
            message: message.into(),
            entities: Vec::new(),
            location: Some(location),
            measured_value: None,
            remediation: None,
        });
    }

    /// Add a warning diagnostic with a measured geometric value (§13.1).
    pub fn warn_with_measure(
        &mut self,
        phase: Stage,
        message: impl Into<String>,
        measured_value: f64,
    ) {
        self.entries.push(DiagnosticEntry {
            severity: Severity::Warning,
            phase,
            message: message.into(),
            entities: Vec::new(),
            location: None,
            measured_value: Some(measured_value),
            remediation: None,
        });
    }

    /// Returns true if any entry has `Error` or `Critical` severity.
    pub fn has_errors(&self) -> bool {
        self.entries
            .iter()
            .any(|e| matches!(e.severity, Severity::Error | Severity::Critical))
    }

    /// Populate missing `Remediation` fields on every `Error`- or
    /// `Critical`-severity entry using the top-level `error`'s
    /// [`BoolError::remediation_hint`] default.
    ///
    /// Entries that already have a custom remediation are left
    /// untouched. Info / Warning entries are skipped because they
    /// don't represent failure conditions. Used by the error-exit
    /// paths in the pipeline orchestrator to guarantee that every
    /// returned `BoolDiagnostics` carries actionable remediation
    /// data (§13.4).
    pub fn fill_missing_remediations(&mut self) {
        let hint = match &self.error {
            Some(e) => e.remediation_hint(),
            None => return,
        };
        for entry in &mut self.entries {
            if entry.remediation.is_none()
                && matches!(entry.severity, Severity::Error | Severity::Critical)
            {
                entry.remediation = Some(hint.clone());
            }
        }
    }

    /// Build a diagnostic from a top-level [`BoolError`], auto-
    /// populating a single error entry with a remediation hint
    /// derived from the error variant (§13.4).
    ///
    /// Use this from failure sites that don't have enough local
    /// context to build richer diagnostics — the resulting
    /// `BoolDiagnostics` has one entry carrying the error message
    /// and the default-mapped [`Remediation`]. Call sites with
    /// richer context (measured values, 3D locations, specific
    /// entity refs) should still hand-build their entries; this
    /// builder is the floor, not the ceiling.
    pub fn from_error(error: BoolError, phase: Stage) -> Self {
        let hint = error.remediation_hint();
        let message = error.to_string();
        Self {
            error: Some(error),
            intersection_curves: Vec::new(),
            face_classifications: Vec::new(),
            failed_stage: Some(phase),
            entries: vec![DiagnosticEntry {
                severity: Severity::Error,
                phase,
                message,
                entities: Vec::new(),
                location: None,
                measured_value: None,
                remediation: Some(hint),
            }],
            phase_metrics: Vec::new(),
        }
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
    /// The requested operation is not applicable (e.g. attempting a
    /// boolean on disjoint bodies with no intersection). No action
    /// is needed — the caller should recognise the no-op.
    NoAction,
    /// Split the operation into smaller parts, or run the caller's
    /// own pre-check to identify problematic face pairs before
    /// retrying.
    Precheck,
}

impl Remediation {
    /// Construct a remediation with an action and rationale and no
    /// tolerance-parameter metadata.
    pub fn new(action: RemediationAction, rationale: impl Into<String>) -> Self {
        Self {
            action,
            parameter: None,
            value: None,
            rationale: rationale.into(),
        }
    }

    /// Construct a tolerance-adjusting remediation that suggests
    /// `new_value` for `parameter`.
    pub fn adjust_tolerance(
        action: RemediationAction,
        parameter: impl Into<String>,
        new_value: f64,
        rationale: impl Into<String>,
    ) -> Self {
        Self {
            action,
            parameter: Some(parameter.into()),
            value: Some(new_value),
            rationale: rationale.into(),
        }
    }
}

impl BoolError {
    /// Best-fit remediation hint for this error variant.
    ///
    /// Returns a structured [`Remediation`] suggestion derived
    /// purely from the error kind — no access to tolerances or
    /// input geometry is needed. Callers can use this as a default
    /// when the failure site didn't populate its own entry-level
    /// `Remediation`, or to build user-facing error messages.
    ///
    /// §13.4 acceptance criterion: every `BoolError` variant
    /// returns a non-`None` hint, so `remediation_hint()` is a
    /// total function.
    pub fn remediation_hint(&self) -> Remediation {
        match self {
            BoolError::NoIntersection => Remediation::new(
                RemediationAction::NoAction,
                "The two solids do not intersect. No boolean result exists \
                 for this operation. If you expected an intersection, verify \
                 that the bodies are positioned as intended — a coordinate \
                 error in one operand is the most common cause.",
            ),
            BoolError::IntersectionFailed(_) => Remediation::new(
                RemediationAction::RetryWithPerturbation,
                "Surface intersection solver failed. Try \
                 `boolean_with_perturbation()` to slightly offset the \
                 geometry and retry, or increase `linear` tolerance to \
                 widen the solver's convergence envelope.",
            ),
            BoolError::SplitFailed { .. } => Remediation::new(
                RemediationAction::IncreaseTolerance,
                "Face splitting failed — typically caused by an \
                 intersection curve that lies too close to an existing \
                 edge or vertex. Try increasing `linear` tolerance, or \
                 heal the input to snap near-coincident features before \
                 retrying.",
            ),
            BoolError::ClassificationFailed { .. } => Remediation::new(
                RemediationAction::IncreaseTolerance,
                "Point-in-solid classification failed for a face \
                 fragment. Try increasing `coincidence` tolerance so the \
                 multi-ray sampler treats near-boundary points \
                 deterministically, or heal the input to eliminate \
                 near-tangent face pairs.",
            ),
            BoolError::MergeFailed { .. } => Remediation::new(
                RemediationAction::HealInput,
                "Assembling the result solid from selected face \
                 fragments failed. This usually indicates non-manifold \
                 topology in the intersection graph. Heal the input \
                 geometry (close gaps, remove duplicate vertices) before \
                 retrying.",
            ),
            BoolError::ValidationFailed(_) => Remediation::new(
                RemediationAction::HealInput,
                "Post-merge validation flagged the result solid as \
                 invalid. The boolean pipeline may have produced slivers, \
                 non-manifold edges, or open loops. Try healing the \
                 result, or retry with tighter tolerances so the merge \
                 phase produces cleaner topology.",
            ),
            BoolError::DegenerateInput(_) => Remediation::new(
                RemediationAction::HealInput,
                "One or both inputs are geometrically degenerate (zero \
                 volume, duplicate vertices, or zero-area faces). Run \
                 the healing pipeline on the inputs before retrying the \
                 boolean.",
            ),
            BoolError::CoincidentFaces { .. } => Remediation::new(
                RemediationAction::UseCoincidenceMode,
                "Two faces of the input solids are geometrically \
                 coincident, creating an ambiguous boundary. Set the \
                 `BoolTolerance::mode` to `AssumedFull` or \
                 `ToleranceEscalation` so the coincidence engine \
                 handles the overlap explicitly instead of the \
                 transversal solver.",
            ),
            BoolError::CoreError(_) => Remediation::new(
                RemediationAction::HealInput,
                "A core topology or geometry operation failed. This is \
                 usually a sign of invalid input; running the healing \
                 pipeline and re-running validation before retrying the \
                 boolean is the safest first step.",
            ),
            BoolError::Cancelled => Remediation::new(
                RemediationAction::NoAction,
                "The operation was cancelled by the caller via a \
                 cancellation token. No action needed — the input \
                 model is unchanged.",
            ),
            BoolError::OpenShellNotSupported { .. } => Remediation::new(
                RemediationAction::Precheck,
                "At least one input is an open shell / sheet body. \
                 Route through `boolean_with_sheets()` instead of the \
                 standard `boolean()` entry point — the sheet-aware \
                 dispatcher handles the operation matrix for open \
                 shells.",
            ),
            BoolError::InterferenceFailed(_) => Remediation::new(
                RemediationAction::RetryWithPerturbation,
                "Interference detection could not determine whether the \
                 bodies overlap. This is usually a near-tangent or \
                 near-coincident configuration. Try \
                 `boolean_with_perturbation()` or increase `linear` \
                 tolerance.",
            ),
            BoolError::SystemError(_) => Remediation::new(
                RemediationAction::Precheck,
                "A system-level resource failure occurred (out of \
                 memory, per-stage iteration limit exceeded). Simplify \
                 the input (fewer faces, decimated NURBS surfaces, \
                 smaller models) and retry.",
            ),
        }
    }
}
