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
}
