//! Error types for popkern-core.

use crate::topology::{EdgeHandle, FaceHandle};
use thiserror::Error;

mod validation_error;
pub use validation_error::{RegionDefect, ValidationError};

/// Errors from healing analysis and repair operations.
#[derive(Debug, Error)]
pub enum HealingError {
    /// Gap repair failed on a specific edge.
    #[error("gap repair failed on edge {edge:?}: {reason}")]
    GapRepairFailed {
        /// The edge where repair failed.
        edge: EdgeHandle,
        /// Why the repair failed.
        reason: String,
    },

    /// Tolerance hierarchy enforcement failed.
    #[error("tolerance enforcement failed: {0}")]
    ToleranceEnforcementFailed(String),

    /// Surface simplification failed on a specific face.
    #[error("surface simplification failed on face {face:?}: {reason}")]
    SimplificationFailed {
        /// The face where simplification failed.
        face: FaceHandle,
        /// Why the simplification failed.
        reason: String,
    },

    /// Repair pipeline failed at a specific step.
    #[error("repair pipeline failed at step {step}: {reason}")]
    PipelineFailed {
        /// The pipeline step that failed.
        step: String,
        /// Why the step failed.
        reason: String,
    },

    /// A repair would displace geometry beyond the configured limit.
    #[error("max displacement {displacement:.6} exceeded limit {limit:.6} on edge {edge:?}")]
    DisplacementExceeded {
        /// The edge where displacement was too large.
        edge: EdgeHandle,
        /// The actual displacement.
        displacement: f64,
        /// The configured limit.
        limit: f64,
    },

    /// Repair iteration made the model worse (more defects after than before).
    #[error("repair made model worse: {before} defects before, {after} after")]
    RepairRegression {
        /// Defect count before the repair pass.
        before: usize,
        /// Defect count after the repair pass.
        after: usize,
    },
}

/// Top-level error type for popkern-core operations.
#[derive(Debug, Error)]
pub enum CoreError {
    /// A validation check failed.
    #[error(transparent)]
    Validation(#[from] ValidationError),

    /// A healing operation failed.
    #[error(transparent)]
    Healing(#[from] HealingError),

    /// A handle does not exist in the BRep.
    #[error("invalid handle: {0}")]
    InvalidHandle(String),

    /// A topological precondition was not met.
    #[error("precondition failed: {0}")]
    PreconditionFailed(String),

    /// A geometric operation violates tolerance.
    #[error("tolerance violation: {0}")]
    ToleranceViolation(String),

    /// Curve is degenerate (zero length or collapsed).
    #[error("degenerate curve: {0}")]
    DegenerateCurve(String),

    /// Surface is degenerate (zero area or collapsed).
    #[error("degenerate surface: {0}")]
    DegenerateSurface(String),

    /// Invalid NURBS data (bad knot vector, wrong control point count, etc).
    #[error("invalid NURBS: {0}")]
    InvalidNurbs(String),

    /// An operation that needs a manifold edge met a fin whose edge does
    /// not have exactly two fins (a free edge, a non-manifold edge, or a
    /// fin not linked into a ring).
    #[error("edge {edge:?} is not a manifold edge: its radial ring does not hold exactly two fins")]
    NoMate {
        /// The edge of the fin asked for its mate.
        edge: crate::topology::EdgeHandle,
    },

    /// A non-manifold Euler operator declined; the body is unchanged.
    #[error(transparent)]
    Euler(#[from] crate::euler::EulerError),

    /// A body is not a kind the operation accepts (a sheet given to an
    /// operation on solids, a general body to one on manifold bodies).
    #[error("{operation} needs {}; the body is a {}", crate::topology::BodyKind::list_label(accepted), found.label())]
    UnsupportedBodyKind {
        /// The operation.
        operation: &'static str,
        /// The body's kind.
        found: crate::topology::BodyKind,
        /// The kinds the operation accepts.
        accepted: &'static [crate::topology::BodyKind],
    },
}

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
