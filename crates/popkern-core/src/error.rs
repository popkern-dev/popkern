//! Error types for popkern-core.

use crate::topology::{
    EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, SolidHandle, VertexHandle,
};
use thiserror::Error;

/// Errors from geometric/topological validation checks.
#[derive(Debug, Error)]
pub enum ValidationError {
    /// Euler-Poincaré invariant violated.
    #[error("Euler-Poincaré violation on solid {solid:?}: expected {expected}, got {actual}")]
    EulerPoincare {
        /// The solid that failed.
        solid: SolidHandle,
        /// Expected value.
        expected: i64,
        /// Actual computed value.
        actual: i64,
    },

    /// Half-edge has no valid twin.
    #[error("Half-edge {half_edge:?} has no valid twin")]
    BrokenTwin {
        /// The half-edge with a broken twin.
        half_edge: HalfEdgeHandle,
    },

    /// Half-edge cycle is broken (following next doesn't return to start).
    #[error("Half-edge cycle broken at {half_edge:?} in loop {loop_:?}")]
    BrokenCycle {
        /// The half-edge where the cycle breaks.
        half_edge: HalfEdgeHandle,
        /// The loop containing the broken cycle.
        loop_: LoopHandle,
    },

    /// Vertex's outgoing half-edge doesn't originate from it.
    #[error("Vertex {vertex:?} outgoing half-edge {half_edge:?} does not originate from it")]
    InconsistentVertex {
        /// The vertex with inconsistent outgoing half-edge.
        vertex: VertexHandle,
        /// The half-edge that should originate from the vertex.
        half_edge: HalfEdgeHandle,
    },

    /// Edge is non-manifold (not exactly 2 half-edges).
    #[error("Edge {edge:?} is non-manifold")]
    NonManifoldEdge {
        /// The non-manifold edge.
        edge: EdgeHandle,
    },

    /// Vertex has valence less than 2.
    #[error("Vertex {vertex:?} has valence {valence} (expected >= 2)")]
    LowValence {
        /// The low-valence vertex.
        vertex: VertexHandle,
        /// The actual valence.
        valence: usize,
    },

    /// Face normal orientation is inconsistent.
    #[error("Face {face:?} has inconsistent normal orientation")]
    InconsistentOrientation {
        /// The face with inconsistent orientation.
        face: FaceHandle,
    },

    /// Loop boundary is not geometrically closed.
    #[error("Loop {loop_:?} is not closed (gap: {gap})")]
    OpenLoop {
        /// The open loop.
        loop_: LoopHandle,
        /// The geometric gap distance.
        gap: f64,
    },

    /// Edge curve deviates from adjacent surface.
    #[error("Edge {edge:?} deviates from surface by {max_deviation}")]
    CurveSurfaceDeviation {
        /// The deviating edge.
        edge: EdgeHandle,
        /// Maximum deviation distance.
        max_deviation: f64,
    },

    /// Gap between adjacent faces at a shared edge.
    #[error("Shell gap at edge {edge:?} (gap: {gap})")]
    ShellGap {
        /// The edge with a gap.
        edge: EdgeHandle,
        /// The gap distance.
        gap: f64,
    },

    /// Two non-adjacent faces intersect each other.
    #[error("Self-intersection between face {face_a:?} and face {face_b:?}")]
    SelfIntersection {
        /// First face involved in intersection.
        face_a: FaceHandle,
        /// Second face involved in intersection.
        face_b: FaceHandle,
    },

    /// A face's surface normal reverses direction (surface folds over itself).
    #[error("Twisted surface on face {face:?}")]
    TwistedSurface {
        /// The face with a twisted surface.
        face: FaceHandle,
    },

    /// Precision constraint violated.
    #[error("precision: {}", .diagnostic.remediation)]
    Precision {
        /// The diagnostic describing the precision violation.
        diagnostic: crate::precision::PrecisionDiagnostic,
    },
}

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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validation_error_display() {
        let err = ValidationError::LowValence {
            vertex: Default::default(),
            valence: 1,
        };
        let msg = err.to_string();
        assert!(msg.contains("valence"));
    }

    #[test]
    fn core_error_from_validation() {
        let val_err = ValidationError::LowValence {
            vertex: Default::default(),
            valence: 1,
        };
        let core_err: CoreError = val_err.into();
        assert!(matches!(core_err, CoreError::Validation(_)));
    }

    #[test]
    fn core_error_invalid_handle() {
        let err = CoreError::InvalidHandle("vertex 42 not found".into());
        assert!(err.to_string().contains("vertex 42"));
    }

    #[test]
    fn core_error_precondition_failed() {
        let err = CoreError::PreconditionFailed("vertices must be in same loop".into());
        assert!(err.to_string().contains("same loop"));
    }

    #[test]
    fn core_error_tolerance_violation() {
        let err = CoreError::ToleranceViolation("points too far apart".into());
        assert!(err.to_string().contains("too far apart"));
    }

    #[test]
    fn core_error_invalid_nurbs() {
        let err = CoreError::InvalidNurbs("knot vector not non-decreasing".into());
        assert!(err.to_string().contains("knot vector"));
    }

    #[test]
    fn validation_error_euler_poincare() {
        let err = ValidationError::EulerPoincare {
            solid: Default::default(),
            expected: 0,
            actual: 2,
        };
        assert!(err.to_string().contains("expected 0"));
        assert!(err.to_string().contains("got 2"));
    }

    #[test]
    fn self_intersection_error_display() {
        let err = ValidationError::SelfIntersection {
            face_a: Default::default(),
            face_b: Default::default(),
        };
        assert!(err.to_string().contains("Self-intersection"));
    }
}
