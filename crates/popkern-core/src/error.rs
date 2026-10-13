//! Error types for popkern-core.

use crate::topology::{
    EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, SolidHandle, VertexHandle,
};
use thiserror::Error;

/// Errors from geometric/topological validation checks.
#[derive(Debug, Error)]
pub enum ValidationError {
    /// A topology entity references another that no longer exists or
    /// does not own it (loop→face, half-edge→loop).
    #[error("dangling topology reference: {what}")]
    DanglingReference {
        /// Description of the broken reference.
        what: String,
    },
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

    /// A closed edge's canonical half-edge traverses the stored
    /// curve backward (the PK-303 canonical-forward convention,
    /// audited by `audit_canonical_forward` — PK-478).
    #[error("Edge {edge:?}: canonical half-edge runs the stored curve backward")]
    CanonicalBackward {
        /// The violating closed edge.
        edge: EdgeHandle,
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

    /// A face's stored orientation flag contradicts its outer loop's
    /// UV winding (interior-left invariant, the adjacent-face
    /// consistency half of check 7). Consumers folding the flag —
    /// edge convexity, mesher triangle winding, STEP `same_sense`,
    /// measured volume sign — mirror against neighbours on such a
    /// face.
    #[error(
        "Winding-dishonest face {face:?}: outer loop signed UV area {signed_uv_area:.3e} \
         contradicts orientation flag {orientation}"
    )]
    WindingDishonest {
        /// The face whose flag lies about its winding.
        face: FaceHandle,
        /// Signed UV area of the outer loop in topological order
        /// (positive = counter-clockwise).
        signed_uv_area: f64,
        /// The stored orientation flag.
        orientation: bool,
    },

    /// A face with surface geometry has loop half-edges without
    /// pcurves. Trim domains, loop-winding audits, and tessellation
    /// all read pcurves; a pcurve-less face silently degrades every
    /// one of them (measured: thicken served pcurve-less bodies for
    /// an unknown period — the winding audit deferred to its
    /// tangency witness and the mesher fell back to untrimmed
    /// full-grid faces; PK-331).
    #[error(
        "Face {face:?} carries surface geometry but {missing} of its {total} loop \
         half-edges have no pcurve"
    )]
    MissingPcurves {
        /// The face whose trim is unreadable.
        face: FaceHandle,
        /// Loop half-edges without a pcurve.
        missing: usize,
        /// Total loop half-edges on the face.
        total: usize,
    },

    /// A face's loops leave a net winding around a periodic chart with
    /// no collapsed row, so they bound no region of it: some edge runs
    /// the long way round the period (PK-792).
    #[error(
        "Face {face:?}: its loops wind {winding_u} times round the u period and \
         {winding_v} times round the v period instead of cancelling"
    )]
    LoopPeriodWinding {
        /// The face whose trim bounds no region.
        face: FaceHandle,
        /// Net winding of all the face's loops round the u period.
        winding_u: i64,
        /// Net winding of all the face's loops round the v period.
        winding_v: i64,
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
