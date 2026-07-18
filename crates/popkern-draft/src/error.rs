//! Error and warning types for draft/taper operations.

use popkern_core::error::CoreError;
use popkern_core::topology::FaceHandle;
use thiserror::Error;

/// Errors that can occur during draft operations.
#[derive(Debug, Error)]
pub enum DraftError {
    /// No draft specifications were provided.
    #[error("no draft specifications provided")]
    EmptySpecs,

    /// The specified face handle is invalid or not present in the BRep.
    #[error("invalid face: {0}")]
    InvalidFace(String),

    /// The draft angle exceeds the maximum allowed value.
    #[error("angle too large: {angle} rad exceeds max {max} rad")]
    AngleTooLarge {
        /// The requested angle in radians.
        angle: f64,
        /// The maximum allowed angle in radians.
        max: f64,
    },

    /// The neutral reference could not be resolved on the face.
    #[error("invalid neutral reference: {0}")]
    InvalidNeutral(String),

    /// The parting line could not be computed or is degenerate.
    #[error("invalid parting line: {0}")]
    InvalidPartingLine(String),

    /// An edge or surface extension operation failed.
    #[error("extension failed: {0}")]
    ExtensionFailed(String),

    /// A topological inconsistency was encountered.
    #[error("topology error: {0}")]
    Topology(String),

    /// An error propagated from popkern-core.
    #[error("core: {0}")]
    Core(#[from] CoreError),
}

/// Non-fatal warnings produced during draft operations.
#[derive(Debug, Clone)]
pub enum DraftWarning {
    /// Adjacent face mode fell back to transition surface generation.
    AdjacentFallbackToTransition {
        /// The face where fallback occurred.
        face: FaceHandle,
    },
    /// A variable draft angle was clamped to the maximum allowed value.
    VariableAngleClamped {
        /// The face where clamping occurred.
        face: FaceHandle,
        /// The original requested angle in radians.
        requested: f64,
        /// The clamped angle in radians.
        clamped: f64,
    },
    /// An edge adjacent to a drafted face kept its ORIGINAL curve
    /// because the adjacent surface pair has no recompute path yet —
    /// the curve may no longer lie on the drafted surface.
    EdgeCurveStale {
        /// The edge whose curve was kept.
        edge: popkern_core::topology::EdgeHandle,
    },
    /// A vertex adjacent to recomputed edge curves snapped to
    /// positions that disagree beyond tolerance between curves.
    VertexSnapDisagreement {
        /// The vertex in question.
        vertex: popkern_core::topology::VertexHandle,
        /// Largest disagreement between per-curve snap targets.
        deviation: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn draft_error_display_empty_specs() {
        let err = DraftError::EmptySpecs;
        assert_eq!(err.to_string(), "no draft specifications provided");
    }

    #[test]
    fn draft_error_display_invalid_face() {
        let err = DraftError::InvalidFace("face not found".into());
        assert!(err.to_string().contains("face not found"));
    }

    #[test]
    fn draft_error_display_angle_too_large() {
        let err = DraftError::AngleTooLarge {
            angle: 1.0,
            max: 0.785,
        };
        assert!(err.to_string().contains("1"));
        assert!(err.to_string().contains("0.785"));
    }

    #[test]
    fn draft_error_display_invalid_neutral() {
        let err = DraftError::InvalidNeutral("no intersection".into());
        assert!(err.to_string().contains("no intersection"));
    }

    #[test]
    fn draft_error_display_invalid_parting_line() {
        let err = DraftError::InvalidPartingLine("degenerate".into());
        assert!(err.to_string().contains("degenerate"));
    }

    #[test]
    fn draft_error_display_extension_failed() {
        let err = DraftError::ExtensionFailed("trim failed".into());
        assert!(err.to_string().contains("trim failed"));
    }

    #[test]
    fn draft_error_display_topology() {
        let err = DraftError::Topology("broken edge".into());
        assert!(err.to_string().contains("broken edge"));
    }

    #[test]
    fn draft_error_from_core() {
        let core_err = CoreError::InvalidHandle("test".into());
        let err: DraftError = core_err.into();
        assert!(matches!(err, DraftError::Core(_)));
    }

    #[test]
    fn draft_warning_variants() {
        let w1 = DraftWarning::AdjacentFallbackToTransition {
            face: Default::default(),
        };
        let w2 = DraftWarning::VariableAngleClamped {
            face: Default::default(),
            requested: 1.0,
            clamped: 0.785,
        };
        let _ = format!("{w1:?}");
        let _ = format!("{w2:?}");
    }
}
