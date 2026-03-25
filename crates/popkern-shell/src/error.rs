//! Error types for shell and offset operations.

use popkern_core::error::CoreError;
use popkern_core::topology::FaceHandle;
use thiserror::Error;

/// Errors that can occur during shell/offset operations.
#[derive(Debug, Error)]
pub enum ShellError {
    /// No faces available to offset.
    #[error("no faces to offset")]
    NoFaces,

    /// Every face failed during offset computation.
    #[error("all faces failed")]
    AllFacesFailed {
        /// Per-face failure details.
        failures: Vec<(FaceHandle, String)>,
    },

    /// The requested offset distance exceeds what the geometry can support.
    #[error("offset too large: {0}")]
    OffsetTooLarge(String),

    /// A topological inconsistency was encountered.
    #[error("topology: {0}")]
    Topology(String),

    /// An error propagated from popkern-core.
    #[error("core: {0}")]
    Core(#[from] CoreError),
}

/// Non-fatal warnings produced during shell operations.
#[derive(Debug, Clone)]
pub enum ShellWarning {
    /// A face's offset surface was approximated (e.g. NURBS fallback).
    FaceOffsetApproximated {
        /// The face whose offset was approximated.
        face: FaceHandle,
    },
    /// A face was automatically removed because its offset failed.
    FaceAutoRemoved {
        /// The face that was removed.
        face: FaceHandle,
        /// Reason for removal.
        reason: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shell_error_display() {
        let err = ShellError::NoFaces;
        assert_eq!(err.to_string(), "no faces to offset");
    }

    #[test]
    fn shell_error_all_faces_failed() {
        let err = ShellError::AllFacesFailed { failures: vec![] };
        assert!(err.to_string().contains("all faces failed"));
    }

    #[test]
    fn shell_error_offset_too_large() {
        let err = ShellError::OffsetTooLarge("radius exceeded".into());
        assert!(err.to_string().contains("radius exceeded"));
    }

    #[test]
    fn shell_error_topology() {
        let err = ShellError::Topology("broken edge".into());
        assert!(err.to_string().contains("broken edge"));
    }

    #[test]
    fn shell_error_from_core() {
        let core_err = CoreError::InvalidHandle("test".into());
        let err: ShellError = core_err.into();
        assert!(matches!(err, ShellError::Core(_)));
    }

    #[test]
    fn shell_warning_variants() {
        let w1 = ShellWarning::FaceOffsetApproximated {
            face: Default::default(),
        };
        let w2 = ShellWarning::FaceAutoRemoved {
            face: Default::default(),
            reason: "test".into(),
        };
        // Just verify they can be constructed and debug-printed.
        let _ = format!("{w1:?}");
        let _ = format!("{w2:?}");
    }
}
