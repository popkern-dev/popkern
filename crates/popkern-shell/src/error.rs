//! Error types for shell and offset operations.

use popkern_core::error::{CoreError, ValidationError};
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
    #[error(
        "offset {requested} exceeds the maximum {max_safe:.6} supported by the \
         geometry ({binding})"
    )]
    OffsetTooLarge {
        /// The offset/thickness the caller asked for.
        requested: f64,
        /// Estimated largest magnitude this geometry supports.
        max_safe: f64,
        /// Human-readable description of the binding entity (the face or
        /// edge that constrains the offset).
        binding: String,
    },

    /// A face's offset could only be approximated and the approximation
    /// missed the requested tolerance even at the finest sampling grid.
    #[error(
        "offset of face {face:?} could not be approximated within tolerance: \
         achieved {achieved:.3e}, requested {requested:.3e}"
    )]
    ApproximationOutOfTolerance {
        /// The face whose offset approximation missed tolerance.
        face: FaceHandle,
        /// Measured maximum deviation of the approximation.
        achieved: f64,
        /// The tolerance the approximation had to meet.
        requested: f64,
    },

    /// The assembled result failed structural validation; serving it
    /// would hand the caller a broken body, so the operation declines
    /// instead (the input model is untouched).
    #[error(
        "result failed structural validation with {} error(s); first: {}",
        errors.len(),
        errors.first().map(|e| e.to_string()).unwrap_or_default()
    )]
    ValidationFailed {
        /// The structural validation errors found on the assembled result.
        errors: Vec<ValidationError>,
    },

    /// The material side of a component could not be measured, so
    /// offset directions have no meaning for it.
    #[error(
        "cannot resolve which side of the body is material (component \
         containing face {face:?}); offset direction is undefined"
    )]
    MaterialSenseUnresolved {
        /// A seed face of the unresolvable component.
        face: FaceHandle,
    },

    /// The input sheet has free boundary half-edges and the operation
    /// cannot yet close them with side walls. Temporary containment:
    /// open-rim closure lands with the topology-changing join stage
    /// (docs/architecture/solid-offset-shell-design.md §5, PK-325).
    #[error(
        "input sheet has {free_half_edges} free boundary half-edge(s); \
         closing open rims into side walls is not yet supported"
    )]
    OpenBoundary {
        /// Number of half-edges with no owning loop.
        free_half_edges: usize,
    },

    /// A topological inconsistency was encountered.
    #[error("topology: {0}")]
    Topology(String),

    /// An error propagated from popkern-core.
    #[error("core: {0}")]
    Core(#[from] CoreError),

    /// A cancellation token fired (PK-276, §17-style clean
    /// termination): the operation stopped between stages and the
    /// caller's model is unchanged — a cancelled shell is
    /// indistinguishable from one that failed before writing.
    #[error("shell operation cancelled by token")]
    Cancelled,
}

/// Non-fatal warnings produced during shell operations.
#[derive(Debug, Clone)]
pub enum ShellWarning {
    /// A face's offset surface was approximated (e.g. NURBS fallback).
    FaceOffsetApproximated {
        /// The face whose offset was approximated.
        face: FaceHandle,
        /// Measured maximum deviation of the approximation from the
        /// exact offset (mid-sample check at the final sampling grid).
        achieved_deviation: f64,
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
        let err = ShellError::OffsetTooLarge {
            requested: 5.0,
            max_safe: 3.2,
            binding: "edge between Face#3 and Face#7".into(),
        };
        let msg = err.to_string();
        assert!(msg.contains('5') && msg.contains("3.2"), "{msg}");
        assert!(msg.contains("Face#3"), "{msg}");
    }

    #[test]
    fn shell_error_validation_failed_reports_first() {
        let err = ShellError::ValidationFailed {
            errors: vec![ValidationError::DanglingReference {
                what: "test ref".into(),
            }],
        };
        let msg = err.to_string();
        assert!(
            msg.contains("1 error(s)") && msg.contains("test ref"),
            "{msg}"
        );
    }

    #[test]
    fn shell_error_open_boundary() {
        let err = ShellError::OpenBoundary { free_half_edges: 4 };
        assert!(err.to_string().contains("4 free boundary"));
    }

    #[test]
    fn shell_error_approximation_out_of_tolerance() {
        let err = ShellError::ApproximationOutOfTolerance {
            face: Default::default(),
            achieved: 1e-3,
            requested: 1e-6,
        };
        let msg = err.to_string();
        assert!(
            msg.contains("1.000e-3") && msg.contains("1.000e-6"),
            "{msg}"
        );
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
            achieved_deviation: 1e-7,
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
