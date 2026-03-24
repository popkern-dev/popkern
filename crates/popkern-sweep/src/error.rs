//! Error types for sweep and loft operations.

use popkern_core::error::CoreError;
use thiserror::Error;

/// Fatal error during sweep/loft.
#[derive(Debug, Error)]
pub enum SweepError {
    /// The profile is invalid for this operation.
    #[error("invalid profile: {0}")]
    InvalidProfile(String),
    /// The path curve is invalid for this operation.
    #[error("invalid path: {0}")]
    InvalidPath(String),
    /// The sweep surface self-intersects at the given parameter.
    #[error("self-intersection at t={t:.4}")]
    SelfIntersection {
        /// Parameter along the path where the intersection occurs.
        t: f64,
    },
    /// Frame computation failed at the given parameter.
    #[error("frame failed at t={t:.4}: {reason}")]
    FrameFailed {
        /// Parameter along the path where the frame failed.
        t: f64,
        /// Description of the failure.
        reason: String,
    },
    /// Loft sections are incompatible (different topology or vertex counts).
    #[error("incompatible sections: {0}")]
    IncompatibleSections(String),
    /// A guide curve deviates too far from the sweep surface.
    #[error("guide deviation: {deviation:.6}")]
    GuideDeviation {
        /// Maximum deviation distance.
        deviation: f64,
    },
    /// A topology construction error.
    #[error("topology: {0}")]
    Topology(String),
    /// A boolean operation error.
    #[error("boolean: {0}")]
    Boolean(String),
    /// An error from popkern-core.
    #[error("core: {0}")]
    Core(#[from] CoreError),
}

/// Non-fatal warning emitted during sweep/loft.
#[derive(Debug, Clone)]
pub enum SweepWarning {
    /// Self-intersection was detected but allowed by options.
    SelfIntersectionAllowed {
        /// Parameter along the path where the intersection occurs.
        t: f64,
    },
    /// Frenet frame was undefined; fell back to Bishop at this parameter.
    FrenetFallback {
        /// Parameter where fallback occurred.
        t: f64,
    },
    /// Helix was approximated (not exact NURBS representation).
    HelixApproximated,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sweep_error_display() {
        let err = SweepError::InvalidProfile("empty wire".into());
        assert!(err.to_string().contains("empty wire"));
    }

    #[test]
    fn self_intersection_display() {
        let err = SweepError::SelfIntersection { t: 0.5 };
        assert!(err.to_string().contains("0.5000"));
    }

    #[test]
    fn frame_failed_display() {
        let err = SweepError::FrameFailed {
            t: 0.25,
            reason: "zero curvature".into(),
        };
        let msg = err.to_string();
        assert!(msg.contains("0.2500"));
        assert!(msg.contains("zero curvature"));
    }

    #[test]
    fn guide_deviation_display() {
        let err = SweepError::GuideDeviation { deviation: 0.001 };
        assert!(err.to_string().contains("0.001000"));
    }

    #[test]
    fn core_error_conversion() {
        let core_err = CoreError::PreconditionFailed("test".into());
        let sweep_err: SweepError = core_err.into();
        assert!(matches!(sweep_err, SweepError::Core(_)));
    }

    #[test]
    fn sweep_warning_clone() {
        let w = SweepWarning::FrenetFallback { t: 0.5 };
        let w2 = w.clone();
        assert!(matches!(w2, SweepWarning::FrenetFallback { t } if (t - 0.5).abs() < f64::EPSILON));
    }
}
