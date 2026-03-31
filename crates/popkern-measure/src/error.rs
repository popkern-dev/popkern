//! Error types for geometric measurement queries.

use thiserror::Error;

/// Errors that can occur during geometric measurement.
#[derive(Debug, Error)]
pub enum MeasureError {
    /// Entity handle is invalid (stale after topology change).
    #[error("invalid handle: {0}")]
    InvalidHandle(String),
    /// Measurement requires a specific surface type (e.g., radius on a plane).
    #[error("incompatible surface: expected {expected}, got {got}")]
    IncompatibleSurface {
        /// The surface type required by this measurement.
        expected: &'static str,
        /// The surface type actually found.
        got: String,
    },
    /// Degenerate geometry that cannot be measured meaningfully.
    #[error("degenerate geometry: {0}")]
    DegenerateGeometry(String),
    /// Tessellation failed during approximate computation.
    #[error("tessellation failed: {0}")]
    TessellationFailed(String),
    /// Face has no surface geometry assigned.
    #[error("face has no surface geometry")]
    MissingSurface,
    /// Edge has no curve geometry assigned.
    #[error("edge has no curve geometry")]
    MissingCurve,
}
