//! Error types for popkern-viz.

use thiserror::Error;

/// Top-level error for the visualization crate.
#[derive(Debug, Error)]
pub enum VizError {
    /// Failed to load a model file.
    #[error("load error: {0}")]
    Load(String),
    /// GPU / rendering error.
    #[error("render error: {0}")]
    Render(String),
    /// IO error.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// JSON serialization/deserialization error.
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
    /// Software renderer error.
    #[error("software render error: {0}")]
    SoftwareRender(String),
}
