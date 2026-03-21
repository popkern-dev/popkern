//! Error types for popkern-core.

use thiserror::Error;

/// Errors that can occur in popkern-core operations.
#[derive(Debug, Error)]
pub enum CoreError {
    /// A validation check failed.
    #[error("validation failed: {0}")]
    ValidationError(String),
}
