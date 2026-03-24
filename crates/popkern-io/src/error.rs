//! Error types for popkern-io.

use thiserror::Error;

/// Top-level error type for IO operations.
#[derive(Debug, Error)]
pub enum IoError {
    /// Filesystem I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Parse error at a specific location.
    #[error("Parse error at line {line}, col {col}: {message}")]
    Parse {
        /// Line number (1-based).
        line: usize,
        /// Column number (1-based).
        col: usize,
        /// Description of the error.
        message: String,
    },

    /// The file's schema is not supported.
    #[error("Unsupported schema: {0}")]
    UnsupportedSchema(String),

    /// Import completed but had errors (see diagnostics).
    #[error("Import failed with {0} errors (see diagnostics)")]
    ImportFailed(usize),

    /// Export failed.
    #[error("Export error: {0}")]
    Export(String),

    /// Validation of the imported model failed.
    #[error("Validation failed: {0}")]
    Validation(#[from] popkern_core::error::CoreError),
}
