//! Error types for tessellation and mesh export.

use thiserror::Error;

/// Errors that can occur during tessellation or mesh export.
#[derive(Debug, Error)]
pub enum MeshError {
    /// Tessellation of a face or solid failed.
    #[error("tessellation failed: {0}")]
    TessellationFailed(String),
    /// A face has no surface geometry and cannot be tessellated.
    #[error("face has no surface geometry")]
    NoSurfaceGeometry,
    /// Serialisation of the exported mesh format failed.
    #[error("export failed: {0}")]
    ExportFailed(String),
    /// An I/O error occurred while writing mesh output.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// A cancellation token fired (PK-276, §17-style clean
    /// termination): tessellation stopped between rounds/stages;
    /// nothing was produced and the input BRep is untouched.
    #[error("tessellation cancelled by token")]
    Cancelled,
}
