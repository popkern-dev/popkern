//! Surface offset computation for shell operations.

use std::collections::HashMap;

use crate::classify::FaceClassification;
use crate::error::{ShellError, ShellWarning};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of computing offset surfaces for classified faces.
#[derive(Debug, Clone)]
pub struct OffsetResult {
    /// Successfully offset surfaces, keyed by outer face handle.
    pub surfaces: HashMap<FaceHandle, SurfaceGeometry>,
    /// Faces that failed to offset.
    pub failed: Vec<(FaceHandle, String)>,
    /// Warnings generated during offset computation.
    pub warnings: Vec<ShellWarning>,
    /// Faces auto-removed due to offset failure (when allow_partial is true).
    pub auto_removed: Vec<FaceHandle>,
}

/// Compute offset surfaces for all classified faces.
///
/// For each face in `classification.offset_faces`, offsets the surface
/// geometry by the specified distance. Failed faces are either collected
/// as errors or auto-removed depending on `allow_partial`.
///
/// # Errors
///
/// Returns [`ShellError::AllFacesFailed`] if every face fails to offset.
pub fn compute_offset_surfaces(
    _brep: &BRep,
    _classification: &FaceClassification,
    _allow_partial: bool,
    _tolerance: &Tolerance,
) -> Result<OffsetResult, ShellError> {
    Err(ShellError::NoFaces)
}
