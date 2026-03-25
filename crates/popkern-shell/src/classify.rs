//! Face classification for shell operations.

use crate::error::ShellError;
use crate::types::ShellDirection;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of classifying faces for a shell operation.
#[derive(Debug, Clone)]
pub struct FaceClassification {
    /// Faces to be offset: (face_handle, signed_distance, inward_flag).
    pub offset_faces: Vec<(FaceHandle, f64, bool)>,
    /// Faces to be removed (open faces of the shell).
    pub removed_faces: Vec<FaceHandle>,
}

/// Classify faces for a shell operation.
///
/// Each face is either marked for removal or assigned an offset distance
/// and direction. The `face_thicknesses` map overrides `default_thickness`
/// for specific faces.
///
/// # Errors
///
/// Returns [`ShellError::NoFaces`] if the BRep has no faces.
pub fn classify_faces(
    _brep: &BRep,
    _default_thickness: f64,
    _face_thicknesses: &[(FaceHandle, f64)],
    _direction: ShellDirection,
    _removed: &[FaceHandle],
) -> Result<FaceClassification, ShellError> {
    Err(ShellError::NoFaces)
}
