//! Inner skin construction: clone, map, and geometry replacement.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// Clone the outer BRep to create an inner skin, building a handle map.
///
/// Faces listed in `removed_faces` are deleted from the inner skin
/// (along with their exclusive loops). Shared edges and vertices are kept.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if the clone or removal encounters
/// inconsistent topology.
pub fn build_inner_skin(
    _outer: &BRep,
    _removed_faces: &[FaceHandle],
) -> Result<(BRep, SkinMap), ShellError> {
    Err(ShellError::Topology("not yet implemented".into()))
}

/// Replace inner skin face geometry with offset surfaces and flip normals.
///
/// For each inner face, looks up the corresponding outer face via the
/// reverse skin map, retrieves the offset surface, and replaces the
/// face's surface geometry. The orientation is flipped so the inner
/// skin normal points inward.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if a face mapping is missing.
pub fn apply_offset_geometry(
    _inner: &mut BRep,
    _offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    _skin_map: &SkinMap,
) -> Result<(), ShellError> {
    Err(ShellError::Topology("not yet implemented".into()))
}
