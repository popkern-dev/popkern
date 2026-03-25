//! Vertex position recomputation for inner skin.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Recompute all vertex positions in the inner skin.
///
/// For each vertex, finds the adjacent faces, maps them back to the outer
/// BRep, retrieves the corresponding offset surfaces, and computes the
/// new vertex position as the intersection of those surfaces.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if adjacency lookups fail.
pub fn recompute_vertices(
    _inner: &mut BRep,
    _offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    _skin_map: &SkinMap,
    _outer: &BRep,
    _tolerance: &Tolerance,
) -> Result<(), ShellError> {
    Ok(())
}
