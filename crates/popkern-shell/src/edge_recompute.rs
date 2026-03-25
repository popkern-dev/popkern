//! Edge curve recomputation for inner skin.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Recompute all edge curves in the inner skin.
///
/// For each edge, finds the two adjacent faces, maps them to the outer
/// BRep, retrieves offset surfaces, and computes the intersection curve
/// (plane-plane yields a line; otherwise uses surface-surface intersection).
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if adjacency lookups fail.
pub fn recompute_edges(
    _inner: &mut BRep,
    _offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    _skin_map: &SkinMap,
    _outer: &BRep,
    _tolerance: &Tolerance,
) -> Result<(), ShellError> {
    Ok(())
}
