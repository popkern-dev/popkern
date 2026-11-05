//! Shared helpers of the bridge operators (two-loop bridge, face-family
//! bridge): the entity counts a bridge's Euler delta is checked against,
//! and the open path of a bridge interval.

use crate::error::CoreError;
use crate::topology::{BRep, HalfEdgeHandle};

/// Vertex, edge, face and loop counts, captured before a bridge so its
/// commit can be checked against the planned delta.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct TopologyCounts {
    pub(super) vertices: i64,
    pub(super) edges: i64,
    pub(super) faces: i64,
    pub(super) loops: i64,
}

impl TopologyCounts {
    pub(super) fn capture(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len() as i64,
            edges: brep.edges.len() as i64,
            faces: brep.faces.len() as i64,
            loops: brep.loops.len() as i64,
        }
    }

    pub(super) fn delta_from(self, before: Self) -> Self {
        Self {
            vertices: self.vertices - before.vertices,
            edges: self.edges - before.edges,
            faces: self.faces - before.faces,
            loops: self.loops - before.loops,
        }
    }
}

/// An interval's half-edges without its last one (whose origin is the
/// far endpoint). `bridge` names the operator in the decline.
pub(super) fn half_open_path<'a>(
    interval: &'a [HalfEdgeHandle],
    bridge: &str,
) -> Result<&'a [HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(format!(
            "{bridge} bridge path needs distinct endpoint origins"
        )));
    }
    Ok(&interval[..interval.len() - 1])
}

/// A path's first or last half-edge, declining an empty path. `path`
/// names which path, `bridge` the operator, `end` the end.
pub(super) fn required_path_handle(
    handle: Option<&HalfEdgeHandle>,
    path: &str,
    bridge: &str,
    end: &str,
) -> Result<HalfEdgeHandle, CoreError> {
    handle.copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!("{path} {bridge} path {end} is empty"))
    })
}
