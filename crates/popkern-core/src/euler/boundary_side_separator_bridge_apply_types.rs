//! Helper types for boundary-side separator bridge mutation.

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

use super::BoundarySideSeparatorBridgePlan;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SeparatorHalfEdges {
    pub(super) edge: EdgeHandle,
    pub(super) forward: HalfEdgeHandle,
    pub(super) reverse: HalfEdgeHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CycleSide {
    First,
    Second,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SeparatorCounts {
    pub(super) vertices: i64,
    pub(super) edges: i64,
    pub(super) faces: i64,
    pub(super) loops: i64,
}

impl SeparatorCounts {
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

pub(super) fn separator_half_open_path(
    interval: &[HalfEdgeHandle],
) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "boundary-side separator bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}

pub(super) fn required_separator_path_handle(
    handle: Option<&HalfEdgeHandle>,
    path: &'static str,
    end: &'static str,
) -> Result<HalfEdgeHandle, CoreError> {
    handle.copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "{path} boundary-side separator bridge path {end} is empty"
        ))
    })
}

pub(super) fn wire_separator_path(brep: &mut BRep, path: &[HalfEdgeHandle]) {
    for pair in path.windows(2) {
        brep.half_edges[pair[0]].next = pair[1];
        brep.half_edges[pair[1]].prev = pair[0];
    }
}

/// The one same-loop cycle of `plan` whose two intervals both run
/// through the separator source (the cycle the separator duplicates).
pub(super) fn duplicate_cycle_index(
    plan: &BoundarySideSeparatorBridgePlan,
) -> Result<usize, CoreError> {
    let matches = plan
        .cycles
        .iter()
        .enumerate()
        .filter(|(_, cycle)| {
            cycle.first_loop == cycle.second_loop
                && interval_uses_separator_source(plan, &cycle.first_interval).unwrap_or(false)
                && interval_uses_separator_source(plan, &cycle.second_interval).unwrap_or(false)
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if matches.len() == 1 {
        return Ok(matches[0]);
    }
    Err(CoreError::PreconditionFailed(format!(
        "boundary-side separator bridge expected one duplicate same-loop cycle, got {}",
        matches.len()
    )))
}

/// Whether `interval` runs through the separator source half-edge.
pub(super) fn interval_uses_separator_source(
    plan: &BoundarySideSeparatorBridgePlan,
    interval: &[HalfEdgeHandle],
) -> Result<bool, CoreError> {
    Ok(separator_half_open_path(interval)?.contains(&plan.separator_source_half_edge))
}
