//! Post-exposure cylinder resource loop partitioning.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle};

use super::super::apply_cylinder_resources::CylinderBridgeResourceSegment;

pub(super) fn partition_remapped_cylinder_loop(
    brep: &BRep,
    loop_: LoopHandle,
    membership: &BTreeMap<HalfEdgeHandle, BTreeSet<usize>>,
) -> Result<Vec<CylinderBridgeResourceSegment>, CoreError> {
    let start = brep
        .loops
        .get(loop_)
        .map(|loop_data| loop_data.half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "remapped cylinder resource loop {loop_:?} is stale"
            ))
        })?;
    let mut segments = Vec::new();
    let mut current_half_edges = Vec::new();
    let mut current_edges = Vec::new();
    let mut current_request_indices = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let data = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "remapped cylinder resource partition references stale half-edge {cursor:?}"
            ))
        })?;
        if data.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(format!(
                "remapped cylinder resource traversal left loop {loop_:?} at {cursor:?}"
            )));
        }
        match membership.get(&cursor) {
            Some(request_indices) => {
                let request_indices = request_indices.iter().copied().collect::<Vec<_>>();
                if current_half_edges.is_empty() || current_request_indices == request_indices {
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                } else {
                    push_remapped_cylinder_resource_segment(
                        &mut segments,
                        &mut current_half_edges,
                        &mut current_edges,
                        &mut current_request_indices,
                    );
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                }
            }
            None => push_remapped_cylinder_resource_segment(
                &mut segments,
                &mut current_half_edges,
                &mut current_edges,
                &mut current_request_indices,
            ),
        }
        cursor = data.next;
        if cursor == start {
            break;
        }
    }
    if cursor != start {
        return Err(CoreError::PreconditionFailed(format!(
            "remapped cylinder resource partition could not close loop {loop_:?}"
        )));
    }
    push_remapped_cylinder_resource_segment(
        &mut segments,
        &mut current_half_edges,
        &mut current_edges,
        &mut current_request_indices,
    );
    Ok(segments)
}

fn push_remapped_cylinder_resource_segment(
    segments: &mut Vec<CylinderBridgeResourceSegment>,
    current_half_edges: &mut Vec<HalfEdgeHandle>,
    current_edges: &mut Vec<EdgeHandle>,
    current_request_indices: &mut Vec<usize>,
) {
    if current_half_edges.is_empty() {
        return;
    }
    segments.push(CylinderBridgeResourceSegment {
        half_edges: std::mem::take(current_half_edges),
        edges: std::mem::take(current_edges),
        request_indices: std::mem::take(current_request_indices),
    });
}
