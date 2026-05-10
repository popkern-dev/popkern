//! Loop-walking helpers for aggregate bridge boundary exposure.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::FaceFamilyBridgePrimitivePlan;

pub(super) fn aggregate_bridge_half_edges(
    primitive: &FaceFamilyBridgePrimitivePlan,
) -> BTreeSet<HalfEdgeHandle> {
    primitive
        .regions()
        .iter()
        .flat_map(|region| {
            region.boundary_cycles().iter().flat_map(|cycle| {
                cycle
                    .support_interval()
                    .iter()
                    .chain(cycle.outside_interval().iter())
            })
        })
        .copied()
        .collect()
}

pub(super) fn bridge_half_edge_is_boundary_ready(
    candidate: &BRep,
    half_edge: HalfEdgeHandle,
) -> bool {
    let Some(data) = candidate.half_edges.get(half_edge) else {
        return false;
    };
    candidate.boundary_edges.contains(&data.edge) && !candidate.half_edges.contains_key(data.twin)
}

pub(super) fn live_face_outer_half_edges(
    candidate: &BRep,
    face: FaceHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let face_data = candidate.faces.get(face)?;
    let loop_data = candidate.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut half_edges = Vec::new();
    for _ in 0..=candidate.half_edges.len() {
        let half_edge = candidate.half_edges.get(cursor)?;
        if half_edge.loop_ != face_data.outer_loop {
            return None;
        }
        half_edges.push(cursor);
        cursor = half_edge.next;
        if cursor == start {
            return Some(half_edges);
        }
    }
    None
}

pub(super) fn live_contiguous_fragments(
    face_half_edges: &[HalfEdgeHandle],
    blocked: &BTreeSet<HalfEdgeHandle>,
) -> Vec<(HalfEdgeHandle, HalfEdgeHandle)> {
    if blocked.is_empty() || blocked.len() == face_half_edges.len() {
        return Vec::new();
    }
    let is_blocked = face_half_edges
        .iter()
        .map(|half_edge| blocked.contains(half_edge))
        .collect::<Vec<_>>();
    let mut starts = Vec::new();
    for i in 0..is_blocked.len() {
        let prev = if i == 0 { is_blocked.len() - 1 } else { i - 1 };
        if is_blocked[i] && !is_blocked[prev] {
            starts.push(i);
        }
    }
    let mut fragments = Vec::new();
    for first_index in starts {
        let mut after_last_index = (first_index + 1) % is_blocked.len();
        while is_blocked[after_last_index] {
            after_last_index = (after_last_index + 1) % is_blocked.len();
            if after_last_index == first_index {
                return Vec::new();
            }
        }
        fragments.push((
            face_half_edges[first_index],
            face_half_edges[after_last_index],
        ));
    }
    fragments
}

pub(super) fn is_non_blocking_exposure_error(message: &str) -> bool {
    message.contains("already contains a boundary edge")
        || message.contains("after-last half-edge is invalid")
        || (message.contains("Half-edge cycle broken") && message.contains("LoopHandle(null)"))
}
