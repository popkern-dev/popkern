use popkern_core::topology::{BRep, HalfEdgeHandle, LoopHandle};

use super::FaceFamilyBridgeLoopIntervals;

pub(super) fn loop_intervals_between_origins(
    brep: &BRep,
    loop_handle: LoopHandle,
    first: HalfEdgeHandle,
    second: HalfEdgeHandle,
) -> Option<FaceFamilyBridgeLoopIntervals> {
    let loop_len = loop_half_edge_count(brep, loop_handle, first)?;
    Some(FaceFamilyBridgeLoopIntervals {
        first_to_second: loop_path_between_origins(brep, loop_handle, first, second)?,
        second_to_first: loop_path_between_origins(brep, loop_handle, second, first)?,
        loop_len,
    })
}

fn loop_half_edge_count(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
) -> Option<usize> {
    let mut cursor = start;
    let mut count = 0;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        count += 1;
        cursor = half_edge.next;
        if cursor == start {
            return Some(count);
        }
    }
    None
}

fn loop_path_between_origins(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
    end: HalfEdgeHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let mut cursor = start;
    let mut path = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        path.push(cursor);
        if cursor == end {
            return Some(path);
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    None
}

pub(super) fn loops_share_edges(brep: &BRep, first: LoopHandle, second: LoopHandle) -> bool {
    let mut first_edges = std::collections::BTreeSet::new();
    let mut cursor = brep.loops[first].half_edge;
    loop {
        first_edges.insert(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[first].half_edge {
            break;
        }
    }
    let mut cursor = brep.loops[second].half_edge;
    loop {
        if first_edges.contains(&brep.half_edges[cursor].edge) {
            return true;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[second].half_edge {
            break;
        }
    }
    false
}
