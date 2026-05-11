//! Shared-start bridge cycle shape and preflight diagnostics.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::euler::{FaceFamilyBridgeCycleRequest, SameLoopFaceBridgeRequest};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SharedStartCycleShape {
    pub(super) cycles: usize,
    pub(super) cross_loop_cycles: usize,
    pub(super) same_loop_cycles: usize,
    pub(super) same_loop_overlap_cycles: usize,
    pub(super) same_loop_overlapping_edges: usize,
    pub(super) same_loop_shared_start_overlaps: usize,
    pub(super) same_loop_shared_end_overlaps: usize,
    pub(super) same_loop_crossed_endpoint_overlaps: usize,
    pub(super) same_loop_interior_overlaps: usize,
}

impl SharedStartCycleShape {
    pub(super) fn empty() -> Self {
        Self {
            cycles: 0,
            cross_loop_cycles: 0,
            same_loop_cycles: 0,
            same_loop_overlap_cycles: 0,
            same_loop_overlapping_edges: 0,
            same_loop_shared_start_overlaps: 0,
            same_loop_shared_end_overlaps: 0,
            same_loop_crossed_endpoint_overlaps: 0,
            same_loop_interior_overlaps: 0,
        }
    }

    pub(super) fn from_cycles(cycles: &[FaceFamilyBridgeCycleRequest], brep: &BRep) -> Self {
        let mut same_loop_cycles = 0;
        let mut same_loop_overlap_cycles = 0;
        let mut same_loop_overlapping_edges = 0;
        let mut same_loop_shared_start_overlaps = 0;
        let mut same_loop_shared_end_overlaps = 0;
        let mut same_loop_crossed_endpoint_overlaps = 0;
        let mut same_loop_interior_overlaps = 0;
        for cycle in cycles
            .iter()
            .filter(|cycle| cycle.first_loop == cycle.second_loop)
        {
            same_loop_cycles += 1;
            let overlap = same_loop_overlap_summary(brep, cycle);
            let overlap_edges = overlap.edges;
            if overlap_edges > 0 {
                same_loop_overlap_cycles += 1;
                same_loop_overlapping_edges += overlap_edges;
                match overlap.position {
                    SameLoopOverlapPosition::SharedStart => same_loop_shared_start_overlaps += 1,
                    SameLoopOverlapPosition::SharedEnd => same_loop_shared_end_overlaps += 1,
                    SameLoopOverlapPosition::CrossedEndpoint => {
                        same_loop_crossed_endpoint_overlaps += 1
                    }
                    SameLoopOverlapPosition::Interior => same_loop_interior_overlaps += 1,
                    SameLoopOverlapPosition::None => {}
                }
            }
        }
        Self {
            cycles: cycles.len(),
            cross_loop_cycles: cycles.len() - same_loop_cycles,
            same_loop_cycles,
            same_loop_overlap_cycles,
            same_loop_overlapping_edges,
            same_loop_shared_start_overlaps,
            same_loop_shared_end_overlaps,
            same_loop_crossed_endpoint_overlaps,
            same_loop_interior_overlaps,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SharedStartCycleMixSummary {
    pub(super) cycles: usize,
    pub(super) cross_loop_cycles: usize,
    pub(super) same_loop_cycles: usize,
    pub(super) same_loop_overlap_cycles: usize,
    pub(super) same_loop_overlapping_edges: usize,
    pub(super) same_loop_shared_start_overlaps: usize,
    pub(super) same_loop_shared_end_overlaps: usize,
    pub(super) same_loop_crossed_endpoint_overlaps: usize,
    pub(super) same_loop_interior_overlaps: usize,
    pub(super) expected_edges: i64,
    pub(super) expected_faces: i64,
    pub(super) expected_loops: i64,
}

pub(super) fn preflight_shared_start_cycle_mix(
    brep: &BRep,
    cycles: Vec<FaceFamilyBridgeCycleRequest>,
) -> Result<SharedStartCycleMixSummary, CoreError> {
    let cycle_shape = SharedStartCycleShape::from_cycles(&cycles, brep);
    let mut cross_loop_cycles = Vec::new();
    let mut same_loop_cycles = Vec::new();
    for cycle in cycles {
        if cycle.first_loop == cycle.second_loop {
            same_loop_cycles.push(cycle);
        } else {
            cross_loop_cycles.push(cycle);
        }
    }
    let mut expected_edges = 0;
    let mut expected_faces = 0;
    let mut expected_loops = 0;
    let cross_loop_count = cross_loop_cycles.len();
    let same_loop_count = same_loop_cycles.len();
    if !cross_loop_cycles.is_empty() {
        let plan = brep.preflight_face_family_bridge(cross_loop_cycles)?;
        expected_edges += plan.expected_delta.edges;
        expected_faces += plan.expected_delta.faces;
        expected_loops += plan.expected_delta.loops;
    }
    for cycle in &same_loop_cycles {
        let plan = brep.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_: cycle.first_loop,
            first_interval: cycle.first_interval.clone(),
            second_interval: cycle.second_interval.clone(),
        })?;
        expected_edges += plan.expected_delta.edges;
        expected_faces += plan.expected_delta.faces;
        expected_loops += plan.expected_delta.loops;
    }
    Ok(SharedStartCycleMixSummary {
        cycles: cross_loop_count + same_loop_count,
        cross_loop_cycles: cross_loop_count,
        same_loop_cycles: same_loop_count,
        same_loop_overlap_cycles: cycle_shape.same_loop_overlap_cycles,
        same_loop_overlapping_edges: cycle_shape.same_loop_overlapping_edges,
        same_loop_shared_start_overlaps: cycle_shape.same_loop_shared_start_overlaps,
        same_loop_shared_end_overlaps: cycle_shape.same_loop_shared_end_overlaps,
        same_loop_crossed_endpoint_overlaps: cycle_shape.same_loop_crossed_endpoint_overlaps,
        same_loop_interior_overlaps: cycle_shape.same_loop_interior_overlaps,
        expected_edges,
        expected_faces,
        expected_loops,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SameLoopOverlapSummary {
    edges: usize,
    position: SameLoopOverlapPosition,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SameLoopOverlapPosition {
    None,
    SharedStart,
    SharedEnd,
    CrossedEndpoint,
    Interior,
}

fn same_loop_overlap_summary(
    brep: &BRep,
    cycle: &FaceFamilyBridgeCycleRequest,
) -> SameLoopOverlapSummary {
    let first_edges = match bridge_path_edges(brep, &cycle.first_interval) {
        Ok(edges) => edges,
        Err(_) => {
            return SameLoopOverlapSummary {
                edges: 0,
                position: SameLoopOverlapPosition::None,
            };
        }
    };
    let second_edges = match bridge_path_edges(brep, &cycle.second_interval) {
        Ok(edges) => edges,
        Err(_) => {
            return SameLoopOverlapSummary {
                edges: 0,
                position: SameLoopOverlapPosition::None,
            };
        }
    };
    let edges = first_edges.set.intersection(&second_edges.set).count();
    let position = if edges == 0 {
        SameLoopOverlapPosition::None
    } else if first_edges.first == second_edges.first {
        SameLoopOverlapPosition::SharedStart
    } else if first_edges.last == second_edges.last {
        SameLoopOverlapPosition::SharedEnd
    } else if first_edges.first == second_edges.last || first_edges.last == second_edges.first {
        SameLoopOverlapPosition::CrossedEndpoint
    } else {
        SameLoopOverlapPosition::Interior
    };
    SameLoopOverlapSummary { edges, position }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BridgePathEdges {
    set: BTreeSet<EdgeHandle>,
    first: Option<EdgeHandle>,
    last: Option<EdgeHandle>,
}

fn bridge_path_edges(
    brep: &BRep,
    interval: &[HalfEdgeHandle],
) -> Result<BridgePathEdges, CoreError> {
    let mut edges = BTreeSet::new();
    let mut first = None;
    let mut last = None;
    for &half_edge in half_open_bridge_path(interval)? {
        let edge = brep
            .half_edges
            .get(half_edge)
            .map(|data| data.edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle("shared-start bridge interval has stale half-edge".into())
            })?;
        if first.is_none() {
            first = Some(edge);
        }
        last = Some(edge);
        edges.insert(edge);
    }
    Ok(BridgePathEdges {
        set: edges,
        first,
        last,
    })
}

fn half_open_bridge_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "shared-start bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}
