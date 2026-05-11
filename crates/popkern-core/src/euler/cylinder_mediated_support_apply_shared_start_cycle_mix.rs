//! Shared-start bridge cycle shape and preflight diagnostics.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::euler::{FaceFamilyBridgeCycleRequest, SameLoopFaceBridgeRequest};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SharedStartCycleShape {
    pub(super) cycles: usize,
    pub(super) cross_loop_cycles: usize,
    pub(super) same_loop_cycles: usize,
    pub(super) one_edge_interior_overlap_candidates: usize,
    pub(super) repeated_one_edge_overlap_occurrences: usize,
    pub(super) one_edge_interior_overlap_plans: usize,
    pub(super) one_edge_interior_overlap_boundary_edges: usize,
}

impl SharedStartCycleShape {
    pub(super) fn empty() -> Self {
        Self {
            cycles: 0,
            cross_loop_cycles: 0,
            same_loop_cycles: 0,
            one_edge_interior_overlap_candidates: 0,
            repeated_one_edge_overlap_occurrences: 0,
            one_edge_interior_overlap_plans: 0,
            one_edge_interior_overlap_boundary_edges: 0,
        }
    }

    pub(super) fn from_cycles(cycles: &[FaceFamilyBridgeCycleRequest], brep: &BRep) -> Self {
        let mut same_loop_cycles = 0;
        let mut one_edge_interior_overlap_candidates = 0;
        let mut repeated_one_edge_overlap_occurrences = 0;
        let mut one_edge_interior_overlap_plans = 0;
        let mut one_edge_interior_overlap_boundary_edges = 0;
        for cycle in cycles
            .iter()
            .filter(|cycle| cycle.first_loop == cycle.second_loop)
        {
            same_loop_cycles += 1;
            let overlap = same_loop_overlap_summary(brep, cycle);
            let overlap_edges = overlap.edges;
            if overlap_edges == 1
                && overlap.position == SameLoopOverlapPosition::Interior
                && let Ok(decomposition) = one_edge_interior_overlap_decomposition(brep, cycle)
            {
                if decomposition.unique_overlap_edge {
                    one_edge_interior_overlap_candidates += 1;
                }
                if decomposition.repeated_overlap_occurrences {
                    repeated_one_edge_overlap_occurrences += 1;
                }
                if decomposition.shared_edge_already_boundary {
                    one_edge_interior_overlap_boundary_edges += 1;
                }
                if decomposition.plan.is_some() {
                    one_edge_interior_overlap_plans += 1;
                }
            }
        }
        Self {
            cycles: cycles.len(),
            cross_loop_cycles: cycles.len() - same_loop_cycles,
            same_loop_cycles,
            one_edge_interior_overlap_candidates,
            repeated_one_edge_overlap_occurrences,
            one_edge_interior_overlap_plans,
            one_edge_interior_overlap_boundary_edges,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SharedStartCycleMixSummary {
    pub(super) cycles: usize,
    pub(super) cross_loop_cycles: usize,
    pub(super) same_loop_cycles: usize,
    pub(super) one_edge_interior_overlap_candidates: usize,
    pub(super) repeated_one_edge_overlap_occurrences: usize,
    pub(super) one_edge_interior_overlap_plans: usize,
    pub(super) one_edge_interior_overlap_boundary_edges: usize,
    pub(super) one_edge_interior_overlap_normalizations: usize,
    pub(super) expected_edges: i64,
    pub(super) expected_faces: i64,
    pub(super) expected_loops: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SameLoopInteriorOverlapPlan {
    pub(super) loop_: LoopHandle,
    pub(super) shared_edge: EdgeHandle,
    pub(super) first_shared_half_edge: HalfEdgeHandle,
    pub(super) second_shared_half_edge: HalfEdgeHandle,
    pub(super) first_before_shared: Vec<HalfEdgeHandle>,
    pub(super) first_after_shared: Vec<HalfEdgeHandle>,
    pub(super) second_before_shared: Vec<HalfEdgeHandle>,
    pub(super) second_after_shared: Vec<HalfEdgeHandle>,
}

pub(super) fn preflight_shared_start_cycle_mix(
    brep: &BRep,
    cycles: Vec<FaceFamilyBridgeCycleRequest>,
) -> Result<SharedStartCycleMixSummary, CoreError> {
    let mut cross_loop_cycles = Vec::new();
    let mut same_loop_cycles = Vec::new();
    let mut interior_overlap_plans = Vec::new();
    let mut one_edge_interior_overlap_candidates = 0;
    let mut repeated_one_edge_overlap_occurrences = 0;
    let mut one_edge_interior_overlap_boundary_edges = 0;
    for cycle in cycles {
        if cycle.first_loop == cycle.second_loop {
            let decomposition = one_edge_interior_overlap_decomposition(brep, &cycle)?;
            if decomposition.unique_overlap_edge {
                one_edge_interior_overlap_candidates += 1;
            }
            if decomposition.repeated_overlap_occurrences {
                repeated_one_edge_overlap_occurrences += 1;
            }
            if decomposition.shared_edge_already_boundary {
                one_edge_interior_overlap_boundary_edges += 1;
            }
            if let Some(plan) = decomposition.plan {
                interior_overlap_plans.push(plan);
            }
            same_loop_cycles.push(cycle.clone());
        } else {
            cross_loop_cycles.push(cycle);
        }
    }
    let normalized_same_loop =
        preflight_normalized_same_loop_cycles(brep, &interior_overlap_plans, &same_loop_cycles)?;
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
    expected_edges += normalized_same_loop.expected_edges;
    expected_faces += normalized_same_loop.expected_faces;
    expected_loops += normalized_same_loop.expected_loops;
    Ok(SharedStartCycleMixSummary {
        cycles: cross_loop_count + same_loop_count,
        cross_loop_cycles: cross_loop_count,
        same_loop_cycles: same_loop_count,
        one_edge_interior_overlap_candidates,
        repeated_one_edge_overlap_occurrences,
        one_edge_interior_overlap_plans: interior_overlap_plans.len(),
        one_edge_interior_overlap_boundary_edges,
        one_edge_interior_overlap_normalizations: normalized_same_loop.normalizations,
        expected_edges,
        expected_faces,
        expected_loops,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct NormalizedSameLoopInteriorOverlaps {
    normalizations: usize,
    expected_edges: i64,
    expected_faces: i64,
    expected_loops: i64,
}

fn preflight_normalized_same_loop_cycles(
    brep: &BRep,
    plans: &[SameLoopInteriorOverlapPlan],
    cycles: &[FaceFamilyBridgeCycleRequest],
) -> Result<NormalizedSameLoopInteriorOverlaps, CoreError> {
    let mut expected_edges = 0;
    let mut expected_faces = 0;
    let mut expected_loops = 0;
    if plans.is_empty() {
        for cycle in cycles {
            let plan = brep.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
                loop_: cycle.first_loop,
                first_interval: cycle.first_interval.clone(),
                second_interval: cycle.second_interval.clone(),
            })?;
            expected_edges += plan.expected_delta.edges;
            expected_faces += plan.expected_delta.faces;
            expected_loops += plan.expected_delta.loops;
        }
        return Ok(NormalizedSameLoopInteriorOverlaps {
            normalizations: 0,
            expected_edges,
            expected_faces,
            expected_loops,
        });
    }
    let mut candidate = brep.clone();
    for plan in plans {
        if candidate.half_edges[plan.first_shared_half_edge].edge != plan.shared_edge
            || candidate.half_edges[plan.second_shared_half_edge].edge != plan.shared_edge
        {
            return Err(CoreError::PreconditionFailed(
                "same-loop overlap plan is stale".into(),
            ));
        }
        candidate.expose_same_loop_edge_as_boundary(plan.first_shared_half_edge)?;
        expected_edges += 1;
    }
    for cycle in cycles {
        let plan = candidate.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_: cycle.first_loop,
            first_interval: cycle.first_interval.clone(),
            second_interval: cycle.second_interval.clone(),
        })?;
        expected_edges += plan.expected_delta.edges;
        expected_faces += plan.expected_delta.faces;
        expected_loops += plan.expected_delta.loops;
    }
    Ok(NormalizedSameLoopInteriorOverlaps {
        normalizations: plans.len(),
        expected_edges,
        expected_faces,
        expected_loops,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OneEdgeInteriorOverlapDecomposition {
    unique_overlap_edge: bool,
    repeated_overlap_occurrences: bool,
    shared_edge_already_boundary: bool,
    plan: Option<SameLoopInteriorOverlapPlan>,
}

fn one_edge_interior_overlap_decomposition(
    brep: &BRep,
    cycle: &FaceFamilyBridgeCycleRequest,
) -> Result<OneEdgeInteriorOverlapDecomposition, CoreError> {
    let first_path = indexed_bridge_path_edges(brep, &cycle.first_interval)?;
    let second_path = indexed_bridge_path_edges(brep, &cycle.second_interval)?;
    let first_edges = first_path
        .iter()
        .map(|entry| entry.edge)
        .collect::<BTreeSet<_>>();
    let second_edges = second_path
        .iter()
        .map(|entry| entry.edge)
        .collect::<BTreeSet<_>>();
    let overlap_edges = first_edges
        .intersection(&second_edges)
        .copied()
        .collect::<Vec<_>>();
    if overlap_edges.len() != 1 {
        return Ok(OneEdgeInteriorOverlapDecomposition {
            unique_overlap_edge: false,
            repeated_overlap_occurrences: false,
            shared_edge_already_boundary: false,
            plan: None,
        });
    }
    let shared_edge = overlap_edges[0];
    let shared_edge_already_boundary = brep.boundary_edges.contains(&shared_edge);
    let shared = first_path
        .iter()
        .flat_map(|first| {
            second_path
                .iter()
                .filter(move |second| first.edge == shared_edge && second.edge == shared_edge)
                .map(move |second| (first, second))
        })
        .collect::<Vec<_>>();
    if shared.len() != 1 {
        return Ok(OneEdgeInteriorOverlapDecomposition {
            unique_overlap_edge: true,
            repeated_overlap_occurrences: true,
            shared_edge_already_boundary,
            plan: None,
        });
    }
    let (first_shared, second_shared) = shared[0];
    Ok(OneEdgeInteriorOverlapDecomposition {
        unique_overlap_edge: true,
        repeated_overlap_occurrences: false,
        shared_edge_already_boundary,
        plan: Some(SameLoopInteriorOverlapPlan {
            loop_: cycle.first_loop,
            shared_edge: first_shared.edge,
            first_shared_half_edge: first_shared.half_edge,
            second_shared_half_edge: second_shared.half_edge,
            first_before_shared: first_path[..first_shared.index]
                .iter()
                .map(|entry| entry.half_edge)
                .collect(),
            first_after_shared: first_path[first_shared.index + 1..]
                .iter()
                .map(|entry| entry.half_edge)
                .collect(),
            second_before_shared: second_path[..second_shared.index]
                .iter()
                .map(|entry| entry.half_edge)
                .collect(),
            second_after_shared: second_path[second_shared.index + 1..]
                .iter()
                .map(|entry| entry.half_edge)
                .collect(),
        }),
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IndexedBridgePathEdge {
    index: usize,
    half_edge: HalfEdgeHandle,
    edge: EdgeHandle,
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

fn indexed_bridge_path_edges(
    brep: &BRep,
    interval: &[HalfEdgeHandle],
) -> Result<Vec<IndexedBridgePathEdge>, CoreError> {
    half_open_bridge_path(interval)?
        .iter()
        .copied()
        .enumerate()
        .map(|(index, half_edge)| {
            let edge = brep
                .half_edges
                .get(half_edge)
                .map(|data| data.edge)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(
                        "shared-start bridge interval has stale half-edge".into(),
                    )
                })?;
            Ok(IndexedBridgePathEdge {
                index,
                half_edge,
                edge,
            })
        })
        .collect()
}

fn half_open_bridge_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "shared-start bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}
