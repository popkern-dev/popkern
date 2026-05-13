//! Core operand derivation for outside/outside support continuation.

use popkern_core::euler::{
    CylinderMediatedSupportContinuationDelta, CylinderMediatedSupportContinuationPlan,
    CylinderMediatedSupportContinuationRequest,
};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_outside::{
    CylinderIntervalOperand, CylinderMediatedCorePreflight, CylinderMediatedOutsideOutsideClass,
    LoopAdjacencyOperand,
};
use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportContinuationPlan, OutsideOutsideSupportTrimPlan,
};
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;

pub(super) fn outside_cylinder_adjacencies(
    knowledge: &TopologyKnowledgeGraph,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    cylinder_face: FaceHandle,
) -> Vec<LoopAdjacencyOperand> {
    knowledge
        .loop_adjacencies()
        .filter_map(|adjacency| {
            if adjacency.face() == outside_face
                && adjacency.loop_handle() == outside_loop
                && adjacency.adjacent_face() == cylinder_face
            {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.loop_handle(),
                    source_half_edge: adjacency.half_edge(),
                    cylinder_loop: adjacency.adjacent_loop(),
                    cylinder_half_edge: adjacency.adjacent_half_edge(),
                });
            }
            if adjacency.adjacent_face() == outside_face
                && adjacency.adjacent_loop() == outside_loop
                && adjacency.face() == cylinder_face
            {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.adjacent_loop(),
                    source_half_edge: adjacency.adjacent_half_edge(),
                    cylinder_loop: adjacency.loop_handle(),
                    cylinder_half_edge: adjacency.half_edge(),
                });
            }
            None
        })
        .collect()
}

pub(super) fn support_cylinder_adjacencies(
    knowledge: &TopologyKnowledgeGraph,
    support_face: FaceHandle,
    cylinder_face: FaceHandle,
) -> Vec<LoopAdjacencyOperand> {
    knowledge
        .loop_adjacencies()
        .filter_map(|adjacency| {
            if adjacency.face() == support_face && adjacency.adjacent_face() == cylinder_face {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.loop_handle(),
                    source_half_edge: adjacency.half_edge(),
                    cylinder_loop: adjacency.adjacent_loop(),
                    cylinder_half_edge: adjacency.adjacent_half_edge(),
                });
            }
            if adjacency.adjacent_face() == support_face && adjacency.face() == cylinder_face {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.adjacent_loop(),
                    source_half_edge: adjacency.adjacent_half_edge(),
                    cylinder_loop: adjacency.loop_handle(),
                    cylinder_half_edge: adjacency.half_edge(),
                });
            }
            None
        })
        .collect()
}

pub(super) fn cylinder_interval_candidates(
    brep: &BRep,
    outside_adjacencies: &[LoopAdjacencyOperand],
    support_adjacencies: &[LoopAdjacencyOperand],
) -> Vec<CylinderIntervalOperand> {
    let mut intervals = outside_adjacencies
        .iter()
        .flat_map(|outside| {
            support_adjacencies
                .iter()
                .filter(move |support| outside.cylinder_loop == support.cylinder_loop)
                .filter_map(move |support| {
                    directed_loop_path(
                        brep,
                        outside.cylinder_loop,
                        outside.cylinder_half_edge,
                        support.cylinder_half_edge,
                    )
                    .map(|interval| CylinderIntervalOperand {
                        outside_loop: outside.source_loop,
                        outside_cylinder_source_half_edge: outside.source_half_edge,
                        outside_cylinder_half_edge: outside.cylinder_half_edge,
                        support_loop: support.source_loop,
                        support_cylinder_source_half_edge: support.source_half_edge,
                        support_cylinder_half_edge: support.cylinder_half_edge,
                        cylinder_loop: outside.cylinder_loop,
                        interval,
                    })
                })
        })
        .collect::<Vec<_>>();
    intervals.sort_by_key(|interval| {
        (
            interval.interval.len(),
            interval.outside_cylinder_source_half_edge,
            interval.support_cylinder_source_half_edge,
        )
    });
    intervals
}

pub(super) fn validate_core_preflights(
    region: &OutsideOutsideSupportContinuationPlan,
    topology_class: &CylinderMediatedOutsideOutsideClass,
    knowledge: &TopologyKnowledgeGraph,
    brep: &BRep,
) -> Vec<CylinderMediatedCorePreflight> {
    region
        .trims()
        .iter()
        .map(|trim| validate_trim_core_preflight(trim, topology_class, knowledge, brep))
        .collect()
}

fn validate_trim_core_preflight(
    trim: &OutsideOutsideSupportTrimPlan,
    topology_class: &CylinderMediatedOutsideOutsideClass,
    knowledge: &TopologyKnowledgeGraph,
    brep: &BRep,
) -> CylinderMediatedCorePreflight {
    if trim.start_face() != trim.end_face() || trim.start_loop() != trim.end_loop() {
        return CylinderMediatedCorePreflight::CrossOutsideOwners {
            fillet_index: trim.fillet_index(),
            start_face: trim.start_face(),
            start_loop: trim.start_loop(),
            end_face: trim.end_face(),
            end_loop: trim.end_loop(),
            start_half_edge: trim.start_half_edge(),
            end_half_edge: trim.end_half_edge(),
        };
    }
    let Some(interval) = topology_class
        .cylinder_intervals
        .iter()
        .find(|interval| interval.outside_loop == trim.start_loop())
    else {
        return CylinderMediatedCorePreflight::Rejected {
            fillet_index: trim.fillet_index(),
            reason: format!(
                "missing cylinder interval for outside loop {:?}",
                trim.start_loop()
            ),
        };
    };
    let request = CylinderMediatedSupportContinuationRequest {
        support_face: trim.support_face(),
        outside_face: trim.start_face(),
        cylinder_face: trim.cylinder_face(),
        start: trim.start(),
        end: trim.end(),
        outside_loop: trim.start_loop(),
        outside_start_half_edge: trim.start_half_edge(),
        outside_end_half_edge: trim.end_half_edge(),
        support_loop: interval.support_loop,
        support_cylinder_half_edge: interval.support_cylinder_source_half_edge,
        outside_cylinder_half_edge: interval.outside_cylinder_source_half_edge,
        cylinder_interval: interval.interval.clone(),
        context_owner_faces: knowledge.context_faces().to_vec(),
    };
    match brep.preflight_cylinder_mediated_support_continuation(request) {
        Ok(plan) => CylinderMediatedCorePreflight::Accepted {
            fillet_index: trim.fillet_index(),
            outside_loop: trim.start_loop(),
            support_loop: interval.support_loop,
            outside_interval_len: plan.outside_interval.len(),
            outside_complement_interval_len: plan.outside_complement_interval.len(),
            reason: preflight_delta_reason(&plan),
            request: plan.request,
        },
        Err(err) => CylinderMediatedCorePreflight::Rejected {
            fillet_index: trim.fillet_index(),
            reason: err.to_string(),
        },
    }
}

fn directed_loop_path(
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

fn preflight_delta_reason(plan: &CylinderMediatedSupportContinuationPlan) -> &'static str {
    match plan.expected_delta {
        CylinderMediatedSupportContinuationDelta::Unknown { reason } => reason,
    }
}
