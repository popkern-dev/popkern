//! Diagnostics for cylinder-mediated support continuation authoring.

use std::collections::BTreeMap;

use crate::euler::CylinderMediatedSupportContinuationAggregateDelta;
use crate::topology::FaceHandle;

use super::apply_bridge_cycles::SegmentBridgeCycleOwnershipStatus;
use super::apply_cylinder_resources::BridgeResourceReconciliationStatus;
use super::{
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement,
};

pub(super) fn bridge_refusal_reason(
    ownership_status: SegmentBridgeCycleOwnershipStatus,
    resource_status: BridgeResourceReconciliationStatus,
) -> &'static str {
    if ownership_status == SegmentBridgeCycleOwnershipStatus::MissingSupportOwner {
        return "cylinder-mediated support continuation segment bridge has no support owner";
    }
    match resource_status {
        BridgeResourceReconciliationStatus::ExactResourcePairing => {
            "cylinder-mediated support continuation exact resource bridge construction not implemented"
        }
        BridgeResourceReconciliationStatus::SharedOutsideTransitionNeedsResource => {
            "cylinder-mediated support continuation shared outside transition bridge construction not implemented"
        }
        BridgeResourceReconciliationStatus::SingleSupportCoveredBySharedCylinderResource => {
            "cylinder-mediated support continuation shared cylinder resource bridge construction not implemented"
        }
        BridgeResourceReconciliationStatus::SharedStartSplitCylinderResourcesReady => {
            "cylinder-mediated support continuation shared-start split resource bridge construction not implemented"
        }
        BridgeResourceReconciliationStatus::SplitCylinderResourcesNeedConnectorAuthoring => {
            "cylinder-mediated support continuation split cylinder resource connector authoring failed"
        }
        BridgeResourceReconciliationStatus::SplitCylinderResourcesNeedMerge => {
            "cylinder-mediated support continuation split cylinder resource merge construction not implemented"
        }
        BridgeResourceReconciliationStatus::IncompatibleResourceSets => {
            "cylinder-mediated support continuation bridge resources are incompatible"
        }
    }
}

pub(super) fn segment_authoring_summary(
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
) -> String {
    let boundary_compatible = authoring
        .iter()
        .filter(|segment| segment.boundary_bridge_half_edges.len() == segment.half_edges.len())
        .count();
    let cylinder_owned = authoring
        .iter()
        .map(|segment| segment.cylinder_owned_half_edges.len())
        .sum::<usize>();
    let outside_owned = authoring
        .iter()
        .map(|segment| segment.outside_owned_half_edges.len())
        .sum::<usize>();
    let support_owned = authoring
        .iter()
        .map(|segment| segment.support_owned_half_edges.len())
        .sum::<usize>();
    let unassigned_twins = authoring
        .iter()
        .map(|segment| segment.unassigned_twin_half_edges.len())
        .sum::<usize>();
    let unsupported_owned = authoring
        .iter()
        .map(|segment| segment.unsupported_owned_half_edges.len())
        .sum::<usize>();
    let context_owned = authoring
        .iter()
        .map(|segment| segment.context_owned_half_edges.len())
        .sum::<usize>();
    let owner_repair_spans = authoring
        .iter()
        .map(|segment| segment.owner_repair_spans.len())
        .sum::<usize>();
    let owner_repair_half_edges = authoring
        .iter()
        .flat_map(|segment| segment.owner_repair_spans.iter())
        .map(|span| span.owner_half_edges.len())
        .sum::<usize>();
    let unsupported_owner_faces = unsupported_owner_face_summary(authoring);
    let context_owner_faces = context_owner_face_summary(authoring);
    let owned_reclassification_segments = authoring
        .iter()
        .filter(|segment| {
            segment.requirement
                == CylinderMediatedSupportContinuationSegmentRequirement::OwnedEdgeReclassification
        })
        .count();
    let unassigned_repair_segments = authoring
        .iter()
        .filter(|segment| {
            segment.requirement
                == CylinderMediatedSupportContinuationSegmentRequirement::UnassignedTwinRepair
        })
        .count();
    let unsupported_segments = authoring
        .iter()
        .filter(|segment| {
            segment.requirement
                == CylinderMediatedSupportContinuationSegmentRequirement::UnsupportedOwnedEdge
        })
        .count();
    format!(
        "segment_authoring=segments={} boundary_compatible={} owned_reclassification_segments={} unassigned_repair_segments={} unsupported_segments={} cylinder_owned_half_edges={} outside_owned_half_edges={} support_owned_half_edges={} context_owned_half_edges={} unassigned_twin_half_edges={} unsupported_owned_half_edges={} owner_repair_spans={} owner_repair_half_edges={} context_owner_faces=[{}] unsupported_owner_faces=[{}]",
        authoring.len(),
        boundary_compatible,
        owned_reclassification_segments,
        unassigned_repair_segments,
        unsupported_segments,
        cylinder_owned,
        outside_owned,
        support_owned,
        context_owned,
        unassigned_twins,
        unsupported_owned,
        owner_repair_spans,
        owner_repair_half_edges,
        context_owner_faces,
        unsupported_owner_faces
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AuthoredTopologyDelta {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

pub(super) fn authored_delta_label(owner_repair_spans: usize, bridge_segments: usize) -> String {
    let delta = authored_topology_delta(owner_repair_spans, bridge_segments);
    format!(
        "dV={} dE={} dF={} dL={}",
        delta.vertices, delta.edges, delta.faces, delta.loops
    )
}

pub(super) fn execution_delta_label(
    owner_repair_spans: usize,
    blocker_exposure_spans: usize,
    bridge_segments: usize,
) -> String {
    let delta =
        execution_topology_delta(owner_repair_spans, blocker_exposure_spans, bridge_segments);
    format!(
        "dV={} dE={} dF={} dL={}",
        delta.vertices, delta.edges, delta.faces, delta.loops
    )
}

pub(super) fn topology_delta_gap_label(
    expected: &CylinderMediatedSupportContinuationAggregateDelta,
    owner_repair_spans: usize,
    blocker_exposure_spans: usize,
    bridge_segments: usize,
) -> String {
    let execution =
        execution_topology_delta(owner_repair_spans, blocker_exposure_spans, bridge_segments);
    match expected {
        CylinderMediatedSupportContinuationAggregateDelta::Known {
            vertices,
            edges,
            faces,
            loops,
        } => format!(
            "topology_delta_gap=dV={} dE={} dF={} dL={}",
            execution.vertices - vertices,
            execution.edges - edges,
            execution.faces - faces,
            execution.loops - loops
        ),
        CylinderMediatedSupportContinuationAggregateDelta::Unknown { reason } => {
            format!("topology_delta_gap=unknown expected_reason={reason}")
        }
    }
}

fn execution_topology_delta(
    owner_repair_spans: usize,
    blocker_exposure_spans: usize,
    bridge_segments: usize,
) -> AuthoredTopologyDelta {
    AuthoredTopologyDelta {
        vertices: 0,
        edges: (owner_repair_spans + blocker_exposure_spans + bridge_segments * 2) as i64,
        faces: bridge_segments as i64,
        loops: bridge_segments as i64,
    }
}

fn authored_topology_delta(
    owner_repair_spans: usize,
    bridge_segments: usize,
) -> AuthoredTopologyDelta {
    AuthoredTopologyDelta {
        vertices: 0,
        edges: (owner_repair_spans + bridge_segments * 2) as i64,
        faces: bridge_segments as i64,
        loops: bridge_segments as i64,
    }
}

fn unsupported_owner_face_summary(
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
) -> String {
    let mut counts = BTreeMap::<FaceHandle, usize>::new();
    for face in authoring
        .iter()
        .flat_map(|segment| segment.unsupported_owner_faces.iter())
    {
        *counts.entry(face.face).or_default() += face.half_edges;
    }
    counts
        .into_iter()
        .map(|(face, half_edges)| format!("{face:?}:{half_edges}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn context_owner_face_summary(
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
) -> String {
    let mut counts = BTreeMap::<FaceHandle, usize>::new();
    for face in authoring
        .iter()
        .flat_map(|segment| segment.context_owner_faces.iter())
    {
        *counts.entry(face.face).or_default() += face.half_edges;
    }
    counts
        .into_iter()
        .map(|(face, half_edges)| format!("{face:?}:{half_edges}"))
        .collect::<Vec<_>>()
        .join(",")
}
