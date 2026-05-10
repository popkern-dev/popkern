//! Diagnostics for cylinder-mediated support continuation authoring.

use std::collections::BTreeMap;

use crate::topology::FaceHandle;

use super::{
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement,
};

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
