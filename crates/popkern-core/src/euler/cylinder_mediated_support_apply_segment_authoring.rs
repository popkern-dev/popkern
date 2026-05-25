//! Segment authoring classification for cylinder-mediated continuation.

use std::collections::BTreeMap;

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::apply_owner_spans::owner_repair_spans;
use super::{
    CylinderMediatedContextOwnerFace, CylinderMediatedSupportContinuationBatchPlan,
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement, CylinderMediatedUnsupportedOwnerFace,
    segment_source_fact_ids,
};

pub(super) fn segment_authoring_plan(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> Result<Vec<CylinderMediatedSupportContinuationSegmentAuthoringPlan>, CoreError> {
    let mut authoring = Vec::new();
    for loop_partition in &plan.outside_interval_partition.loops {
        for (segment_index, segment) in loop_partition.segments.iter().enumerate() {
            let mut segment_plan = CylinderMediatedSupportContinuationSegmentAuthoringPlan {
                outside_loop: loop_partition.outside_loop,
                segment_index,
                request_indices: segment.request_indices.clone(),
                half_edges: segment.half_edges.clone(),
                boundary_bridge_half_edges: Vec::new(),
                cylinder_owned_half_edges: Vec::new(),
                outside_owned_half_edges: Vec::new(),
                support_owned_half_edges: Vec::new(),
                unassigned_twin_half_edges: Vec::new(),
                unsupported_owned_half_edges: Vec::new(),
                context_owned_half_edges: Vec::new(),
                unsupported_owner_faces: Vec::new(),
                context_owner_faces: Vec::new(),
                owner_repair_spans: Vec::new(),
                requirement: CylinderMediatedSupportContinuationSegmentRequirement::BoundaryBridge,
                source_fact_ids: segment_source_fact_ids(brep, &segment.half_edges),
            };
            for &half_edge in &segment.half_edges {
                classify_segment_half_edge(brep, plan, &mut segment_plan, half_edge)?;
            }
            segment_plan.unsupported_owner_faces = unsupported_owner_faces(brep, &segment_plan)?;
            segment_plan.context_owner_faces = context_owner_faces(brep, &segment_plan)?;
            segment_plan.owner_repair_spans = owner_repair_spans(brep, &segment_plan)?;
            segment_plan.requirement = segment_requirement(&segment_plan);
            authoring.push(segment_plan);
        }
    }
    Ok(authoring)
}

fn classify_segment_half_edge(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment_plan: &mut CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    half_edge: HalfEdgeHandle,
) -> Result<(), CoreError> {
    let data = brep.half_edges.get(half_edge).ok_or_else(|| {
        CoreError::InvalidHandle(format!(
            "cylinder-mediated segment references stale half-edge {half_edge:?}"
        ))
    })?;
    if brep.boundary_edges.contains(&data.edge) && !brep.half_edges.contains_key(data.twin) {
        segment_plan.boundary_bridge_half_edges.push(half_edge);
        return Ok(());
    }
    let Some(twin) = brep.half_edges.get(data.twin) else {
        segment_plan.unsupported_owned_half_edges.push(half_edge);
        return Ok(());
    };
    let Some(twin_loop) = brep.loops.get(twin.loop_) else {
        segment_plan.unassigned_twin_half_edges.push(half_edge);
        return Ok(());
    };
    if selected_owner_face(
        plan,
        &segment_plan.request_indices,
        twin_loop.face,
        OwnerRole::Cylinder,
    ) {
        segment_plan.cylinder_owned_half_edges.push(half_edge);
        return Ok(());
    }
    if selected_owner_face(
        plan,
        &segment_plan.request_indices,
        twin_loop.face,
        OwnerRole::Outside,
    ) {
        segment_plan.outside_owned_half_edges.push(half_edge);
        return Ok(());
    }
    if selected_owner_face(
        plan,
        &segment_plan.request_indices,
        twin_loop.face,
        OwnerRole::Support,
    ) {
        segment_plan.support_owned_half_edges.push(half_edge);
        return Ok(());
    }
    if selected_owner_face(
        plan,
        &segment_plan.request_indices,
        twin_loop.face,
        OwnerRole::Context,
    ) {
        segment_plan.context_owned_half_edges.push(half_edge);
    } else {
        segment_plan.unsupported_owned_half_edges.push(half_edge);
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OwnerRole {
    Cylinder,
    Outside,
    Support,
    Context,
}

fn selected_owner_face(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    request_indices: &[usize],
    face: FaceHandle,
    role: OwnerRole,
) -> bool {
    request_indices.iter().any(|&request_index| {
        plan.plans
            .get(request_index)
            .map(|continuation| match role {
                OwnerRole::Cylinder => continuation.request.cylinder_face == face,
                OwnerRole::Outside => continuation.request.outside_face == face,
                OwnerRole::Support => continuation.request.support_face == face,
                OwnerRole::Context => continuation.request.context_owner_faces.contains(&face),
            })
            .unwrap_or(false)
    })
}

fn unsupported_owner_faces(
    brep: &BRep,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<CylinderMediatedUnsupportedOwnerFace>, CoreError> {
    let counts = owner_face_counts(brep, &segment.unsupported_owned_half_edges)?;
    Ok(counts
        .into_iter()
        .map(|(face, half_edges)| CylinderMediatedUnsupportedOwnerFace { face, half_edges })
        .collect())
}

fn context_owner_faces(
    brep: &BRep,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<CylinderMediatedContextOwnerFace>, CoreError> {
    let counts = owner_face_counts(brep, &segment.context_owned_half_edges)?;
    Ok(counts
        .into_iter()
        .map(|(face, half_edges)| CylinderMediatedContextOwnerFace { face, half_edges })
        .collect())
}

fn owner_face_counts(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Result<BTreeMap<FaceHandle, usize>, CoreError> {
    let mut counts = BTreeMap::<FaceHandle, usize>::new();
    for &half_edge in half_edges {
        let twin = brep
            .half_edges
            .get(half_edge)
            .and_then(|data| brep.half_edges.get(data.twin))
            .ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "owner half-edge {half_edge:?} has no live twin"
                ))
            })?;
        let face = brep
            .loops
            .get(twin.loop_)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "owner half-edge {half_edge:?} has no live twin loop"
                ))
            })?;
        *counts.entry(face).or_default() += 1;
    }
    Ok(counts)
}

fn segment_requirement(
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> CylinderMediatedSupportContinuationSegmentRequirement {
    if !segment.unsupported_owned_half_edges.is_empty() {
        return CylinderMediatedSupportContinuationSegmentRequirement::UnsupportedOwnedEdge;
    }
    if !segment.unassigned_twin_half_edges.is_empty() {
        return CylinderMediatedSupportContinuationSegmentRequirement::UnassignedTwinRepair;
    }
    if !segment.cylinder_owned_half_edges.is_empty()
        || !segment.outside_owned_half_edges.is_empty()
        || !segment.support_owned_half_edges.is_empty()
        || !segment.context_owned_half_edges.is_empty()
    {
        return CylinderMediatedSupportContinuationSegmentRequirement::OwnedEdgeReclassification;
    }
    CylinderMediatedSupportContinuationSegmentRequirement::BoundaryBridge
}
