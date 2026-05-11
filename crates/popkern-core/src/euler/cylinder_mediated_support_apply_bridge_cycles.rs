//! Bridge-cycle operand derivation after owner-side exposure.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::{
    CylinderMediatedSupportContinuationBatchPlan,
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
};
use crate::euler::CylinderMediatedSupportContinuationRequest;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SegmentBridgeCycleOperandSet {
    operands: Vec<SegmentBridgeCycleOperand>,
}

impl SegmentBridgeCycleOperandSet {
    pub(super) fn segment_count(&self) -> usize {
        self.operands.len()
    }

    pub(super) fn request_index_sets(&self) -> impl Iterator<Item = &[usize]> {
        self.operands
            .iter()
            .map(|operand| operand.request_indices.as_slice())
    }

    pub(super) fn boundary_endpoint_pair(
        &self,
        segment_index: usize,
    ) -> Option<(VertexHandle, VertexHandle)> {
        self.operands
            .get(segment_index)
            .map(|operand| (operand.start, operand.end))
    }

    pub(super) fn boundary_half_edge_count(&self, segment_index: usize) -> Option<usize> {
        self.operands
            .get(segment_index)
            .map(|operand| operand.boundary_half_edges.len())
    }

    pub(super) fn ownership_status(&self) -> SegmentBridgeCycleOwnershipStatus {
        if self
            .operands
            .iter()
            .any(|operand| operand.support_faces.is_empty())
        {
            return SegmentBridgeCycleOwnershipStatus::MissingSupportOwner;
        }
        if self.has_shared_support_segments() {
            return SegmentBridgeCycleOwnershipStatus::SharedSupportTransitionReady;
        }
        SegmentBridgeCycleOwnershipStatus::SingleSupportReady
    }

    pub(super) fn has_shared_support_segments(&self) -> bool {
        self.operands
            .iter()
            .any(|operand| operand.support_faces.len() > 1)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let shared_support_segments = self
            .operands
            .iter()
            .filter(|operand| operand.support_faces.len() > 1)
            .count();
        let missing_support_segments = self
            .operands
            .iter()
            .filter(|operand| operand.support_faces.is_empty())
            .count();
        let max_support_faces = self
            .operands
            .iter()
            .map(|operand| operand.support_faces.len())
            .max()
            .unwrap_or_default();
        let total_half_edges = self
            .operands
            .iter()
            .map(|operand| operand.boundary_half_edges.len())
            .sum::<usize>();
        let degenerate_endpoint_segments = self
            .operands
            .iter()
            .filter(|operand| operand.start == operand.end)
            .count();
        let max_cylinder_faces = self
            .operands
            .iter()
            .map(|operand| operand.cylinder_faces.len())
            .max()
            .unwrap_or_default();
        let missing_opposite_path_segments = self
            .operands
            .iter()
            .filter(|operand| operand.opposite_cylinder_half_edges.is_empty())
            .count();
        let opposite_cylinder_half_edges = self
            .operands
            .iter()
            .map(|operand| operand.opposite_cylinder_half_edges.len())
            .sum::<usize>();
        let duplicate_opposite_cylinder_half_edges = duplicate_half_edge_count(
            self.operands
                .iter()
                .flat_map(|operand| operand.opposite_cylinder_half_edges.iter().copied()),
        );
        format!(
            "segment_bridge_candidate=segments={} boundary_half_edges={} ownership={} shared_support_segments={} missing_support_segments={} max_support_faces={} degenerate_endpoint_segments={} max_cylinder_faces={} opposite_cylinder_half_edges={} duplicate_opposite_cylinder_half_edges={} missing_opposite_path_segments={}",
            self.operands.len(),
            total_half_edges,
            self.ownership_status().diagnostic_label(),
            shared_support_segments,
            missing_support_segments,
            max_support_faces,
            degenerate_endpoint_segments,
            max_cylinder_faces,
            opposite_cylinder_half_edges,
            duplicate_opposite_cylinder_half_edges,
            missing_opposite_path_segments
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SegmentBridgeCycleOwnershipStatus {
    SingleSupportReady,
    SharedSupportTransitionReady,
    MissingSupportOwner,
}

impl SegmentBridgeCycleOwnershipStatus {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::SingleSupportReady => "single-support-ready",
            Self::SharedSupportTransitionReady => "shared-support-transition-ready",
            Self::MissingSupportOwner => "missing-support-owner",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SegmentBridgeCycleOperand {
    request_indices: Vec<usize>,
    boundary_half_edges: Vec<HalfEdgeHandle>,
    start: VertexHandle,
    end: VertexHandle,
    support_faces: Vec<FaceHandle>,
    cylinder_faces: Vec<FaceHandle>,
    opposite_cylinder_half_edges: Vec<HalfEdgeHandle>,
}

pub(super) fn derive_segment_bridge_cycle_operands(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
) -> Result<SegmentBridgeCycleOperandSet, CoreError> {
    let operands = authoring
        .iter()
        .map(|segment| segment_bridge_cycle_operand(brep, plan, segment))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(SegmentBridgeCycleOperandSet { operands })
}

fn segment_bridge_cycle_operand(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<SegmentBridgeCycleOperand, CoreError> {
    let first = segment.half_edges.first().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "segment bridge cycle has empty segment {:?}",
            segment.segment_index
        ))
    })?;
    let last = segment.half_edges.last().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "segment bridge cycle has empty segment {:?}",
            segment.segment_index
        ))
    })?;
    for &half_edge in &segment.half_edges {
        require_boundary_ready_half_edge(brep, half_edge)?;
    }
    let start = brep.half_edges[first].origin;
    let after_last = brep.half_edges[last].next;
    let end = brep
        .half_edges
        .get(after_last)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "segment bridge cycle after-last half-edge {after_last:?} is stale"
            ))
        })?;
    Ok(SegmentBridgeCycleOperand {
        request_indices: segment.request_indices.clone(),
        boundary_half_edges: segment.half_edges.clone(),
        start,
        end,
        support_faces: selected_support_faces(plan, segment)?,
        cylinder_faces: selected_cylinder_faces(plan, segment)?,
        opposite_cylinder_half_edges: common_cylinder_path(plan, segment)?,
    })
}

fn require_boundary_ready_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Result<(), CoreError> {
    let data = brep.half_edges.get(half_edge).ok_or_else(|| {
        CoreError::InvalidHandle(format!("segment bridge half-edge {half_edge:?} is stale"))
    })?;
    if !brep.boundary_edges.contains(&data.edge) || brep.half_edges.contains_key(data.twin) {
        return Err(CoreError::PreconditionFailed(format!(
            "segment bridge half-edge {half_edge:?} is not boundary-ready after owner exposure"
        )));
    }
    Ok(())
}

fn selected_support_faces(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<FaceHandle>, CoreError> {
    selected_faces(plan, segment, |request| request.support_face)
}

fn selected_cylinder_faces(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<FaceHandle>, CoreError> {
    selected_faces(plan, segment, |request| request.cylinder_face)
}

fn common_cylinder_path(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let first_request_index = *segment.request_indices.first().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "segment bridge {:?} has no request owner",
            segment.segment_index
        ))
    })?;
    let first = plan.plans.get(first_request_index).ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "segment bridge references missing request index {first_request_index}"
        ))
    })?;
    let other_sets = segment
        .request_indices
        .iter()
        .skip(1)
        .map(|&request_index| {
            plan.plans
                .get(request_index)
                .map(|plan| plan.request.cylinder_interval.iter().copied().collect())
                .ok_or_else(|| {
                    CoreError::PreconditionFailed(format!(
                        "segment bridge references missing request index {request_index}"
                    ))
                })
        })
        .collect::<Result<Vec<BTreeSet<_>>, _>>()?;
    Ok(first
        .request
        .cylinder_interval
        .iter()
        .copied()
        .filter(|half_edge| other_sets.iter().all(|set| set.contains(half_edge)))
        .collect())
}

fn selected_faces(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    face: impl Fn(&CylinderMediatedSupportContinuationRequest) -> FaceHandle,
) -> Result<Vec<FaceHandle>, CoreError> {
    let mut faces = Vec::new();
    for &request_index in &segment.request_indices {
        let request = plan.plans.get(request_index).ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "segment bridge references missing request index {request_index}"
            ))
        })?;
        let selected = face(&request.request);
        if !faces.contains(&selected) {
            faces.push(selected);
        }
    }
    Ok(faces)
}

fn duplicate_half_edge_count(half_edges: impl Iterator<Item = HalfEdgeHandle>) -> usize {
    let mut seen = BTreeSet::new();
    let mut duplicates = BTreeSet::new();
    for half_edge in half_edges {
        if !seen.insert(half_edge) {
            duplicates.insert(half_edge);
        }
    }
    duplicates.len()
}
