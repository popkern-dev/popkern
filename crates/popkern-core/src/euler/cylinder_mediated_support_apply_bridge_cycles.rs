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
    pub(super) fn diagnostic_label(&self) -> String {
        let ambiguous_support_segments = self
            .operands
            .iter()
            .filter(|operand| operand.support_faces.len() != 1)
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
        format!(
            "segment_bridge_candidate=segments={} boundary_half_edges={} ambiguous_support_segments={} max_support_faces={} degenerate_endpoint_segments={} max_cylinder_faces={}",
            self.operands.len(),
            total_half_edges,
            ambiguous_support_segments,
            max_support_faces,
            degenerate_endpoint_segments,
            max_cylinder_faces
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SegmentBridgeCycleOperand {
    boundary_half_edges: Vec<HalfEdgeHandle>,
    start: VertexHandle,
    end: VertexHandle,
    support_faces: Vec<FaceHandle>,
    cylinder_faces: Vec<FaceHandle>,
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
        boundary_half_edges: segment.half_edges.clone(),
        start,
        end,
        support_faces: selected_support_faces(plan, segment)?,
        cylinder_faces: selected_cylinder_faces(plan, segment)?,
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

fn selected_faces(
    plan: &CylinderMediatedSupportContinuationBatchPlan,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    face: impl Fn(&CylinderMediatedSupportContinuationRequest) -> FaceHandle,
) -> Result<Vec<FaceHandle>, CoreError> {
    let mut faces = BTreeSet::new();
    for &request_index in &segment.request_indices {
        let request = plan.plans.get(request_index).ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "segment bridge references missing request index {request_index}"
            ))
        })?;
        faces.insert(face(&request.request));
    }
    Ok(faces.into_iter().collect())
}
