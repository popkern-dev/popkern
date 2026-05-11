//! Cylinder-side resource partitioning for mediated support continuation.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle};

use super::CylinderMediatedSupportContinuationBatchPlan;

#[path = "cylinder_mediated_support_apply_cylinder_resource_reconciliation.rs"]
mod reconciliation;

pub(super) use reconciliation::BridgeResourceReconciliationStatus;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourcePartition {
    segments: Vec<CylinderBridgeResourceSegment>,
}

impl CylinderBridgeResourcePartition {
    pub(super) fn from_remapped_segments(
        segments: &[CylinderBridgeResourceSegment],
    ) -> CylinderBridgeResourcePartition {
        CylinderBridgeResourcePartition {
            segments: segments.to_vec(),
        }
    }

    pub(super) fn segments(&self) -> impl Iterator<Item = &CylinderBridgeResourceSegment> {
        self.segments.iter()
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let shared_segments = self
            .segments
            .iter()
            .filter(|segment| segment.request_indices.len() > 1)
            .count();
        let max_request_owners = self
            .segments
            .iter()
            .map(|segment| segment.request_indices.len())
            .max()
            .unwrap_or_default();
        let half_edges = self
            .segments
            .iter()
            .map(|segment| segment.half_edges.len())
            .sum::<usize>();
        format!(
            "cylinder_resource_partition=segments={} half_edges={} shared_segments={} max_request_owners={}",
            self.segments.len(),
            half_edges,
            shared_segments,
            max_request_owners
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourceSegment {
    pub(super) half_edges: Vec<HalfEdgeHandle>,
    pub(super) edges: Vec<EdgeHandle>,
    pub(super) request_indices: Vec<usize>,
}

pub(super) fn derive_cylinder_bridge_resource_partition(
    brep: &BRep,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> Result<CylinderBridgeResourcePartition, CoreError> {
    let mut selected_by_loop = BTreeMap::<_, BTreeMap<_, BTreeSet<usize>>>::new();
    for (request_index, request_plan) in plan.plans.iter().enumerate() {
        for &half_edge in &request_plan.request.cylinder_interval {
            let loop_ = brep
                .half_edges
                .get(half_edge)
                .map(|data| data.loop_)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(format!(
                        "cylinder resource partition references stale half-edge {half_edge:?}"
                    ))
                })?;
            selected_by_loop
                .entry(loop_)
                .or_default()
                .entry(half_edge)
                .or_default()
                .insert(request_index);
        }
    }
    let mut segments = Vec::new();
    for (loop_, membership) in selected_by_loop {
        segments.extend(partition_cylinder_loop(brep, loop_, &membership)?);
    }
    Ok(CylinderBridgeResourcePartition { segments })
}

fn partition_cylinder_loop(
    brep: &BRep,
    loop_: LoopHandle,
    membership: &BTreeMap<HalfEdgeHandle, BTreeSet<usize>>,
) -> Result<Vec<CylinderBridgeResourceSegment>, CoreError> {
    let start = brep
        .loops
        .get(loop_)
        .map(|loop_data| loop_data.half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!("cylinder resource loop {loop_:?} is stale"))
        })?;
    let mut segments = Vec::new();
    let mut current_half_edges = Vec::new();
    let mut current_edges = Vec::new();
    let mut current_request_indices = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let data = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "cylinder resource partition references stale loop half-edge {cursor:?}"
            ))
        })?;
        if data.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder resource partition traversal left loop {loop_:?} at {cursor:?}"
            )));
        }
        match membership.get(&cursor) {
            Some(request_indices) => {
                let request_indices = request_indices.iter().copied().collect::<Vec<_>>();
                if current_half_edges.is_empty() || current_request_indices == request_indices {
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                } else {
                    push_cylinder_resource_segment(
                        &mut segments,
                        &mut current_half_edges,
                        &mut current_edges,
                        &mut current_request_indices,
                    );
                    current_half_edges.push(cursor);
                    current_edges.push(data.edge);
                    current_request_indices = request_indices;
                }
            }
            None => push_cylinder_resource_segment(
                &mut segments,
                &mut current_half_edges,
                &mut current_edges,
                &mut current_request_indices,
            ),
        }
        cursor = data.next;
        if cursor == start {
            break;
        }
    }
    if cursor != start {
        return Err(CoreError::PreconditionFailed(format!(
            "cylinder resource partition could not close loop {loop_:?}"
        )));
    }
    push_cylinder_resource_segment(
        &mut segments,
        &mut current_half_edges,
        &mut current_edges,
        &mut current_request_indices,
    );
    Ok(segments)
}

fn push_cylinder_resource_segment(
    segments: &mut Vec<CylinderBridgeResourceSegment>,
    current_half_edges: &mut Vec<HalfEdgeHandle>,
    current_edges: &mut Vec<EdgeHandle>,
    current_request_indices: &mut Vec<usize>,
) {
    if current_half_edges.is_empty() {
        return;
    }
    segments.push(CylinderBridgeResourceSegment {
        half_edges: std::mem::take(current_half_edges),
        edges: std::mem::take(current_edges),
        request_indices: std::mem::take(current_request_indices),
    });
}
