//! Connector-authoring mutation dispatch for reconciled cylinder resources.

use crate::error::CoreError;
use crate::euler::{
    BoundarySideSeparatorBridgeRequest, FaceFamilyBridgeCycleRequest, SameLoopFaceBridgeRequest,
};
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::super::super::apply_bridge_cycles::SegmentBridgeCycleOperandSet;
use super::super::CylinderBridgeResourceSegment;
use super::split_merge::apply_multi_run_bridge_face_plans;
use super::{
    BridgeResourceConstructionActionKind, BridgeResourceReconciliation,
    BridgeResourceReconciliationStatus,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(in crate::euler::cylinder_mediated_support_apply) struct ConnectorAuthoringMutationResult {
    pub(in crate::euler::cylinder_mediated_support_apply) bridge_faces: Vec<FaceHandle>,
    pub(in crate::euler::cylinder_mediated_support_apply) bridge_loops: Vec<LoopHandle>,
}

impl BridgeResourceReconciliation {
    pub(in crate::euler::cylinder_mediated_support_apply) fn apply_connector_authoring_mutation(
        &self,
        brep: &mut BRep,
        bridge_cycles: &SegmentBridgeCycleOperandSet,
        cylinder_segments: &[CylinderBridgeResourceSegment],
    ) -> Result<Option<ConnectorAuthoringMutationResult>, CoreError> {
        if self.status()
            != BridgeResourceReconciliationStatus::SplitCylinderResourcesNeedConnectorAuthoring
        {
            return Ok(None);
        }
        let mut bridge_faces = Vec::new();
        let mut bridge_loops = Vec::new();
        let exact_cycles = self.exact_pair_bridge_cycles(brep, bridge_cycles, cylinder_segments)?;
        if !exact_cycles.is_empty() {
            apply_exact_bridge_cycles(brep, exact_cycles, &mut bridge_faces, &mut bridge_loops)?;
        }
        if let Some(result) =
            apply_multi_run_bridge_face_plans(brep, &self.split_resource_merge_operands)?
        {
            bridge_faces.extend(result.bridge_faces);
            bridge_loops.extend(result.bridge_loops);
        }
        if bridge_faces.is_empty() {
            return Ok(None);
        }
        Ok(Some(ConnectorAuthoringMutationResult {
            bridge_faces,
            bridge_loops,
        }))
    }

    fn exact_pair_bridge_cycles(
        &self,
        brep: &BRep,
        bridge_cycles: &SegmentBridgeCycleOperandSet,
        cylinder_segments: &[CylinderBridgeResourceSegment],
    ) -> Result<Vec<FaceFamilyBridgeCycleRequest>, CoreError> {
        self.actions
            .iter()
            .filter(|action| action.kind == BridgeResourceConstructionActionKind::ExactPair)
            .map(|action| {
                let cylinder_index = action.cylinder_segment_index.ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "exact bridge resource action is missing a cylinder segment".into(),
                    )
                })?;
                let cylinder_segment = cylinder_segments.get(cylinder_index).ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "exact bridge resource action references a missing cylinder segment".into(),
                    )
                })?;
                Ok(FaceFamilyBridgeCycleRequest {
                    first_loop: bridge_cycles
                        .boundary_loop(action.outside_segment_index)
                        .ok_or_else(|| {
                            CoreError::PreconditionFailed(
                                "exact bridge resource action has no outside loop".into(),
                            )
                        })?,
                    first_interval: bridge_cycles
                        .boundary_interval(action.outside_segment_index)
                        .ok_or_else(|| {
                            CoreError::PreconditionFailed(
                                "exact bridge resource action has no outside interval".into(),
                            )
                        })?,
                    second_loop: cylinder_segment_loop(brep, cylinder_segment)?,
                    second_interval: cylinder_segment_interval(brep, cylinder_segment)?,
                })
            })
            .collect()
    }
}

fn apply_exact_bridge_cycles(
    brep: &mut BRep,
    exact_cycles: Vec<FaceFamilyBridgeCycleRequest>,
    bridge_faces: &mut Vec<FaceHandle>,
    bridge_loops: &mut Vec<LoopHandle>,
) -> Result<(), CoreError> {
    if let Ok(plan) = brep.preflight_face_family_bridge(exact_cycles.clone()) {
        let result = brep.face_family_bridge(&plan)?;
        bridge_faces.extend(result.cycles.iter().map(|cycle| cycle.bridge_face));
        bridge_loops.extend(result.cycles.iter().map(|cycle| cycle.bridge_loop));
        return Ok(());
    }
    if let Ok(plan) =
        brep.preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: exact_cycles.clone(),
        })
    {
        let result = brep.boundary_side_separator_bridge(&plan)?;
        bridge_faces.extend(result.bridge_faces);
        bridge_loops.extend(result.bridge_loops);
        return Ok(());
    }
    for cycle in exact_cycles
        .into_iter()
        .filter(|cycle| cycle.first_loop == cycle.second_loop)
    {
        let plan = brep.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_: cycle.first_loop,
            first_interval: cycle.first_interval,
            second_interval: cycle.second_interval,
        })?;
        let result = brep.same_loop_face_bridge(&plan)?;
        bridge_faces.push(result.bridge_face);
        bridge_loops.push(result.bridge_loop);
    }
    Ok(())
}

fn cylinder_segment_loop(
    brep: &BRep,
    segment: &CylinderBridgeResourceSegment,
) -> Result<LoopHandle, CoreError> {
    let first = segment.half_edges.first().copied().ok_or_else(|| {
        CoreError::PreconditionFailed("cylinder bridge resource segment is empty".into())
    })?;
    brep.half_edges
        .get(first)
        .map(|half_edge| half_edge.loop_)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "cylinder bridge resource segment has stale first half-edge".into(),
            )
        })
}

fn cylinder_segment_interval(
    brep: &BRep,
    segment: &CylinderBridgeResourceSegment,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let last = segment.half_edges.last().copied().ok_or_else(|| {
        CoreError::PreconditionFailed("cylinder bridge resource segment is empty".into())
    })?;
    let after_last = brep
        .half_edges
        .get(last)
        .map(|half_edge| half_edge.next)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "cylinder bridge resource segment has stale last half-edge".into(),
            )
        })?;
    let mut interval = segment.half_edges.clone();
    interval.push(after_last);
    Ok(interval)
}
