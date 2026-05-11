//! Transactional boundary for cylinder-mediated support continuation.

use std::collections::BTreeMap;

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::contract::EulerTopologyCounts;
use super::cylinder_mediated_support::{
    CylinderMediatedSupportContinuationAggregateDelta, CylinderMediatedSupportContinuationBatchPlan,
};
use super::cylinder_mediated_support_partition::CylinderMediatedOutsideIntervalPartition;

#[path = "cylinder_mediated_support_apply_bridge_cycles.rs"]
mod apply_bridge_cycles;
#[path = "cylinder_mediated_support_apply_cylinder_resource_blockers.rs"]
mod apply_cylinder_resource_blockers;
#[path = "cylinder_mediated_support_apply_cylinder_resource_remap.rs"]
mod apply_cylinder_resource_remap;
#[path = "cylinder_mediated_support_apply_cylinder_resources.rs"]
mod apply_cylinder_resources;
#[path = "cylinder_mediated_support_apply_diagnostics.rs"]
mod apply_diagnostics;
#[path = "cylinder_mediated_support_apply_owner_spans.rs"]
mod apply_owner_spans;

use apply_bridge_cycles::derive_segment_bridge_cycle_operands;
use apply_cylinder_resource_blockers::expose_cylinder_resource_blockers;
use apply_cylinder_resource_remap::remap_cylinder_bridge_resources;
use apply_cylinder_resources::{
    BridgeResourceReconciliationStatus, derive_cylinder_bridge_resource_partition,
};
use apply_diagnostics::segment_authoring_summary;
use apply_owner_spans::{expose_owner_repair_spans, owner_repair_spans};

#[cfg(test)]
#[path = "cylinder_mediated_support_apply_tests.rs"]
mod tests;

/// Per-segment topology authoring requirements for a cylinder-mediated batch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedSupportContinuationSegmentAuthoringPlan {
    /// Outside loop containing the segment.
    pub outside_loop: LoopHandle,
    /// Segment index within the loop partition.
    pub segment_index: usize,
    /// Batch request indices that select this segment.
    pub request_indices: Vec<usize>,
    /// Half-open outside/cap half-edges in the segment.
    pub half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges already compatible with a boundary bridge.
    pub boundary_bridge_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges whose twin side is already owned by a selected cylinder face.
    pub cylinder_owned_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges whose twin side is already owned by a selected outside/cap face.
    pub outside_owned_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges whose twin side is already owned by a selected support face.
    pub support_owned_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges whose twin exists but is not assigned to a live loop.
    pub unassigned_twin_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges owned by another face outside the accepted cylinder/outside set.
    pub unsupported_owned_half_edges: Vec<HalfEdgeHandle>,
    /// Half-edges owned by local context faces that are not selected operands.
    pub context_owned_half_edges: Vec<HalfEdgeHandle>,
    /// Unsupported owner faces and how many segment half-edges they own.
    pub unsupported_owner_faces: Vec<CylinderMediatedUnsupportedOwnerFace>,
    /// Context owner faces and how many segment half-edges they own.
    pub context_owner_faces: Vec<CylinderMediatedContextOwnerFace>,
    /// Validated owner-side spans that must be exposed or reclassified before
    /// this segment can be consumed by a bridge-face cycle.
    pub owner_repair_spans: Vec<CylinderMediatedOwnerRepairSpan>,
    /// Required topology authoring strategy for this segment.
    pub requirement: CylinderMediatedSupportContinuationSegmentRequirement,
}

/// Unsupported owner face observed on a partition segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CylinderMediatedUnsupportedOwnerFace {
    /// Unsupported twin-side face.
    pub face: FaceHandle,
    /// Number of segment half-edges whose twin side belongs to this face.
    pub half_edges: usize,
}

/// Local context owner face observed on a partition segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CylinderMediatedContextOwnerFace {
    /// Context twin-side face.
    pub face: FaceHandle,
    /// Number of segment half-edges whose twin side belongs to this face.
    pub half_edges: usize,
}

/// Owner-side path paired with a segment of outside/cap half-edges.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedOwnerRepairSpan {
    /// Face currently owning the opposite side of the outside/cap segment span.
    pub owner_face: FaceHandle,
    /// Loop currently owning the opposite side of the outside/cap segment span.
    pub owner_loop: LoopHandle,
    /// Outside/cap half-edges in outside-loop order.
    pub outside_half_edges: Vec<HalfEdgeHandle>,
    /// Opposite half-edges in owner-loop order.
    pub owner_half_edges: Vec<HalfEdgeHandle>,
}

/// Required topology strategy for one cylinder-mediated partition segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderMediatedSupportContinuationSegmentRequirement {
    /// Segment can be consumed by ordinary boundary bridge authoring.
    BoundaryBridge,
    /// Segment needs outside/cap or cylinder owned-edge reclassification.
    OwnedEdgeReclassification,
    /// Segment has unassigned twin half-edges that must be assigned or repaired first.
    UnassignedTwinRepair,
    /// Segment contains a face owner outside the accepted outside/cylinder set.
    UnsupportedOwnedEdge,
}

/// Result of the transactional cylinder-mediated continuation operator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedSupportContinuationBatchResult {
    /// Expected topology delta consumed by the operator.
    pub expected_delta: CylinderMediatedSupportContinuationAggregateDelta,
    /// Segment authoring plan consumed by the transaction.
    pub segment_authoring: Vec<CylinderMediatedSupportContinuationSegmentAuthoringPlan>,
    /// Postcondition status for the attempted transaction.
    pub postconditions: CylinderMediatedSupportContinuationPostconditions,
}

/// Postcondition status for a cylinder-mediated continuation transaction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CylinderMediatedSupportContinuationPostconditions {
    /// The topology was committed and all postconditions passed.
    Satisfied,
    /// The transaction was refused before mutation.
    Refused {
        /// Stable diagnostic reason.
        reason: &'static str,
    },
}

impl BRep {
    /// Transactionally apply a cylinder-mediated support continuation batch.
    ///
    /// This operator boundary is intentionally conservative: it revalidates the
    /// preflighted operands and aggregate partition against the current BRep,
    /// verifies the aggregate delta is known, then refuses before mutation until
    /// the topology rewiring and postconditions are implemented.
    pub fn cylinder_mediated_support_continuation_batch(
        &mut self,
        plan: &CylinderMediatedSupportContinuationBatchPlan,
    ) -> Result<CylinderMediatedSupportContinuationBatchResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let current = self.preflight_cylinder_mediated_support_continuation_batch(
            plan.plans
                .iter()
                .map(|continuation| continuation.request.clone())
                .collect(),
        )?;
        if &current != plan {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated support continuation plan is stale".into(),
            ));
        }
        require_known_aggregate_delta(&plan.expected_delta)?;
        require_non_empty_partition(&plan.outside_interval_partition)?;
        let segment_authoring = self.segment_authoring_plan(plan)?;
        let cylinder_resources = derive_cylinder_bridge_resource_partition(self, plan)?;
        let repair_result = {
            let mut candidate = self.clone();
            let repair_result = expose_owner_repair_spans(&mut candidate, &segment_authoring)?;
            let initial_cylinder_resource_remap =
                remap_cylinder_bridge_resources(&candidate, &cylinder_resources, plan)?;
            let blocker_exposure = expose_cylinder_resource_blockers(
                &mut candidate,
                &initial_cylinder_resource_remap,
            )?;
            let cylinder_resource_remap =
                remap_cylinder_bridge_resources(&candidate, &cylinder_resources, plan)?;
            let bridge_cycles =
                derive_segment_bridge_cycle_operands(&candidate, plan, &segment_authoring)?;
            let reconciliation = cylinder_resources.reconcile(&bridge_cycles);
            (
                repair_result,
                blocker_exposure,
                cylinder_resource_remap,
                bridge_cycles,
                reconciliation,
            )
        };
        let after = EulerTopologyCounts::from_brep(self);
        if before != after {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated support continuation refusal mutated topology".into(),
            ));
        }
        let bridge_refusal =
            bridge_refusal_reason(repair_result.3.ownership_status(), repair_result.4.status());
        Err(CoreError::PreconditionFailed(format!(
            "{bridge_refusal}: expected_delta={} authored_delta={} {} {} {} {} {} {} {}",
            plan.expected_delta.diagnostic_label(),
            authored_delta_label(
                repair_result.0.exposed_spans,
                repair_result.3.segment_count()
            ),
            segment_authoring_summary(&segment_authoring),
            repair_result.0.diagnostic_label(),
            repair_result.3.diagnostic_label(),
            cylinder_resources.diagnostic_label(),
            repair_result.1.diagnostic_label(),
            repair_result.2.diagnostic_label(),
            repair_result.4.diagnostic_label()
        )))
    }

    fn segment_authoring_plan(
        &self,
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
                    requirement:
                        CylinderMediatedSupportContinuationSegmentRequirement::BoundaryBridge,
                };
                for &half_edge in &segment.half_edges {
                    self.classify_segment_half_edge(plan, &mut segment_plan, half_edge)?;
                }
                segment_plan.unsupported_owner_faces =
                    unsupported_owner_faces(self, &segment_plan)?;
                segment_plan.context_owner_faces = context_owner_faces(self, &segment_plan)?;
                segment_plan.owner_repair_spans = owner_repair_spans(self, &segment_plan)?;
                segment_plan.requirement = segment_requirement(&segment_plan);
                authoring.push(segment_plan);
            }
        }
        Ok(authoring)
    }

    fn classify_segment_half_edge(
        &self,
        plan: &CylinderMediatedSupportContinuationBatchPlan,
        segment_plan: &mut CylinderMediatedSupportContinuationSegmentAuthoringPlan,
        half_edge: HalfEdgeHandle,
    ) -> Result<(), CoreError> {
        let data = self.half_edges.get(half_edge).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "cylinder-mediated segment references stale half-edge {half_edge:?}"
            ))
        })?;
        if self.boundary_edges.contains(&data.edge) && !self.half_edges.contains_key(data.twin) {
            segment_plan.boundary_bridge_half_edges.push(half_edge);
            return Ok(());
        }
        let Some(twin) = self.half_edges.get(data.twin) else {
            segment_plan.unsupported_owned_half_edges.push(half_edge);
            return Ok(());
        };
        let Some(twin_loop) = self.loops.get(twin.loop_) else {
            segment_plan.unassigned_twin_half_edges.push(half_edge);
            return Ok(());
        };
        let selected_cylinder = segment_plan.request_indices.iter().any(|&request_index| {
            plan.plans
                .get(request_index)
                .map(|continuation| continuation.request.cylinder_face == twin_loop.face)
                .unwrap_or(false)
        });
        if selected_cylinder {
            segment_plan.cylinder_owned_half_edges.push(half_edge);
            return Ok(());
        }
        let selected_outside = segment_plan.request_indices.iter().any(|&request_index| {
            plan.plans
                .get(request_index)
                .map(|continuation| continuation.request.outside_face == twin_loop.face)
                .unwrap_or(false)
        });
        if selected_outside {
            segment_plan.outside_owned_half_edges.push(half_edge);
            return Ok(());
        }
        let selected_support = segment_plan.request_indices.iter().any(|&request_index| {
            plan.plans
                .get(request_index)
                .map(|continuation| continuation.request.support_face == twin_loop.face)
                .unwrap_or(false)
        });
        if selected_support {
            segment_plan.support_owned_half_edges.push(half_edge);
            return Ok(());
        }
        let context_owner = segment_plan.request_indices.iter().any(|&request_index| {
            plan.plans
                .get(request_index)
                .map(|continuation| {
                    continuation
                        .request
                        .context_owner_faces
                        .contains(&twin_loop.face)
                })
                .unwrap_or(false)
        });
        if context_owner {
            segment_plan.context_owned_half_edges.push(half_edge);
        } else {
            segment_plan.unsupported_owned_half_edges.push(half_edge);
        }
        Ok(())
    }
}

fn unsupported_owner_faces(
    brep: &BRep,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<CylinderMediatedUnsupportedOwnerFace>, CoreError> {
    let mut counts = BTreeMap::<FaceHandle, usize>::new();
    for &half_edge in &segment.unsupported_owned_half_edges {
        let twin = brep
            .half_edges
            .get(half_edge)
            .and_then(|data| brep.half_edges.get(data.twin))
            .ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "unsupported owner half-edge {half_edge:?} has no live twin"
                ))
            })?;
        let face = brep
            .loops
            .get(twin.loop_)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(format!(
                    "unsupported owner half-edge {half_edge:?} has no live twin loop"
                ))
            })?;
        *counts.entry(face).or_default() += 1;
    }
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

fn bridge_refusal_reason(
    ownership_status: apply_bridge_cycles::SegmentBridgeCycleOwnershipStatus,
    resource_status: BridgeResourceReconciliationStatus,
) -> &'static str {
    if ownership_status
        == apply_bridge_cycles::SegmentBridgeCycleOwnershipStatus::MissingSupportOwner
    {
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
        BridgeResourceReconciliationStatus::IncompatibleResourceSets => {
            "cylinder-mediated support continuation bridge resources are incompatible"
        }
    }
}

fn require_known_aggregate_delta(
    delta: &CylinderMediatedSupportContinuationAggregateDelta,
) -> Result<(), CoreError> {
    match delta {
        CylinderMediatedSupportContinuationAggregateDelta::Known { .. } => Ok(()),
        CylinderMediatedSupportContinuationAggregateDelta::Unknown { reason } => {
            Err(CoreError::PreconditionFailed(format!(
                "cylinder-mediated support continuation aggregate delta is unknown: {reason}"
            )))
        }
    }
}

fn require_non_empty_partition(
    partition: &CylinderMediatedOutsideIntervalPartition,
) -> Result<(), CoreError> {
    if partition.loops.is_empty()
        || partition
            .loops
            .iter()
            .all(|loop_partition| loop_partition.segments.is_empty())
    {
        return Err(CoreError::PreconditionFailed(
            "cylinder-mediated support continuation requires a non-empty outside interval partition"
                .into(),
        ));
    }
    Ok(())
}

fn authored_delta_label(owner_repair_spans: usize, bridge_segments: usize) -> String {
    format!(
        "dV=0 dE={} dF={} dL={}",
        owner_repair_spans + bridge_segments * 2,
        bridge_segments,
        bridge_segments
    )
}
