//! Transactional boundary for cylinder-mediated support continuation.

use crate::error::CoreError;
use crate::topology::{BRep, BrepRepairFactId, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::contract::EulerTopologyCounts;
use super::cylinder_mediated_support::{
    CylinderMediatedSupportContinuationAggregateDelta, CylinderMediatedSupportContinuationBatchPlan,
};

#[path = "cylinder_mediated_support_apply_bridge_cycles.rs"]
mod apply_bridge_cycles;
#[path = "cylinder_mediated_support_apply_candidate_coverage.rs"]
mod apply_candidate_coverage;
#[path = "cylinder_mediated_support_apply_cylinder_resource_blockers.rs"]
mod apply_cylinder_resource_blockers;
#[path = "cylinder_mediated_support_apply_cylinder_resource_remap.rs"]
mod apply_cylinder_resource_remap;
#[path = "cylinder_mediated_support_apply_cylinder_resources.rs"]
mod apply_cylinder_resources;
#[path = "cylinder_mediated_support_apply_diagnostics.rs"]
mod apply_diagnostics;
#[path = "cylinder_mediated_support_apply_exposure_candidates.rs"]
mod apply_exposure_candidates;
#[path = "cylinder_mediated_support_apply_owner_spans.rs"]
mod apply_owner_spans;
#[path = "cylinder_mediated_support_apply_preconditions.rs"]
mod apply_preconditions;
#[path = "cylinder_mediated_support_apply_resource_candidates.rs"]
mod apply_resource_candidates;
#[path = "cylinder_mediated_support_apply_segment_authoring.rs"]
mod apply_segment_authoring;
#[path = "cylinder_mediated_support_apply_shared_start_bridge.rs"]
mod apply_shared_start_bridge;
#[path = "cylinder_mediated_support_apply_shared_start_bridge_mutation.rs"]
mod apply_shared_start_bridge_mutation;
#[path = "cylinder_mediated_support_apply_shared_start_cycle_mix.rs"]
mod apply_shared_start_cycle_mix;

use super::cylinder_mediated_support_facts::segment_source_fact_ids;
use apply_bridge_cycles::derive_segment_bridge_cycle_operands;
use apply_candidate_coverage::CylinderMediatedCandidateCoverage;
use apply_cylinder_resource_blockers::expose_cylinder_resource_blockers;
use apply_cylinder_resource_remap::remap_cylinder_bridge_resources;
use apply_cylinder_resources::{
    CylinderBridgeResourcePartition, derive_cylinder_bridge_resource_partition,
};
use apply_diagnostics::{
    authored_delta_label, bridge_refusal_reason, execution_delta_label, segment_authoring_summary,
    topology_delta_gap_label,
};
use apply_owner_spans::expose_owner_repair_spans;
use apply_preconditions::{require_known_aggregate_delta, require_non_empty_partition};
use apply_shared_start_bridge::{SharedStartBridgePreflight, preflight_shared_start_bridge};
use apply_shared_start_bridge_mutation::{
    SharedStartBridgeMutationStatus, apply_shared_start_bridge_mutation,
};

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
    /// Shared repair knowledge facts used to classify this segment.
    pub source_fact_ids: Vec<BrepRepairFactId>,
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
    /// Bridge faces produced by a committed shared-start continuation.
    pub bridge_faces: Vec<FaceHandle>,
    /// Bridge loops produced by a committed shared-start continuation.
    pub bridge_loops: Vec<LoopHandle>,
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
        let mut candidate = self.clone();
        let owner_exposure = expose_owner_repair_spans(&mut candidate, &segment_authoring)?;
        let initial_cylinder_resource_remap =
            remap_cylinder_bridge_resources(&candidate, &cylinder_resources, plan)?;
        let blocker_exposure =
            expose_cylinder_resource_blockers(&mut candidate, &initial_cylinder_resource_remap)?;
        let cylinder_resource_remap =
            remap_cylinder_bridge_resources(&candidate, &cylinder_resources, plan)?;
        let bridge_cycles =
            derive_segment_bridge_cycle_operands(&candidate, plan, &segment_authoring)?;
        let reconciliation = CylinderBridgeResourcePartition::from_remapped_segments(
            cylinder_resource_remap.segments(),
        )
        .reconcile(&candidate, &bridge_cycles)?;
        let shared_start_pairings = reconciliation.shared_start_bridge_pairings();
        let shared_start_bridge_preflight = if let Some(pairings) = shared_start_pairings.as_ref() {
            preflight_shared_start_bridge(
                &candidate,
                pairings,
                &segment_authoring,
                cylinder_resource_remap.segments(),
            )
        } else {
            SharedStartBridgePreflight::not_applicable()
        };
        let after = EulerTopologyCounts::from_brep(self);
        if before != after {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated support continuation refusal mutated topology".into(),
            ));
        }
        let candidate_coverage = CylinderMediatedCandidateCoverage::from_pipeline(
            owner_exposure,
            blocker_exposure,
            reconciliation.status(),
            &shared_start_bridge_preflight,
        )
        .require_complete()?;
        match apply_shared_start_bridge_mutation(
            &mut candidate,
            shared_start_pairings.as_deref(),
            &segment_authoring,
            cylinder_resource_remap.segments(),
        )? {
            SharedStartBridgeMutationStatus::Committed(result) => {
                candidate
                    .validate_half_edge_pairing()
                    .map_err(CoreError::from)?;
                candidate
                    .validate_half_edge_cycles()
                    .map_err(CoreError::from)?;
                candidate
                    .validate_manifold_edges()
                    .map_err(CoreError::from)?;
                *self = candidate;
                return Ok(CylinderMediatedSupportContinuationBatchResult {
                    expected_delta: plan.expected_delta.clone(),
                    segment_authoring,
                    bridge_faces: result.bridge_faces.to_vec(),
                    bridge_loops: result.bridge_loops.to_vec(),
                    postconditions: CylinderMediatedSupportContinuationPostconditions::Satisfied,
                });
            }
            SharedStartBridgeMutationStatus::NotApplicable => {}
        }
        if let Some(result) = reconciliation.apply_connector_authoring_mutation(
            &mut candidate,
            &bridge_cycles,
            cylinder_resource_remap.segments(),
        )? {
            candidate
                .validate_half_edge_pairing()
                .map_err(CoreError::from)?;
            candidate
                .validate_half_edge_cycles()
                .map_err(CoreError::from)?;
            candidate
                .validate_manifold_edges()
                .map_err(CoreError::from)?;
            *self = candidate;
            return Ok(CylinderMediatedSupportContinuationBatchResult {
                expected_delta: plan.expected_delta.clone(),
                segment_authoring,
                bridge_faces: result.bridge_faces,
                bridge_loops: result.bridge_loops,
                postconditions: CylinderMediatedSupportContinuationPostconditions::Satisfied,
            });
        }
        let bridge_refusal =
            bridge_refusal_reason(bridge_cycles.ownership_status(), reconciliation.status());
        Err(CoreError::PreconditionFailed(format!(
            "{bridge_refusal}: expected_delta={} authored_delta={} execution_delta={} {} {} {} {} {} {} {} {} {} {}",
            plan.expected_delta.diagnostic_label(),
            authored_delta_label(owner_exposure.exposed_spans, bridge_cycles.segment_count()),
            execution_delta_label(
                owner_exposure.exposed_spans,
                blocker_exposure.exposed_spans(),
                bridge_cycles.segment_count()
            ),
            topology_delta_gap_label(
                &plan.expected_delta,
                owner_exposure.exposed_spans,
                blocker_exposure.exposed_spans(),
                bridge_cycles.segment_count()
            ),
            segment_authoring_summary(&segment_authoring),
            owner_exposure.diagnostic_label(),
            bridge_cycles.diagnostic_label(),
            cylinder_resources.diagnostic_label(),
            blocker_exposure.diagnostic_label(),
            cylinder_resource_remap.diagnostic_label(),
            reconciliation.diagnostic_label(),
            shared_start_bridge_preflight.diagnostic_label(),
            candidate_coverage.diagnostic_label()
        )))
    }

    fn segment_authoring_plan(
        &self,
        plan: &CylinderMediatedSupportContinuationBatchPlan,
    ) -> Result<Vec<CylinderMediatedSupportContinuationSegmentAuthoringPlan>, CoreError> {
        apply_segment_authoring::segment_authoring_plan(self, plan)
    }
}
