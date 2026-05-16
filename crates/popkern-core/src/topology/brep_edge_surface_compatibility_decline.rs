//! Decline constructors for edge/surface compatibility repair.

use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityRepairDecision,
    BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::{BrepRepairFactId, EdgeHandle};

pub(super) fn declined(
    edge: EdgeHandle,
    reason: BrepEdgeSurfaceCompatibilityDeclineReason,
    evidence: BrepEdgeSurfaceRepairEvidence,
    shared_fact_ids: Vec<BrepRepairFactId>,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(reason),
        endpoint_gap: 0.0,
        residual_deviation: 0.0,
        evidence,
        shared_fact_ids,
    }
}

pub(super) fn declined_without_evidence(
    edge: EdgeHandle,
    reason: BrepEdgeSurfaceCompatibilityDeclineReason,
    shared_fact_ids: Vec<BrepRepairFactId>,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(reason),
        endpoint_gap: 0.0,
        residual_deviation: 0.0,
        evidence: BrepEdgeSurfaceRepairEvidence::unavailable(),
        shared_fact_ids,
    }
}
