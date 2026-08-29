//! Durable surface-boundary audit event recording for fillet construction.

use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::topology::{
    BRep, HalfEdgeHandle, SharedTrimEdgeFailureReason, SharedTrimEdgeValidation,
    SurfaceBoundaryAuditEvent, SurfaceBoundaryAuditOutcome, SurfaceBoundaryHintDeclineReason,
    SurfaceBoundarySharedTrimFailureReason, SurfaceBoundarySharedTrimValidation,
};

pub(crate) struct BoundaryAuditRecordInput {
    pub(crate) half_edge: HalfEdgeHandle,
    pub(crate) outcome: SurfaceBoundaryAuditOutcome,
    pub(crate) boundary: Option<BoundarySide>,
    pub(crate) decline_reason: Option<SurfaceBoundaryHintDeclineReason>,
    pub(crate) shared_trim_failure: Option<SurfaceBoundarySharedTrimFailureReason>,
    pub(crate) shared_trim_validation: Option<SurfaceBoundarySharedTrimValidation>,
    pub(crate) tolerance: Option<f64>,
}

impl Default for BoundaryAuditRecordInput {
    fn default() -> Self {
        Self {
            half_edge: HalfEdgeHandle::default(),
            outcome: SurfaceBoundaryAuditOutcome::ClassifiedDecline,
            boundary: None,
            decline_reason: None,
            shared_trim_failure: None,
            shared_trim_validation: None,
            tolerance: None,
        }
    }
}

pub(crate) fn record_boundary_audit(brep: &mut BRep, input: BoundaryAuditRecordInput) {
    let half_edge = input.half_edge;
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return;
    };
    let edge = half_edge_data.edge;
    let face = brep
        .loops
        .get(half_edge_data.loop_)
        .map(|loop_data| loop_data.face);
    let twin_face = brep
        .half_edges
        .get(half_edge_data.twin)
        .and_then(|twin| brep.loops.get(twin.loop_))
        .map(|loop_data| loop_data.face);
    brep.record_surface_boundary_audit_event(SurfaceBoundaryAuditEvent {
        half_edge,
        edge: Some(edge),
        face,
        twin_face,
        outcome: input.outcome,
        boundary: input.boundary,
        decline_reason: input.decline_reason,
        shared_trim_failure: input.shared_trim_failure,
        shared_trim_validation: input.shared_trim_validation,
        tolerance: input.tolerance,
        provenance: "gregory_analytic_boundary_audit".to_string(),
    });
}

pub(crate) fn surface_boundary_shared_trim_failure_reason(
    reason: SharedTrimEdgeFailureReason,
) -> SurfaceBoundarySharedTrimFailureReason {
    match reason {
        SharedTrimEdgeFailureReason::MissingTopology => {
            SurfaceBoundarySharedTrimFailureReason::MissingTopology
        }
        SharedTrimEdgeFailureReason::HalfEdgeEdgeMismatch => {
            SurfaceBoundarySharedTrimFailureReason::HalfEdgeEdgeMismatch
        }
        SharedTrimEdgeFailureReason::HalfEdgesNotTwins => {
            SurfaceBoundarySharedTrimFailureReason::HalfEdgesNotTwins
        }
        SharedTrimEdgeFailureReason::EndpointMismatch => {
            SurfaceBoundarySharedTrimFailureReason::EndpointMismatch
        }
        SharedTrimEdgeFailureReason::CurveSurfaceMismatch => {
            SurfaceBoundarySharedTrimFailureReason::CurveSurfaceMismatch
        }
        SharedTrimEdgeFailureReason::SideAgreementMismatch => {
            SurfaceBoundarySharedTrimFailureReason::SideAgreementMismatch
        }
        SharedTrimEdgeFailureReason::ValidationFailed => {
            SurfaceBoundarySharedTrimFailureReason::ValidationFailed
        }
    }
}

pub(crate) fn surface_boundary_shared_trim_validation(
    validation: SharedTrimEdgeValidation,
) -> SurfaceBoundarySharedTrimValidation {
    SurfaceBoundarySharedTrimValidation {
        max_curve_surface_deviation: validation.max_curve_surface_deviation,
        max_side_agreement_deviation: validation.max_side_agreement_deviation,
        max_endpoint_gap: validation.max_endpoint_gap,
    }
}
