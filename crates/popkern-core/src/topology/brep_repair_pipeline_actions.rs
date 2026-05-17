//! Action-report population helpers for the generic BRep repair pipeline.

use crate::topology::{
    BrepCurveTrimReconstructionFailureReason, BrepEdgeCurveCompletionFailureReason,
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityRepairFailureReason,
    BrepEdgeSurfaceCompatibilityRepairOutcome, BrepEdgeSurfaceCompatibilityRepairReport,
    BrepRepairActionDecision, BrepRepairActionOutcome, BrepRepairActionReport, BrepRepairFactId,
    BrepRepairPipelineStage, BrepRepairTopologyHandle, BrepSewingDeclineReason,
    BrepSnapSewingDeclineReason, BrepSnapSewingFailureReason, BrepSnapSewingReport,
    BrepTrimRepairFailureReason,
};

pub(super) fn append_snap_sewing_actions(
    actions: &mut BrepRepairActionReport,
    report: &BrepSnapSewingReport,
) {
    if let Some(clean) = report.clean_report.as_ref() {
        for &(first, second) in &clean.accepted_pair_edges {
            actions.push(BrepRepairActionDecision {
                stage: BrepRepairPipelineStage::SnapSewing,
                outcome: BrepRepairActionOutcome::Authored,
                subject: "sewn_boundary_pair".to_string(),
                reason: None,
                max_motion: report.max_vertex_motion,
                max_deviation: clean.max_surface_gap.max(clean.max_endpoint_gap),
                affected_handles: edge_pair_handles(first, second),
                shared_fact_ids: clean.shared_fact_ids.clone(),
            });
        }
        for decline in &clean.declined_pairs {
            actions.push(BrepRepairActionDecision {
                stage: BrepRepairPipelineStage::SnapSewing,
                outcome: BrepRepairActionOutcome::Declined,
                subject: format!(
                    "{}|{}",
                    decline.pair.provenance.first, decline.pair.provenance.second
                ),
                reason: Some(sewing_decline_reason_label(&decline.reason).to_string()),
                max_motion: report.max_vertex_motion,
                max_deviation: clean.max_surface_gap.max(clean.max_endpoint_gap),
                affected_handles: edge_pair_handles(decline.pair.first, decline.pair.second),
                shared_fact_ids: decline.shared_fact_ids.clone(),
            });
        }
    }
    for decline in &report.declined_pairs {
        actions.push(BrepRepairActionDecision {
            stage: BrepRepairPipelineStage::SnapSewing,
            outcome: BrepRepairActionOutcome::Declined,
            subject: format!(
                "{}|{}",
                decline.pair.provenance.first, decline.pair.provenance.second
            ),
            reason: Some(snap_decline_reason_label(&decline.reason)),
            max_motion: report.max_vertex_motion,
            max_deviation: 0.0,
            affected_handles: edge_pair_handles(decline.pair.first, decline.pair.second),
            shared_fact_ids: decline.shared_fact_ids.clone(),
        });
    }
}

pub(super) fn append_edge_surface_actions(
    actions: &mut BrepRepairActionReport,
    report: &BrepEdgeSurfaceCompatibilityRepairReport,
) {
    for decision in &report.decisions {
        let (outcome, reason) = match decision.outcome {
            BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(reason) => (
                BrepRepairActionOutcome::Declined,
                Some(edge_surface_decline_reason_label(reason).to_string()),
            ),
            _ => (BrepRepairActionOutcome::Authored, None),
        };
        actions.push(BrepRepairActionDecision {
            stage: BrepRepairPipelineStage::EdgeSurfaceCompatibility,
            outcome,
            subject: format!("{:?}", decision.edge),
            reason,
            max_motion: decision.endpoint_gap,
            max_deviation: decision.residual_deviation,
            affected_handles: vec![BrepRepairTopologyHandle::Edge(decision.edge)],
            shared_fact_ids: decision.shared_fact_ids.clone(),
        });
    }
}

pub(super) fn sewing_decline_reason_label(reason: &BrepSewingDeclineReason) -> &'static str {
    match reason {
        BrepSewingDeclineReason::SameEdge => "same_edge",
        BrepSewingDeclineReason::MissingEdge => "missing_edge",
        BrepSewingDeclineReason::EdgeAlreadyUsed => "edge_already_used",
        BrepSewingDeclineReason::NotBoundaryEdge => "not_boundary_edge",
        BrepSewingDeclineReason::NonManifoldConflict => "nonmanifold_conflict",
        BrepSewingDeclineReason::AlreadyJoined => "already_joined",
        BrepSewingDeclineReason::AmbiguousTopology => "ambiguous_topology",
        BrepSewingDeclineReason::EndpointGapExceeded { .. } => "endpoint_gap_exceeded",
        BrepSewingDeclineReason::MissingSurfaceSampleData => "missing_surface_sample_data",
        BrepSewingDeclineReason::SurfaceGapExceeded { .. } => "surface_gap_exceeded",
    }
}

pub(super) fn snap_decline_reason_label(reason: &BrepSnapSewingDeclineReason) -> String {
    match reason {
        BrepSnapSewingDeclineReason::RejectedClassification(inner) => {
            format!(
                "rejected_classification:{}",
                sewing_decline_reason_label(inner)
            )
        }
        BrepSnapSewingDeclineReason::MissingEdge => "missing_edge".to_string(),
        BrepSnapSewingDeclineReason::VertexMotionExceeded { .. } => {
            "vertex_motion_exceeded".to_string()
        }
    }
}

pub(super) fn edge_surface_decline_reason_label(
    reason: BrepEdgeSurfaceCompatibilityDeclineReason,
) -> &'static str {
    match reason {
        BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology => "missing_topology",
        BrepEdgeSurfaceCompatibilityDeclineReason::MissingAdjacentSurfaces => {
            "missing_adjacent_surfaces"
        }
        BrepEdgeSurfaceCompatibilityDeclineReason::DegenerateEdge => "degenerate_edge",
        BrepEdgeSurfaceCompatibilityDeclineReason::UnsupportedSurfacePair => {
            "unsupported_surface_pair"
        }
        BrepEdgeSurfaceCompatibilityDeclineReason::SurfaceBoundaryReconstructionFailed => {
            "surface_boundary_reconstruction_failed"
        }
        BrepEdgeSurfaceCompatibilityDeclineReason::IntersectionFailed => "intersection_failed",
        BrepEdgeSurfaceCompatibilityDeclineReason::EndpointMismatch => "endpoint_mismatch",
        BrepEdgeSurfaceCompatibilityDeclineReason::AnchorMotionExceeded => "anchor_motion_exceeded",
        BrepEdgeSurfaceCompatibilityDeclineReason::IncidentEndpointConflict => {
            "incident_endpoint_conflict"
        }
        BrepEdgeSurfaceCompatibilityDeclineReason::ResidualDeviation => "residual_deviation",
    }
}

pub(super) fn snap_sewing_failure_reason_label(
    reason: &BrepSnapSewingFailureReason,
) -> &'static str {
    match reason {
        BrepSnapSewingFailureReason::CleanSewingFailed(_) => "clean_sewing_failed",
    }
}

pub(super) fn edge_curve_completion_failure_reason_label(
    reason: &BrepEdgeCurveCompletionFailureReason,
) -> &'static str {
    match reason {
        BrepEdgeCurveCompletionFailureReason::UnresolvedMissingEdgeCurves => {
            "unresolved_missing_edge_curves"
        }
        BrepEdgeCurveCompletionFailureReason::ValidationFailed(_) => "validation_failed",
    }
}

pub(super) fn trim_repair_failure_reason_label(
    reason: &BrepTrimRepairFailureReason,
) -> &'static str {
    match reason {
        BrepTrimRepairFailureReason::UnresolvedTrimFindings => "unresolved_trim_findings",
        BrepTrimRepairFailureReason::ValidationFailed(_) => "validation_failed",
    }
}

pub(super) fn curve_trim_reconstruction_failure_reason_label(
    reason: &BrepCurveTrimReconstructionFailureReason,
) -> &'static str {
    match reason {
        BrepCurveTrimReconstructionFailureReason::UnresolvedTrimReconstruction => {
            "unresolved_curve_trim_reconstruction"
        }
        BrepCurveTrimReconstructionFailureReason::ValidationFailed(_) => "validation_failed",
    }
}

pub(super) fn edge_surface_compatibility_failure_reason_label(
    reason: &BrepEdgeSurfaceCompatibilityRepairFailureReason,
) -> &'static str {
    match reason {
        BrepEdgeSurfaceCompatibilityRepairFailureReason::UnresolvedIncompatibleEdges => {
            "unresolved_edge_surface_compatibility"
        }
        BrepEdgeSurfaceCompatibilityRepairFailureReason::ValidationFailed(_) => "validation_failed",
    }
}

pub(super) fn append_counted_stage_action(
    actions: &mut BrepRepairActionReport,
    stage: BrepRepairPipelineStage,
    authored: usize,
    declined: usize,
    subject: &str,
    shared_fact_ids: &[BrepRepairFactId],
) {
    for _ in 0..authored {
        actions.push(BrepRepairActionDecision {
            stage,
            outcome: BrepRepairActionOutcome::Authored,
            subject: subject.to_string(),
            reason: None,
            max_motion: 0.0,
            max_deviation: 0.0,
            affected_handles: Vec::new(),
            shared_fact_ids: shared_fact_ids.to_vec(),
        });
    }
    for _ in 0..declined {
        actions.push(BrepRepairActionDecision {
            stage,
            outcome: BrepRepairActionOutcome::Declined,
            subject: subject.to_string(),
            reason: Some("declined".to_string()),
            max_motion: 0.0,
            max_deviation: 0.0,
            affected_handles: Vec::new(),
            shared_fact_ids: shared_fact_ids.to_vec(),
        });
    }
}

fn edge_pair_handles(
    first: crate::topology::EdgeHandle,
    second: crate::topology::EdgeHandle,
) -> Vec<BrepRepairTopologyHandle> {
    vec![
        BrepRepairTopologyHandle::Edge(first),
        BrepRepairTopologyHandle::Edge(second),
    ]
}
