use popkern_core::topology::{
    BrepRepairActionArtifact, BrepRepairActionRecord, BrepRepairActionSummary,
};

use super::super::types_exact_boundary::ExactBoundaryRouteDecisionOutcome;
use super::{
    ApexConstructionDiagnostic, BlendRepairReport, support_seam_decline_reason_label,
    support_seam_rollback_reason_label,
};

/// Build a serialization-ready repair-action artifact from blend repair diagnostics.
pub fn build_blend_repair_action_artifact(report: &BlendRepairReport) -> BrepRepairActionArtifact {
    let seam = &report.support_seam_relocation;
    let rolled_back = usize::from(seam.rollback_reason.is_some());
    let mut actions: Vec<_> = seam
        .decisions
        .iter()
        .map(|decision| BrepRepairActionRecord {
            stage: "support_seam_relocation".to_string(),
            outcome: if decision.accepted {
                "authored".to_string()
            } else {
                "declined".to_string()
            },
            subject: format!(
                "cap={:?};half_edge={:?}",
                decision.cap_face, decision.half_edge
            ),
            reason: decision
                .reason
                .map(support_seam_decline_reason_label)
                .map(str::to_string),
            max_motion: decision.max_endpoint_motion,
            max_deviation: decision.max_sampled_disagreement,
            affected_handles: Vec::new(),
            shared_fact_ids: Vec::new(),
        })
        .collect();
    if let Some(reason) = seam.rollback_reason {
        actions.push(BrepRepairActionRecord {
            stage: "support_seam_relocation".to_string(),
            outcome: "rolled_back".to_string(),
            subject: "support_seam_relocation_transaction".to_string(),
            reason: Some(support_seam_rollback_reason_label(reason).to_string()),
            max_motion: seam.max_endpoint_motion,
            max_deviation: seam.max_sampled_disagreement,
            affected_handles: Vec::new(),
            shared_fact_ids: Vec::new(),
        });
    }
    if report.edge_curve_endpoint_coherence_repairs > 0 {
        actions.push(BrepRepairActionRecord {
            stage: "edge_curve_endpoint_coherence".to_string(),
            outcome: "authored".to_string(),
            subject: "edge_curve_endpoints".to_string(),
            reason: Some("stale_curve_endpoint".to_string()),
            max_motion: 0.0,
            max_deviation: 0.0,
            affected_handles: Vec::new(),
            shared_fact_ids: Vec::new(),
        });
    }
    let mut apex_attempted = 0usize;
    let mut apex_authored = 0usize;
    let mut apex_declined = 0usize;
    let mut apex_max_deviation = 0.0f64;
    append_route_decision_actions(report, &mut actions);
    append_apex_diagnostic_actions(
        report,
        &mut actions,
        &mut apex_attempted,
        &mut apex_authored,
        &mut apex_declined,
        &mut apex_max_deviation,
    );
    BrepRepairActionArtifact {
        summary: BrepRepairActionSummary {
            attempted_actions: seam.candidates
                + rolled_back
                + apex_attempted
                + report.edge_curve_endpoint_coherence_repairs,
            authored_actions: seam.authored
                + apex_authored
                + report.edge_curve_endpoint_coherence_repairs,
            declined_actions: seam.declined + apex_declined,
            rolled_back_actions: rolled_back,
            max_motion: seam.max_endpoint_motion,
            max_deviation: seam.max_sampled_disagreement.max(apex_max_deviation),
        },
        actions,
    }
}

fn append_route_decision_actions(
    report: &BlendRepairReport,
    actions: &mut Vec<BrepRepairActionRecord>,
) {
    for route in &report.apex_routes {
        for (index, step) in route
            .exact_boundary_report
            .decision_chain
            .iter()
            .enumerate()
        {
            actions.push(BrepRepairActionRecord {
                stage: "exact_boundary_cap.route_decision".to_string(),
                outcome: route_decision_outcome_label(step.outcome).to_string(),
                subject: format!(
                    "vertex={:?};decision_index={};route={}",
                    route.vertex, index, step.route_label
                ),
                reason: step.diagnostic_code.clone(),
                max_motion: 0.0,
                max_deviation: 0.0,
                affected_handles: route
                    .incident_edges
                    .iter()
                    .map(|edge| format!("edge:{edge:?}"))
                    .collect(),
                shared_fact_ids: route.exact_boundary_report.source_fact_ids.clone(),
            });
        }
    }
}

fn append_apex_diagnostic_actions(
    report: &BlendRepairReport,
    actions: &mut Vec<BrepRepairActionRecord>,
    apex_attempted: &mut usize,
    apex_authored: &mut usize,
    apex_declined: &mut usize,
    apex_max_deviation: &mut f64,
) {
    for diagnostic in &report.apex_diagnostics {
        match diagnostic {
            ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) => {
                *apex_attempted += 1;
                *apex_declined += 1;
                let deviation = diagnostic.residual.unwrap_or(0.0);
                *apex_max_deviation = apex_max_deviation.max(deviation);
                actions.push(BrepRepairActionRecord {
                    stage: diagnostic.stage.clone(),
                    outcome: "declined".to_string(),
                    subject: format!(
                        "sector={};edge={:?};support_half_edge={:?}",
                        diagnostic.sector_name, diagnostic.edge, diagnostic.support_half_edge
                    ),
                    reason: Some(diagnostic.code.clone()),
                    max_motion: 0.0,
                    max_deviation: deviation,
                    affected_handles: vec![
                        format!("edge:{:?}", diagnostic.edge),
                        format!("support_half_edge:{:?}", diagnostic.support_half_edge),
                        format!("support_face:{:?}", diagnostic.support_face),
                    ],
                    shared_fact_ids: diagnostic.shared_fact_ids.iter().map(|id| id.0).collect(),
                });
            }
            ApexConstructionDiagnostic::CylinderTrimMaterialization(_) => {}
            ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => {
                *apex_attempted += 1;
                *apex_declined += 1;
                let worst_g0_gap = diagnostic.continuity.worst_g0_gap.unwrap_or(0.0);
                let worst_deviation = worst_g0_gap.max(diagnostic.worst_g1_coplanarity);
                *apex_max_deviation = apex_max_deviation.max(worst_deviation);
                actions.push(BrepRepairActionRecord {
                    stage: "apex_cap_support_continuity".to_string(),
                    outcome: "declined".to_string(),
                    subject: format!(
                        "scope={};cap={:?};support={:?};edge={:?};sample={:?};cap_uv={:?};support_uv={:?};cap_curve_gap={:?};support_curve_gap={:?}",
                        diagnostic.measurement_scope,
                        diagnostic.cap_face,
                        diagnostic.support_face,
                        diagnostic.edge,
                        diagnostic.sample_parameter,
                        diagnostic.cap_uv,
                        diagnostic.support_uv,
                        diagnostic.cap_curve_gap,
                        diagnostic.support_curve_gap
                    ),
                    reason: diagnostic
                        .continuity
                        .decline_code
                        .clone()
                        .or_else(|| Some("cap_outer_g1_exceeds_gate".to_string())),
                    max_motion: 0.0,
                    max_deviation: worst_deviation,
                    affected_handles: vec![
                        format!("face:{:?}", diagnostic.cap_face),
                        format!("face:{:?}", diagnostic.support_face),
                        format!("edge:{:?}", diagnostic.edge),
                    ],
                    shared_fact_ids: Vec::new(),
                });
            }
            ApexConstructionDiagnostic::SideSupportContinuation(diagnostic) => {
                *apex_attempted += 1;
                if diagnostic.applied {
                    *apex_authored += 1;
                } else {
                    *apex_declined += 1;
                }
                actions.push(BrepRepairActionRecord {
                    stage: diagnostic.stage.clone(),
                    outcome: if diagnostic.applied {
                        "authored".to_string()
                    } else {
                        "declined".to_string()
                    },
                    subject: format!(
                        "cap={:?};primitive={};bridge_faces={};bridge_loops={}",
                        diagnostic.cap_face,
                        diagnostic.primitive,
                        diagnostic.bridge_face_count,
                        diagnostic.bridge_loop_count
                    ),
                    reason: Some(diagnostic.code.clone()),
                    max_motion: 0.0,
                    max_deviation: 0.0,
                    affected_handles: vec![format!("face:{:?}", diagnostic.cap_face)],
                    shared_fact_ids: Vec::new(),
                });
            }
            ApexConstructionDiagnostic::SupportCellSplit(_) => {}
            ApexConstructionDiagnostic::CoherenceValidation(_) => {}
        }
    }
}

fn route_decision_outcome_label(outcome: ExactBoundaryRouteDecisionOutcome) -> &'static str {
    match outcome {
        ExactBoundaryRouteDecisionOutcome::NotAttempted => "not_attempted",
        ExactBoundaryRouteDecisionOutcome::Accepted => "accepted",
        ExactBoundaryRouteDecisionOutcome::CommitEligible => "commit_eligible",
        ExactBoundaryRouteDecisionOutcome::Declined => "declined",
    }
}
