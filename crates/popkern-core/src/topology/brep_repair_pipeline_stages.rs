//! Stage runners for the generic BRep repair pipeline.

use crate::topology::{
    BRep, BrepClassifiedEdgePair, BrepRepairActionDecision, BrepRepairActionOutcome,
    BrepRepairPipelineFailure, BrepRepairPipelineFailureReason, BrepRepairPipelineOptions,
    BrepRepairPipelineReport, BrepRepairPipelineStage, BrepSnapSewingReport,
    complete_missing_edge_curves_transactionally, reconstruct_curve_backed_trims_transactionally,
    repair_brep_edge_surface_compatibility_transactionally,
    repair_brep_face_use_trim_authority_transactionally, repair_brep_trims_transactionally,
    snap_and_sew_classified_edge_pairs,
};

use super::pipeline_actions;

pub(super) fn run_snap_sewing(
    brep: &BRep,
    classified_pairs: &[BrepClassifiedEdgePair],
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let sewn = snap_and_sew_classified_edge_pairs(brep, classified_pairs, options.snap_sewing)
        .map_err(|err| {
            let reason = BrepRepairPipelineFailureReason::SnapSewingFailed(
                pipeline_actions::snap_sewing_failure_reason_label(&err.reason).to_string(),
            );
            let mut failure_report = report.clone();
            pipeline_actions::append_snap_sewing_actions(&mut failure_report.actions, &err.report);
            failure_report.actions.push(BrepRepairActionDecision {
                stage: BrepRepairPipelineStage::SnapSewing,
                outcome: BrepRepairActionOutcome::RolledBack,
                subject: "snap_sewing_transaction".to_string(),
                reason: Some(
                    pipeline_actions::snap_sewing_failure_reason_label(&err.reason).to_string(),
                ),
                max_motion: err.report.max_vertex_motion,
                max_deviation: err
                    .report
                    .clean_report
                    .as_ref()
                    .map(|clean| clean.max_surface_gap.max(clean.max_endpoint_gap))
                    .unwrap_or(0.0),
                affected_handles: Vec::new(),
                shared_fact_ids: err.report.shared_fact_ids.clone(),
            });
            failure_report.snap_sewing = Some(err.report);
            failure(BrepRepairPipelineStage::SnapSewing, reason, failure_report)
        })?;
    pipeline_actions::append_snap_sewing_actions(&mut report.actions, &sewn.report);
    if options.require_all_pairs_repaired && snap_report_has_declines(&sewn.report) {
        report.snap_sewing = Some(sewn.report);
        return Err(failure(
            BrepRepairPipelineStage::SnapSewing,
            BrepRepairPipelineFailureReason::SnapSewingDeclined,
            report.clone(),
        ));
    }
    report.snap_sewing = Some(sewn.report);
    Ok(sewn.brep)
}

pub(super) fn run_edge_curve_completion(
    brep: &BRep,
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let repaired = complete_missing_edge_curves_transactionally(
        brep,
        options.edge_curve_completion,
    )
    .map_err(|err| {
        let reason = BrepRepairPipelineFailureReason::EdgeCurveCompletionFailed(
            pipeline_actions::edge_curve_completion_failure_reason_label(&err.reason).to_string(),
        );
        let mut failure_report = report.clone();
        failure_report.actions.push(BrepRepairActionDecision {
            stage: BrepRepairPipelineStage::EdgeCurveCompletion,
            outcome: BrepRepairActionOutcome::RolledBack,
            subject: "edge_curve_completion_transaction".to_string(),
            reason: Some(
                pipeline_actions::edge_curve_completion_failure_reason_label(&err.reason)
                    .to_string(),
            ),
            max_motion: 0.0,
            max_deviation: 0.0,
            affected_handles: Vec::new(),
            shared_fact_ids: err.report.shared_fact_ids.clone(),
        });
        failure_report.edge_curve_completion = Some(err.report);
        failure(
            BrepRepairPipelineStage::EdgeCurveCompletion,
            reason,
            failure_report,
        )
    })?;
    pipeline_actions::append_edge_curve_completion_actions(&mut report.actions, &repaired.report);
    report.edge_curve_completion = Some(repaired.report);
    Ok(repaired.brep)
}

pub(super) fn run_trim_repair(
    brep: &BRep,
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let repaired = repair_brep_trims_transactionally(brep, options.trim_repair).map_err(|err| {
        let reason = BrepRepairPipelineFailureReason::TrimRepairFailed(
            pipeline_actions::trim_repair_failure_reason_label(&err.reason).to_string(),
        );
        let mut failure_report = report.clone();
        failure_report.actions.push(BrepRepairActionDecision {
            stage: BrepRepairPipelineStage::TrimRepair,
            outcome: BrepRepairActionOutcome::RolledBack,
            subject: "trim_repair_transaction".to_string(),
            reason: Some(
                pipeline_actions::trim_repair_failure_reason_label(&err.reason).to_string(),
            ),
            max_motion: 0.0,
            max_deviation: 0.0,
            affected_handles: Vec::new(),
            shared_fact_ids: err.report.shared_fact_ids.clone(),
        });
        failure_report.trim_repair = Some(err.report);
        failure(BrepRepairPipelineStage::TrimRepair, reason, failure_report)
    })?;
    pipeline_actions::append_trim_repair_actions(&mut report.actions, &repaired.report);
    report.trim_repair = Some(repaired.report);
    Ok(repaired.brep)
}

pub(super) fn run_curve_trim_reconstruction(
    brep: &BRep,
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let repaired =
        reconstruct_curve_backed_trims_transactionally(brep, options.curve_trim_reconstruction)
            .map_err(|err| {
                let reason = BrepRepairPipelineFailureReason::CurveTrimReconstructionFailed(
                    pipeline_actions::curve_trim_reconstruction_failure_reason_label(&err.reason)
                        .to_string(),
                );
                let mut failure_report = report.clone();
                failure_report.actions.push(BrepRepairActionDecision {
                    stage: BrepRepairPipelineStage::CurveTrimReconstruction,
                    outcome: BrepRepairActionOutcome::RolledBack,
                    subject: "curve_trim_reconstruction_transaction".to_string(),
                    reason: Some(
                        pipeline_actions::curve_trim_reconstruction_failure_reason_label(
                            &err.reason,
                        )
                        .to_string(),
                    ),
                    max_motion: 0.0,
                    max_deviation: err.report.max_sampled_disagreement,
                    affected_handles: Vec::new(),
                    shared_fact_ids: err.report.shared_fact_ids.clone(),
                });
                failure_report.curve_trim_reconstruction = Some(err.report);
                failure(
                    BrepRepairPipelineStage::CurveTrimReconstruction,
                    reason,
                    failure_report,
                )
            })?;
    pipeline_actions::append_curve_trim_reconstruction_actions(
        &mut report.actions,
        &repaired.report,
    );
    report.curve_trim_reconstruction = Some(repaired.report);
    Ok(repaired.brep)
}

pub(super) fn run_edge_surface_compatibility(
    brep: &BRep,
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let repaired = repair_brep_edge_surface_compatibility_transactionally(
        brep,
        options.edge_surface_compatibility,
    )
    .map_err(|err| {
        let reason = BrepRepairPipelineFailureReason::EdgeSurfaceCompatibilityFailed(
            pipeline_actions::edge_surface_compatibility_failure_reason_label(&err.reason)
                .to_string(),
        );
        let mut failure_report = report.clone();
        pipeline_actions::append_edge_surface_actions(&mut failure_report.actions, &err.report);
        failure_report.actions.push(BrepRepairActionDecision {
            stage: BrepRepairPipelineStage::EdgeSurfaceCompatibility,
            outcome: BrepRepairActionOutcome::RolledBack,
            subject: "edge_surface_compatibility_transaction".to_string(),
            reason: Some(
                pipeline_actions::edge_surface_compatibility_failure_reason_label(&err.reason)
                    .to_string(),
            ),
            max_motion: err.report.max_endpoint_gap,
            max_deviation: err.report.after.max_deviation,
            affected_handles: Vec::new(),
            shared_fact_ids: err.report.shared_fact_ids.clone(),
        });
        failure_report.edge_surface_compatibility = Some(err.report);
        failure(
            BrepRepairPipelineStage::EdgeSurfaceCompatibility,
            reason,
            failure_report,
        )
    })?;
    pipeline_actions::append_edge_surface_actions(&mut report.actions, &repaired.report);
    report.edge_surface_compatibility = Some(repaired.report);
    Ok(repaired.brep)
}

pub(super) fn run_face_use_trim_authority(
    brep: &BRep,
    options: BrepRepairPipelineOptions,
    report: &mut BrepRepairPipelineReport,
) -> Result<BRep, Box<BrepRepairPipelineFailure>> {
    let repaired =
        repair_brep_face_use_trim_authority_transactionally(brep, options.face_use_trim_authority)
            .map_err(|err| {
                let reason = BrepRepairPipelineFailureReason::FaceUseTrimAuthorityFailed(
                    "unresolved_face_use_trim_authority".to_string(),
                );
                let mut failure_report = report.clone();
                pipeline_actions::append_face_use_trim_authority_actions(
                    &mut failure_report.actions,
                    &err.report,
                );
                failure_report.actions.push(BrepRepairActionDecision {
                    stage: BrepRepairPipelineStage::FaceUseTrimAuthority,
                    outcome: BrepRepairActionOutcome::RolledBack,
                    subject: "face_use_trim_authority_transaction".to_string(),
                    reason: Some("unresolved_face_use_trim_authority".to_string()),
                    max_motion: 0.0,
                    max_deviation: err
                        .report
                        .after
                        .max_endpoint_gap
                        .max(err.report.after.max_sampled_gap),
                    affected_handles: Vec::new(),
                    shared_fact_ids: Vec::new(),
                });
                failure_report.face_use_trim_authority = Some(err.report);
                failure(
                    BrepRepairPipelineStage::FaceUseTrimAuthority,
                    reason,
                    failure_report,
                )
            })?;
    pipeline_actions::append_face_use_trim_authority_actions(&mut report.actions, &repaired.report);
    report.face_use_trim_authority = Some(repaired.report);
    Ok(repaired.brep)
}

fn snap_report_has_declines(report: &BrepSnapSewingReport) -> bool {
    !report.declined_pairs.is_empty()
        || report
            .clean_report
            .as_ref()
            .is_some_and(|clean| !clean.declined_pairs.is_empty())
}

fn failure(
    stage: BrepRepairPipelineStage,
    reason: BrepRepairPipelineFailureReason,
    report: BrepRepairPipelineReport,
) -> Box<BrepRepairPipelineFailure> {
    Box::new(BrepRepairPipelineFailure {
        stage,
        reason,
        report,
    })
}
