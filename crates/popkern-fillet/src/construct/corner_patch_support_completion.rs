//! Result-level completion for Gregory cap/support seam ownership.
#![allow(clippy::duplicate_mod)]

use std::{borrow::Cow, collections::BTreeSet};

use popkern_core::error::ValidationError;
use popkern_core::topology::BRep;
use popkern_core::topology::FaceHandle;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};

#[path = "corner_patch_support_pcurve_classification.rs"]
mod cap_support_pcurve_classification;
use cap_support_pcurve_classification::{
    ExternalCapSupportPcurveClassification, classify_external_cap_support_pcurves,
};
#[path = "kp20_polar_cap_support_boundary_solve.rs"]
mod cap_support_boundary_solve;
#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_completion_apply.rs"]
mod support_cell_completion_apply;
#[path = "kp20_polar_support_cell_completion_dispatch.rs"]
mod support_cell_completion_dispatch;
#[path = "kp20_polar_support_cell_completion_rank.rs"]
mod support_cell_completion_rank;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

#[path = "kp20_polar_support_cell_completion_candidates.rs"]
mod support_cell_completion_candidates;
use support_cell_completion_candidates::SupportCellCompletionCandidate;

#[derive(Clone)]
struct SupportCellCompletionSelection {
    brep: BRep,
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone)]
struct SupportCellCompletionCandidateScore {
    candidate: SupportCellCompletionCandidate,
    status: &'static str,
    operation_status: Cow<'static, str>,
    rank_reason: &'static str,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone)]
pub(super) struct SupportCellCompletionApplication {
    operator_succeeded: bool,
    operation_status: Cow<'static, str>,
}

pub(crate) fn complete_corner_patch_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    complete_support_cell_pcurves_transactionally(brep, corner_patches, diagnostics);
    complete_external_cap_support_pcurves(brep, corner_patches, diagnostics);
    emit_coherence_validation_diagnostic(brep, corner_patches, diagnostics);
}

/// Stage B1 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
/// coherence-invariants.md`: at the end of corner-patch support
/// completion, scope `validate_coherence` to the cap faces and emit
/// the result as an `ApexConstructionDiagnostic::CoherenceValidation`.
/// Stage B1 is observational — the report is recorded, no behavior
/// gating yet.
fn emit_coherence_validation_diagnostic(
    brep: &BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let cap_faces = surfaced_cap_face_set(brep, corner_patches);
    let report = brep.validate_coherence(&popkern_core::topology::CoherenceValidationScope::Faces(
        cap_faces,
    ));
    diagnostics.push(ApexConstructionDiagnostic::CoherenceValidation(
        crate::types::CoherenceValidationDiagnostic {
            stage: "corner_patch_support_completion".to_string(),
            report,
        },
    ));
}

pub(crate) fn complete_external_cap_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let cap_faces = surfaced_cap_face_set(brep, corner_patches);
    let mut eligible_edges = 0usize;
    let mut authored_cap_boundary_pcurves = 0usize;
    let mut repaired_cap_boundary_pcurves = 0usize;
    let mut inferred_cap_boundary_pcurves = 0usize;
    let mut declined_cap_boundary_pcurves = 0usize;
    let mut authored_support_pcurves = 0usize;
    let mut repaired_support_pcurves = 0usize;
    for cap_face in cap_faces.iter().copied() {
        let before = brep.clone();
        let completion_report =
            super::three_quad_support_pcurve_completion::complete_missing_external_cap_support_pcurves(
                brep, cap_face, &cap_faces,
            );
        eligible_edges += completion_report.eligible_edges;
        authored_cap_boundary_pcurves += completion_report.authored_cap_boundary_pcurves;
        repaired_cap_boundary_pcurves += completion_report.repaired_cap_boundary_pcurves;
        inferred_cap_boundary_pcurves += completion_report.inferred_cap_boundary_pcurves;
        declined_cap_boundary_pcurves += completion_report.declined_cap_boundary_pcurves;
        authored_support_pcurves += completion_report.authored_support_pcurves;
        repaired_support_pcurves += completion_report.repaired_support_pcurves;
        if let Err(error) = brep.validate() {
            *brep = before;
            diagnostics.push(cap_support_completion_rollback_diagnostic(cap_face, &error));
        }
    }
    let cap_pcurve_classification = classify_external_cap_support_pcurves(brep, &cap_faces);
    let before_boundary_solve = brep.clone();
    match cap_support_boundary_solve::solve_external_cap_support_boundaries(brep, &cap_faces) {
        Ok(report) => {
            let validation_error = brep.validate().err();
            let status_detail = validation_error
                .as_ref()
                .map(|error| format!("validation_error={error:?}"));
            let status = if validation_error.is_some() {
                *brep = before_boundary_solve;
                "rolled_back_invalid_brep"
            } else if let Some(status) = cap_support_boundary_solve_decline_status(report) {
                *brep = before_boundary_solve;
                status
            } else {
                "applied"
            };
            diagnostics.push(cap_support_boundary_solve_fact(
                report,
                status,
                status_detail.as_deref(),
            ));
        }
        Err(error) => {
            *brep = before_boundary_solve;
            diagnostics.push(cap_support_boundary_solve_error_fact(&error));
        }
    }
    diagnostics.push(cap_support_completion_fact(
        cap_faces.len(),
        eligible_edges,
        authored_cap_boundary_pcurves,
        repaired_cap_boundary_pcurves,
        inferred_cap_boundary_pcurves,
        declined_cap_boundary_pcurves,
        authored_support_pcurves,
        repaired_support_pcurves,
        cap_pcurve_classification,
    ));
}

fn surfaced_cap_face_set(brep: &BRep, corner_patches: &[CornerPatchInfo]) -> BTreeSet<FaceHandle> {
    corner_patches
        .iter()
        .filter_map(|patch| face_has_surface(brep, patch.face).then_some(patch.face))
        .collect()
}

fn face_has_surface(brep: &BRep, face: FaceHandle) -> bool {
    brep.faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .is_some()
}

fn cap_support_completion_rollback_diagnostic(
    cap_face: FaceHandle,
    error: &ValidationError,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(crate::types::ExactBoundaryCapDiagnostic {
        stage: "cap_support_completion".to_string(),
        code: "cap_support_completion_rolled_back".to_string(),
        detail_code: Some(format!("{error:?}")),
        sector_name: format!("{cap_face:?}"),
        edge: None,
        support_half_edge: None,
        support_face: Some(cap_face),
        residual: validation_gap(error),
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

fn cap_support_completion_fact(
    cap_face_count: usize,
    eligible_edges: usize,
    authored_cap_boundary_pcurves: usize,
    repaired_cap_boundary_pcurves: usize,
    inferred_cap_boundary_pcurves: usize,
    declined_cap_boundary_pcurves: usize,
    authored_support_pcurves: usize,
    repaired_support_pcurves: usize,
    cap_pcurve_classification: ExternalCapSupportPcurveClassification,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(crate::types::ExactBoundaryCapDiagnostic {
        stage: "cap_support_completion".to_string(),
        code: "cap_support_completion.fact".to_string(),
        detail_code: Some(format!(
            "cap_face_count={cap_face_count};eligible_edges={eligible_edges};authored_cap_boundary_pcurves={authored_cap_boundary_pcurves};repaired_cap_boundary_pcurves={repaired_cap_boundary_pcurves};inferred_cap_boundary_pcurves={inferred_cap_boundary_pcurves};declined_cap_boundary_pcurves={declined_cap_boundary_pcurves};authored_support_pcurves={authored_support_pcurves};repaired_support_pcurves={repaired_support_pcurves};unrecognized_cap_boundary_pcurve_edges={};unrecognized_cap_boundary_hint_edges={};recognized_cap_boundary_hint_edges={}",
            cap_pcurve_classification.unrecognized_cap_boundary_pcurve_edges,
            cap_pcurve_classification.unrecognized_cap_boundary_hint_edges,
            cap_pcurve_classification.recognized_cap_boundary_hint_edges,
        )),
        sector_name: "external_cap_support_pcurves".to_string(),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

fn cap_support_boundary_solve_decline_status(
    report: cap_support_boundary_solve::CapSupportBoundarySolveReport,
) -> Option<&'static str> {
    if report.sampled_edges == 0 {
        return None;
    }
    if report.displacement_unbounded_rows > 0
        || report.max_boundary_displacement
            > cap_support_boundary_solve::MAX_BOUNDARY_ROW_DISPLACEMENT
    {
        return Some("rolled_back_displacement_exceeded");
    }
    if report.max_after_gap > 5.0e-3 {
        return Some("rolled_back_postcondition_failed");
    }
    if report.max_before_gap > 5.0e-3 && report.max_after_gap + 1.0e-12 >= report.max_before_gap {
        return Some("rolled_back_quality_regressed");
    }
    None
}

fn cap_support_boundary_solve_fact(
    report: cap_support_boundary_solve::CapSupportBoundarySolveReport,
    status: &str,
    status_detail: Option<&str>,
) -> ApexConstructionDiagnostic {
    let status_detail = status_detail
        .map(|detail| format!(";{detail}"))
        .unwrap_or_default();
    ApexConstructionDiagnostic::ExactBoundaryCap(crate::types::ExactBoundaryCapDiagnostic {
        stage: "cap_support_boundary_solve".to_string(),
        code: "cap_support_boundary_solve.fact".to_string(),
        detail_code: Some(format!(
            "status={status}{status_detail};eligible_edges={};sampled_edges={};fit_rows={};row_fit_candidates={};constrained_row_fit_candidates={};constrained_row_singular_candidates={};max_hard_satisfying_constrained_displacement={:.3e};min_hard_satisfying_constrained_displacement={:.3e};hard_satisfying_rows_over_cap={};worst_hard_satisfying_face={:?};worst_hard_satisfying_side={:?};max_best_unscaled_constrained_residual={:.3e};hard_constraint_rows={};hard_constraint_samples={};hard_constraint_conflicts={};hard_constraint_overfull_rows={};zero_regularization_rows={};displacement_scaled_candidates={};displacement_rejected_candidates={};displacement_bounded_rows={};displacement_unbounded_rows={};multi_edge_rows={};max_samples_per_row={};reauthored_edge_curves={};missing_cap_pcurve_edges={};missing_support_pcurve_edges={};unrecognized_cap_pcurve_edges={};missing_support_surface_edges={};non_rectangular_cap_edges={};max_row_parameter_span={:.3e};max_row_before_sample_gap={:.3e};max_row_after_sample_gap={:.3e};max_row_after_sample_face={:?};max_row_after_sample_side={:?};max_row_after_sample_parameter={:.3e};max_row_after_sample_edge={:?};max_row_after_sample_cap_half_edge={:?};max_row_after_sample_support_face={:?};max_row_after_sample_support_parameter={:.3e};max_row_after_sample_target_source={:?};max_row_after_sample_before_gap={:.3e};max_row_after_sample_target_edge_curve_gap={:.3e};max_hard_constraint_residual={:.3e};max_row_target_edge_curve_gap={:.3e};max_row_target_edge_curve_edge={:?};max_row_target_edge_curve_face={:?};max_row_target_edge_curve_side={:?};max_row_target_edge_curve_parameter={:.3e};max_row_after_edge_curve_gap={:.3e};max_row_after_edge_curve_edge={:?};max_row_after_edge_curve_face={:?};max_row_after_edge_curve_side={:?};max_row_after_edge_curve_parameter={:.3e};max_before_gap={:.3e};max_after_gap={:.3e};max_boundary_displacement={:.3e}",
            report.eligible_edges,
            report.sampled_edges,
            report.fit_rows,
            report.row_fit_candidates,
            report.constrained_row_fit_candidates,
            report.constrained_row_singular_candidates,
            report.max_hard_satisfying_constrained_displacement,
            report.min_hard_satisfying_constrained_displacement,
            report.hard_satisfying_rows_over_cap,
            report.worst_hard_satisfying_face,
            report.worst_hard_satisfying_side,
            report.max_best_unscaled_constrained_residual,
            report.hard_constraint_rows,
            report.hard_constraint_samples,
            report.hard_constraint_conflicts,
            report.hard_constraint_overfull_rows,
            report.zero_regularization_rows,
            report.displacement_scaled_candidates,
            report.displacement_rejected_candidates,
            report.displacement_bounded_rows,
            report.displacement_unbounded_rows,
            report.multi_edge_rows,
            report.max_samples_per_row,
            report.reauthored_edge_curves,
            report.missing_cap_pcurve_edges,
            report.missing_support_pcurve_edges,
            report.unrecognized_cap_pcurve_edges,
            report.missing_support_surface_edges,
            report.non_rectangular_cap_edges,
            report.max_row_parameter_span,
            report.max_row_before_sample_gap,
            report.max_row_after_sample_gap,
            report.max_row_after_sample_face,
            report.max_row_after_sample_side,
            report.max_row_after_sample_parameter,
            report.max_row_after_sample_edge,
            report.max_row_after_sample_cap_half_edge,
            report.max_row_after_sample_support_face,
            report.max_row_after_sample_support_parameter,
            report.max_row_after_sample_target_source,
            report.max_row_after_sample_before_gap,
            report.max_row_after_sample_target_edge_curve_gap,
            report.max_hard_constraint_residual,
            report.max_row_target_edge_curve_gap,
            report.max_row_target_edge_curve_edge,
            report.max_row_target_edge_curve_face,
            report.max_row_target_edge_curve_side,
            report.max_row_target_edge_curve_parameter,
            report.max_row_after_edge_curve_gap,
            report.max_row_after_edge_curve_edge,
            report.max_row_after_edge_curve_face,
            report.max_row_after_edge_curve_side,
            report.max_row_after_edge_curve_parameter,
            report.max_before_gap,
            report.max_after_gap,
            report.max_boundary_displacement,
        )),
        sector_name: "external_cap_support_boundary_solve".to_string(),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: Some(report.max_after_gap),
        outer_side: None,
        sample_parameter: None,
        tolerance: Some(5.0e-3),
        shared_fact_ids: Vec::new(),
    })
}

fn cap_support_boundary_solve_error_fact(
    error: &popkern_core::error::CoreError,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(crate::types::ExactBoundaryCapDiagnostic {
        stage: "cap_support_boundary_solve".to_string(),
        code: "cap_support_boundary_solve.error".to_string(),
        detail_code: Some(format!("{error:?}")),
        sector_name: "external_cap_support_boundary_solve".to_string(),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance: Some(5.0e-3),
        shared_fact_ids: Vec::new(),
    })
}

fn validation_gap(error: &ValidationError) -> Option<f64> {
    match error {
        ValidationError::OpenLoop { gap, .. } => Some(*gap),
        ValidationError::ShellGap { gap, .. } => Some(*gap),
        _ => None,
    }
}

fn complete_support_cell_pcurves_transactionally(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let selected = shared_trim_quality::support_cell_selected_label(diagnostics);
    for group in shared_trim_quality::support_cell_face_groups(corner_patches) {
        let before =
            shared_trim_quality::measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        let (selection, scores) =
            select_support_cell_completion_candidate(brep, &group.faces, before);
        let (accepted, selected_candidate, after) = if let Some(selection) = selection {
            let candidate = selection.candidate.label();
            let quality = selection.quality;
            *brep = selection.brep;
            (true, candidate, quality)
        } else {
            (false, "none", before)
        };
        for score in scores {
            diagnostics.push(
                shared_trim_quality::support_cell_shared_trim_candidate_fact(
                    group.n_sides,
                    selected,
                    score.candidate.label(),
                    score.status,
                    score.operation_status.as_ref(),
                    score.rank_reason,
                    before,
                    score.quality,
                ),
            );
        }
        diagnostics.push(shared_trim_quality::support_cell_shared_trim_repair_fact(
            group.n_sides,
            selected,
            selected_candidate,
            accepted,
            before,
            after,
        ));
    }
}

fn select_support_cell_completion_candidate(
    brep: &BRep,
    faces: &[popkern_core::topology::FaceHandle],
    before: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> (
    Option<SupportCellCompletionSelection>,
    Vec<SupportCellCompletionCandidateScore>,
) {
    let mut best: Option<SupportCellCompletionSelection> = None;
    let mut scores = Vec::new();
    for candidate in SupportCellCompletionCandidate::ALL {
        let mut candidate_brep = brep.clone();

        let application = support_cell_completion_dispatch::apply_support_cell_completion_candidate(
            &mut candidate_brep,
            faces,
            candidate,
        );
        crate::construct::post_pipeline::determinism_trace(
            &format!("post_candidate_{}", candidate.label()),
            &candidate_brep,
        );
        if !application.operator_succeeded {
            scores.push(SupportCellCompletionCandidateScore {
                candidate,
                status: "operator_failed",
                operation_status: application.operation_status,
                rank_reason: "not_ranked",
                quality: before,
            });
            continue;
        }
        let quality = shared_trim_quality::measure_support_cell_shared_trim_quality(
            &candidate_brep,
            faces,
            9,
        );
        let status = if candidate_brep.validate().is_err() {
            "invalid_brep"
        } else if !quality.is_no_worse_than(before) {
            "quality_regressed"
        } else {
            "accepted_not_selected"
        };
        scores.push(SupportCellCompletionCandidateScore {
            candidate,
            status,
            operation_status: application.operation_status,
            rank_reason: "not_ranked",
            quality,
        });
        if status != "accepted_not_selected" {
            continue;
        }
        let selection = SupportCellCompletionSelection {
            brep: candidate_brep,
            candidate,
            quality,
        };
        let rank_decision = support_cell_completion_rank::selection_decision(
            support_cell_rank(selection.candidate, selection.quality),
            best.as_ref()
                .map(|selection| support_cell_rank(selection.candidate, selection.quality)),
        );
        if rank_decision.candidate_is_better {
            best = Some(selection);
        }
    }
    if let Some(selection) = best.as_ref() {
        for score in &mut scores {
            if score.candidate == selection.candidate {
                score.status = "accepted_selected";
                score.rank_reason = "selected_best_rank";
            } else if score.status == "accepted_not_selected" {
                score.rank_reason = support_cell_completion_rank::selection_decision(
                    support_cell_rank(score.candidate, score.quality),
                    Some(support_cell_rank(selection.candidate, selection.quality)),
                )
                .reason;
            }
        }
    }
    (best, scores)
}

fn support_cell_rank(
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> support_cell_completion_rank::SupportCellCompletionRank {
    support_cell_completion_rank::SupportCellCompletionRank {
        boundary_safe: candidate.is_boundary_safe(),
        quality,
    }
}
