//! Result-level completion for Gregory cap/support seam ownership.
#![allow(clippy::duplicate_mod)]

use std::borrow::Cow;

use popkern_core::error::ValidationError;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;
use popkern_core::topology::FaceHandle;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

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
    let mut completed_cap_faces = Vec::new();
    for patch in corner_patches {
        if patch
            .sub_patch_provenance
            .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        {
            continue;
        }
        if completed_cap_faces.contains(&patch.face) {
            continue;
        }
        if !face_is_gregory_cap(brep, patch.face) {
            continue;
        }
        let before = brep.clone();
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, patch.face,
        );
        if let Err(error) = brep.validate() {
            *brep = before;
            diagnostics.push(cap_support_completion_rollback_diagnostic(
                patch.face, &error,
            ));
        }
        completed_cap_faces.push(patch.face);
    }
}

fn face_is_gregory_cap(brep: &BRep, face: FaceHandle) -> bool {
    brep.faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Gregory(_)))
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
        if support_cell_completion_rank::selection_is_better(
            support_cell_rank(selection.candidate, selection.quality),
            best.as_ref()
                .map(|selection| support_cell_rank(selection.candidate, selection.quality)),
        ) {
            best = Some(selection);
        }
    }
    if let Some(selection) = best.as_ref() {
        for score in &mut scores {
            if score.candidate == selection.candidate {
                score.status = "accepted_selected";
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
