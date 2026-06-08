//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SupportCellCompletionCandidate {
    PcurveOnly,
    BoundaryOnly,
    PcurveThenBoundary,
    PcurveThenBoundaryIfG0Safe,
}

impl SupportCellCompletionCandidate {
    const ALL: [Self; 4] = [
        Self::PcurveOnly,
        Self::BoundaryOnly,
        Self::PcurveThenBoundary,
        Self::PcurveThenBoundaryIfG0Safe,
    ];

    fn label(self) -> &'static str {
        match self {
            Self::PcurveOnly => "pcurve_only",
            Self::BoundaryOnly => "boundary_only",
            Self::PcurveThenBoundary => "pcurve_then_boundary",
            Self::PcurveThenBoundaryIfG0Safe => "pcurve_then_boundary_if_g0_safe",
        }
    }
}

#[derive(Clone)]
struct SupportCellCompletionSelection {
    brep: BRep,
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone, Copy)]
struct SupportCellCompletionCandidateScore {
    candidate: SupportCellCompletionCandidate,
    status: &'static str,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
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
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, patch.face,
        );
        completed_cap_faces.push(patch.face);
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
    let mut best = None;
    let mut scores = Vec::new();
    for candidate in SupportCellCompletionCandidate::ALL {
        let mut candidate_brep = brep.clone();
        if !apply_support_cell_completion_candidate(&mut candidate_brep, faces, candidate) {
            scores.push(SupportCellCompletionCandidateScore {
                candidate,
                status: "operator_failed",
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
        if support_cell_completion_selection_is_better(&selection, best.as_ref()) {
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

fn apply_support_cell_completion_candidate(
    brep: &mut BRep,
    faces: &[popkern_core::topology::FaceHandle],
    candidate: SupportCellCompletionCandidate,
) -> bool {
    match candidate {
        SupportCellCompletionCandidate::PcurveOnly => {
            support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_ok()
        }
        SupportCellCompletionCandidate::BoundaryOnly => {
            support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                brep, faces,
            )
            .is_ok()
        }
        SupportCellCompletionCandidate::PcurveThenBoundary => {
            support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_ok()
                && support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                    brep, faces,
                )
                .is_ok()
        }
        SupportCellCompletionCandidate::PcurveThenBoundaryIfG0Safe => {
            apply_boundary_after_pcurve_if_quality_safe(brep, faces)
        }
    }
}

fn apply_boundary_after_pcurve_if_quality_safe(
    brep: &mut BRep,
    faces: &[popkern_core::topology::FaceHandle],
) -> bool {
    if support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_err() {
        return false;
    }
    let after_pcurve =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let before_boundary = brep.clone();
    if support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(brep, faces)
        .is_err()
    {
        *brep = before_boundary;
        return true;
    }
    let after_boundary =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    if brep.validate().is_ok() && after_boundary.is_no_worse_than(after_pcurve) {
        true
    } else {
        *brep = before_boundary;
        true
    }
}

fn support_cell_completion_selection_is_better(
    candidate: &SupportCellCompletionSelection,
    current: Option<&SupportCellCompletionSelection>,
) -> bool {
    let Some(current) = current else {
        return true;
    };
    let a = candidate.quality;
    let b = current.quality;
    (a.pcurve_out_of_unit_count, a.degenerate_pcurve_count)
        .cmp(&(b.pcurve_out_of_unit_count, b.degenerate_pcurve_count))
        == std::cmp::Ordering::Less
        || (a.pcurve_out_of_unit_count == b.pcurve_out_of_unit_count
            && a.degenerate_pcurve_count == b.degenerate_pcurve_count
            && a.worst_g0_gap + 1.0e-12 < b.worst_g0_gap)
}
