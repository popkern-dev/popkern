//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

#[derive(Debug, Clone, Copy)]
enum SupportCellCompletionCandidate {
    PcurveOnly,
    BoundaryOnly,
    PcurveThenBoundary,
}

impl SupportCellCompletionCandidate {
    const ALL: [Self; 3] = [
        Self::PcurveOnly,
        Self::BoundaryOnly,
        Self::PcurveThenBoundary,
    ];

    fn label(self) -> &'static str {
        match self {
            Self::PcurveOnly => "pcurve_only",
            Self::BoundaryOnly => "boundary_only",
            Self::PcurveThenBoundary => "pcurve_then_boundary",
        }
    }
}

#[derive(Clone)]
struct SupportCellCompletionSelection {
    brep: BRep,
    candidate: SupportCellCompletionCandidate,
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
        let selection = select_support_cell_completion_candidate(brep, &group.faces, before);
        let (accepted, selected_candidate, after) = if let Some(selection) = selection {
            let candidate = selection.candidate.label();
            let quality = selection.quality;
            *brep = selection.brep;
            (true, candidate, quality)
        } else {
            (false, "none", before)
        };
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
) -> Option<SupportCellCompletionSelection> {
    let mut best = None;
    for candidate in SupportCellCompletionCandidate::ALL {
        let mut candidate_brep = brep.clone();
        if !apply_support_cell_completion_candidate(&mut candidate_brep, faces, candidate) {
            continue;
        }
        let quality = shared_trim_quality::measure_support_cell_shared_trim_quality(
            &candidate_brep,
            faces,
            9,
        );
        if candidate_brep.validate().is_err() || !quality.is_no_worse_than(before) {
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
    best
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
