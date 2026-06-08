//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

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
        let before_brep = brep.clone();
        let before =
            shared_trim_quality::measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        let _ =
            support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, &group.faces);
        let _ = support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
            brep,
            &group.faces,
        );
        let after =
            shared_trim_quality::measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        let accepted = brep.validate().is_ok() && after.is_no_worse_than(before);
        if !accepted {
            *brep = before_brep;
        }
        diagnostics.push(shared_trim_quality::support_cell_shared_trim_repair_fact(
            group.n_sides,
            selected,
            accepted,
            before,
            after,
        ));
    }
}
