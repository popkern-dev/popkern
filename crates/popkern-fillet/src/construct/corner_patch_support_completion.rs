//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use crate::types::{CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

pub(crate) fn complete_corner_patch_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
) {
    let support_cell_faces = corner_patches
        .iter()
        .filter(|patch| {
            patch
                .sub_patch_provenance
                .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        })
        .map(|patch| patch.face)
        .collect::<Vec<_>>();
    let _ =
        support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, &support_cell_faces);
    let _ = support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
        brep,
        &support_cell_faces,
    );
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
