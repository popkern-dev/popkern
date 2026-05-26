//! Result-level completion for Gregory cap/support seam ownership.

use popkern_core::topology::BRep;

use crate::types::CornerPatchInfo;

pub(crate) fn complete_corner_patch_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
) {
    let mut completed_cap_faces = Vec::new();
    for patch in corner_patches {
        if completed_cap_faces.contains(&patch.face) {
            continue;
        }
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, patch.face,
        );
        completed_cap_faces.push(patch.face);
    }
}
