//! Acceptance rule for support-cell corner patch multiplicity.

use crate::types::{CornerPatchInfo, CornerSubPatchRole};

pub(super) fn corner_patch_count_is_acceptable(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> bool {
    patches.len() == n_sides || support_cell_patch_count_is_acceptable(patches, n_sides)
}

fn support_cell_patch_count_is_acceptable(patches: &[CornerPatchInfo], n_sides: usize) -> bool {
    patches.len() == 4 * n_sides
        && patches.iter().all(|patch| {
            patch.n_sides == n_sides
                && patch
                    .sub_patch_provenance
                    .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{ContinuityLevel, CornerMethod, CornerSubPatchProvenance};
    use crate::vertex_classifier::VertexBlendClassification;
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_core::topology::{FaceHandle, VertexHandle};

    #[test]
    fn accepts_four_support_cells_per_side() {
        let patches = (0..24)
            .map(|index| support_cell_patch(index))
            .collect::<Vec<_>>();

        assert!(corner_patch_count_is_acceptable(&patches, 6));
    }

    #[test]
    fn rejects_unprovenanced_patch_count_mismatch() {
        let patches = (0..24)
            .map(|index| {
                let mut patch = support_cell_patch(index);
                patch.sub_patch_provenance = None;
                patch
            })
            .collect::<Vec<_>>();

        assert!(!corner_patch_count_is_acceptable(&patches, 6));
    }

    fn support_cell_patch(index: usize) -> CornerPatchInfo {
        let (vertex, face) = handles();
        CornerPatchInfo {
            vertex,
            n_sides: 6,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification: VertexBlendClassification::SameSideConvex,
            face,
            bridge_fallbacks: Vec::new(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            sub_patch_provenance: Some(CornerSubPatchProvenance::kp20_support_cell(
                index / 4,
                6,
                index,
                24,
            )),
        }
    }

    fn handles() -> (VertexHandle, FaceHandle) {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
        (mvfs.vertex, mvfs.face)
    }
}
