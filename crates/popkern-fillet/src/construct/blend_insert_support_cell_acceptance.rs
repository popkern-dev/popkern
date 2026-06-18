//! Acceptance rule for support-cell corner patch multiplicity.

use std::collections::BTreeMap;

use crate::types::{CornerPatchInfo, CornerSubPatchRole};

pub(super) fn corner_patch_count_is_acceptable(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> bool {
    patches.len() == n_sides || support_cell_patch_count_is_acceptable(patches, n_sides)
}

pub(super) fn support_cell_patch_count_is_acceptable(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> bool {
    patches.len() == 4 * n_sides
        && patches.iter().all(|patch| {
            patch.n_sides == n_sides
                && patch
                    .sub_patch_provenance
                    .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        })
}

pub(super) fn partial_support_cell_coverage_detail(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> Option<String> {
    if patches.is_empty() {
        return None;
    }
    let mut sector_counts = BTreeMap::<usize, usize>::new();
    for patch in patches {
        let provenance = patch.sub_patch_provenance?;
        if patch.n_sides != n_sides
            || provenance.sector_count != n_sides
            || provenance.role != CornerSubPatchRole::SupportCell
        {
            return None;
        }
        *sector_counts.entry(provenance.sector_index).or_default() += 1;
    }
    let target_sectors = sector_counts
        .keys()
        .map(|sector| sector.to_string())
        .collect::<Vec<_>>()
        .join("|");
    let cells_per_target_sector = common_target_sector_cell_count(&sector_counts)
        .map(|count| count.to_string())
        .unwrap_or_else(|| "mixed".to_string());
    Some(format!(
        "expected_full_support_cell_patch_count={};standard_patch_count={n_sides};actual_patch_count={};target_count={};target_sectors={target_sectors};cells_per_target_sector={cells_per_target_sector}",
        4 * n_sides,
        patches.len(),
        sector_counts.len(),
    ))
}

fn common_target_sector_cell_count(sector_counts: &BTreeMap<usize, usize>) -> Option<usize> {
    let first = *sector_counts.values().next()?;
    sector_counts
        .values()
        .all(|&count| count == first)
        .then_some(first)
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

    #[test]
    fn summarizes_partial_support_cell_coverage() {
        let patches = [0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5]
            .into_iter()
            .enumerate()
            .map(|(index, sector)| support_cell_patch_in_sector(index, sector, 12))
            .collect::<Vec<_>>();

        let detail =
            partial_support_cell_coverage_detail(&patches, 6).expect("partial support-cell detail");

        assert!(detail.contains("expected_full_support_cell_patch_count=24"));
        assert!(detail.contains("actual_patch_count=12"));
        assert!(detail.contains("target_sectors=0|2|5"));
        assert!(detail.contains("cells_per_target_sector=4"));
    }

    fn support_cell_patch(index: usize) -> CornerPatchInfo {
        support_cell_patch_in_sector(index, index / 4, 24)
    }

    fn support_cell_patch_in_sector(
        index: usize,
        sector: usize,
        patch_count: usize,
    ) -> CornerPatchInfo {
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
                sector,
                6,
                index,
                patch_count,
            )),
        }
    }

    fn handles() -> (VertexHandle, FaceHandle) {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
        (mvfs.vertex, mvfs.face)
    }
}
