//! Acceptance rule for support-cell corner patch multiplicity.

use std::collections::BTreeMap;

use crate::types::{CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchRole};

pub(super) fn corner_patch_count_is_acceptable(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> bool {
    patches.len() == n_sides
        || support_cell_patch_count_is_acceptable(patches, n_sides)
        || mixed_support_cell_patch_count_is_acceptable(patches, n_sides)
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
                    .is_some_and(|provenance| provenance.role.is_support_cell())
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
            || !provenance.role.is_support_cell()
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

pub(super) fn mixed_support_cell_patch_count_is_acceptable(
    patches: &[CornerPatchInfo],
    n_sides: usize,
) -> bool {
    let mut sectors = vec![SectorCoverage::default(); n_sides];
    for patch in patches {
        let Some(provenance) = patch.sub_patch_provenance else {
            return false;
        };
        if patch.n_sides != n_sides
            || provenance.family != CornerSubPatchFamily::Kp20Bi6
            || provenance.sector_count != n_sides
            || provenance.sector_index >= n_sides
        {
            return false;
        }
        if !sectors[provenance.sector_index].record(provenance.role) {
            return false;
        }
    }
    sectors.iter().any(SectorCoverage::has_support_cells)
        && sectors.iter().all(SectorCoverage::is_complete)
}

fn common_target_sector_cell_count(sector_counts: &BTreeMap<usize, usize>) -> Option<usize> {
    let first = *sector_counts.values().next()?;
    sector_counts
        .values()
        .all(|&count| count == first)
        .then_some(first)
}

#[derive(Clone, Default)]
struct SectorCoverage {
    macro_roles: Vec<CornerSubPatchRole>,
    generic_support_cell_count: usize,
    support_cell_owner_counts: Vec<(CornerSubPatchRole, usize)>,
}

impl SectorCoverage {
    fn record(&mut self, role: CornerSubPatchRole) -> bool {
        if role.is_support_cell() {
            if let Some(owner_role) = role.support_cell_owner_role() {
                if let Some((_, count)) = self
                    .support_cell_owner_counts
                    .iter_mut()
                    .find(|(owner, _)| *owner == owner_role)
                {
                    *count += 1;
                } else {
                    self.support_cell_owner_counts.push((owner_role, 1));
                }
            } else {
                self.generic_support_cell_count += 1;
            }
            return true;
        }
        if !macro_support_cell_coverage_role(role) || self.macro_roles.contains(&role) {
            return false;
        }
        self.macro_roles.push(role);
        true
    }

    fn is_complete(&self) -> bool {
        if self.generic_support_cell_count == 0 && self.support_cell_owner_counts.is_empty() {
            return self.macro_roles.len() == 4;
        }
        if self.generic_support_cell_count > 0 {
            return self.support_cell_owner_counts.is_empty()
                && matches!(self.generic_support_cell_count, 4 | 8 | 16)
                && self.macro_roles.len() == 3
                && self
                    .macro_roles
                    .contains(&CornerSubPatchRole::NearOriginLeft)
                && self
                    .macro_roles
                    .contains(&CornerSubPatchRole::NearOriginRight)
                && (self
                    .macro_roles
                    .contains(&CornerSubPatchRole::NearMidpointLeft)
                    ^ self
                        .macro_roles
                        .contains(&CornerSubPatchRole::NearMidpointRight));
        }
        self.support_cell_owner_counts
            .iter()
            .all(|(_, count)| matches!(*count, 4 | 8 | 16))
            && self.macro_roles.len() + self.support_cell_owner_counts.len() == 4
            && macro_support_cell_owner_roles().iter().all(|role| {
                let replaced = self
                    .support_cell_owner_counts
                    .iter()
                    .any(|(owner, _)| owner == role);
                self.macro_roles.contains(role) != replaced
            })
    }

    fn support_cell_count(&self) -> usize {
        self.generic_support_cell_count
            + self
                .support_cell_owner_counts
                .iter()
                .map(|(_, count)| *count)
                .sum::<usize>()
    }

    fn has_support_cells(&self) -> bool {
        self.support_cell_count() > 0
    }
}

fn macro_support_cell_owner_roles() -> [CornerSubPatchRole; 4] {
    [
        CornerSubPatchRole::NearOriginLeft,
        CornerSubPatchRole::NearOriginRight,
        CornerSubPatchRole::NearMidpointRight,
        CornerSubPatchRole::NearMidpointLeft,
    ]
}

fn macro_support_cell_coverage_role(role: CornerSubPatchRole) -> bool {
    matches!(
        role,
        CornerSubPatchRole::NearOriginLeft
            | CornerSubPatchRole::NearOriginRight
            | CornerSubPatchRole::NearMidpointRight
            | CornerSubPatchRole::NearMidpointLeft
    )
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

    #[test]
    fn accepts_complete_mixed_support_cell_coverage() {
        let patches = mixed_support_cell_patches(3);

        assert!(mixed_support_cell_patch_count_is_acceptable(&patches, 6));
        assert!(corner_patch_count_is_acceptable(&patches, 6));
    }

    #[test]
    fn accepts_complete_mixed_support_cell_coverage_for_right_owner() {
        let patches = mixed_support_cell_patches(2);

        assert!(mixed_support_cell_patch_count_is_acceptable(&patches, 6));
        assert!(corner_patch_count_is_acceptable(&patches, 6));
    }

    #[test]
    fn accepts_complete_mixed_support_cell_coverage_for_two_owner_targets() {
        let mut patches = Vec::new();
        for sector in 0..6 {
            for piece in 0..4 {
                let owner = match (sector, piece) {
                    (0, 2) => Some(CornerSubPatchRole::NearMidpointRight),
                    (0, 3) => Some(CornerSubPatchRole::NearMidpointLeft),
                    _ => None,
                };
                if let Some(owner) = owner {
                    for index in 0..4 {
                        patches.push(patch_with_provenance(
                            CornerSubPatchProvenance::kp20_support_cell_for_owner_role(
                                sector, 6, index, 8, owner,
                            ),
                        ));
                    }
                } else {
                    patches.push(patch_with_provenance(CornerSubPatchProvenance::macro_cap(
                        CornerSubPatchFamily::Kp20Bi6,
                        sector,
                        6,
                        piece,
                    )));
                }
            }
        }

        assert!(mixed_support_cell_patch_count_is_acceptable(&patches, 6));
        assert!(corner_patch_count_is_acceptable(&patches, 6));
    }

    fn support_cell_patch(index: usize) -> CornerPatchInfo {
        support_cell_patch_in_sector(index, index / 4, 24)
    }

    fn support_cell_patch_in_sector(
        index: usize,
        sector: usize,
        patch_count: usize,
    ) -> CornerPatchInfo {
        patch_with_provenance(CornerSubPatchProvenance::kp20_support_cell(
            sector,
            6,
            index,
            patch_count,
        ))
    }

    fn mixed_support_cell_patches(replaced_piece: usize) -> Vec<CornerPatchInfo> {
        let mut patches = Vec::new();
        for sector in 0..6 {
            for piece in 0..4 {
                if piece == replaced_piece && matches!(sector, 0 | 2 | 5) {
                    for index in 0..4 {
                        patches.push(patch_with_provenance(
                            CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12),
                        ));
                    }
                } else {
                    patches.push(patch_with_provenance(CornerSubPatchProvenance::macro_cap(
                        CornerSubPatchFamily::Kp20Bi6,
                        sector,
                        6,
                        piece,
                    )));
                }
            }
        }
        patches
    }

    fn patch_with_provenance(provenance: CornerSubPatchProvenance) -> CornerPatchInfo {
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
            sub_patch_provenance: Some(provenance),
        }
    }

    fn handles() -> (VertexHandle, FaceHandle) {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
        (mvfs.vertex, mvfs.face)
    }
}
