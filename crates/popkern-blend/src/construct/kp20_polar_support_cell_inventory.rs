use popkern_core::geometry::gregory::Kp20PieceId;

use crate::types::{CornerSubPatchProvenance, CornerSubPatchRole};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Kp20SupportCellInventory {
    pub(crate) sector_count: usize,
    pub(crate) target_sectors: Vec<usize>,
    pub(crate) cells_per_target_sector: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20SupportCellOwnerTarget {
    pub(crate) sector: usize,
    pub(crate) target_piece: Kp20PieceId,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Kp20SupportCellOwnerInventory {
    pub(crate) sector_count: usize,
    pub(crate) targets: Vec<Kp20SupportCellOwnerTarget>,
    pub(crate) cells_per_target_sector: usize,
}

impl Kp20SupportCellInventory {
    pub(crate) fn full_ring(sector_count: usize, surface_count: usize) -> Option<Self> {
        if sector_count < 3 || surface_count == 0 || !surface_count.is_multiple_of(sector_count) {
            return None;
        }
        let cells_per_target_sector = surface_count / sector_count;
        support_cell_cells_per_target_sector_is_supported(cells_per_target_sector).then(|| Self {
            sector_count,
            target_sectors: (0..sector_count).collect(),
            cells_per_target_sector,
        })
    }
}

impl Kp20SupportCellOwnerInventory {
    pub(crate) fn target_range(
        &self,
        sector: usize,
        target_piece: Kp20PieceId,
    ) -> Option<std::ops::Range<usize>> {
        self.targets
            .iter()
            .position(|target| target.sector == sector && target.target_piece == target_piece)
            .map(|ordinal| {
                let start = ordinal * self.cells_per_target_sector;
                start..start + self.cells_per_target_sector
            })
    }

    pub(crate) fn support_cell_surface_count(&self) -> usize {
        self.targets.len() * self.cells_per_target_sector
    }
}

pub(crate) fn support_cell_inventory_from_provenance(
    sector_count: usize,
    surface_count: usize,
    provenance: &[CornerSubPatchProvenance],
) -> Result<Kp20SupportCellInventory, String> {
    if sector_count < 3 {
        return Err(format!(
            "kp20 support-cell inventory expected at least 3 sectors, got {sector_count}"
        ));
    }
    if surface_count == 0 || provenance.len() != surface_count {
        return Err(format!(
            "kp20 support-cell inventory expected one provenance item per surface, got surfaces={surface_count} provenance={}",
            provenance.len()
        ));
    }

    let mut target_sectors = Vec::new();
    let mut group_counts = Vec::new();
    let mut closed_sectors = Vec::new();
    let mut current_sector = None;
    let mut owner_piece = None;
    let mut saw_generic_support_cell = false;
    for (source_index, item) in provenance.iter().enumerate() {
        validate_support_cell_provenance(sector_count, surface_count, source_index, item)?;
        match support_cell_owner_piece(item.role) {
            Some(item_owner) => {
                if saw_generic_support_cell {
                    return Err(
                        "kp20 support-cell inventory cannot mix generic and owner-specific support-cell provenance"
                            .to_string(),
                    );
                }
                if let Some(existing) = owner_piece {
                    if existing != item_owner {
                        return Err(format!(
                            "kp20 support-cell inventory requires a single owner piece, got {existing:?} and {item_owner:?}"
                        ));
                    }
                } else {
                    owner_piece = Some(item_owner);
                }
            }
            None => {
                if owner_piece.is_some() {
                    return Err(
                        "kp20 support-cell inventory cannot mix generic and owner-specific support-cell provenance"
                            .to_string(),
                    );
                }
                saw_generic_support_cell = true;
            }
        }
        if current_sector == Some(item.sector_index) {
            if let Some(count) = group_counts.last_mut() {
                *count += 1;
            }
            continue;
        }
        if let Some(previous) = current_sector {
            closed_sectors.push(previous);
        }
        if closed_sectors.contains(&item.sector_index)
            || target_sectors.contains(&item.sector_index)
        {
            return Err(format!(
                "kp20 support-cell inventory sector {} is not stored as one contiguous source group",
                item.sector_index
            ));
        }
        current_sector = Some(item.sector_index);
        target_sectors.push(item.sector_index);
        group_counts.push(1);
    }

    let Some(&cells_per_target_sector) = group_counts.first() else {
        return Err("kp20 support-cell inventory has no target sectors".to_string());
    };
    if !support_cell_cells_per_target_sector_is_supported(cells_per_target_sector) {
        return Err(format!(
            "kp20 support-cell inventory does not support {cells_per_target_sector} cells per target sector"
        ));
    }
    if group_counts
        .iter()
        .any(|&count| count != cells_per_target_sector)
    {
        return Err(format!(
            "kp20 support-cell inventory expected {cells_per_target_sector} cells for every target sector, got {group_counts:?}"
        ));
    }
    if target_sectors.len() * cells_per_target_sector != surface_count {
        return Err(format!(
            "kp20 support-cell inventory target sector count {} and cells per sector {cells_per_target_sector} do not cover {surface_count} surfaces",
            target_sectors.len()
        ));
    }

    Ok(Kp20SupportCellInventory {
        sector_count,
        target_sectors,
        cells_per_target_sector,
    })
}

pub(crate) fn support_cell_owner_inventory_from_provenance(
    sector_count: usize,
    surface_count: usize,
    provenance: &[CornerSubPatchProvenance],
) -> Result<Kp20SupportCellOwnerInventory, String> {
    if sector_count < 3 {
        return Err(format!(
            "kp20 owner support-cell inventory expected at least 3 sectors, got {sector_count}"
        ));
    }
    if surface_count == 0 || provenance.len() != surface_count {
        return Err(format!(
            "kp20 owner support-cell inventory expected one provenance item per surface, got surfaces={surface_count} provenance={}",
            provenance.len()
        ));
    }

    let mut targets = Vec::new();
    let mut group_counts = Vec::new();
    let mut closed_targets = Vec::new();
    let mut current_target = None;
    for (source_index, item) in provenance.iter().enumerate() {
        validate_owner_support_cell_provenance(sector_count, surface_count, source_index, item)?;
        let target = Kp20SupportCellOwnerTarget {
            sector: item.sector_index,
            target_piece: support_cell_owner_piece(item.role).ok_or_else(|| {
                format!(
                    "kp20 owner support-cell inventory source {source_index} lacks an owner-specific support-cell role"
                )
            })?,
        };
        if current_target == Some(target) {
            if let Some(count) = group_counts.last_mut() {
                *count += 1;
            }
            continue;
        }
        if let Some(previous) = current_target {
            closed_targets.push(previous);
        }
        if closed_targets.contains(&target) || targets.contains(&target) {
            return Err(format!(
                "kp20 owner support-cell inventory target sector {} piece {:?} is not stored as one contiguous source group",
                target.sector, target.target_piece
            ));
        }
        current_target = Some(target);
        targets.push(target);
        group_counts.push(1);
    }

    let Some(&cells_per_target_sector) = group_counts.first() else {
        return Err("kp20 owner support-cell inventory has no target groups".to_string());
    };
    if !support_cell_cells_per_target_sector_is_supported(cells_per_target_sector) {
        return Err(format!(
            "kp20 owner support-cell inventory does not support {cells_per_target_sector} cells per target"
        ));
    }
    if group_counts
        .iter()
        .any(|&count| count != cells_per_target_sector)
    {
        return Err(format!(
            "kp20 owner support-cell inventory expected {cells_per_target_sector} cells for every target, got {group_counts:?}"
        ));
    }
    if targets.len() * cells_per_target_sector != surface_count {
        return Err(format!(
            "kp20 owner support-cell inventory target count {} and cells per target {cells_per_target_sector} do not cover {surface_count} surfaces",
            targets.len()
        ));
    }

    Ok(Kp20SupportCellOwnerInventory {
        sector_count,
        targets,
        cells_per_target_sector,
    })
}

pub(crate) fn support_cell_inventory_from_patchwork_provenance(
    sector_count: usize,
    provenance: &[CornerSubPatchProvenance],
) -> Result<Kp20SupportCellInventory, String> {
    let support_cell_provenance = provenance
        .iter()
        .copied()
        .filter(|item| item.role.is_support_cell())
        .collect::<Vec<_>>();
    support_cell_inventory_from_provenance(
        sector_count,
        support_cell_provenance.len(),
        &support_cell_provenance,
    )
}

pub(crate) fn support_cell_owner_inventory_from_patchwork_provenance(
    sector_count: usize,
    provenance: &[CornerSubPatchProvenance],
) -> Result<Kp20SupportCellOwnerInventory, String> {
    let support_cell_provenance = provenance
        .iter()
        .copied()
        .filter(|item| support_cell_owner_piece(item.role).is_some())
        .collect::<Vec<_>>();
    support_cell_owner_inventory_from_provenance(
        sector_count,
        support_cell_provenance.len(),
        &support_cell_provenance,
    )
}

pub(crate) fn support_cell_cells_per_target_sector_is_supported(cells: usize) -> bool {
    matches!(cells, 4 | 8 | 16)
}

pub(crate) fn support_cell_owner_piece(role: CornerSubPatchRole) -> Option<Kp20PieceId> {
    match role {
        CornerSubPatchRole::SupportCellNearMidpointRight => Some(Kp20PieceId::NearMidpointRight),
        CornerSubPatchRole::SupportCellNearMidpointLeft => Some(Kp20PieceId::NearMidpointLeft),
        _ => None,
    }
}

fn validate_support_cell_provenance(
    sector_count: usize,
    surface_count: usize,
    source_index: usize,
    item: &CornerSubPatchProvenance,
) -> Result<(), String> {
    if !item.role.is_support_cell() {
        return Err(format!(
            "kp20 support-cell inventory source {source_index} has non-support-cell role {:?}",
            item.role
        ));
    }
    if item.sector_count != sector_count || item.sector_index >= sector_count {
        return Err(format!(
            "kp20 support-cell inventory source {source_index} has sector {}/{} for ring sector count {sector_count}",
            item.sector_index, item.sector_count
        ));
    }
    if item.piece_count != surface_count {
        return Err(format!(
            "kp20 support-cell inventory source {source_index} has piece_count {}, expected {surface_count}",
            item.piece_count
        ));
    }
    Ok(())
}

fn validate_owner_support_cell_provenance(
    sector_count: usize,
    surface_count: usize,
    source_index: usize,
    item: &CornerSubPatchProvenance,
) -> Result<(), String> {
    validate_support_cell_provenance(sector_count, surface_count, source_index, item)?;
    if support_cell_owner_piece(item.role).is_none() {
        return Err(format!(
            "kp20 owner support-cell inventory source {source_index} has generic support-cell role"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::CornerSubPatchProvenance;

    #[test]
    fn validates_sparse_target_sector_groups() {
        let provenance = [0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5]
            .into_iter()
            .enumerate()
            .map(|(index, sector)| {
                CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12)
            })
            .collect::<Vec<_>>();

        let inventory = support_cell_inventory_from_provenance(6, provenance.len(), &provenance)
            .expect("valid sparse support-cell inventory");

        assert_eq!(inventory.sector_count, 6);
        assert_eq!(inventory.target_sectors, vec![0, 2, 5]);
        assert_eq!(inventory.cells_per_target_sector, 4);
    }

    #[test]
    fn rejects_noncontiguous_sector_sources() {
        let provenance = [0, 0, 2, 2, 0, 0, 5, 5, 5, 5, 2, 2]
            .into_iter()
            .enumerate()
            .map(|(index, sector)| {
                CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12)
            })
            .collect::<Vec<_>>();

        assert!(support_cell_inventory_from_provenance(6, provenance.len(), &provenance).is_err());
    }

    #[test]
    fn extracts_support_cell_inventory_from_mixed_patchwork() {
        let mut provenance = Vec::new();
        for sector in 0..6 {
            for piece in 0..3 {
                provenance.push(CornerSubPatchProvenance::macro_cap(
                    crate::types::CornerSubPatchFamily::Kp20Bi6,
                    sector,
                    6,
                    piece,
                ));
            }
            if matches!(sector, 0 | 2 | 5) {
                for index in 0..4 {
                    provenance.push(CornerSubPatchProvenance::kp20_support_cell(
                        sector, 6, index, 12,
                    ));
                }
            } else {
                provenance.push(CornerSubPatchProvenance::macro_cap(
                    crate::types::CornerSubPatchFamily::Kp20Bi6,
                    sector,
                    6,
                    3,
                ));
            }
        }

        let inventory = support_cell_inventory_from_patchwork_provenance(6, &provenance)
            .expect("support-cell inventory from mixed patchwork");

        assert_eq!(inventory.target_sectors, vec![0, 2, 5]);
        assert_eq!(inventory.cells_per_target_sector, 4);
    }

    #[test]
    fn validates_owner_specific_target_groups() {
        let provenance = [
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointRight),
            (0, CornerSubPatchRole::NearMidpointRight),
            (0, CornerSubPatchRole::NearMidpointRight),
            (0, CornerSubPatchRole::NearMidpointRight),
            (2, CornerSubPatchRole::NearMidpointLeft),
            (2, CornerSubPatchRole::NearMidpointLeft),
            (2, CornerSubPatchRole::NearMidpointLeft),
            (2, CornerSubPatchRole::NearMidpointLeft),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (sector, owner))| {
            CornerSubPatchProvenance::kp20_support_cell_for_owner_role(sector, 6, index, 12, owner)
        })
        .collect::<Vec<_>>();

        let inventory =
            support_cell_owner_inventory_from_provenance(6, provenance.len(), &provenance)
                .expect("valid owner-specific support-cell inventory");

        assert_eq!(inventory.sector_count, 6);
        assert_eq!(inventory.cells_per_target_sector, 4);
        assert_eq!(
            inventory.targets,
            vec![
                Kp20SupportCellOwnerTarget {
                    sector: 0,
                    target_piece: Kp20PieceId::NearMidpointLeft,
                },
                Kp20SupportCellOwnerTarget {
                    sector: 0,
                    target_piece: Kp20PieceId::NearMidpointRight,
                },
                Kp20SupportCellOwnerTarget {
                    sector: 2,
                    target_piece: Kp20PieceId::NearMidpointLeft,
                },
            ]
        );
        assert_eq!(
            inventory.target_range(0, Kp20PieceId::NearMidpointRight),
            Some(4..8)
        );
    }

    #[test]
    fn single_owner_inventory_rejects_mixed_owner_provenance() {
        let provenance = [
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (0, CornerSubPatchRole::NearMidpointLeft),
            (2, CornerSubPatchRole::NearMidpointRight),
            (2, CornerSubPatchRole::NearMidpointRight),
            (2, CornerSubPatchRole::NearMidpointRight),
            (2, CornerSubPatchRole::NearMidpointRight),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (sector, owner))| {
            CornerSubPatchProvenance::kp20_support_cell_for_owner_role(sector, 6, index, 8, owner)
        })
        .collect::<Vec<_>>();

        assert!(support_cell_inventory_from_provenance(6, provenance.len(), &provenance).is_err());
    }
}
