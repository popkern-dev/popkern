use popkern_core::geometry::gregory::GregoryPatch;

use crate::construct::kp20_polar_support_cell_inventory::{
    Kp20SupportCellInventory, Kp20SupportCellOwnerInventory,
    support_cell_inventory_from_patchwork_provenance,
    support_cell_owner_inventory_from_patchwork_provenance,
};
use crate::types::CornerSubPatchProvenance;

pub(super) fn kp20_support_cell_inventory(
    support_surfaces: &[GregoryPatch],
    n_sides: usize,
    provenance: Option<&[CornerSubPatchProvenance]>,
) -> Option<Kp20SupportCellResolvedInventory> {
    if let Some(provenance) = provenance {
        if let Ok(inventory) =
            support_cell_owner_inventory_from_patchwork_provenance(n_sides, provenance)
        {
            return Some(Kp20SupportCellResolvedInventory::Owner(inventory));
        }
        return support_cell_inventory_from_patchwork_provenance(n_sides, provenance)
            .ok()
            .map(Kp20SupportCellResolvedInventory::Legacy);
    }
    Kp20SupportCellInventory::full_ring(n_sides, support_surfaces.len())
        .map(Kp20SupportCellResolvedInventory::Legacy)
}

#[derive(Clone)]
pub(super) enum Kp20SupportCellResolvedInventory {
    Legacy(Kp20SupportCellInventory),
    Owner(Kp20SupportCellOwnerInventory),
}

impl Kp20SupportCellResolvedInventory {
    pub(super) fn cells_per_target_sector(&self) -> usize {
        match self {
            Self::Legacy(inventory) => inventory.cells_per_target_sector,
            Self::Owner(inventory) => inventory.cells_per_target_sector,
        }
    }

    pub(super) fn support_cell_surface_count(&self) -> usize {
        match self {
            Self::Legacy(inventory) => {
                inventory.target_sectors.len() * inventory.cells_per_target_sector
            }
            Self::Owner(inventory) => inventory.support_cell_surface_count(),
        }
    }
}
