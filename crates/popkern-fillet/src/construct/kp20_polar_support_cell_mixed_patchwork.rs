use popkern_core::geometry::gregory::{GregoryPatch, Kp20PieceId};

use crate::construct::kp20_polar_support_cell_inventory::{
    Kp20SupportCellInventory, support_cell_inventory_from_provenance,
};
use crate::types::{CornerSubPatchFamily, CornerSubPatchProvenance};

pub(super) struct Kp20MixedSupportCellPatchwork {
    pub(super) surfaces: Vec<GregoryPatch>,
    pub(super) provenance: Vec<CornerSubPatchProvenance>,
}

pub(super) fn build_mixed_support_cell_patchwork(
    sector_count: usize,
    macro_surfaces: &[GregoryPatch],
    support_cell_surfaces: &[GregoryPatch],
    support_cell_provenance: &[CornerSubPatchProvenance],
    target_piece: Kp20PieceId,
) -> Result<Kp20MixedSupportCellPatchwork, String> {
    if macro_surfaces.len() != 4 * sector_count {
        return Err(format!(
            "mixed support-cell patchwork expected {} macro surfaces for {sector_count} sectors, got {}",
            4 * sector_count,
            macro_surfaces.len(),
        ));
    }
    let inventory = support_cell_inventory_from_provenance(
        sector_count,
        support_cell_surfaces.len(),
        support_cell_provenance,
    )?;
    let expected =
        4 * sector_count + inventory.target_sectors.len() * (inventory.cells_per_target_sector - 1);
    let mut surfaces = Vec::with_capacity(expected);
    let mut provenance = Vec::with_capacity(expected);
    for sector in 0..sector_count {
        for piece in Kp20PieceId::all() {
            if piece == target_piece
                && let Some(range) = support_cell_source_range(&inventory, sector)
            {
                for source in range {
                    surfaces.push(support_cell_surfaces[source].clone());
                    provenance.push(support_cell_provenance[source]);
                }
                continue;
            }
            let piece_index = piece.index();
            surfaces.push(macro_surfaces[4 * sector + piece_index].clone());
            provenance.push(CornerSubPatchProvenance::macro_cap(
                CornerSubPatchFamily::Kp20Bi6,
                sector,
                sector_count,
                piece_index,
            ));
        }
    }
    Ok(Kp20MixedSupportCellPatchwork {
        surfaces,
        provenance,
    })
}

pub(super) fn trace_support_cell_patchwork_if_requested(surfaces: &[GregoryPatch]) {
    if std::env::var_os("POPKERN_DETERMINISM_TRACE").is_none() {
        return;
    }
    let mut hash = 0xcbf29ce484222325u64;
    for patch in surfaces {
        if let GregoryPatch::RectangularBiNonic(rectangular) = patch {
            for row in &rectangular.grid {
                for point in row {
                    for value in [point.x, point.y, point.z] {
                        hash ^= value.to_bits();
                        hash = hash.wrapping_mul(0x100000001b3);
                    }
                }
            }
        }
    }
    eprintln!("dtrace: lift_support_cell_patches {hash:016x}");
}

fn support_cell_source_range(
    inventory: &Kp20SupportCellInventory,
    sector: usize,
) -> Option<std::ops::Range<usize>> {
    inventory
        .target_sectors
        .iter()
        .position(|&target| target == sector)
        .map(|ordinal| {
            let start = ordinal * inventory.cells_per_target_sector;
            start..start + inventory.cells_per_target_sector
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

    #[test]
    fn mixed_patchwork_replaces_target_support_piece_by_cells() {
        let macro_surfaces = (0..24).map(dummy_surface).collect::<Vec<_>>();
        let support_surfaces = (0..12).map(dummy_surface).collect::<Vec<_>>();
        let support_provenance = [0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5]
            .into_iter()
            .enumerate()
            .map(|(index, sector)| {
                CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12)
            })
            .collect::<Vec<_>>();

        let patchwork = build_mixed_support_cell_patchwork(
            6,
            &macro_surfaces,
            &support_surfaces,
            &support_provenance,
            popkern_core::geometry::gregory::Kp20PieceId::NearMidpointLeft,
        )
        .expect("mixed patchwork");

        assert_eq!(patchwork.surfaces.len(), 33);
        assert_eq!(patchwork.provenance.len(), 33);
        assert_eq!(patchwork.provenance[3].sector_index, 0);
        assert_eq!(
            patchwork.provenance[3].role,
            crate::types::CornerSubPatchRole::SupportCell
        );
        assert_eq!(patchwork.provenance[14].sector_index, 2);
        assert_eq!(patchwork.provenance[32].sector_index, 5);
    }

    #[test]
    fn mixed_patchwork_can_replace_right_support_piece_by_cells() {
        let macro_surfaces = (0..24).map(dummy_surface).collect::<Vec<_>>();
        let support_surfaces = (0..12).map(dummy_surface).collect::<Vec<_>>();
        let support_provenance = [0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5]
            .into_iter()
            .enumerate()
            .map(|(index, sector)| {
                CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12)
            })
            .collect::<Vec<_>>();

        let patchwork = build_mixed_support_cell_patchwork(
            6,
            &macro_surfaces,
            &support_surfaces,
            &support_provenance,
            Kp20PieceId::NearMidpointRight,
        )
        .expect("mixed patchwork");

        assert_eq!(patchwork.surfaces.len(), 33);
        assert_eq!(
            patchwork.provenance[2].role,
            crate::types::CornerSubPatchRole::SupportCell
        );
        assert_eq!(
            patchwork.provenance[6].role,
            crate::types::CornerSubPatchRole::NearMidpointLeft
        );
    }

    fn dummy_surface(index: usize) -> GregoryPatch {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for row in &mut grid {
            for point in row {
                point.z = index as f64;
            }
        }
        GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
    }
}
