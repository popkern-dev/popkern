use popkern_core::geometry::gregory::Kp20PieceId;

use crate::types::{CornerSubPatchProvenance, CornerSubPatchRole};

pub(super) fn support_cell_provenance(
    sector_index: usize,
    sector_count: usize,
    patch_index: usize,
    patch_count: usize,
    target_piece: Kp20PieceId,
) -> CornerSubPatchProvenance {
    match target_piece {
        Kp20PieceId::NearMidpointRight => {
            CornerSubPatchProvenance::kp20_support_cell_for_owner_role(
                sector_index,
                sector_count,
                patch_index,
                patch_count,
                CornerSubPatchRole::NearMidpointRight,
            )
        }
        Kp20PieceId::NearMidpointLeft => {
            CornerSubPatchProvenance::kp20_support_cell_for_owner_role(
                sector_index,
                sector_count,
                patch_index,
                patch_count,
                CornerSubPatchRole::NearMidpointLeft,
            )
        }
        _ => CornerSubPatchProvenance::kp20_support_cell(
            sector_index,
            sector_count,
            patch_index,
            patch_count,
        ),
    }
}
