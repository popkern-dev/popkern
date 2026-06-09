pub(super) fn support_cell_piece_count(
    patches: &[popkern_fillet::CornerPatchInfo],
) -> Option<usize> {
    patches.iter().try_fold(None, |current, patch| {
        let provenance = patch.sub_patch_provenance?;
        if !matches!(
            provenance.role,
            popkern_fillet::CornerSubPatchRole::SupportCell
        ) {
            return None;
        }
        match current {
            Some(piece_count) if piece_count != provenance.piece_count => None,
            Some(piece_count) => Some(Some(piece_count)),
            None => Some(Some(provenance.piece_count)),
        }
    })?
}

pub(super) fn detail_contains_support_cell_inventory(
    detail: &str,
    piece_count: usize,
    cells_per_sector: usize,
) -> bool {
    detail.contains(&format!("support_cell_primitive_patch_count={piece_count}"))
        && detail.contains("support_cell_primitive_selection_policy=")
        && match cells_per_sector {
            4 => support_cell_grid_detail_matches(detail, &[(1, 4), (4, 1), (2, 2)]),
            8 => support_cell_grid_detail_matches(detail, &[(2, 4), (4, 2)]),
            _ => false,
        }
}

pub(super) fn support_cell_dispatch_details(
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
) -> Vec<String> {
    diagnostics
        .iter()
        .filter_map(|diagnostic| {
            let popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic
            else {
                return None;
            };
            (cap.code == "kp20_polar_macro.dispatch.fact" && cap.stage.contains("support_cell"))
                .then(|| cap.detail_code.clone())
                .flatten()
        })
        .collect()
}

fn support_cell_grid_detail_matches(detail: &str, grids: &[(usize, usize)]) -> bool {
    grids.iter().any(|&(cells_u, cells_v)| {
        detail.contains(&format!("support_cell_primitive_cells_u={cells_u}"))
            && detail.contains(&format!("support_cell_primitive_cells_v={cells_v}"))
    })
}
