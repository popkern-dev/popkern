use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::{
    Kp20InnerConstructionMode, apply_side_aligned_inner_patches, changed_inner_count,
    origin_piece_for_strip,
};
use crate::construct::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use crate::construct::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use crate::construct::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport;
use crate::construct::kp20_polar_support_surface_split::propagate_support_pole_interfaces;

pub(super) fn apply_exact_worst_coupled_support_inner(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    support_changed: Vec<[bool; 2]>,
    mode: Kp20InnerConstructionMode,
    interior_scale: f64,
) -> usize {
    let Some(site) = accepted_summary.worst_fold_site() else {
        return 0;
    };
    let support_piece = support_piece_for_site(site.piece);
    if site.sector >= macro_patches.len() {
        return 0;
    }
    let before = changed_patches.to_vec();
    apply_side_aligned_inner_patches(macro_patches, changed_patches, interfaces, inner);
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    reauthor_coupled_inner_shoulder(
        macro_patches,
        changed_patches,
        site.sector,
        support_piece,
        site.v,
        mode == Kp20InnerConstructionMode::ExactWorstCoupledSupportInnerExpand,
        interior_scale,
    );
    changed_inner_count(changed_patches, &before)
}

fn support_piece_for_site(piece: Kp20PieceId) -> Kp20PieceId {
    match piece {
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearMidpointLeft => {
            Kp20PieceId::NearMidpointLeft
        }
        Kp20PieceId::NearOriginRight | Kp20PieceId::NearMidpointRight => {
            Kp20PieceId::NearMidpointRight
        }
    }
}

fn reauthor_coupled_inner_shoulder(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    sector: usize,
    support_piece: Kp20PieceId,
    site_v: f64,
    expand: bool,
    interior_scale: f64,
) {
    let origin_piece = origin_piece_for_strip(support_piece);
    let support = *macro_patches[sector].piece(support_piece);
    let origin = macro_patches[sector].piece_mut(origin_piece);
    let target_column = if site_v <= 0.5 { 1_usize } else { 5_usize };
    for column in 0..7 {
        let weight = coupled_shoulder_column_weight(column, target_column);
        if weight <= 0.0 {
            continue;
        }
        let boundary = origin.cps[0][column];
        let support_adjacent = support.cps[5][column];
        let tangent_target = boundary + (boundary - support_adjacent);
        let current = origin.cps[1][column];
        let directional_target = if expand {
            current + (tangent_target - current) * weight
        } else {
            current + (current - tangent_target) * (0.5 * weight)
        };
        origin.cps[1][column] =
            current + (directional_target - current) * interior_scale.clamp(0.0, 1.0);
    }
    changed_patches[sector][origin_piece.index()] = true;
}

fn coupled_shoulder_column_weight(column: usize, target_column: usize) -> f64 {
    match column.abs_diff(target_column) {
        0 => 1.0,
        1 => 0.45,
        _ => 0.0,
    }
}
