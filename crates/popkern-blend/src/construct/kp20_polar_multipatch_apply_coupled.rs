use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::{
    Kp20InnerConstructionMode, apply_side_aligned_inner_patches, changed_inner_count,
    origin_piece_for_strip,
};
use crate::construct::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use crate::construct::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use crate::construct::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport;
use crate::construct::kp20_polar_support_surface_split::propagate_support_pole_interfaces;

#[allow(clippy::too_many_arguments)]
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
    if matches!(
        mode,
        Kp20InnerConstructionMode::ExactFoldSetSharedSupportInnerExpand
            | Kp20InnerConstructionMode::ExactFoldSetSharedSupportInnerContract
    ) {
        return apply_exact_fold_set_shared_support_inner(
            macro_patches,
            changed_patches,
            accepted_summary,
            interfaces,
            inner,
            support_changed,
            mode,
            interior_scale,
        );
    }
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
        CoupledShoulderMode::from_inner_mode(mode),
        interior_scale,
    );
    changed_inner_count(changed_patches, &before)
}

#[allow(clippy::too_many_arguments)]
fn apply_exact_fold_set_shared_support_inner(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    support_changed: Vec<[bool; 2]>,
    mode: Kp20InnerConstructionMode,
    interior_scale: f64,
) -> usize {
    let before = changed_patches.to_vec();
    apply_side_aligned_inner_patches(macro_patches, changed_patches, interfaces, inner);
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    let expand = matches!(
        mode,
        Kp20InnerConstructionMode::ExactFoldSetSharedSupportInnerExpand
    );
    for target in fold_set_targets(accepted_summary, macro_patches.len()) {
        reauthor_coupled_inner_shoulder(
            macro_patches,
            changed_patches,
            target.sector,
            target.support_piece,
            target.site_v,
            CoupledShoulderMode {
                expand,
                shared_frame: true,
            },
            interior_scale,
        );
    }
    changed_inner_count(changed_patches, &before)
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct FoldSetTarget {
    sector: usize,
    support_piece: Kp20PieceId,
    target_column: usize,
    site_v: f64,
}

fn fold_set_targets(summary: &Kp20SamplerFoldSummary, sector_count: usize) -> Vec<FoldSetTarget> {
    let mut targets = Vec::new();
    for site in &summary.fold_sites {
        if site.sector >= sector_count {
            continue;
        }
        let support_piece = support_piece_for_site(site.piece);
        let target_column = if site.v <= 0.5 { 1_usize } else { 5_usize };
        if targets.iter().any(|target: &FoldSetTarget| {
            target.sector == site.sector
                && target.support_piece == support_piece
                && target.target_column == target_column
        }) {
            continue;
        }
        targets.push(FoldSetTarget {
            sector: site.sector,
            support_piece,
            target_column,
            site_v: if target_column == 1 { 0.0 } else { 1.0 },
        });
    }
    targets
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
    mode: CoupledShoulderMode,
    interior_scale: f64,
) {
    let origin_piece = origin_piece_for_strip(support_piece);
    let target_column = if site_v <= 0.5 { 1_usize } else { 5_usize };
    if mode.shared_frame {
        reauthor_shared_support_inner_frame(
            macro_patches,
            sector,
            support_piece,
            origin_piece,
            target_column,
            mode.expand,
            interior_scale,
        );
    } else {
        let support = *macro_patches[sector].piece(support_piece);
        let origin = macro_patches[sector].piece_mut(origin_piece);
        reauthor_inner_shoulder_from_support(
            origin,
            &support,
            target_column,
            mode.expand,
            interior_scale,
        );
    }
    if mode.shared_frame {
        changed_patches[sector][support_piece.index()] = true;
    }
    changed_patches[sector][origin_piece.index()] = true;
}

#[derive(Debug, Clone, Copy)]
struct CoupledShoulderMode {
    expand: bool,
    shared_frame: bool,
}

impl CoupledShoulderMode {
    fn from_inner_mode(mode: Kp20InnerConstructionMode) -> Self {
        Self {
            expand: matches!(
                mode,
                Kp20InnerConstructionMode::ExactWorstCoupledSupportInnerExpand
                    | Kp20InnerConstructionMode::ExactWorstSharedSupportInnerExpand
            ),
            shared_frame: matches!(
                mode,
                Kp20InnerConstructionMode::ExactWorstSharedSupportInnerExpand
                    | Kp20InnerConstructionMode::ExactWorstSharedSupportInnerContract
            ),
        }
    }
}

fn reauthor_inner_shoulder_from_support(
    origin: &mut popkern_core::geometry::gregory::Kp20Bi6Patch,
    support: &popkern_core::geometry::gregory::Kp20Bi6Patch,
    target_column: usize,
    expand: bool,
    interior_scale: f64,
) {
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
}

fn reauthor_shared_support_inner_frame(
    macro_patches: &mut [Kp20MacroPatch],
    sector: usize,
    support_piece: Kp20PieceId,
    origin_piece: Kp20PieceId,
    target_column: usize,
    expand: bool,
    interior_scale: f64,
) {
    for column in 0..7 {
        let weight = coupled_shoulder_column_weight(column, target_column);
        if weight <= 0.0 {
            continue;
        }
        let support_current = macro_patches[sector].piece(support_piece).cps[5][column];
        let origin_current = macro_patches[sector].piece(origin_piece).cps[1][column];
        let interface = macro_patches[sector].piece(support_piece).cps[6][column];
        let half_step = (origin_current - support_current) * 0.5;
        let scale = if expand { 1.15 } else { 0.75 };
        let support_target = interface - half_step * scale;
        let origin_target = interface + half_step * scale;
        let factor = (weight * interior_scale).clamp(0.0, 1.0);
        macro_patches[sector].piece_mut(support_piece).cps[5][column] =
            support_current + (support_target - support_current) * factor;
        macro_patches[sector].piece_mut(origin_piece).cps[1][column] =
            origin_current + (origin_target - origin_current) * factor;
    }
}

fn coupled_shoulder_column_weight(column: usize, target_column: usize) -> f64 {
    match column.abs_diff(target_column) {
        0 => 1.0,
        1 => 0.45,
        _ => 0.0,
    }
}
