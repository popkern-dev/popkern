//! Candidate application helpers for KP20 multipatch trials.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_multipatch_inner::{
    Kp20MultiPatchInnerBuildReport, Kp20MultiPatchInnerPatch,
};
use super::kp20_polar_multipatch_inner_shoulder::{Kp20OriginLiftMode, exact_worst_origin_surface};
use super::kp20_polar_multipatch_interface::{
    Kp20MultiPatchInterfaceBuildReport, Kp20MultiPatchInterfaceCurve,
};
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;
use super::kp20_polar_multipatch_strip_optimize::{
    Kp20StripInteriorOptimizerMode, optimized_strip_surface,
};
use super::kp20_polar_support_surface_split::{
    propagate_support_pole_interfaces, support_side_index,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kp20InnerConstructionMode {
    SupportPropagation,
    SideAlignedInnerPatch,
    ExactWorstInnerCorePositive,
    ExactWorstInnerCoreNegative,
    ExactWorstInnerShoulderPositive,
    ExactWorstInnerShoulderNegative,
    ExactWorstPropagatedShoulderPositive,
    ExactWorstPropagatedShoulderNegative,
}

impl Kp20InnerConstructionMode {
    pub(super) const fn all() -> [Self; 2] {
        [Self::SupportPropagation, Self::SideAlignedInnerPatch]
    }

    pub(super) const fn exact_guided_all() -> [Self; 8] {
        [
            Self::SupportPropagation,
            Self::SideAlignedInnerPatch,
            Self::ExactWorstInnerCorePositive,
            Self::ExactWorstInnerCoreNegative,
            Self::ExactWorstInnerShoulderPositive,
            Self::ExactWorstInnerShoulderNegative,
            Self::ExactWorstPropagatedShoulderPositive,
            Self::ExactWorstPropagatedShoulderNegative,
        ]
    }

    pub(super) const fn label(self) -> &'static str {
        match self {
            Self::SupportPropagation => "support_propagation",
            Self::SideAlignedInnerPatch => "side_aligned_inner_patch",
            Self::ExactWorstInnerCorePositive => "exact_worst_inner_core_positive",
            Self::ExactWorstInnerCoreNegative => "exact_worst_inner_core_negative",
            Self::ExactWorstInnerShoulderPositive => "exact_worst_inner_shoulder_positive",
            Self::ExactWorstInnerShoulderNegative => "exact_worst_inner_shoulder_negative",
            Self::ExactWorstPropagatedShoulderPositive => {
                "exact_worst_propagated_shoulder_positive"
            }
            Self::ExactWorstPropagatedShoulderNegative => {
                "exact_worst_propagated_shoulder_negative"
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct Kp20AppliedMultiPatchCandidate {
    pub(super) candidate: Vec<Kp20MacroPatch>,
    pub(super) changed_patches: Vec<[bool; 4]>,
    pub(super) changed_strip_piece_count: usize,
    pub(super) changed_inner_piece_count: usize,
    pub(super) changed_piece_count: usize,
}

pub(super) fn apply_multipatch_candidate(
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    strips: &Kp20MultiPatchStripBuildReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    strip_mode: Kp20StripInteriorOptimizerMode,
    inner_mode: Kp20InnerConstructionMode,
    interior_scale: f64,
) -> Kp20AppliedMultiPatchCandidate {
    let mut candidate = accepted_patches.to_vec();
    let mut changed_patches = vec![[false; 4]; candidate.len()];
    let strip_apply = apply_strip_patches(
        accepted_patches,
        &mut candidate,
        &mut changed_patches,
        accepted_summary,
        strips,
        strip_mode,
        interior_scale,
    );
    let changed_inner_piece_count = match inner_mode {
        Kp20InnerConstructionMode::SupportPropagation => propagate_strip_owned_inner_pieces(
            &mut candidate,
            &mut changed_patches,
            strip_apply.support_changed,
        ),
        Kp20InnerConstructionMode::SideAlignedInnerPatch => apply_side_aligned_inner_patches(
            &mut candidate,
            &mut changed_patches,
            interfaces,
            inner,
        ),
        Kp20InnerConstructionMode::ExactWorstInnerCorePositive
        | Kp20InnerConstructionMode::ExactWorstInnerCoreNegative
        | Kp20InnerConstructionMode::ExactWorstInnerShoulderPositive
        | Kp20InnerConstructionMode::ExactWorstInnerShoulderNegative => {
            apply_exact_worst_inner_core(
                accepted_patches,
                &mut candidate,
                &mut changed_patches,
                accepted_summary,
                interfaces,
                inner,
                inner_mode,
                interior_scale,
            )
        }
        Kp20InnerConstructionMode::ExactWorstPropagatedShoulderPositive
        | Kp20InnerConstructionMode::ExactWorstPropagatedShoulderNegative => {
            apply_exact_worst_propagated_shoulder(
                &mut candidate,
                &mut changed_patches,
                accepted_summary,
                strip_apply.support_changed,
                inner_mode,
                interior_scale,
            )
        }
    };
    let changed_piece_count = changed_piece_count(&changed_patches);
    Kp20AppliedMultiPatchCandidate {
        candidate,
        changed_patches,
        changed_strip_piece_count: strip_apply.changed_strip_piece_count,
        changed_inner_piece_count,
        changed_piece_count,
    }
}

struct Kp20StripApplyResult {
    changed_strip_piece_count: usize,
    support_changed: Vec<[bool; 2]>,
}

fn apply_strip_patches(
    accepted_patches: &[Kp20MacroPatch],
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    strips: &Kp20MultiPatchStripBuildReport,
    mode: Kp20StripInteriorOptimizerMode,
    interior_scale: f64,
) -> Kp20StripApplyResult {
    let mut changed_count = 0_usize;
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    for strip in &strips.strip_surfaces {
        if strip.sector_index >= macro_patches.len() {
            continue;
        }
        *macro_patches[strip.sector_index].piece_mut(strip.strip_piece) = optimized_strip_surface(
            accepted_patches,
            accepted_summary,
            strip,
            mode,
            interior_scale,
        );
        if !changed_patches[strip.sector_index][strip.strip_piece.index()] {
            changed_count += 1;
        }
        changed_patches[strip.sector_index][strip.strip_piece.index()] = true;
        support_changed[strip.sector_index][support_side_index(strip.strip_piece)] = true;
    }
    Kp20StripApplyResult {
        changed_strip_piece_count: changed_count,
        support_changed,
    }
}

fn propagate_strip_owned_inner_pieces(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    support_changed: Vec<[bool; 2]>,
) -> usize {
    let before = changed_patches.to_vec();
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    changed_patches
        .iter()
        .zip(before)
        .map(|(after, before)| {
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight]
                .into_iter()
                .filter(|piece| after[piece.index()] && !before[piece.index()])
                .count()
        })
        .sum()
}

fn apply_exact_worst_propagated_shoulder(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    support_changed: Vec<[bool; 2]>,
    mode: Kp20InnerConstructionMode,
    interior_scale: f64,
) -> usize {
    let Some(site) = accepted_summary.worst_fold_site() else {
        return 0;
    };
    if !matches!(
        site.piece,
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight
    ) {
        return 0;
    }
    let before = changed_patches.to_vec();
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    if site.sector >= macro_patches.len() {
        return changed_inner_count(changed_patches, &before);
    }
    let base = *macro_patches[site.sector].piece(site.piece);
    *macro_patches[site.sector].piece_mut(site.piece) =
        exact_worst_origin_surface(&base, base, site, origin_lift_mode(mode), interior_scale);
    changed_patches[site.sector][site.piece.index()] = true;
    changed_inner_count(changed_patches, &before)
}

fn apply_exact_worst_inner_core(
    accepted_patches: &[Kp20MacroPatch],
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    mode: Kp20InnerConstructionMode,
    interior_scale: f64,
) -> usize {
    let Some(site) = accepted_summary.worst_fold_site() else {
        return 0;
    };
    if !matches!(
        site.piece,
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight
    ) {
        return 0;
    }
    let before = changed_patches.to_vec();
    apply_side_aligned_inner_patches(macro_patches, changed_patches, interfaces, inner);
    if site.sector >= macro_patches.len() {
        return changed_inner_count(changed_patches, &before);
    }
    let Some(accepted_piece) = accepted_patches
        .get(site.sector)
        .map(|patch| *patch.piece(site.piece))
    else {
        return changed_inner_count(changed_patches, &before);
    };
    let lifted = exact_worst_origin_surface(
        &accepted_piece,
        *macro_patches[site.sector].piece(site.piece),
        site,
        origin_lift_mode(mode),
        interior_scale,
    );
    *macro_patches[site.sector].piece_mut(site.piece) = lifted;
    changed_patches[site.sector][site.piece.index()] = true;
    changed_inner_count(changed_patches, &before)
}

fn apply_side_aligned_inner_patches(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
) -> usize {
    let mut changed_count = 0_usize;
    for patch in &inner.inner_patches {
        if patch.sector_index >= macro_patches.len() {
            continue;
        }
        for (boundary_index, interface_name) in patch.boundary_interface_names.iter().enumerate() {
            let Some(interface) = interface_by_name(interfaces, interface_name) else {
                continue;
            };
            let origin_piece = origin_piece_for_strip(interface.strip_piece);
            *macro_patches[patch.sector_index].piece_mut(origin_piece) =
                side_aligned_inner_surface(patch, boundary_index);
            if !changed_patches[patch.sector_index][origin_piece.index()] {
                changed_count += 1;
            }
            changed_patches[patch.sector_index][origin_piece.index()] = true;
        }
    }
    changed_count
}

fn changed_inner_count(changed_patches: &[[bool; 4]], before: &[[bool; 4]]) -> usize {
    changed_patches
        .iter()
        .zip(before)
        .map(|(after, before)| {
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight]
                .into_iter()
                .filter(|piece| after[piece.index()] && !before[piece.index()])
                .count()
        })
        .sum()
}

fn side_aligned_inner_surface(
    patch: &Kp20MultiPatchInnerPatch,
    boundary_index: usize,
) -> Kp20Bi6Patch {
    match boundary_index {
        0 => patch.surface,
        1 => flip_u(patch.surface),
        2 => transpose(patch.surface),
        _ => flip_u(transpose(patch.surface)),
    }
}

fn interface_by_name<'a>(
    interfaces: &'a Kp20MultiPatchInterfaceBuildReport,
    name: &str,
) -> Option<&'a Kp20MultiPatchInterfaceCurve> {
    interfaces
        .interface_curves
        .iter()
        .find(|interface| interface.interface_name == name)
}

fn origin_piece_for_strip(strip_piece: Kp20PieceId) -> Kp20PieceId {
    match strip_piece {
        Kp20PieceId::NearMidpointLeft => Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearMidpointRight => Kp20PieceId::NearOriginRight,
        _ => unreachable!("multipatch interfaces are support-strip owned"),
    }
}

fn changed_piece_count(changed_patches: &[[bool; 4]]) -> usize {
    changed_patches
        .iter()
        .flat_map(|pieces| pieces.iter())
        .filter(|changed| **changed)
        .count()
}

fn origin_lift_mode(mode: Kp20InnerConstructionMode) -> Kp20OriginLiftMode {
    match mode {
        Kp20InnerConstructionMode::ExactWorstInnerCoreNegative => Kp20OriginLiftMode::CoreNegative,
        Kp20InnerConstructionMode::ExactWorstInnerShoulderPositive => {
            Kp20OriginLiftMode::ShoulderPositive
        }
        Kp20InnerConstructionMode::ExactWorstInnerShoulderNegative => {
            Kp20OriginLiftMode::ShoulderNegative
        }
        Kp20InnerConstructionMode::ExactWorstPropagatedShoulderPositive => {
            Kp20OriginLiftMode::ShoulderPositive
        }
        Kp20InnerConstructionMode::ExactWorstPropagatedShoulderNegative => {
            Kp20OriginLiftMode::ShoulderNegative
        }
        _ => Kp20OriginLiftMode::CorePositive,
    }
}

fn flip_u(surface: Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut flipped = Kp20Bi6Patch::zero();
    for u in 0..7 {
        flipped.cps[u] = surface.cps[6 - u];
    }
    flipped
}

fn transpose(surface: Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut transposed = Kp20Bi6Patch::zero();
    for u in 0..7 {
        for v in 0..7 {
            transposed.cps[u][v] = surface.cps[v][u];
        }
    }
    transposed
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inner::build_multipatch_inner_patches;
    use crate::construct::kp20_polar_multipatch_interface::build_multipatch_interfaces;
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
    use crate::construct::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;

    #[test]
    fn support_propagation_marks_matching_near_origin_piece() {
        let inventory = multipatch_inventory_from_fold_summary(&single_left_fold_summary());
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
        let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch]);
        let interfaces = build_multipatch_interfaces(&strips);
        let inner = build_multipatch_inner_patches(&inventory, &interfaces);

        let applied = apply_multipatch_candidate(
            &[Kp20MacroPatch::zero()],
            &empty_fold_summary(),
            &strips,
            &interfaces,
            &inner,
            Kp20StripInteriorOptimizerMode::FairBlend,
            Kp20InnerConstructionMode::SupportPropagation,
            1.0,
        );

        assert_eq!(applied.changed_strip_piece_count, 1);
        assert_eq!(applied.changed_inner_piece_count, 1);
        assert!(applied.changed_patches[0][Kp20PieceId::NearMidpointLeft.index()]);
        assert!(applied.changed_patches[0][Kp20PieceId::NearOriginLeft.index()]);
        assert!(!applied.changed_patches[0][Kp20PieceId::NearOriginRight.index()]);
        assert_eq!(applied.changed_piece_count, 2);
        assert_eq!(
            applied.candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[0][3],
            Point3::new(1.0, 2.0, 3.0)
        );
    }

    #[test]
    fn side_aligned_inner_patch_targets_matching_origin_side() {
        let inventory = multipatch_inventory_from_fold_summary(&single_left_fold_summary());
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][3] =
            Point3::new(0.25, 0.5, 0.75);
        let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch]);
        let interfaces = build_multipatch_interfaces(&strips);
        let inner = build_multipatch_inner_patches(&inventory, &interfaces);

        let applied = apply_multipatch_candidate(
            &[Kp20MacroPatch::zero()],
            &empty_fold_summary(),
            &strips,
            &interfaces,
            &inner,
            Kp20StripInteriorOptimizerMode::FairBlend,
            Kp20InnerConstructionMode::SideAlignedInnerPatch,
            1.0,
        );

        assert_eq!(applied.changed_inner_piece_count, 1);
        assert!(applied.changed_patches[0][Kp20PieceId::NearOriginLeft.index()]);
        assert!(!applied.changed_patches[0][Kp20PieceId::NearOriginRight.index()]);
        assert_eq!(
            applied.candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[0][3],
            Point3::new(1.0, 2.0, 3.0)
        );
        assert_ne!(
            applied.candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[1][3],
            applied.candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[0][3]
        );
    }

    fn single_left_fold_summary() -> Kp20SamplerFoldSummary {
        Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.5,
                v: 0.5,
                min_jacobian_ratio: -0.1,
                min_jacobian_normal_cos: -0.1,
            }],
        }
    }

    fn empty_fold_summary() -> Kp20SamplerFoldSummary {
        Kp20SamplerFoldSummary {
            fold_count: 0,
            min_jacobian_normal_cos: 1.0,
            piece_fold_counts: [0; 4],
            fold_details: Vec::new(),
            fold_sites: Vec::new(),
        }
    }
}
