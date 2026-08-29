//! Support-lane trial patchwork for expanded KP20 topology search.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20PieceId};

use super::split_bi6_patch_u;
use super::{
    Kp20ExpandedSupportMode, Kp20ExpandedSupportPatchwork, Kp20FoldSite, push_support_lane,
    reauthored_support_lanes,
};

pub(super) fn expanded_support_lane_patchwork(
    macro_patches: &[popkern_core::geometry::gregory::Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    support_piece: Kp20PieceId,
    split_parameter: f64,
    reauthor: Option<(Kp20ExpandedSupportMode, f64)>,
) -> Kp20ExpandedSupportPatchwork {
    let patch_count = targets.len() * 2;
    let mut patches = Vec::with_capacity(patch_count);
    let mut provenance = Vec::with_capacity(patch_count);
    for target in targets {
        let base = *macro_patches[target.sector].piece(support_piece);
        let lanes = support_lanes(&base, *target, split_parameter, reauthor);
        let sector_count = macro_patches.len();
        push_support_lane(
            &mut patches,
            &mut provenance,
            lanes.0,
            *target,
            sector_count,
            patch_count,
            true,
        );
        push_support_lane(
            &mut patches,
            &mut provenance,
            lanes.1,
            *target,
            sector_count,
            patch_count,
            false,
        );
    }
    Kp20ExpandedSupportPatchwork {
        patches,
        provenance,
    }
}

fn support_lanes(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    split_parameter: f64,
    reauthor: Option<(Kp20ExpandedSupportMode, f64)>,
) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    match reauthor {
        Some((mode, scale)) => reauthored_support_lanes(base, target, split_parameter, mode, scale),
        None => split_bi6_patch_u(base, split_parameter),
    }
}
