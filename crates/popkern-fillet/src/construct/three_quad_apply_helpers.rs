//! Small helpers for three-quad cap application.

use crate::construct::promotion_ladder::PromotionResult;
use crate::types::{CornerSubPatchFamily, CornerSubPatchProvenance};

pub(super) fn sub_patch_provenance(
    macro_family: Option<CornerSubPatchFamily>,
    n_sides: usize,
    sub_patch_index: usize,
) -> Option<CornerSubPatchProvenance> {
    macro_family.map(|family| {
        CornerSubPatchProvenance::macro_cap(
            family,
            sub_patch_index / 4,
            n_sides,
            sub_patch_index % 4,
        )
    })
}

pub(super) fn macro_family_from_promotion(
    promotion_result: PromotionResult,
) -> Option<CornerSubPatchFamily> {
    match promotion_result {
        PromotionResult::Standard => None,
        PromotionResult::Macro4n(family) => Some(family),
    }
}

pub(super) fn decomp_self_intersects(
    decomp: &crate::blend::gregory_three_quad::NQuadDecomposition,
    tolerance: f64,
) -> bool {
    decomp
        .sub_patches
        .iter()
        .any(|patch| patch.self_intersects_with_tolerance(tolerance))
}
