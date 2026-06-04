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
    let summary = single_left_fold_summary();
    let inventory = multipatch_inventory_from_fold_summary(&summary);
    let mut macro_patch = Kp20MacroPatch::zero();
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
    let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);
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
    let summary = single_left_fold_summary();
    let inventory = multipatch_inventory_from_fold_summary(&summary);
    let mut macro_patch = Kp20MacroPatch::zero();
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][3] = Point3::new(0.25, 0.5, 0.75);
    let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);
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

#[test]
fn side_aligned_propagated_inner_patch_preserves_constructed_interface_g0() {
    let summary = single_left_fold_summary();
    let inventory = multipatch_inventory_from_fold_summary(&summary);
    let mut macro_patch = Kp20MacroPatch::zero();
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][3] = Point3::new(0.25, 0.5, 0.75);
    let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);
    let interfaces = build_multipatch_interfaces(&strips);
    let inner = build_multipatch_inner_patches(&inventory, &interfaces);

    let applied = apply_multipatch_candidate(
        &[Kp20MacroPatch::zero()],
        &empty_fold_summary(),
        &strips,
        &interfaces,
        &inner,
        Kp20StripInteriorOptimizerMode::FairBlend,
        Kp20InnerConstructionMode::SideAlignedPropagatedInnerPatch,
        1.0,
    );

    let origin = applied.candidate[0].piece(Kp20PieceId::NearOriginLeft);
    assert_eq!(applied.changed_inner_piece_count, 1);
    assert_eq!(origin.cps[0][3], Point3::new(1.0, 2.0, 3.0));
    assert_ne!(origin.cps[1][3], origin.cps[0][3]);
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
