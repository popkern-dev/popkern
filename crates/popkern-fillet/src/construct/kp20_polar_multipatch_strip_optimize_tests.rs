use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};

use super::*;
use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};

#[test]
fn strip_damping_preserves_tensor_boundary_and_applies_interface_g0_row() {
    let mut accepted = Kp20MacroPatch::zero();
    let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
    piece.cps[0][2] = Point3::new(0.0, 0.0, 0.0);
    piece.cps[1][2] = Point3::new(1.0, 0.0, 0.0);
    piece.cps[2][2] = Point3::new(2.0, 0.0, 0.0);
    piece.cps[3][2] = Point3::new(3.0, 0.0, 0.0);
    piece.cps[5][2] = Point3::new(5.0, 0.0, 0.0);
    piece.cps[6][2] = Point3::new(6.0, 0.0, 0.0);
    let mut strip_surface = *accepted.piece(Kp20PieceId::NearMidpointLeft);
    strip_surface.cps[2][2] = Point3::new(10.0, 0.0, 0.0);
    strip_surface.cps[3][2] = Point3::new(11.0, 0.0, 0.0);
    strip_surface.cps[5][2] = Point3::new(9.0, 0.0, 0.0);
    strip_surface.cps[6][2] = Point3::new(12.0, 0.0, 0.0);
    let strip = test_strip(strip_surface);

    let scaled = optimized_strip_surface(
        &[accepted],
        &empty_fold_summary(),
        &strip,
        Kp20StripInteriorOptimizerMode::FairBlend,
        0.25,
    );

    assert_eq!(scaled.cps[0][2], Point3::new(0.0, 0.0, 0.0));
    assert_eq!(scaled.cps[1][2], Point3::new(1.0, 0.0, 0.0));
    assert_eq!(scaled.cps[2][2], Point3::new(2.0, 0.0, 0.0));
    assert_eq!(scaled.cps[3][2], Point3::new(5.0, 0.0, 0.0));
    assert_eq!(scaled.cps[5][2], Point3::new(6.0, 0.0, 0.0));
    assert_eq!(scaled.cps[6][2], Point3::new(12.0, 0.0, 0.0));
}

#[test]
fn fold_normal_mode_moves_unconstrained_rows_and_keeps_interface_g0_constructed() {
    let mut accepted = Kp20MacroPatch::zero();
    let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
    for row in 0..7 {
        for column in 0..7 {
            piece.cps[row][column] = Point3::new(row as f64, column as f64, 0.0);
        }
    }
    let strip = test_strip(*accepted.piece(Kp20PieceId::NearMidpointLeft));

    let shifted = optimized_strip_surface(
        &[accepted],
        &single_fold_summary(),
        &strip,
        Kp20StripInteriorOptimizerMode::FoldNormalPositive,
        1.0,
    );

    assert_eq!(shifted.cps[2][3], Point3::new(2.0, 3.0, 0.0));
    assert_eq!(shifted.cps[1][3], Point3::new(1.0, 3.0, 0.0));
    assert!(shifted.cps[3][3].z > 0.0);
    assert!(shifted.cps[5][3].z > 0.0);
    assert_eq!(shifted.cps[6][3], Point3::new(6.0, 3.0, 0.0));
}

#[test]
fn shoulder_row_target_is_bounded_to_constructed_row() {
    let mut accepted = Kp20Bi6Patch::zero();
    let mut constructed = Kp20Bi6Patch::zero();
    accepted.cps[5][4] = Point3::new(1.0, 0.0, 0.0);
    constructed.cps[5][4] = Point3::new(5.0, 0.0, 0.0);

    assert_eq!(
        shoulder_row_target(
            &accepted,
            &constructed,
            &empty_fold_summary(),
            &test_strip(constructed),
            Kp20StripInteriorOptimizerMode::FairBlend,
            4,
            1.0,
            0.25
        ),
        Point3::new(2.0, 0.0, 0.0)
    );
    assert_eq!(
        shoulder_row_target(
            &accepted,
            &constructed,
            &empty_fold_summary(),
            &test_strip(constructed),
            Kp20StripInteriorOptimizerMode::FairBlend,
            4,
            1.0,
            3.0
        ),
        Point3::new(5.0, 0.0, 0.0)
    );
}

#[test]
fn worst_edge_column_displacement_is_local_to_worst_support_edge() {
    let accepted = flat_strip_surface();
    let strip = test_strip(accepted);
    let summary = support_edge_fold_summary();

    let local = worst_edge_column_displacement(
        &accepted,
        &summary,
        &strip,
        Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive,
        4,
        1,
    );
    let distant = worst_edge_column_displacement(
        &accepted,
        &summary,
        &strip,
        Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive,
        4,
        4,
    );

    assert!(local.length() > 0.0);
    assert_eq!(distant, Vec3::zero());
}

#[test]
fn edge_parameter_reauthoring_expands_or_contracts_around_worst_site() {
    let expanded = edge_reauthor_row_parameters(0.52, true, 1.0);
    let contracted = edge_reauthor_row_parameters(0.52, false, 1.0);

    assert!(expanded[0] < 0.25);
    assert!(expanded[1] < 0.5);
    assert!(expanded[2] > 0.75);
    assert!(contracted[0] > 0.25);
    assert!(contracted[1] > 0.5);
    assert!(contracted[2] < 0.75);
    assert!(expanded[0] < expanded[1] && expanded[1] < expanded[2]);
    assert!(contracted[0] < contracted[1] && contracted[1] < contracted[2]);
}

#[test]
fn edge_parameter_mode_preserves_interface_g0_and_moves_target_column() {
    let accepted = flat_strip_surface();
    let strip = test_strip(accepted);
    let shifted = optimized_strip_surface(
        &[macro_patch_with_strip(accepted)],
        &support_edge_fold_summary(),
        &strip,
        Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand,
        1.0,
    );

    assert_eq!(shifted.cps[6][1], Point3::new(6.0, 1.0, 0.0));
    assert_ne!(shifted.cps[4][1], Point3::new(4.0, 1.0, 0.0));
    assert_eq!(shifted.cps[4][4], Point3::new(4.0, 4.0, 0.0));
}

#[test]
fn fold_site_weight_is_local_and_severity_bounded() {
    let center_weight = fold_site_weight(0.5, 0.5, 0.5, 0.5, -1.0);
    let adjacent_weight = fold_site_weight(0.5, 0.5, 0.5, 0.65, -1.0);
    let distant_weight = fold_site_weight(0.5, 0.5, 0.5, 1.0, -1.0);
    let mild_weight = fold_site_weight(0.5, 0.5, 0.5, 0.5, 0.75);

    assert_eq!(center_weight, 1.0);
    assert!(adjacent_weight > 0.0);
    assert!(adjacent_weight < center_weight);
    assert_eq!(distant_weight, 0.0);
    assert_eq!(mild_weight, 0.25);
}

#[test]
fn fold_site_displacement_is_capped() {
    let capped = capped_displacement(Vec3::new(1.0, 0.0, 0.0), 0.08);
    let uncapped = capped_displacement(Vec3::new(0.01, 0.02, 0.0), 0.08);

    assert!((capped.length() - 0.08).abs() < 1.0e-12);
    assert_eq!(uncapped, Vec3::new(0.01, 0.02, 0.0));
}

fn flat_strip_surface() -> Kp20Bi6Patch {
    let mut accepted = Kp20Bi6Patch::zero();
    for row in 0..7 {
        for column in 0..7 {
            accepted.cps[row][column] = Point3::new(row as f64, column as f64, 0.0);
        }
    }
    accepted
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

fn single_fold_summary() -> Kp20SamplerFoldSummary {
    Kp20SamplerFoldSummary {
        fold_count: 1,
        min_jacobian_normal_cos: -0.5,
        piece_fold_counts: [1, 0, 0, 0],
        fold_details: Vec::new(),
        fold_sites: vec![Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearOriginLeft,
            u: 0.5,
            v: 0.5,
            min_jacobian_ratio: -0.5,
            min_jacobian_normal_cos: -0.5,
        }],
    }
}

fn support_edge_fold_summary() -> Kp20SamplerFoldSummary {
    Kp20SamplerFoldSummary {
        fold_count: 1,
        min_jacobian_normal_cos: -1.0,
        piece_fold_counts: [0, 0, 0, 1],
        fold_details: Vec::new(),
        fold_sites: vec![Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearMidpointLeft,
            u: 0.68,
            v: 0.06,
            min_jacobian_ratio: -1.0,
            min_jacobian_normal_cos: -1.0,
        }],
    }
}

fn test_strip(surface: Kp20Bi6Patch) -> Kp20MultiPatchStripSurface {
    Kp20MultiPatchStripSurface {
        patch_name: "strip".to_string(),
        sector_index: 0,
        strip_piece: Kp20PieceId::NearMidpointLeft,
        surface,
    }
}

fn macro_patch_with_strip(surface: Kp20Bi6Patch) -> Kp20MacroPatch {
    let mut patch = Kp20MacroPatch::zero();
    *patch.piece_mut(Kp20PieceId::NearMidpointLeft) = surface;
    patch
}
