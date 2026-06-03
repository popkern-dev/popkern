use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::*;
use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;

#[test]
fn strip_builder_reauthors_unconstrained_interior_rows() {
    let summary = Kp20SamplerFoldSummary {
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
    };
    let inventory = multipatch_inventory_from_fold_summary(&summary);
    let mut macro_patch = Kp20MacroPatch::zero();
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[2][2] = Point3::new(1.0, 2.0, 3.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[1][2] = Point3::new(0.0, 0.0, 0.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][2] = Point3::new(0.0, 0.0, 4.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][2] = Point3::new(0.0, 0.0, 5.0);

    let report = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);

    assert_eq!(report.attempted_strip_patch_count, 1);
    assert_eq!(report.built_strip_patch_count, 1);
    assert_eq!(report.reauthored_strip_patch_count, 1);
    assert!(report.max_reauthored_interior_displacement > 0.0);
    assert_eq!(report.finite_control_point_count, 49);
    assert_eq!(report.invalid_control_point_count, 0);
    assert!(report.max_row_spacing_ratio.is_finite());
    assert!(report.max_row_curvature_ratio.is_finite());
    assert_eq!(report.decline_reason, None);
    assert_eq!(report.strip_surfaces[0].sector_index, 0);
    assert_eq!(
        report.strip_surfaces[0].strip_piece,
        Kp20PieceId::NearMidpointLeft
    );
    assert!(report.strip_surfaces[0].surface.cps[2][2].z > 1.0);
    assert!(report.strip_surfaces[0].surface.cps[2][2].z < 2.0);
    assert_eq!(
        report.strip_surfaces[0].surface.cps[5][2],
        Point3::new(0.0, 0.0, 4.0)
    );
    assert_eq!(
        report.strip_surfaces[0].surface.cps[6][2],
        Point3::new(0.0, 0.0, 5.0)
    );
}

#[test]
fn strip_parameter_map_is_monotone_and_fold_site_aware() {
    let summary = Kp20SamplerFoldSummary {
        fold_count: 1,
        min_jacobian_normal_cos: -0.9,
        piece_fold_counts: [1, 0, 0, 0],
        fold_details: Vec::new(),
        fold_sites: vec![Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearOriginLeft,
            u: 0.85,
            v: 0.5,
            min_jacobian_ratio: -0.1,
            min_jacobian_normal_cos: -0.9,
        }],
    };

    let map = strip_column_parameter_map(
        &summary,
        0,
        Kp20PieceId::NearMidpointLeft,
        3,
        balanced_column_geometry(),
    );

    assert!(map.row_parameters[0] > UNIFORM_ROW_PARAMETERS[0]);
    assert!(map.min_step() >= MIN_PARAMETER_STEP - 1.0e-12);
    assert!(map.max_skew() > 0.0);
}

#[test]
fn strip_parameter_map_uses_uniform_rows_without_owned_fold_sites() {
    let summary = empty_fold_summary();

    let map = strip_column_parameter_map(
        &summary,
        0,
        Kp20PieceId::NearMidpointLeft,
        3,
        balanced_column_geometry(),
    );

    assert_eq!(map.row_parameters, UNIFORM_ROW_PARAMETERS);
    assert_eq!(map.max_skew(), 0.0);
}

#[test]
fn strip_parameter_map_applies_bounded_tangent_balance_bias() {
    let geometry = Kp20StripColumnGeometry {
        support_boundary: Point3::new(0.0, 0.0, 0.0),
        support_adjacent: Point3::new(2.0, 0.0, 0.0),
        interface_adjacent: Point3::new(4.5, 0.0, 0.0),
        interface_boundary: Point3::new(5.0, 0.0, 0.0),
    };

    let map = strip_column_parameter_map(
        &empty_fold_summary(),
        0,
        Kp20PieceId::NearMidpointLeft,
        3,
        geometry,
    );

    assert!(map.row_parameters[0] > UNIFORM_ROW_PARAMETERS[0]);
    assert!(map.min_step() >= MIN_PARAMETER_STEP - 1.0e-12);
    assert!(geometry.tangent_step_imbalance() > 0.0);
    assert!(geometry.tangent_alignment() > 0.0);
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

fn balanced_column_geometry() -> Kp20StripColumnGeometry {
    Kp20StripColumnGeometry {
        support_boundary: Point3::new(0.0, 0.0, 0.0),
        support_adjacent: Point3::new(1.0, 0.0, 0.0),
        interface_adjacent: Point3::new(4.0, 0.0, 0.0),
        interface_boundary: Point3::new(5.0, 0.0, 0.0),
    }
}
