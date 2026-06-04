use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::*;
use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
use crate::construct::kp20_polar_multipatch_interface::build_multipatch_interfaces;
use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
use crate::construct::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;

#[test]
fn inner_builder_copies_named_interface_boundary_controls() {
    let summary = single_left_fold_summary();
    let inventory = multipatch_inventory_from_fold_summary(&summary);
    let mut macro_patch = Kp20MacroPatch::zero();
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
    macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][3] = Point3::new(0.25, 0.5, 0.75);
    let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);
    let interfaces = build_multipatch_interfaces(&strips);

    let report = build_multipatch_inner_patches(&inventory, &interfaces);

    assert_eq!(report.attempted_inner_patch_count, 1);
    assert_eq!(report.built_inner_patch_count, 1);
    assert_eq!(report.finite_control_point_count, 49);
    assert_eq!(report.invalid_control_point_count, 0);
    assert_eq!(report.boundary_constraint_count, 7);
    assert_eq!(report.adjacent_row_constraint_count, 7);
    assert_eq!(report.interior_constraint_count, 9);
    assert_eq!(report.exact_boundary_owned_endpoint_count, 0);
    assert_eq!(report.max_boundary_g0_residual, 0.0);
    assert_eq!(report.max_adjacent_row_residual, 0.0);
    assert_eq!(report.max_interior_residual, 0.0);
    assert!(report.min_adjacent_row_step.is_finite());
    assert!(report.max_coupled_core_displacement > 0.0);
    assert_eq!(report.decline_reason, None);
    assert_eq!(report.inner_patches[0].sector_index, 0);
    assert_eq!(
        report.inner_patches[0].patch_name,
        "kp20_polar_macro.sector_0.inner_pole_patch"
    );
    assert_eq!(
        report.inner_patches[0].boundary_interface_names,
        vec!["kp20_polar_macro.sector_0.left_support_interface"]
    );
    assert_eq!(
        report.inner_patches[0].surface.cps[0][3],
        Point3::new(1.0, 2.0, 3.0)
    );
    assert_eq!(
        report.inner_patches[0].surface.cps[1][3],
        Point3::new(1.75, 3.5, 5.25)
    );
    assert_ne!(
        report.inner_patches[0].surface.cps[3][3],
        Point3::new(0.0, 0.0, 0.0)
    );
}

#[test]
fn coupled_inner_core_target_uses_adjacent_frame_axes() {
    let mut surface = Kp20Bi6Patch::zero();
    surface.cps[3][1] = Point3::new(0.0, 2.0, 0.0);
    surface.cps[3][5] = Point3::new(6.0, 2.0, 0.0);
    surface.cps[1][3] = Point3::new(2.0, 0.0, 0.0);
    surface.cps[5][3] = Point3::new(2.0, 6.0, 0.0);

    let target = coupled_inner_core_target(&surface, 3, 3);

    assert_eq!(target, Point3::new(2.5, 2.5, 0.0));
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
