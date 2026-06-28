use super::super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::super::kp20_bi6_cc_ring::kp20_catmull_clark_characteristic_ring;
use super::super::kp20_bi6_cnet::kp20_tensor_border_from_regular_cnet_rows;
use super::super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauF3ScalarSolveInput, kp20_planar_tau_diagonal_symmetry_constraints,
    kp20_planar_tau_internal_c2_constraints, kp20_planar_tau_rotated_sector_g1_g2_constraints,
    kp20_planar_tau_tensor_border_constraints, solve_kp20_planar_tau_f3_with_scalar_constraints,
};
use super::super::kp20_bi6_tau_solve_assemble::{
    Kp20PlanarTauPaperRowFamily, Kp20PlanarTauPaperSolveInput,
    kp20_planar_tau_paper_family_residuals, solve_kp20_planar_tau_paper_sector_f3,
    solve_kp20_planar_tau_rotational_paper_sector_f3,
};
use super::*;

fn assert_close(a: Point3, b: Point3, tol: f64) {
    assert!(
        a.approx_eq(&b, tol),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

fn affine_grid(offset: f64) -> [[Point3; 4]; 4] {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = Point3::new(
                offset + i as f64,
                10.0 * offset + j as f64,
                i as f64 - j as f64,
            );
        }
    }
    grid
}

#[test]
fn stam_tile_dispatch_matches_equations_7_to_9() {
    let lower = kp20_stam_tile_for(0.75, 0.25).unwrap();
    assert_eq!(lower.patch, Kp20StamPatch::LowerRight);
    assert_eq!(lower.level, 1);
    assert_eq!((lower.u, lower.v), (0.5, 0.5));

    let center = kp20_stam_tile_for(0.375, 0.375).unwrap();
    assert_eq!(center.patch, Kp20StamPatch::CenterRight);
    assert_eq!(center.level, 2);
    assert_eq!((center.u, center.v), (0.5, 0.5));

    let upper = kp20_stam_tile_for(0.125, 0.375).unwrap();
    assert_eq!(upper.patch, Kp20StamPatch::UpperLeft);
    assert_eq!(upper.level, 2);
    assert_eq!((upper.u, upper.v), (0.5, 0.5));
}

#[test]
fn stam_characteristic_eval_scales_by_lambda_per_tile_level() {
    let grids = [affine_grid(1.0), affine_grid(2.0), affine_grid(3.0)];
    let lambda = 0.42;

    let level_one = kp20_stam_eval_characteristic(&grids, lambda, 0.75, 0.25).unwrap();
    let direct_one = eval_cubic_bspline_grid(&grids[0], 0.5, 0.5).unwrap();
    assert_close(level_one, direct_one, 1e-12);

    let level_two = kp20_stam_eval_characteristic(&grids, lambda, 0.375, 0.375).unwrap();
    let direct_two = eval_cubic_bspline_grid(&grids[1], 0.5, 0.5).unwrap();
    assert_close(level_two, scale_point(direct_two, lambda), 1e-12);
}

#[test]
fn cubic_bspline_eval_preserves_constant_grid() {
    let grid = [[Point3::new(2.0, -3.0, 5.0); 4]; 4];
    for (u, v) in [(0.0, 0.0), (0.25, 0.5), (1.0, 1.0)] {
        assert_close(
            eval_cubic_bspline_grid(&grid, u, v).unwrap(),
            Point3::new(2.0, -3.0, 5.0),
            1e-12,
        );
    }
}

#[test]
fn figure6_picking_grids_follow_stam_index_layout() {
    let grids = kp20_stam_figure6_picking_grids(5).unwrap();
    assert_eq!(grids[0][0], [17, 16, 12, 20]);
    assert_eq!(grids[0][3], [8, 7, 15, 23]);
    assert_eq!(grids[1][0], [26, 25, 24, 19]);
    assert_eq!(grids[2][0], [27, 26, 25, 24]);
    assert_eq!(grids[2][3], [2, 1, 6, 14]);
}

#[test]
fn figure6_picking_grids_reject_invalid_valence() {
    assert!(kp20_stam_figure6_picking_grids(2).is_err());
}

#[test]
fn pr08_stam_visual_projection_is_rotationally_covariant() {
    let n = 7;
    let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let theta = core::f64::consts::TAU / n as f64;
    let (s, c) = theta.sin_cos();

    for sector in 0..n {
        let next = (sector + 1) % n;
        for patch in 0..3 {
            for row in 0..4 {
                for col in 0..4 {
                    let current = projected.sectors[sector][patch][row][col];
                    let expected = Point3::new(
                        c * current.x - s * current.y,
                        s * current.x + c * current.y,
                        0.0,
                    );
                    assert_close(projected.sectors[next][patch][row][col], expected, 1e-10);
                }
            }
        }
    }
}

#[test]
fn pr08_middle_patch_normalization_localizes_visual_projection_mismatch() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let mut best = f64::INFINITY;
    let mut best_label = "";

    for (label, grid) in [
        (
            "lower_right",
            &projected.sectors[0][Kp20StamPatch::LowerRight.index()],
        ),
        (
            "center_right",
            &projected.sectors[0][Kp20StamPatch::CenterRight.index()],
        ),
        (
            "upper_left",
            &projected.sectors[0][Kp20StamPatch::UpperLeft.index()],
        ),
    ] {
        for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            let value = eval_cubic_bspline_grid(grid, u, v).unwrap();
            let distance = value.distance_to(&expected);
            if distance < best {
                best = distance;
                best_label = label;
            }
        }
    }

    assert_eq!(best_label, "upper_left");
    assert!(
        best > 0.25 && best < 0.27,
        "expected current visual projection to miss PR08 normalization by about 0.26, got {best}"
    );
}

#[test]
fn pr08_middle_patch_normalization_rejects_transformed_visual_projection() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let candidates = [
        (
            "lower-right-identity",
            Kp20StamPatch::LowerRight,
            Kp20StamGridTransform::Identity,
        ),
        (
            "center-right-identity",
            Kp20StamPatch::CenterRight,
            Kp20StamGridTransform::Identity,
        ),
        (
            "upper-left-identity",
            Kp20StamPatch::UpperLeft,
            Kp20StamGridTransform::Identity,
        ),
        (
            "upper-left-rotate180",
            Kp20StamPatch::UpperLeft,
            Kp20StamGridTransform::Rotate180,
        ),
    ];
    let mut best = f64::INFINITY;
    let mut best_label = "";

    for (label, patch, transform) in candidates {
        let grid =
            kp20_pr08_stam_projected_transformed_grid(&projected, 0, patch, transform).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            let value = eval_cubic_bspline_grid(&grid, u, v).unwrap();
            let distance = value.distance_to(&expected);
            if distance < best {
                best = distance;
                best_label = label;
            }
        }
    }

    assert_eq!(best_label, "upper-left-rotate180");
    assert!(
        best > 0.25 && best < 0.27,
        "expected square-transformed visual projection to miss PR08 normalization by about 0.26, got {best}"
    );
}

#[test]
fn pr08_stam_appendix_tensor_order_still_misses_middle_patch_normalization() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_appendix_tensor_projected_grids(&ring).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let mut best = f64::INFINITY;
    let mut best_label = "";

    for (label, grid) in [
        (
            "lower-right",
            &projected.sectors[0][Kp20StamPatch::LowerRight.index()],
        ),
        (
            "center-right",
            &projected.sectors[0][Kp20StamPatch::CenterRight.index()],
        ),
        (
            "upper-left",
            &projected.sectors[0][Kp20StamPatch::UpperLeft.index()],
        ),
    ] {
        for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            let value = eval_cubic_bspline_grid(grid, u, v).unwrap();
            let distance = value.distance_to(&expected);
            if distance < best {
                best = distance;
                best_label = label;
            }
        }
    }

    assert_eq!(best_label, "upper-left");
    assert!(
        best > 0.25 && best < 0.27,
        "expected Appendix-B tensor-order projection to preserve the known 0.26 normalization miss, got {best}"
    );
}

#[test]
fn pr08_stam_appendix_tensor_order_is_visual_rotate90_basis_order() {
    let ring = kp20_catmull_clark_characteristic_ring(7).unwrap();
    let visual = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let appendix = kp20_pr08_stam_appendix_tensor_projected_grids(&ring).unwrap();

    for sector in 0..ring.sectors.len() {
        for patch in [
            Kp20StamPatch::LowerRight,
            Kp20StamPatch::CenterRight,
            Kp20StamPatch::UpperLeft,
        ] {
            let transformed = kp20_pr08_stam_projected_transformed_grid(
                &visual,
                sector,
                patch,
                Kp20StamGridTransform::Rotate90,
            )
            .unwrap();
            for u in 0..4 {
                for v in 0..4 {
                    assert_close(
                        appendix.sectors[sector][patch.index()][u][v],
                        transformed[u][v],
                        1e-12,
                    );
                }
            }
        }
    }
}

#[test]
fn pr08_middle_normalized_projection_matches_printed_normalization() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_middle_normalized_projected_grids(&ring).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let grid = &projected.sectors[0][Kp20StamPatch::LowerRight.index()];
    assert_close(
        eval_cubic_bspline_grid(grid, 1.0, 1.0).unwrap(),
        expected,
        1.0e-12,
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source visual projection residual diagnostic"]
fn pr08_middle_normalized_projected_sides_are_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_middle_normalized_projected_grids(&ring).unwrap();
    let mut best = f64::INFINITY;
    let mut best_side = Kp20StamGridSide::Top;
    let mut best_reversed = false;
    for side in Kp20StamGridSide::all() {
        for reversed_depth in [false, true] {
            let rows = kp20_pr08_stam_projected_tensor_border_side_rows(
                &projected,
                Kp20StamPatch::LowerRight,
                side,
                reversed_depth,
            )
            .unwrap();
            let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
            let solved = solve_kp20_planar_tau_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: 5,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
            if solved.max_constraint_residual < best {
                best = solved.max_constraint_residual;
                best_side = side;
                best_reversed = reversed_depth;
            }
        }
    }

    assert!(
        best > 107.5 && best < 107.8,
        "PR08 middle-normalized projection should remain rejected near 107.64, got {best} side {} reversed={best_reversed}",
        best_side.name()
    );
    assert_eq!(best_side, Kp20StamGridSide::Right);
    assert!(best_reversed);
}

#[test]
fn pr08_stam_visual_projected_sides_are_rejected_by_paper_tau_constraints() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    for patch in [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ] {
        for side in Kp20StamGridSide::all() {
            for reversed_depth in [false, true] {
                let rows = kp20_pr08_stam_projected_tensor_border_side_rows(
                    &projected,
                    patch,
                    side,
                    reversed_depth,
                )
                .unwrap();
                let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                let solved = solve_kp20_planar_tau_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                    valence: 5,
                    sector: 0,
                    characteristic_tensor_border: &tensor_border,
                    initial: Kp20MacroPatch::zero(),
                    regularization_weight: 1e-10,
                })
                .unwrap();

                assert!(
                    solved.max_constraint_residual > 1.0,
                    "visual PR08-to-Stam projection unexpectedly passed patch {:?} side {} reversed={}: residual {}",
                    patch,
                    side.name(),
                    reversed_depth,
                    solved.max_constraint_residual
                );
            }
        }
    }
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source visual projection family diagnostic"]
fn pr08_stam_visual_best_candidate_residual_is_localized_by_family() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let rows = kp20_pr08_stam_projected_tensor_border_side_rows(
        &projected,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridSide::Bottom,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1e-10,
    })
    .unwrap();
    let residuals =
        kp20_planar_tau_paper_family_residuals(5, 0, &tensor_border, &solved.tau_macro_patch)
            .unwrap();

    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::SectorG1G2NearMidpoint
            && residual.max_residual > 45.0
            && residual.max_residual < 46.0
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::TensorBorderEq14
            && residual.max_residual < 1.0e-6
    }));
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source visual projection residual diagnostic"]
fn pr08_stam_visual_best_candidate_rejects_current_rotational_sector_rows() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let rows = kp20_pr08_stam_projected_tensor_border_side_rows(
        &projected,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridSide::Bottom,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1e-10,
    })
    .unwrap();

    assert!(
        solved.max_constraint_residual > 59.0 && solved.max_constraint_residual < 61.0,
        "expected current best rotational sector-row residual near 59.77, got {}",
        solved.max_constraint_residual
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source visual rotational orientation diagnostic"]
fn pr08_stam_visual_rotational_sector_orientation_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let rows = kp20_pr08_stam_projected_tensor_border_side_rows(
        &projected,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridSide::Bottom,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();

    let mut best_residual = f64::INFINITY;
    let mut best_label = "";
    for (label, previous_o, current_o, previous_m, current_m) in [
        (
            "right-to-left",
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "left-to-right",
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        ),
    ] {
        for (rotation_label, rotation) in [("-angle", -angle), ("+angle", angle)] {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                previous_o,
                current_o,
                kp20_b_coeffs_near_o(c),
                rotation,
            ));
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                previous_m,
                current_m,
                kp20_b_coeffs_near_m(c),
                rotation,
            ));
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints,
                    regularization_weight: 1e-10,
                },
            )
            .unwrap();
            if solved.max_constraint_residual < best_residual {
                best_residual = solved.max_constraint_residual;
                best_label = if label == "right-to-left" && rotation_label == "+angle" {
                    "right-to-left/+angle"
                } else {
                    "other"
                };
            }
        }
    }

    assert_eq!(best_label, "right-to-left/+angle");
    assert!(
        best_residual > 33.0 && best_residual < 34.0,
        "expected best rotational orientation residual near 33.66, got {best_residual}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source transformed visual side diagnostic"]
fn pr08_stam_visual_transformed_sides_do_not_reach_best_row_map_frontier() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projected,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridTransform::Rotate180,
        Kp20StamGridSide::Top,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1e-10,
    })
    .unwrap();

    assert!(
        solved.max_constraint_residual > 61.0 && solved.max_constraint_residual < 63.0,
        "expected transformed visual side residual near 61.95, got {}",
        solved.max_constraint_residual
    );
}
