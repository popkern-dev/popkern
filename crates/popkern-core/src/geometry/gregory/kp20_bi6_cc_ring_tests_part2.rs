#[allow(unused_imports)]
use super::super::super::kp_polar_shared::weighted_sum;
use super::*;

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_scaled_appendix_cubic_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();

    for (label, tau_family) in [
        ("appendix-a", Kp20TauFamily::Paper),
        ("appendix-a-refined", Kp20TauFamily::PaperRefinedApex),
    ] {
        let local_appendix_tau = sample_kp20_macro_patch_from_tau_family_maps(
            &cc_test_local_identity_tri6_guide(),
            valence,
            sector,
            tau_family,
            kp20_paper_quadrant_sector_piece_maps(),
        )
        .expect("Appendix-A tau samples over paper quadrants");
        let appendix_tau =
            cc_test_local_planar_patch_to_global(valence, sector, local_appendix_tau);
        let vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));

        let mut cubic_rows = Vec::new();
        for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            for i in 0..=3 {
                for j in 0..=(3 - i) {
                    let target = appendix_tau.piece(piece_id).cps[i][j];
                    let coeff = Kp20PlanarTauCoeff::new(piece_id, i, j);
                    for (axis, target_axis) in [target.x, target.y].into_iter().enumerate() {
                        cubic_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                            terms: vec![
                                (
                                    Kp20PlanarTauAugmentedScalarCoeff::Tau(
                                        Kp20PlanarTauScalarCoeff::new(coeff, axis),
                                    ),
                                    1.0,
                                ),
                                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_axis),
                            ],
                            target: 0.0,
                        });
                    }
                }
            }
        }

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_scale) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &vector_constraints,
                &scalar_constraints,
                &cubic_rows,
                1,
            );
        let mut final_vector_constraints = vector_constraints;
        final_vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: final_vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: cubic_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 self-diagonal scaled {label} cubic: residual={residual:e} projection_residual={projection_residual:e} source_delta={source_delta:e} projection_scale={projection_scale:e} final_scale={final_scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            projection_residual = projection_residual,
            source_delta = max_source_delta,
            projection_scale = projection_scale[0],
            final_scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            p00x = p00.x,
            p00y = p00.y,
            dux = du[0],
            duy = du[1],
            dvx = dv[0],
            dvy = dv[1],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_left_scaled_appendix_cubic_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();

    for (label, tau_family) in [
        ("appendix-a", Kp20TauFamily::Paper),
        ("appendix-a-refined", Kp20TauFamily::PaperRefinedApex),
    ] {
        let local_appendix_tau = sample_kp20_macro_patch_from_tau_family_maps(
            &cc_test_local_identity_tri6_guide(),
            valence,
            sector,
            tau_family,
            kp20_paper_quadrant_sector_piece_maps(),
        )
        .expect("Appendix-A tau samples over paper quadrants");
        let appendix_tau =
            cc_test_local_planar_patch_to_global(valence, sector, local_appendix_tau);
        let vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));

        let mut cubic_rows = Vec::new();
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                let target = appendix_tau.piece(Kp20PieceId::NearOriginLeft).cps[i][j];
                let coeff = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, i, j);
                for (axis, target_axis) in [target.x, target.y].into_iter().enumerate() {
                    cubic_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                        terms: vec![
                            (
                                Kp20PlanarTauAugmentedScalarCoeff::Tau(
                                    Kp20PlanarTauScalarCoeff::new(coeff, axis),
                                ),
                                1.0,
                            ),
                            (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_axis),
                        ],
                        target: 0.0,
                    });
                }
            }
        }

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_scale) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &vector_constraints,
                &scalar_constraints,
                &cubic_rows,
                1,
            );
        let mut final_vector_constraints = vector_constraints;
        final_vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: final_vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: cubic_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 self-diagonal left scaled {label} cubic: residual={residual:e} projection_residual={projection_residual:e} source_delta={source_delta:e} projection_scale={projection_scale:e} final_scale={final_scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            projection_residual = projection_residual,
            source_delta = max_source_delta,
            projection_scale = projection_scale[0],
            final_scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            p00x = p00.x,
            p00y = p00.y,
            dux = du[0],
            duy = du[1],
            dvx = dv[0],
            dvy = dv[1],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_pair_diagonal_component_variants_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let variants = [
        ("none", None),
        ("swap", Some([local_axes[1], local_axes[0]])),
        ("same", Some([local_axes[0], local_axes[1]])),
        (
            "swap-neg",
            Some([
                cc_test_neg_component(local_axes[1]),
                cc_test_neg_component(local_axes[0]),
            ]),
        ),
        (
            "same-neg",
            Some([
                cc_test_neg_component(local_axes[0]),
                cc_test_neg_component(local_axes[1]),
            ]),
        ),
        (
            "first-neg-swap",
            Some([cc_test_neg_component(local_axes[1]), local_axes[0]]),
        ),
        (
            "second-neg-swap",
            Some([local_axes[1], cc_test_neg_component(local_axes[0])]),
        ),
    ];

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (label, upper_components) in variants {
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        if let Some(upper_components) = upper_components {
            scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows_with_components(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes,
                upper_components,
            ));
        }
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let augmented_scalar_constraints =
            cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 pair variant {label}: residual={residual:e} source_delta={max_source_delta:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if validity.domain.outside_triangle_samples < best.1
            || (validity.domain.outside_triangle_samples == best.1
                && validity.orientation.negative_determinants < best.2)
            || (validity.domain.outside_triangle_samples == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                label,
                validity.domain.outside_triangle_samples,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                max_source_delta,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!("scaled q13/q2 best pair variant diagnostic: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_fig8_diagonal_piece_correspondence_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();
    let maps = kp20_section42_map_candidates(angle);
    let cases = [
        (
            "corrected-current",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "same-side-orientation-winner",
            0_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "corrected-low-neg-rejected",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointRight,
        ),
        (
            "corrected-algebraic-collapse",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ),
    ];

    let mut best = (
        "",
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearMidpointRight,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    let mut current_layout = None;

    for (case_label, map_index, self_a, self_b, lower_right, upper_left) in cases {
        let (map_label, map) = maps[map_index];
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            self_a,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            self_b,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows(
            lower_right,
            upper_left,
            local_axes[0],
            local_axes[1],
        ));

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let augmented_scalar_constraints =
            cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();

        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let score = (
            validity.orientation.negative_determinants,
            square.outside,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            max_source_delta,
            -validity.orientation.min_determinant,
        );
        eprintln!(
            "scaled q13/q2 Fig.8 diagonal correspondence case={case_label} map={map_label} self=({self_a:?},{self_b:?}) pair=({lower_right:?},{upper_left:?}) residual={residual:e} source_delta={source_delta:e} scale={scale:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
            residual = solved.max_constraint_residual,
            source_delta = max_source_delta,
            scale = solved.extra_values[0],
            square_outside = square.outside,
            square_samples = square.samples,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            min_u = square.min_u,
            max_u = square.max_u,
            min_v = square.min_v,
            max_v = square.max_v,
        );

        if case_label == "corrected-current" {
            current_layout = Some(score);
        }
        if score < (best.5, best.6, best.7, best.8, best.9, best.10) {
            best = (
                case_label,
                self_a,
                self_b,
                lower_right,
                upper_left,
                score.0,
                score.1,
                score.2,
                score.3,
                score.4,
                score.5,
            );
        }
    }

    eprintln!(
        "scaled q13/q2 Fig.8 diagonal correspondence best map={} self=({:?},{:?}) pair=({:?},{:?}) neg={} square_outside={} triangle_outside={} residual={:e} source_delta={:e} neg_min_det={:e} corrected_current={current_layout:?}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7, best.8, best.9, best.10
    );
    assert!(current_layout.is_some());
    assert!(best.8.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_pair_diagonal_lemma_owned_rows_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
        (
            "exclude-apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i != 0 || j != 0,
            ),
        ),
        (
            "exclude-first-order-square",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i > 1 || j > 1,
            ),
        ),
        (
            "exclude-total-degree-one",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 1,
            ),
        ),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
    ] {
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(pair_rows);

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let augmented_scalar_constraints =
            cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 pair lemma-owned {label}: residual={residual:e} source_delta={source_delta:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e}",
            residual = solved.max_constraint_residual,
            source_delta = max_source_delta,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn kp20_fig8_paper_quadrant_pair_correspondence_is_derived_from_domains() {
    let lower_right = Kp20AffineGuideMap {
        origin: [0.5, 0.0],
        du: [0.5, 0.0],
        dv: [0.0, 0.5],
    };
    let upper_left_coefficient_chart = Kp20AffineGuideMap {
        origin: [0.5, 1.0],
        du: [-0.5, 0.0],
        dv: [0.0, -0.5],
    };

    for i in 0..=6 {
        for j in 0..=6 {
            let lower_point = cc_test_affine_map_point(lower_right, i as f64 / 6.0, j as f64 / 6.0);
            let diagonal_image = [lower_point[1], lower_point[0]];
            let upper_uv = cc_test_affine_map_inverse(upper_left_coefficient_chart, diagonal_image)
                .expect("paper quadrant upper-left coefficient chart is nonsingular");
            let mapped_i = cc_test_nearest_degree6_index(upper_uv[0]);
            let mapped_j = cc_test_nearest_degree6_index(upper_uv[1]);
            assert_eq!(
                (mapped_i, mapped_j),
                (6 - j, 6 - i),
                "KP20 Fig.8 coefficient-chart diagonal should map NearOriginRight[{i},{j}] to NearMidpointLeft[{},{}]",
                6 - j,
                6 - i,
            );
        }
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_with_exact_pair_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &centered_rows,
            1,
        );

    for (label, include_centered_in_final) in [
        ("projection-only", false),
        ("projection-and-centered-final", true),
    ] {
        let mut vector_constraints = fig8_rows.clone();
        vector_constraints.extend(projected_tensor_rows.clone());
        let (tau_macro_patch, max_constraint_residual, final_scale) = if include_centered_in_final {
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            (
                solved.tau_macro_patch,
                solved.max_constraint_residual,
                Some(solved.extra_values[0]),
            )
        } else {
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            (solved.tau_macro_patch, solved.max_constraint_residual, None)
        };
        let local_tau = cc_test_global_planar_patch_to_local(valence, sector, tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 exact-pair Lemma4 source projection {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            projected_scale = projection_extra_values[0],
            residual = max_constraint_residual,
            final_scale = final_scale.unwrap_or(f64::NAN),
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            p00x = p00.x,
            p00y = p00.y,
            dux = du[0],
            duy = du[1],
            dvx = dv[0],
            dvy = dv[1],
        );
        assert!(max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_cubic_lemma4_projection_with_exact_pair_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let cubic_rows = cc_test_paper_refined_apex_cubic_derivative_rows(
        valence,
        sector,
        Kp20PieceId::NearOriginLeft,
        local_axes,
        0,
    );
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &cubic_rows,
            1,
        );

    let mut vector_constraints = fig8_rows;
    vector_constraints.extend(projected_tensor_rows);
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: cubic_rows,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
    let p00 = piece.cps[0][0];
    let p10 = piece.cps[1][0];
    let p01 = piece.cps[0][1];
    let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
    let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
    let apex_det = du[0] * dv[1] - dv[0] * du[1];
    eprintln!(
        "scaled q13/q2 exact-pair cubic-Lemma4 source projection: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
        projected_scale = projection_extra_values[0],
        final_scale = solved.extra_values[0],
        residual = solved.max_constraint_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
        p00x = p00.x,
        p00y = p00.y,
        dux = du[0],
        duy = du[1],
        dvx = dv[0],
        dvy = dv[1],
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_without_full_diagonal_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    for (label, centered_rows) in [
        (
            "first-order",
            cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0),
        ),
        (
            "cubic-derivative",
            cc_test_paper_refined_apex_cubic_derivative_rows(
                valence,
                sector,
                Kp20PieceId::NearOriginLeft,
                local_axes,
                0,
            ),
        ),
    ] {
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: centered_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 no-full-diagonal Lemma4 source projection {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            p00x = p00.x,
            p00y = p00.y,
            dux = du[0],
            duy = du[1],
            dvx = dv[0],
            dvy = dv[1],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_fig8_average_completions_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (label, fig8_rows) in [
        (
            "paper-single",
            cc_test_fig8_indicated_average_variant_rows(false, false),
        ),
        (
            "diagonal-mirror",
            cc_test_fig8_indicated_average_variant_rows(true, false),
        ),
        (
            "rotated-quadrants",
            cc_test_fig8_indicated_average_variant_rows(false, true),
        ),
        (
            "rotated-quadrants-plus-mirror",
            cc_test_fig8_indicated_average_variant_rows(true, true),
        ),
    ] {
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );

        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows.clone(),
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 Fig8 average completion {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            worst = validity.orientation.worst_piece,
            u = validity.orientation.worst_u,
            v = validity.orientation.worst_v,
        );
        if validity.domain.outside_triangle_samples < best.1
            || (validity.domain.outside_triangle_samples == best.1
                && validity.orientation.negative_determinants < best.2)
            || (validity.domain.outside_triangle_samples == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                label,
                validity.domain.outside_triangle_samples,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                max_source_delta,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!("scaled q13/q2 best Fig8 average completion diagnostic: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_deltas_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &centered_rows,
            1,
        );

    let mut deltas = original_tensor_rows
        .iter()
        .zip(projected_tensor_rows.iter())
        .enumerate()
        .map(|(row_index, (original, projected))| {
            let coeff = original.terms[0].0;
            let dx = projected.target[0] - original.target[0];
            let dy = projected.target[1] - original.target[1];
            let norm = (dx * dx + dy * dy).sqrt();
            (
                norm,
                row_index,
                coeff,
                dx,
                dy,
                original.target,
                projected.target,
            )
        })
        .collect::<Vec<_>>();
    deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    for (rank, (norm, row_index, coeff, dx, dy, original, projected)) in
        deltas.iter().take(12).enumerate()
    {
        eprintln!(
            "scaled q13/q2 Lemma4 eq14 projection delta rank={rank} row={row_index} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
            piece = coeff.piece,
            i = coeff.i,
            j = coeff.j,
            ox = original[0],
            oy = original[1],
            px = projected[0],
            py = projected[1],
        );
    }
    eprintln!(
        "scaled q13/q2 Lemma4 eq14 projection delta summary: projection_residual={projection_residual:e} max_source_delta={max_source_delta:e} projected_scale={projected_scale:e}",
        projected_scale = projection_extra_values[0],
    );
    assert!(projection_residual.is_finite());
    assert!(max_source_delta.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_c_origin_eq14_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (label, normalize_left, normalize_right) in [
        ("none", false, false),
        ("left-c-origin", true, false),
        ("right-c-origin", false, true),
        ("both-c-origin", true, true),
    ] {
        let mut tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let left_c = tensor_rows[6].target;
        let right_c = tensor_rows[21 + 6].target;
        for (row_index, row) in tensor_rows.iter_mut().enumerate() {
            let origin = if row_index < 21 {
                normalize_left.then_some(left_c)
            } else {
                normalize_right.then_some(right_c)
            };
            if let Some(origin) = origin {
                row.target[0] -= origin[0];
                row.target[1] -= origin[1];
            }
        }
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
        let mut vector_constraints = fig8_rows.clone();
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows.clone(),
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 c-origin eq14 normalization {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            worst = validity.orientation.worst_piece,
            u = validity.orientation.worst_u,
            v = validity.orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_per_sector_characteristic_eq14_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let normalized_rows = cc_test_per_sector_characteristic_tensor_rows(&raw_rows, valence, sector);
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (label, rows) in [
        ("raw-global", raw_rows),
        ("per-sector-characteristic", normalized_rows),
    ] {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();

        for include_projection in [false, true] {
            let (source_rows, source_delta, projection_residual, projected_scale) =
                if include_projection {
                    let (
                        projected_tensor_rows,
                        max_source_delta,
                        projection_residual,
                        projection_extra_values,
                    ) = cc_test_project_tensor_rows_for_augmented_scalar_basis(
                        tensor_rows.clone(),
                        &fig8_rows,
                        &scalar_constraints,
                        &centered_rows,
                        1,
                    );
                    (
                        projected_tensor_rows,
                        max_source_delta,
                        projection_residual,
                        projection_extra_values[0],
                    )
                } else {
                    (tensor_rows.clone(), 0.0, 0.0, f64::NAN)
                };
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(source_rows);
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
            let p00 = piece.cps[0][0];
            let p10 = piece.cps[1][0];
            let p01 = piece.cps[0][1];
            let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
            let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
            let apex_det = du[0] * dv[1] - dv[0] * du[1];
            eprintln!(
                "scaled q13/q2 {label} eq14 characteristic normalization projection={include_projection}: projection_residual={projection_residual:e} source_delta={source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} worst={worst:?}@({u:.2},{v:.2})",
                final_scale = solved.extra_values[0],
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
                worst = validity.orientation.worst_piece,
                u = validity.orientation.worst_u,
                v = validity.orientation.worst_v,
            );
            assert!(solved.max_constraint_residual.is_finite());
        }
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_characteristic_dyadic_level_with_lemma4_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    let mut best = (
        "",
        0,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (normalization_label, normalized_base) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
        ),
    ] {
        for level_exponent in -2..=10 {
            let rows = scale_tensor_border_rows(&normalized_base, ring.lambda.powi(level_exponent));
            let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
            let tensor_rows =
                kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
            let (
                projected_tensor_rows,
                max_source_delta,
                projection_residual,
                projection_extra_values,
            ) = cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "scaled q13/q2 characteristic dyadic {normalization_label} level={level_exponent}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                projected_scale = projection_extra_values[0],
                final_scale = solved.extra_values[0],
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            if validity.orientation.negative_determinants < best.3
                || (validity.orientation.negative_determinants == best.3
                    && validity.domain.outside_triangle_samples < best.2)
                || (validity.orientation.negative_determinants == best.3
                    && validity.domain.outside_triangle_samples == best.2
                    && validity.orientation.min_determinant > best.4)
            {
                best = (
                    normalization_label,
                    level_exponent,
                    validity.domain.outside_triangle_samples,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    max_source_delta,
                    solved.max_constraint_residual,
                );
            }
        }
    }
    eprintln!("scaled q13/q2 best characteristic dyadic Lemma4 diagnostic: {best:?}");
    assert!(best.6.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_linear_characteristic_source_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (normalization_label, rows) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
        ),
    ] {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let mut augmented_scalar_constraints = centered_rows.clone();
        augmented_scalar_constraints
            .extend(cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1));
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: fig8_rows.clone(),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: 5,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let a00 = 1.0 + solved.extra_values[1];
        let a01 = solved.extra_values[2];
        let a10 = solved.extra_values[3];
        let a11 = 1.0 + solved.extra_values[4];
        let det_a = a00 * a11 - a01 * a10;
        let frob_delta = (solved.extra_values[1] * solved.extra_values[1]
            + solved.extra_values[2] * solved.extra_values[2]
            + solved.extra_values[3] * solved.extra_values[3]
            + solved.extra_values[4] * solved.extra_values[4])
            .sqrt();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 linear characteristic source normalization {normalization_label}: residual={residual:e} center_scale={center_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} frob_delta={frob_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            worst = validity.orientation.worst_piece,
            u = validity.orientation.worst_u,
            v = validity.orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
        assert!(det_a.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_linear_source_collapse_is_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let rows = cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector);
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair = cc_test_sector_local_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
    );

    for (label, diagonal_blocks) in [
        (
            "section42-no-diagonal",
            Vec::<&[Kp20PlanarTauScalarLinearConstraint]>::new(),
        ),
        ("section42-self-origin", vec![self_origin.as_slice()]),
        ("section42-self-midpoint", vec![self_midpoint.as_slice()]),
        (
            "section42-self-diagonals",
            vec![self_origin.as_slice(), self_midpoint.as_slice()],
        ),
        ("section42-pair-only", vec![pair.as_slice()]),
        (
            "section42-all-diagonal",
            vec![
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair.as_slice(),
            ],
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        for block in diagonal_blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints = centered_rows.clone();
        augmented_scalar_constraints
            .extend(cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1));
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: fig8_rows.clone(),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: 5,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let a00 = 1.0 + solved.extra_values[1];
        let a01 = solved.extra_values[2];
        let a10 = solved.extra_values[3];
        let a11 = 1.0 + solved.extra_values[4];
        let det_a = a00 * a11 - a01 * a10;
        eprintln!(
            "scaled q13/q2 linear source collapse {label}: residual={residual:e} center_scale={center_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e}",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[0],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_anchored_linear_characteristic_source_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );

    for (normalization_label, rows) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
        ),
    ] {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
            let mut augmented_scalar_constraints = cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                0,
            );
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: anchored_scale,
            });
            augmented_scalar_constraints
                .extend(cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1));
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: fig8_rows.clone(),
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints,
                    extra_variable_count: 5,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            let a00 = 1.0 + solved.extra_values[1];
            let a01 = solved.extra_values[2];
            let a10 = solved.extra_values[3];
            let a11 = 1.0 + solved.extra_values[4];
            let det_a = a00 * a11 - a01 * a10;
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "scaled q13/q2 anchored linear characteristic source {normalization_label} scale={anchored_scale:e}: residual={residual:e} final_scale={final_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                final_scale = solved.extra_values[0],
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(solved.max_constraint_residual.is_finite());
        }
    }
}

#[test]
fn stam_appendix_a_regular_subpatch_anchored_characteristic_source_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );

    let mut best = (f64::INFINITY, 0, false, 0, false, false, 0.0, 0.0, 0_usize);
    for patch_index in 0..3 {
        for transposed in [false, true] {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    for reverse_cols in [false, true] {
                        let rows = stam_appendix_b_projected_rows(
                            &ring,
                            patch_index,
                            transposed,
                            row_start,
                            reverse_rows,
                            reverse_cols,
                        )
                        .unwrap();
                        let rows =
                            cc_test_per_sector_characteristic_tensor_rows(&rows, valence, sector);
                        let tensor_border =
                            kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                        let tensor_rows =
                            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector)
                                .unwrap();
                        for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
                            let mut augmented_scalar_constraints =
                                cc_test_first_order_centered_source_rows(
                                    Kp20PieceId::NearOriginLeft,
                                    local_axes,
                                    0,
                                );
                            augmented_scalar_constraints.push(
                                Kp20PlanarTauAugmentedScalarLinearConstraint {
                                    terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                                    target: anchored_scale,
                                },
                            );
                            augmented_scalar_constraints.extend(
                                cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1),
                            );
                            let solved =
                                solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                                    &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                                        initial: Kp20MacroPatch::zero(),
                                        vector_constraints: fig8_rows.clone(),
                                        scalar_constraints: scalar_constraints.clone(),
                                        augmented_scalar_constraints,
                                        extra_variable_count: 5,
                                        regularization_weight: 1.0e-10,
                                        extra_regularization_weight: 1.0e-10,
                                    },
                                )
                                .unwrap();
                            let a00 = 1.0 + solved.extra_values[1];
                            let a01 = solved.extra_values[2];
                            let a10 = solved.extra_values[3];
                            let a11 = 1.0 + solved.extra_values[4];
                            let det_a = a00 * a11 - a01 * a10;
                            let local_tau = cc_test_global_planar_patch_to_local(
                                valence,
                                sector,
                                solved.tau_macro_patch,
                            );
                            let validity = kp20_planar_bi6_tau_validity_report(
                                &local_tau,
                                &[0.0, 0.5, 1.0],
                                1.0e-10,
                                1.0e-10,
                            );
                            if solved.max_constraint_residual < best.0 {
                                best = (
                                    solved.max_constraint_residual,
                                    patch_index,
                                    transposed,
                                    row_start,
                                    reverse_rows,
                                    reverse_cols,
                                    anchored_scale,
                                    det_a,
                                    validity.orientation.negative_determinants,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!("Stam regular subpatch anchored characteristic source best: {best:?}");
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_promoted_section42_lemma4_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let diagonal_rows =
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, 0).unwrap();

    let mut best = (
        f64::INFINITY,
        "",
        0,
        0,
        false,
        0.0,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
    );
    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for patch_index in 0..3 {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
                        let mut augmented_scalar_constraints =
                            stam_appendix_a_unknown_top_vector_tensor_border_constraints(
                                &ring,
                                patch_index,
                                row_start,
                                reverse_rows,
                            );
                        augmented_scalar_constraints.extend(
                            cc_test_first_order_centered_source_rows(
                                Kp20PieceId::NearOriginLeft,
                                local_axes,
                                source_extra_count,
                            ),
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(
                                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                                    1.0,
                                )],
                                target: anchored_scale,
                            },
                        );
                        let mut scalar_constraints =
                            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
                        scalar_constraints.extend(diagonal_rows.clone());
                        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                                initial: Kp20MacroPatch::zero(),
                                vector_constraints: fig8_rows.clone(),
                                scalar_constraints,
                                augmented_scalar_constraints,
                                extra_variable_count: source_extra_count + 1,
                                regularization_weight: 1.0e-10,
                                extra_regularization_weight: 1.0e-12,
                            },
                        )
                        .unwrap();
                        let local_tau = cc_test_global_planar_patch_to_local(
                            valence,
                            0,
                            solved.tau_macro_patch,
                        );
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                        );
                        eprintln!(
                            "unknown top promoted §4.2 Lemma4 map={map_label} patch={patch_index} row_start={row_start} reverse={reverse_rows} scale={anchored_scale:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                            residual = solved.max_constraint_residual,
                            outside = validity.domain.outside_triangle_samples,
                            domain_samples = validity.domain.sample_count,
                            neg = validity.orientation.negative_determinants,
                            orientation_samples = validity.orientation.sample_count,
                            zero = validity.orientation.near_zero_determinants,
                            min_det = validity.orientation.min_determinant,
                            max_det = validity.orientation.max_determinant,
                        );
                        if solved.max_constraint_residual < best.0
                            || (solved.max_constraint_residual <= best.0 + 1.0e-12
                                && validity.orientation.negative_determinants < best.7)
                            || (solved.max_constraint_residual <= best.0 + 1.0e-12
                                && validity.orientation.negative_determinants == best.7
                                && validity.orientation.min_determinant > best.8)
                        {
                            best = (
                                solved.max_constraint_residual,
                                map_label,
                                patch_index,
                                row_start,
                                reverse_rows,
                                anchored_scale,
                                validity.domain.outside_triangle_samples,
                                validity.orientation.negative_determinants,
                                validity.orientation.min_determinant,
                            );
                        }
                    }
                }
            }
        }
    }
    eprintln!("unknown top promoted §4.2 Lemma4 best: {best:?}");
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_lemma4_zero_lock_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);

    for (label, vector_constraints, scalar_constraints) in [
        (
            "eq14-source-only",
            Vec::<Kp20PlanarTauLinearConstraint>::new(),
            Vec::<Kp20PlanarTauScalarLinearConstraint>::new(),
        ),
        (
            "fig8-eq14",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            Vec::<Kp20PlanarTauScalarLinearConstraint>::new(),
        ),
        (
            "section42-no-diagonal",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
        ),
        (
            "section42-diagonal",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            {
                let mut rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
                rows.extend(
                    kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, 0).unwrap(),
                );
                rows
            },
        ),
    ] {
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "unknown top Lemma4 zero-lock {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_unknown_top_vector_diagonal_block_zero_lock_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair = cc_test_sector_local_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
    );

    for (label, diagonal_blocks) in [
        ("none", Vec::<&[Kp20PlanarTauScalarLinearConstraint]>::new()),
        ("self-origin", vec![self_origin.as_slice()]),
        ("self-midpoint", vec![self_midpoint.as_slice()]),
        (
            "self-both",
            vec![self_origin.as_slice(), self_midpoint.as_slice()],
        ),
        ("pair-only", vec![pair.as_slice()]),
        (
            "all",
            vec![
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair.as_slice(),
            ],
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        for block in diagonal_blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "unknown top diagonal block {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_unknown_top_vector_pair_diagonal_component_variants_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);
    let variants = [
        ("none", None),
        ("swap", Some([local_axes[1], local_axes[0]])),
        ("same", Some([local_axes[0], local_axes[1]])),
        (
            "swap-neg",
            Some([
                cc_test_neg_component(local_axes[1]),
                cc_test_neg_component(local_axes[0]),
            ]),
        ),
        (
            "same-neg",
            Some([
                cc_test_neg_component(local_axes[0]),
                cc_test_neg_component(local_axes[1]),
            ]),
        ),
        (
            "first-neg-swap",
            Some([cc_test_neg_component(local_axes[1]), local_axes[0]]),
        ),
        (
            "second-neg-swap",
            Some([local_axes[1], cc_test_neg_component(local_axes[0])]),
        ),
    ];

    let mut best = ("", f64::INFINITY, usize::MAX, f64::NEG_INFINITY);
    for (label, upper_components) in variants {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        if let Some(upper_components) = upper_components {
            scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows_with_components(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes,
                upper_components,
            ));
        }
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "unknown top pair variant {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if solved.max_constraint_residual < best.1
            || (solved.max_constraint_residual <= best.1 + 1.0e-12
                && validity.orientation.negative_determinants < best.2)
            || (solved.max_constraint_residual <= best.1 + 1.0e-12
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                label,
                solved.max_constraint_residual,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
            );
        }
    }
    eprintln!("unknown top best pair variant: {best:?}");
    assert!(best.1.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_pair_diagonal_subset_zero_lock_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 0 && j == 0,
            ),
        ),
        (
            "first-order-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i <= 1 && j <= 1,
            ),
        ),
        (
            "near-origin-block",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i <= 2 && j <= 2,
            ),
        ),
        (
            "domain-diagonal",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == j,
            ),
        ),
        (
            "outer-boundary",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 6 || j == 6,
            ),
        ),
        (
            "full",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(pair_rows);
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "unknown top pair subset {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_unknown_top_vector_pair_diagonal_lemma_owned_rows_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
        (
            "exclude-apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i != 0 || j != 0,
            ),
        ),
        (
            "exclude-first-order-square",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i > 1 || j > 1,
            ),
        ),
        (
            "exclude-total-degree-one",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 1,
            ),
        ),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(pair_rows);
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "unknown top pair lemma-owned {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_unknown_top_vector_lemma_owned_pair_source_modes_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    ));
    scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    ));
    scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    ));
    let mut augmented_scalar_constraints = source_constraints;
    augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        source_extra_count,
    ));
    augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();
    assert!(solved.max_constraint_residual.is_finite());

    let odd = (0..valence)
        .map(|sector| {
            Complex::new(
                solved.extra_values[1 + 2 * sector],
                solved.extra_values[top_count + 1 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let even = (0..valence)
        .map(|sector| {
            Complex::new(
                solved.extra_values[2 + 2 * sector],
                solved.extra_values[top_count + 2 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let odd_residual = fourier_mode_residual(&odd, 1);
    let even_residual = fourier_mode_residual(&even, 1);
    let odd_magnitudes = fourier_mode_magnitudes(&odd);
    let even_magnitudes = fourier_mode_magnitudes(&even);
    eprintln!(
        "unknown top lemma-owned pair source modes: residual={residual:e} center_scale={center_scale:e} odd_mode1_residual={odd_residual:e} even_mode1_residual={even_residual:e} odd_modes={odd_magnitudes:?} even_modes={even_magnitudes:?}",
        residual = solved.max_constraint_residual,
        center_scale = solved.extra_values[source_extra_count],
    );
    assert!(odd_residual.is_finite());
    assert!(even_residual.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_anchored_lemma_owned_pair_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let source_constraints =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(pair_rows);
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(1), 1.0)],
            target: ring.sectors[0][1].x,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let odd = (0..valence)
            .map(|sector| {
                Complex::new(
                    solved.extra_values[1 + 2 * sector],
                    solved.extra_values[top_count + 1 + 2 * sector],
                )
            })
            .collect::<Vec<_>>();
        let even = (0..valence)
            .map(|sector| {
                Complex::new(
                    solved.extra_values[2 + 2 * sector],
                    solved.extra_values[top_count + 2 + 2 * sector],
                )
            })
            .collect::<Vec<_>>();
        eprintln!(
            "unknown top anchored lemma-owned pair {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] odd_mode1_residual={odd_residual:e} even_mode1_residual={even_residual:e}",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            odd_residual = fourier_mode_residual(&odd, 1),
            even_residual = fourier_mode_residual(&even, 1),
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = (
        f64::INFINITY,
        "",
        0_usize,
        0_usize,
        false,
        0.0_f64,
        f64::INFINITY,
        0_usize,
        0_usize,
    );

    for (pair_label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
        for rotation_direction in [-1.0, 1.0] {
            for patch_index in 0..3 {
                for row_start in [0, 1] {
                    for reverse_rows in [false, true] {
                        let mut scalar_constraints =
                            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
                        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
                            Kp20PieceId::NearOriginLeft,
                            local_axes[0],
                            local_axes[1],
                        ));
                        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
                            Kp20PieceId::NearMidpointRight,
                            local_axes[0],
                            local_axes[1],
                        ));
                        scalar_constraints.extend(pair_rows.clone());

                        let mut augmented_scalar_constraints =
                            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                                &ring,
                                patch_index,
                                row_start,
                                reverse_rows,
                                rotation_direction,
                            );
                        augmented_scalar_constraints.extend(
                            cc_test_first_order_centered_source_rows(
                                Kp20PieceId::NearOriginLeft,
                                local_axes,
                                source_extra_count,
                            ),
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                                target: ring.sectors[0][1].x,
                            },
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(
                                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                                    1.0,
                                )],
                                target: anchored_scale,
                            },
                        );
                        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                                initial: Kp20MacroPatch::zero(),
                                vector_constraints: cc_test_fig8_indicated_average_variant_rows(
                                    false, false,
                                ),
                                scalar_constraints,
                                augmented_scalar_constraints,
                                extra_variable_count: source_extra_count + 1,
                                regularization_weight: 1.0e-10,
                                extra_regularization_weight: 1.0e-12,
                            },
                        )
                        .unwrap();
                        let local_tau = cc_test_global_planar_patch_to_local(
                            valence,
                            0,
                            solved.tau_macro_patch,
                        );
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau,
                            &[0.0, 0.5, 1.0],
                            1.0e-10,
                            1.0e-10,
                        );
                        if solved.max_constraint_residual < best.0 {
                            best = (
                                solved.max_constraint_residual,
                                pair_label,
                                patch_index,
                                row_start,
                                reverse_rows,
                                rotation_direction,
                                solved.extra_values[source_extra_count],
                                validity.orientation.negative_determinants,
                                validity.orientation.near_zero_determinants,
                            );
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "rotational top lemma-owned pair best: residual={:e} pair={} patch={} row_start={} reverse={} rotation={} center_scale={:e} neg={} zero={}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7, best.8
    );
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_residual_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let vector_constraints = cc_test_fig8_indicated_average_variant_rows(false, false);
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_rows = cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    );
    let source_rows =
        stam_appendix_a_rotational_top_vector_tensor_border_constraints(&ring, 0, 0, false, -1.0);
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 4);
    let source_anchor = vec![Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: ring.sectors[0][1].x,
    }];
    let scale_anchor = vec![Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    }];

    let mut scalar_constraints = section42_rows.clone();
    scalar_constraints.extend(self_origin.clone());
    scalar_constraints.extend(self_midpoint.clone());
    scalar_constraints.extend(pair_rows.clone());
    let mut augmented_scalar_constraints = source_rows.clone();
    augmented_scalar_constraints.extend(centered_rows.clone());
    augmented_scalar_constraints.extend(source_anchor.clone());
    augmented_scalar_constraints.extend(scale_anchor.clone());
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: vector_constraints.clone(),
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();

    eprintln!(
        "rotational top lemma-owned residual families: total={total:e} fig8={fig8:e} section42={section42:e} self_origin={self_origin_res:e} self_midpoint={self_midpoint_res:e} pair_cut={pair_cut:e} eq14_source={eq14_source:e} centered={centered:e} source_anchor={source_anchor_res:e} scale_anchor={scale_anchor_res:e}",
        total = solved.max_constraint_residual,
        fig8 = cc_test_vector_rows_max_residual(&solved.tau_macro_patch, &vector_constraints),
        section42 = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &section42_rows),
        self_origin_res = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &self_origin),
        self_midpoint_res =
            cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &self_midpoint),
        pair_cut = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &pair_rows),
        eq14_source = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &source_rows
        ),
        centered = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &centered_rows
        ),
        source_anchor_res = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &source_anchor
        ),
        scale_anchor_res = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &scale_anchor
        ),
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_source_anchor_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    ));
    scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    ));
    scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    ));

    for anchor_extra in 0..source_extra_count {
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(anchor_extra), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        eprintln!(
            "rotational top lemma-owned source anchor extra{anchor_extra}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_source_collapse_localizes_by_section42_family() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_cut = cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    );
    let cases: Vec<(&str, Vec<&[Kp20PlanarTauScalarLinearConstraint]>)> = vec![
        ("fig8-eq14-only", vec![]),
        ("c0", vec![c0_rows.as_slice()]),
        ("eq11", vec![eq11_rows.as_slice()]),
        ("eq12", vec![eq12_rows.as_slice()]),
        ("eq13", vec![eq13_rows.as_slice()]),
        (
            "section42",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "section42-self",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
                self_origin.as_slice(),
                self_midpoint.as_slice(),
            ],
        ),
        (
            "section42-self-pair-cut",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair_cut.as_slice(),
            ],
        ),
    ];

    for (label, row_blocks) in cases {
        let mut scalar_constraints = Vec::new();
        for block in row_blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        eprintln!(
            "rotational top source collapse section42 {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_c0_collapse_sweeps_symbol_maps() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for (c0_label, include_plain, include_underlined) in [
            ("plain", true, false),
            ("underlined", false, true),
            ("both", true, true),
        ] {
            let mut scalar_constraints = Vec::new();
            for i in 0..=6 {
                if include_plain {
                    add_kp20_section42_vector_rows(
                        &mut scalar_constraints,
                        map,
                        &[
                            section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                            section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
                        ],
                    );
                }
                if include_underlined {
                    add_kp20_section42_vector_rows(
                        &mut scalar_constraints,
                        map,
                        &[
                            section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                            section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                        ],
                    );
                }
            }
            let mut augmented_scalar_constraints =
                stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                    &ring, 0, 0, false, -1.0,
                );
            augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                source_extra_count,
            ));
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: 1.0,
            });
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                    1.0,
                )],
                target: anchored_scale,
            });
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                    scalar_constraints,
                    augmented_scalar_constraints,
                    extra_variable_count: source_extra_count + 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-12,
                },
            )
            .unwrap();
            eprintln!(
                "rotational top C0 map sweep {map_label}/{c0_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
                residual = solved.max_constraint_residual,
                center_scale = solved.extra_values[source_extra_count],
                s0 = solved.extra_values[0],
                s1 = solved.extra_values[1],
                s2 = solved.extra_values[2],
                s3 = solved.extra_values[3],
            );
            assert!(solved.max_constraint_residual.is_finite());
        }
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_full_section42_sweeps_symbol_maps() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            Kp20PieceId::NearMidpointRight,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
            local_axes[0],
            local_axes[1],
            |i, j| i + j > 3,
        ));
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        eprintln!(
            "rotational top full §4.2 map sweep {map_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) triangle_outside={triangle_outside}/{domain_samples} square_outside={square_outside}/{square_samples} square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            triangle_outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            square_outside = square.outside,
            square_samples = square.samples,
            min_u = square.min_u,
            max_u = square.max_u,
            min_v = square.min_v,
            max_v = square.max_v,
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_source_collapse_localizes_c0_companions() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let cases: Vec<(&str, Vec<&[Kp20PlanarTauScalarLinearConstraint]>)> = vec![
        ("c0", vec![c0_rows.as_slice()]),
        ("c0-eq11", vec![c0_rows.as_slice(), eq11_rows.as_slice()]),
        ("c0-eq12", vec![c0_rows.as_slice(), eq12_rows.as_slice()]),
        ("c0-eq13", vec![c0_rows.as_slice(), eq13_rows.as_slice()]),
        (
            "c0-eq11-eq12",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
            ],
        ),
        (
            "c0-eq11-eq13",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "c0-eq12-eq13",
            vec![
                c0_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "section42",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
    ];

    for (label, blocks) in cases {
        let mut scalar_constraints = Vec::new();
        for block in blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top source collapse C0 companions {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected §4.2 C0 companion set {label} should close with a nonzero rotational source, got residual {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected §4.2 C0 companion set {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_split_rows_are_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);

    let mut cases: Vec<(&str, Vec<Kp20PlanarTauScalarLinearConstraint>)> = Vec::new();
    cases.push(("baseline-c0-eq11", {
        let mut rows = c0_rows.clone();
        rows.extend(eq11_rows.clone());
        rows
    }));
    cases.push(("all-eq12", {
        let mut rows = c0_rows.clone();
        rows.extend(eq11_rows.clone());
        rows.extend(eq12_rows);
        rows
    }));
    for split in Kp20Section42Eq12Split::ALL {
        let mut rows = c0_rows.clone();
        rows.extend(eq11_rows.clone());
        rows.extend(kp20_section42_eq12_split_rows(map, split));
        cases.push((split.label(), rows));
    }

    for (label, scalar_constraints) in cases {
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected eq12 split {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected eq12 split {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected eq12 split {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_endpoint_depth_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);

    let endpoint_halves = [
        (
            "plain",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            "arrow",
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];

    for (half_label, upper, lower) in endpoint_halves {
        for j in 0..=2 {
            let mut scalar_constraints = c0_rows.clone();
            scalar_constraints.extend(eq11_rows.clone());
            scalar_constraints.extend(kp20_section42_eq12_endpoint_row(map, upper, lower, j));

            let mut augmented_scalar_constraints =
                stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                    &ring, 0, 0, false, -1.0,
                );
            augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                source_extra_count,
            ));
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: 1.0,
            });
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                    1.0,
                )],
                target: anchored_scale,
            });
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                    scalar_constraints,
                    augmented_scalar_constraints,
                    extra_variable_count: source_extra_count + 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-12,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
            eprintln!(
                "rotational top corrected eq12 endpoint {half_label} j={j}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                center_scale = solved.extra_values[source_extra_count],
                s0 = solved.extra_values[0],
                s1 = solved.extra_values[1],
                s2 = solved.extra_values[2],
                s3 = solved.extra_values[3],
                neg = validity.orientation.negative_determinants,
                samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(
                solved.max_constraint_residual < 1.0e-10,
                "corrected eq12 endpoint {half_label} j={j} should close algebraically, got {}",
                solved.max_constraint_residual
            );
        }
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_symbol_direction_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let current_pairs = [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];
    let reversed_pairs = [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ];
    let reversed_plain_pairs = [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];
    let reversed_arrow_pairs = [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ];
    let cases = [
        ("current", current_pairs.as_slice()),
        ("reversed-both", reversed_pairs.as_slice()),
        ("reversed-plain", reversed_plain_pairs.as_slice()),
        ("reversed-arrow", reversed_arrow_pairs.as_slice()),
    ];

    for (label, pairs) in cases {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(kp20_section42_eq12_rows_with_symbol_pairs(map, pairs));

        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected eq12 direction {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected eq12 direction {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected eq12 direction {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_with_eq13_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let current_eq12_rows = kp20_section42_only_eq12_rows(map);
    let reversed_eq12_rows = kp20_section42_eq12_rows_with_symbol_pairs(
        map,
        &[
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ],
    );
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);

    for (label, eq12_rows, include_eq13) in [
        ("current-eq12-no-eq13", current_eq12_rows.as_slice(), false),
        (
            "reversed-eq12-no-eq13",
            reversed_eq12_rows.as_slice(),
            false,
        ),
        ("current-eq12-eq13", current_eq12_rows.as_slice(), true),
        ("reversed-eq12-eq13", reversed_eq12_rows.as_slice(), true),
    ] {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(eq12_rows.iter().cloned());
        if include_eq13 {
            scalar_constraints.extend(eq13_rows.clone());
        }

        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected reversed-eq12 full {case_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            case_label = label,
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed-eq12 case {} should close algebraically, got {}",
            label,
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_eq13_bands_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let reversed_eq12_rows = kp20_section42_eq12_rows_with_symbol_pairs(
        map,
        &[
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ],
    );

    let cases = [
        (
            "current-bands",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
        ),
        (
            "swapped-bands",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
        (
            "all-plain-band",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
        ),
        (
            "all-underlined-band",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];

    for (label, strip_p, strip_q, row4_p, row4_q) in cases {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(reversed_eq12_rows.clone());
        scalar_constraints.extend(kp20_section42_eq13_rows_with_bands(
            map,
            2.0 * c,
            strip_p,
            strip_q,
            row4_p,
            row4_q,
        ));

        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected reversed-eq12 eq13 bands {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed-eq12 eq13 band case {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_split_rows_are_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);

    let splits = [
        (
            "reversed-plain-endpoint",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainEndpoint,
        ),
        (
            "reversed-plain-plus-row",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainPlusRow,
        ),
        (
            "reversed-plain-minus-row",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainMinusRow,
        ),
        (
            "reversed-arrow-endpoint",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowEndpoint,
        ),
        (
            "reversed-arrow-plus-row",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowPlusRow,
        ),
        (
            "reversed-arrow-minus-row",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowMinusRow,
        ),
    ];

    for (label, upper, lower, split) in splits {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(kp20_section42_eq12_split_rows_with_symbols(
            map, upper, lower, split,
        ));

        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected reversed eq12 split {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed eq12 split {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_corrected_section42_source_anchor_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);

    for anchor_extra in 0..source_extra_count {
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(anchor_extra), 1.0)],
            target: 1.0,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected §4.2 source anchor extra{anchor_extra}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected §4.2 source anchor extra{anchor_extra} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            validity.orientation.negative_determinants > 0,
            "corrected §4.2 source anchor extra{anchor_extra} unexpectedly passed coarse orientation"
        );
    }
}
