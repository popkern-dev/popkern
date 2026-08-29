use super::*;

#[test]
fn rotational_endpoint_aware_tensor_depth_policy_after_near_origin_fix_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let cases: [(&str, &[usize]); 7] = [
        ("depths-0-1-2", &[0, 1, 2]),
        ("depths-0-2-eq13-owns-1", &[0, 2]),
        ("depths-0-1", &[0, 1]),
        ("depths-1-2", &[1, 2]),
        ("depth-0-only", &[0]),
        ("depth-1-only", &[1]),
        ("depth-2-only", &[2]),
    ];
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, 0.0_f64);
    for (label, kept_depths) in cases {
        let tensor_rows = keep_tensor_chart_depths(&source_rows, source_variant, kept_depths);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("tensor depth policy after near-origin fix solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "tensor depth policy after near-origin fix {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            full_tensor_residual,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best tensor depth policy after near-origin fix: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "tensor depth policy after near-origin fix must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_row_space_seed_regularization_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const WEIGHTS: [f64; 9] = [
        0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for weight in WEIGHTS {
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: weight,
            })
            .expect("corrected row-space seed-regularization solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "corrected row-space seed regularization weight={weight:e}: residual={residual:e} f3={f3:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            f3 = kp20_planar_bi6_tau_f3_energy(&local_tau),
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            weight,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected row-space seed regularization: weight={:e} neg={} outside={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "seed-regularization ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_scale_after_near_origin_fix_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 7] = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        1.0_f64,
        1.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("source cross-depth scale solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            if validity.orientation.negative_determinants <= 18
                || validity.domain.outside_triangle_samples <= 5
            {
                eprintln!(
                    "source cross-depth scale left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                    residual = solved.max_constraint_residual,
                    full_tensor = full_tensor_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
            }
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best source cross-depth scale after near-origin fix: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "source cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_scale_window_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 8] = [0.125, 0.175, 0.225, 0.25, 0.275, 0.325, 0.375, 0.5];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        1.0_f64,
        1.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("source cross-depth scale window solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            if validity.orientation.negative_determinants <= 17
                || validity.domain.outside_triangle_samples <= 1
            {
                eprintln!(
                    "source cross-depth scale window left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                    residual = solved.max_constraint_residual,
                    full_tensor = full_tensor_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
            }
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best source cross-depth scale window: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "source cross-depth scale window ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_spectral_scale_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let scale_cases = [
        ("lambda", ring.lambda),
        ("one-minus-lambda", 1.0 - ring.lambda),
        ("half-one-minus-lambda", 0.5 * (1.0 - ring.lambda)),
        ("lambda-squared", ring.lambda * ring.lambda),
        ("quarter", 0.25),
        ("empirical-window-best", 0.225),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, scale) in scale_cases {
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &sector_rotated_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            scale,
        );
        let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            scale,
        );
        let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral source cross-depth scale solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "source cross-depth spectral scale {label}: scale={scale:e} lambda={lambda:e} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            lambda = ring.lambda,
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best source cross-depth spectral scale: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        (0.5 * (1.0 - ring.lambda) - 0.225).abs() < 1.0e-5,
        "Catmull-Clark spectral complement should explain the previous empirical cross-depth scale numerically"
    );
    assert!(
        best.3.is_finite(),
        "spectral source cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_spectral_source_depth1_projection_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let depth1_cases = [
        ("omit-depth1", None),
        ("same-spectral-depth1", Some(spectral_cross_depth_scale)),
        (
            "half-spectral-depth1",
            Some(0.5 * spectral_cross_depth_scale),
        ),
        (
            "sqrt-spectral-depth1",
            Some(spectral_cross_depth_scale.sqrt()),
        ),
        ("one-minus-lambda-depth1", Some(1.0 - ring.lambda)),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, depth1_scale) in depth1_cases {
        let mut scaled_rows = sector_rotated_rows.clone();
        if let Some(scale) = depth1_scale {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &scaled_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                1,
                scale,
            );
            scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                1,
                scale,
            );
        }
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &scaled_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            spectral_cross_depth_scale,
        );
        scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            spectral_cross_depth_scale,
        );
        let projected_tensor_rows = if depth1_scale.is_some() {
            scaled_rows.clone()
        } else {
            without_tensor_chart_depth(&scaled_rows, source_variant, 1)
        };
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral depth1 source projection solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "spectral source depth1 projection {label}: depth1_scale={depth1:e} depth2_scale={depth2:e} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            depth1 = depth1_scale.unwrap_or(f64::NAN),
            depth2 = spectral_cross_depth_scale,
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best spectral source depth1 projection: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite(),
        "spectral source depth1 projection ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_spectral_source_diagonal_symmetry_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");
    let cases = [
        ("fig8-only", true, false),
        ("diagonal-only", false, true),
        ("fig8-and-diagonal", true, true),
    ];
    for (label, include_fig8, include_diagonal) in cases {
        let mut vector_constraints = Vec::new();
        if include_fig8 {
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        }
        if include_diagonal {
            let mut scalar_constraints_with_diagonal = scalar_constraints.clone();
            scalar_constraints_with_diagonal.extend(diagonal_rows.clone());
            vector_constraints.extend(projected_tensor_rows.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints_with_diagonal,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("spectral source diagonal-symmetry solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "spectral source diagonal symmetry {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                full_tensor = full_tensor_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(
                solved.max_constraint_residual.is_finite(),
                "spectral source diagonal case {label} must be finite"
            );
            continue;
        }
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral source diagonal-symmetry solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "spectral source diagonal symmetry {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "spectral source diagonal case {label} must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_beta_source_under_corrected_chart_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let previous = valence - 1;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let beta_sources = [
        (
            "kp15-beta-split-then-apply",
            kp15_beta_split_then_apply_for_test(&tensor_border).expect("split beta source"),
        ),
        (
            "kp15-beta-full-then-split",
            kp15_beta_full_then_split_for_test(&tensor_border).expect("full beta source"),
        ),
    ];
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let previous_piece = Kp20PieceId::NearMidpointRight;
    let current_piece = Kp20PieceId::NearMidpointLeft;
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        "none",
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (source_label, split) in beta_sources {
        let previous_rows = split[previous].near_m.rows;
        let current_rows = split[sector].near_c.rows;
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            previous_rows,
            current_rows,
            source_variant,
        );
        let rotated_rows = rotate_tensor_chart_piece_targets(&base_rows, previous_piece, -angle);
        for (scale_label, source_rows) in [
            ("raw-depth", rotated_rows.clone()),
            ("spectral-depth2", {
                let previous_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                    &rotated_rows,
                    source_variant,
                    previous_piece,
                    2,
                    spectral_cross_depth_scale,
                );
                scale_tensor_chart_depth_from_boundary_for_piece(
                    &previous_scaled,
                    source_variant,
                    current_piece,
                    2,
                    spectral_cross_depth_scale,
                )
            }),
        ] {
            let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("corrected-chart beta source solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = source_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "corrected-chart beta source {source_label} {scale_label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                full_tensor = full_tensor_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            let row = (
                source_label,
                scale_label,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best corrected-chart beta source: source={} scale={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "corrected-chart beta source ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_row_space_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows.clone());
    let report = kp20_planar_tau_constraint_space_report(&vector_constraints, &scalar_constraints)
        .expect("corrected spectral row-space report");
    eprintln!(
        "corrected spectral source row-space: vars={} dense_rows={} independent_rows={} nullity={}",
        report.scalar_variable_count,
        report.dense_row_count,
        report.independent_row_count,
        report.nullity,
    );

    let mut vector_with_internal_c2 = kp20_planar_tau_internal_c2_constraints();
    vector_with_internal_c2.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_with_internal_c2.extend(projected_tensor_rows);
    let internal_c2_report =
        kp20_planar_tau_constraint_space_report(&vector_with_internal_c2, &scalar_constraints)
            .expect("corrected spectral internal-C2 row-space report");
    eprintln!(
        "corrected spectral source plus internal-C2 row-space: vars={} dense_rows={} independent_rows={} nullity={}",
        internal_c2_report.scalar_variable_count,
        internal_c2_report.dense_row_count,
        internal_c2_report.independent_row_count,
        internal_c2_report.nullity,
    );
    assert_eq!(report.scalar_variable_count, 392);
    assert!(
        report.nullity > 15,
        "corrected spectral source rows should still expose the missing sigma finite-basis reduction"
    );
    assert!(
        internal_c2_report.nullity <= report.nullity,
        "adding internal C2 should not increase the corrected spectral row-space nullity"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_fig8_average_variants_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (label, include_diagonal_mirror, include_rotated_quadrants) in [
        ("baseline", false, false),
        ("diagonal-mirror", true, false),
        ("rotated-quadrants", false, true),
        ("diagonal-mirror+rotated-quadrants", true, true),
    ] {
        let fig8_rows = test_fig8_indicated_average_variant_rows(
            include_diagonal_mirror,
            include_rotated_quadrants,
        );
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(fig8_rows);
        vector_constraints.extend(projected_tensor_rows.clone());
        let report =
            kp20_planar_tau_constraint_space_report(&vector_constraints, &scalar_constraints)
                .expect("Fig. 8 average variant row-space report");
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("Fig. 8 average variant solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "corrected spectral Fig8 average variant {label}: rows={rows} independent_rows={independent} nullity={nullity} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            rows = report.dense_row_count,
            independent = report.independent_row_count,
            nullity = report.nullity,
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            report.nullity,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
            || (row.1 == best.1 && row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected spectral Fig8 average variant: label={} neg={} outside={} nullity={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.4.is_finite(),
        "Fig. 8 average variant ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_internal_c2_solve_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    for (label, include_internal_c2) in [("without-internal-c2", false), ("with-internal-c2", true)]
    {
        let mut vector_constraints = Vec::new();
        if include_internal_c2 {
            vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        }
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("corrected spectral internal-C2 solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "corrected spectral source {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "corrected spectral internal-C2 case {label} must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_internal_c2_source_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let (projected_tensor_rows, projection_report) = kp20_project_rotational_paper_tensor_band(
        projected_source_input.clone(),
        &fixed_vector_constraints,
        &scalar_constraints,
    )
    .expect("corrected spectral internal-C2 source projection succeeds");
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_tensor_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("projected corrected spectral internal-C2 solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "corrected spectral internal-C2 source projection: projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        source_delta = projection_report.max_source_delta,
        adjusted = projection_report.adjusted_scalar_targets,
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    let mut deltas = projected_source_input
        .iter()
        .zip(&projected_tensor_rows)
        .map(|(source, projected)| {
            let coeff = source.terms[0].0;
            let dx = projected.target[0] - source.target[0];
            let dy = projected.target[1] - source.target[1];
            (dx.hypot(dy), coeff, dx, dy, source.target, projected.target)
        })
        .collect::<Vec<_>>();
    deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (rank, (norm, coeff, dx, dy, source, projected)) in deltas.iter().take(8).enumerate() {
        eprintln!(
            "  corrected spectral projection delta rank={rank} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} source=({sx:e},{sy:e}) projected=({px:e},{py:e})",
            piece = coeff.piece,
            i = coeff.i,
            j = coeff.j,
            sx = source[0],
            sy = source[1],
            px = projected[0],
            py = projected[1],
        );
    }
    assert!(
        solved.max_constraint_residual.is_finite()
            && projection_report.max_row_residual.is_finite(),
        "corrected spectral internal-C2 source projection must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_depth2_only_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let free_depth2_coeffs = projected_source_input
        .iter()
        .map(|row| row.terms[0].0)
        .filter(|&coeff| tensor_chart_coeff_depth(coeff, source_variant) == 2)
        .collect::<Vec<_>>();
    let (projected_tensor_rows, projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            projected_source_input.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &free_depth2_coeffs,
        )
        .expect("depth2-only corrected spectral projection succeeds");
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_tensor_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("depth2-only corrected spectral solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "corrected spectral internal-C2 depth2-only projection: projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        source_delta = projection_report.max_source_delta,
        adjusted = projection_report.adjusted_scalar_targets,
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite()
            && projection_report.max_row_residual.is_finite(),
        "depth2-only corrected spectral projection must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_depth2_piece_projection_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let mut best = ("none", f64::INFINITY, usize::MAX, usize::MAX, f64::INFINITY);
    for (label, free_piece) in [
        ("underlined-p-depth2", Kp20PieceId::NearMidpointLeft),
        ("arrow-underlined-p-depth2", Kp20PieceId::NearMidpointRight),
    ] {
        let free_coeffs = projected_source_input
            .iter()
            .map(|row| row.terms[0].0)
            .filter(|&coeff| {
                coeff.piece == free_piece && tensor_chart_coeff_depth(coeff, source_variant) == 2
            })
            .collect::<Vec<_>>();
        let (projected_tensor_rows, projection_report) =
            kp20_project_section42_source_tensor_rows_with_free_coeffs(
                projected_source_input.clone(),
                &fixed_vector_constraints,
                &scalar_constraints,
                &free_coeffs,
            )
            .expect("piece depth2 corrected spectral projection succeeds");
        let mut vector_constraints = fixed_vector_constraints.clone();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("piece depth2 corrected spectral solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "corrected spectral internal-C2 piece depth2 projection {label}: free={} projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            free_coeffs.len(),
            projection_residual = projection_report.max_row_residual,
            source_delta = projection_report.max_source_delta,
            adjusted = projection_report.adjusted_scalar_targets,
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            label,
            solved.max_constraint_residual,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            projection_report.max_source_delta,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected spectral internal-C2 piece depth2 projection: label={} residual={:e} neg={} outside={} source_delta={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(best.1.is_finite(), "piece depth2 ranking must be finite");
}

#[test]
fn rotational_endpoint_aware_pr08_sources_without_projection_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let mut scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
            .expect("sector-local diagonal rows"),
    );

    let mut best = (
        String::from("none"),
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    let mut valid_candidates = 0_usize;
    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
        let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
        let sector_rotated_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &sector_rotated_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            spectral_cross_depth_scale,
        );
        let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            spectral_cross_depth_scale,
        );

        for (depth_label, tensor_rows) in [
            (
                "depth0+2",
                without_tensor_chart_depth(&scaled_rows, source_variant, 1),
            ),
            ("depth0+1+2", scaled_rows.clone()),
        ] {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(tensor_rows.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("unprojected PR08 source solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let is_valid = validity.is_strict_sector_bounded_reparameterization();
            if is_valid {
                valid_candidates += 1;
            }
            eprintln!(
                "unprojected PR08 source {candidate:?} {depth_label}: valid={is_valid} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                full_tensor = full_tensor_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            let row = (
                format!("{candidate:?}/{depth_label}"),
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.1 < best.1
                || (row.1 == best.1 && row.2 < best.2)
                || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
                || (row.1 == best.1 && row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best unprojected PR08 source: label={} valid_candidates={} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, valid_candidates, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "unprojected PR08 source ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_promising_pr08_charts_are_ranked_with_live_and_diagnostic_diagonal() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let base_scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let cases = [
        (
            "best-validity-before-diagonal",
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: true,
                previous_reverse_u: true,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: true,
            },
        ),
        (
            "best-domain-before-diagonal",
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: true,
            },
        ),
        (
            "scaled-top-domain-before-diagonal",
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: false,
            },
        ),
    ];

    for include_diagnostic_diagonal in [false, true] {
        let mut scalar_constraints = base_scalar_constraints.clone();
        if include_diagnostic_diagonal {
            scalar_constraints.extend(
                kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
                    .expect("sector-local diagonal rows"),
            );
        }
        let mode_label = if include_diagnostic_diagonal {
            "diagnostic-diagonal"
        } else {
            "live"
        };
        let mut best = (
            "none",
            usize::MAX,
            usize::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for (label, candidate, source_variant) in cases {
            let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
                .expect("PR08 candidate rows");
            let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
                &raw_rows, valence, sector, true, true,
            )
            .expect("candidate endpoint normalization");
            let tensor_border =
                kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
            let base_rows =
                tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
            let sector_rotated_rows = rotate_tensor_chart_piece_targets(
                &base_rows,
                Kp20PieceId::NearMidpointRight,
                -angle,
            );
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                spectral_cross_depth_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                spectral_cross_depth_scale,
            );
            for (depth_label, tensor_rows) in [
                (
                    "depth0+2",
                    without_tensor_chart_depth(&scaled_rows, source_variant, 1),
                ),
                ("depth0+1+2", scaled_rows.clone()),
            ] {
                let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
                vector_constraints.extend(tensor_rows);
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: planar_tau_seed_macro_patch(valence, sector),
                        vector_constraints,
                        scalar_constraints: scalar_constraints.clone(),
                        regularization_weight: 1.0e-10,
                    },
                )
                .expect("promising PR08 chart solve succeeds");
                let local_tau =
                    global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                let validity =
                    kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
                eprintln!(
                    "promising PR08 chart mode={mode_label} {label} candidate={} chart={} {depth_label}: valid={} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    candidate.name(),
                    source_variant.label(),
                    validity.is_strict_sector_bounded_reparameterization(),
                    residual = solved.max_constraint_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    zero = validity.orientation.near_zero_determinants,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
                let row = (
                    label,
                    validity.orientation.negative_determinants,
                    validity.domain.outside_triangle_samples,
                    solved.max_constraint_residual,
                    validity.orientation.min_determinant,
                );
                if row.1 < best.1
                    || (row.1 == best.1 && row.2 < best.2)
                    || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
                {
                    best = row;
                }
            }
        }
        eprintln!(
            "best promising PR08 chart mode={mode_label}: label={} neg={} outside={} residual={:e} min_det={:e}",
            best.0, best.1, best.2, best.3, best.4,
        );
        assert!(
            best.3.is_finite(),
            "promising PR08 chart ranking must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_production_pr08_inner_q1_source_matches_live_domain_baseline() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(without_tensor_chart_depth(&scaled_rows, source_variant, 1));
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(
                TestSection42ChartMap {
                    p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearOriginLeft,
                        reverse_i: true,
                        reverse_j: false,
                        transpose: false,
                    },
                    arrow_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearOriginRight,
                        reverse_i: false,
                        reverse_j: false,
                        transpose: false,
                    },
                    underlined_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearMidpointLeft,
                        reverse_i: false,
                        reverse_j: false,
                        transpose: false,
                    },
                    arrow_underlined_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearMidpointRight,
                        reverse_i: true,
                        reverse_j: false,
                        transpose: false,
                    },
                },
                2.0 * c,
            ),
            regularization_weight: 1.0e-10,
        })
        .expect("production PR08 inner-Q1 source solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "production PR08 inner-Q1 live chart: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );

    assert_eq!(validity.domain.outside_triangle_samples, 0);
    assert_eq!(validity.orientation.negative_determinants, 37);
    assert!(solved.max_constraint_residual < 7.0e-2);
    assert!(full_tensor_residual < 1.0e-1);
}

#[test]
fn rotational_c2_edge_chart_production_pr08_inner_q1_source_is_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.125,
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(source_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge production PR08 inner-Q1 source solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    let piece_reports = test_piece_domain_orientation_reports(&local_solved, &SAMPLES, 1.0e-10);
    eprintln!(
        "C2-edge production PR08 inner-Q1 source: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    for report in piece_reports {
        eprintln!(
            "C2-edge production piece {:?}: domain outside={}/{} orientation neg={}/{} zero={} det=[{:.6},{:.6}] worst_det={}",
            report.piece,
            report.outside,
            report.domain_samples,
            report.negative,
            report.orientation_samples,
            report.near_zero,
            report.min_det,
            report.max_det,
            report.worst_orientation_sample,
        );
    }

    assert_eq!(validity.domain.outside_triangle_samples, 0);
    assert_eq!(validity.orientation.negative_determinants, 20);
    assert!(solved.max_constraint_residual < 5.5e-2);
    assert!(full_tensor_residual < 9.0e-2);
}

#[test]
fn rotational_c2_edge_eq14_reparameterized_pr08_source_factor_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const FACTORS: [f64; 13] = [
        -1.0, -0.5, -0.25, -0.2, -0.125, 0.0, 0.125, 0.2, 0.25, 0.5, 0.75, 1.0, 1.25,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);

    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for second_order_factor in FACTORS {
        let (previous_rows, current_rows) = split_regular_sector_rows_with_second_order_factor(
            &tensor_border,
            sector,
            source_variant.previous_uses_near_m,
            source_variant.current_uses_near_m,
            second_order_factor,
        );
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            previous_rows,
            current_rows,
            source_variant,
        );
        let source_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let projected_source_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(projected_source_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge eq14 reparameterized PR08 source solve succeeds");
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "C2-edge eq14 reparameterized PR08 source factor={second_order_factor:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            second_order_factor,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge eq14 reparameterized PR08 source factor={:e}: neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge eq14 reparameterized source ranking must remain finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_cross_depth_scale_window_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 12] = [
        -0.25,
        -0.125,
        0.0,
        0.125,
        0.175,
        0.2250058227408514,
        0.275,
        0.325,
        0.375,
        0.5,
        0.75,
        1.0,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );

    let mut best = (
        0.0_f64,
        0.0_f64,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("live cross-depth scale solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.4 < best.4)
                || (row.2 == best.2 && row.4 == best.4 && row.5 < best.5)
            {
                best = row;
            }
            if validity.orientation.negative_determinants <= 30
                || validity.domain.outside_triangle_samples == 0
            {
                eprintln!(
                    "live cross-depth scale left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    residual = solved.max_constraint_residual,
                    full_tensor = full_tensor_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    zero = validity.orientation.near_zero_determinants,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
            }
        }
    }
    eprintln!(
        "best live cross-depth scale: left={:e} right={:e} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7,
    );
    assert!(
        best.5.is_finite(),
        "live cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_cross_depth_signed_refinement_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const LEFT_SCALES: [f64; 11] = [
        -0.35, -0.3, -0.25, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, 0.0,
    ];
    const RIGHT_SCALES: [f64; 9] = [
        0.125,
        0.175,
        0.2250058227408514,
        0.275,
        0.325,
        0.375,
        0.45,
        0.6,
        0.75,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        0.0_f64,
        0.0_f64,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in LEFT_SCALES {
        for right_scale in RIGHT_SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("live signed cross-depth refinement solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.4 < best.4)
                || (row.2 == best.2 && row.4 == best.4 && row.5 < best.5)
            {
                best = row;
            }
            if validity.orientation.negative_determinants <= 26 {
                eprintln!(
                    "live signed cross-depth left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    residual = solved.max_constraint_residual,
                    full_tensor = full_tensor_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    zero = validity.orientation.near_zero_determinants,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
            }
        }
    }
    eprintln!(
        "best live signed cross-depth refinement: left={:e} right={:e} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7,
    );
    assert!(
        best.5.is_finite(),
        "live signed cross-depth refinement must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_cross_depth_best_is_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(
                TestSection42ChartMap {
                    p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearOriginLeft,
                        reverse_i: true,
                        reverse_j: false,
                        transpose: false,
                    },
                    arrow_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearOriginRight,
                        reverse_i: false,
                        reverse_j: false,
                        transpose: false,
                    },
                    underlined_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearMidpointLeft,
                        reverse_i: false,
                        reverse_j: false,
                        transpose: false,
                    },
                    arrow_underlined_p: TestSection42ChartStrip {
                        piece: Kp20PieceId::NearMidpointRight,
                        reverse_i: true,
                        reverse_j: false,
                        transpose: false,
                    },
                },
                2.0 * c,
            ),
            regularization_weight: 1.0e-10,
        })
        .expect("live signed cross-depth best localization solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "live signed cross-depth best: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    for report in piece_reports {
        eprintln!(
            "piece {:?}: domain outside={}/{} u=[{:.6},{:.6}] v=[{:.6},{:.6}] max_excess={:.6} worst_domain={} orientation neg={}/{} zero={} det=[{:.6},{:.6}] worst_det={}",
            report.piece,
            report.outside,
            report.domain_samples,
            report.min_u,
            report.max_u,
            report.min_v,
            report.max_v,
            report.max_excess,
            report.worst_domain_sample,
            report.negative,
            report.orientation_samples,
            report.near_zero,
            report.min_det,
            report.max_det,
            report.worst_orientation_sample,
        );
    }
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "live signed cross-depth best localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_source_with_first_order_frame_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let axes = sector_local_axis_components(valence, sector);
    for (label, include_right) in [("left-origin-frame", false), ("both-origin-frames", true)] {
        let mut augmented_scalar_constraints =
            first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, axes, 0);
        if include_right {
            augmented_scalar_constraints.extend(first_order_centered_source_rows(
                Kp20PieceId::NearOriginRight,
                axes,
                0,
            ));
        }
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints: vector_constraints.clone(),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("live signed source first-order frame solve succeeds");
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        eprintln!(
            "live signed source first-order {label}: residual={residual:e} scale={scale:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "first-order frame solve must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_live_signed_source_first_order_frames_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let axes = sector_local_axis_components(valence, sector);
    let frame_cases = [
        ("identity", [1.0, 0.0], [0.0, 1.0]),
        ("both-negated", [-1.0, 0.0], [0.0, -1.0]),
        ("swap", [0.0, 1.0], [1.0, 0.0]),
        ("swap-negated", [0.0, -1.0], [-1.0, 0.0]),
        ("fig8-left", [0.5, 0.5], [0.0, 1.0]),
        ("fig8-right-positive", [1.0, 0.0], [0.5, 0.5]),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, du, dv) in frame_cases {
        let augmented_scalar_constraints = first_order_source_rows_for_local_frame(
            Kp20PieceId::NearOriginLeft,
            axes,
            [0.0, 0.0],
            du,
            dv,
            0,
        );
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints: vector_constraints.clone(),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("live signed source first-order frame variant solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        eprintln!(
            "live signed source first-order frame {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] left_det={left_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            solved.extra_values[0],
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best live signed source first-order frame: label={} neg={} outside={} residual={:e} min_det={:e} scale={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite(),
        "first-order frame variant ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_source_centered_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let axes = sector_local_axis_components(valence, sector);
    for (label, du, dv) in [
        ("identity", [1.0, 0.0], [0.0, 1.0]),
        ("fig8-left", [0.5, 0.5], [0.0, 1.0]),
    ] {
        let centered_rows = first_order_source_rows_for_local_frame(
            Kp20PieceId::NearOriginLeft,
            axes,
            [0.0, 0.0],
            du,
            dv,
            0,
        );
        let projected_rows = project_tensor_band_with_centered_first_order(
            source_rows.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &centered_rows,
        );
        let max_source_delta = source_rows
            .iter()
            .zip(&projected_rows)
            .map(|(source, projected)| {
                (projected.target[0] - source.target[0])
                    .hypot(projected.target[1] - source.target[1])
            })
            .fold(0.0_f64, f64::max);
        let mut deltas = source_rows
            .iter()
            .zip(&projected_rows)
            .map(|(source, projected)| {
                let coeff = source.terms[0].0;
                let dx = projected.target[0] - source.target[0];
                let dy = projected.target[1] - source.target[1];
                (dx.hypot(dy), coeff, dx, dy, source.target, projected.target)
            })
            .collect::<Vec<_>>();
        deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let mut vector_constraints = fixed_vector_constraints.clone();
        vector_constraints.extend(projected_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("centered projected source solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "live signed centered source projection {label}: residual={residual:e} scale={scale:e} max_source_delta={max_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            max_delta = max_source_delta,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        for (rank, (norm, coeff, dx, dy, source, projected)) in deltas.iter().take(8).enumerate() {
            eprintln!(
                "  centered projection {label} delta rank={rank} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} source=({sx:e},{sy:e}) projected=({px:e},{py:e})",
                piece = coeff.piece,
                i = coeff.i,
                j = coeff.j,
                sx = source[0],
                sy = source[1],
                px = projected[0],
                py = projected[1],
            );
        }
        assert!(
            solved.max_constraint_residual.is_finite() && max_source_delta.is_finite(),
            "centered source projection must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_live_signed_source_appendix_a_refined_pins_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let mut best = (
        "none",
        0_usize,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, tau_family) in [
        ("paper-eq22", Kp20TauFamily::Paper),
        ("appendix-a-refined", Kp20TauFamily::PaperRefinedApex),
    ] {
        let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
            &local_identity_tri6_guide(),
            valence,
            sector,
            tau_family,
            kp20_paper_quadrant_sector_piece_maps(),
        )
        .expect("Appendix A tau samples with identity guide");
        let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
        for (source_label, include_source_rows) in
            [("without-source", false), ("with-source", true)]
        {
            for max_total_degree in [1_usize, 2, 3] {
                let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
                if include_source_rows {
                    vector_constraints.extend(source_rows.clone());
                }
                for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
                    for i in 0..=max_total_degree {
                        for j in 0..=(max_total_degree - i) {
                            let target = paper_tau.piece(piece_id).cps[i][j];
                            vector_constraints.push(Kp20PlanarTauLinearConstraint::pin(
                                Kp20PlanarTauCoeff::new(piece_id, i, j),
                                [target.x, target.y],
                            ));
                        }
                    }
                }
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: planar_tau_seed_macro_patch(valence, sector),
                        vector_constraints,
                        scalar_constraints: scalar_constraints.clone(),
                        regularization_weight: 1.0e-10,
                    },
                )
                .expect("live signed source Appendix A central pin solve succeeds");
                let full_tensor_residual = scaled_rows
                    .iter()
                    .map(|row| {
                        let residual = vector_residual(&solved.tau_macro_patch, row);
                        residual[0].abs().max(residual[1].abs())
                    })
                    .fold(0.0_f64, f64::max);
                let local_solved =
                    global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                let validity =
                    kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
                eprintln!(
                    "live signed Appendix-A pins {label} {source_label} degree={degree}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    degree = max_total_degree,
                    residual = solved.max_constraint_residual,
                    full_tensor = full_tensor_residual,
                    outside = validity.domain.outside_triangle_samples,
                    domain_samples = validity.domain.sample_count,
                    neg = validity.orientation.negative_determinants,
                    orientation_samples = validity.orientation.sample_count,
                    zero = validity.orientation.near_zero_determinants,
                    min_det = validity.orientation.min_determinant,
                    max_det = validity.orientation.max_determinant,
                );
                let row = (
                    label,
                    max_total_degree,
                    validity.orientation.negative_determinants,
                    validity.domain.outside_triangle_samples,
                    solved.max_constraint_residual,
                    validity.orientation.min_determinant,
                    full_tensor_residual,
                );
                if include_source_rows
                    && (row.2 < best.2
                        || (row.2 == best.2 && row.3 < best.3)
                        || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4))
                {
                    best = row;
                }
            }
        }
    }
    eprintln!(
        "best live signed Appendix-A pins: label={} degree={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite() && best.6.is_finite(),
        "Appendix A central pin ranking must remain finite"
    );
}

#[test]
fn rotational_appendix_a_refined_central_pins_row_stack_components_are_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A refined tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mut central_pin_rows = Vec::new();
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        for i in 0..=2 {
            for j in 0..=(2 - i) {
                let target = paper_tau.piece(piece_id).cps[i][j];
                central_pin_rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece_id, i, j),
                    [target.x, target.y],
                ));
            }
        }
    }
    let scalar_rows = test_section42_endpoint_aware_chart_rows(
        TestSection42ChartMap {
            p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginLeft,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
            arrow_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearOriginRight,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointLeft,
                reverse_i: false,
                reverse_j: false,
                transpose: false,
            },
            arrow_underlined_p: TestSection42ChartStrip {
                piece: Kp20PieceId::NearMidpointRight,
                reverse_i: true,
                reverse_j: false,
                transpose: false,
            },
        },
        2.0 * c,
    );
    let cases = [
        ("pins-only", false, false, false),
        ("pins-internal-c2", true, false, false),
        ("pins-fig8", false, true, false),
        ("pins-section42-scalar", false, false, true),
        ("pins-internal-c2-fig8", true, true, false),
        ("pins-internal-c2-section42", true, false, true),
        ("pins-fig8-section42", false, true, true),
        ("pins-full-live-no-source", true, true, true),
    ];
    for (label, include_internal_c2, include_fig8, include_scalar) in cases {
        let mut vector_constraints = central_pin_rows.clone();
        if include_internal_c2 {
            vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        }
        if include_fig8 {
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        }
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: if include_scalar {
                    scalar_rows.clone()
                } else {
                    Vec::new()
                },
                regularization_weight: 1.0e-10,
            })
            .expect("Appendix A refined row component solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "Appendix-A refined central pins row component {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "Appendix A refined row component residual must be finite"
        );
    }

    let map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: true,
            reverse_j: false,
            transpose: false,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: false,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: false,
            transpose: false,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: true,
            reverse_j: false,
            transpose: false,
        },
    };
    let mut grouped_rows = Vec::new();
    let mut c0_rows = Vec::new();
    test_add_section42_chart_c0_rows(&mut c0_rows, map);
    grouped_rows.push(("section42-c0", c0_rows));
    let mut eq11_rows = Vec::new();
    test_add_section42_chart_eq11_rows(&mut eq11_rows, map);
    grouped_rows.push(("eq11-sector-curve", eq11_rows));
    let mut eq12_rows = Vec::new();
    test_add_section42_chart_endpoint_aware_eq12_rows(&mut eq12_rows, map);
    grouped_rows.push(("eq12-strip-join", eq12_rows));
    let mut eq13_rows = Vec::new();
    test_add_section42_chart_eq13_rows(&mut eq13_rows, map, 2.0 * c);
    grouped_rows.push(("eq13-lower-strip", eq13_rows));
    for (label, rows) in grouped_rows {
        let mut vector_constraints = central_pin_rows.clone();
        vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: rows,
                regularization_weight: 1.0e-10,
            })
            .expect("Appendix A refined scalar group solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "Appendix-A refined central pins with internal C2 scalar group {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "Appendix A refined scalar group residual must be finite"
        );
    }

    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let c2_edge_scalar_rows = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut vector_constraints = central_pin_rows;
    vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: c2_edge_scalar_rows,
            regularization_weight: 1.0e-10,
        })
        .expect("Appendix A refined C2-edge chart solve succeeds");
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "Appendix-A refined central pins full rows c2-edge chart: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        residual = solved.max_constraint_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite(),
        "Appendix A refined C2-edge chart residual must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_source_c2_edge_section42_chart_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .expect("PR08 candidate rows");
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
        &raw_rows, valence, sector, true, true,
    )
    .expect("candidate endpoint normalization");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let sector_rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(source_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("live source C2-edge chart solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "live source C2-edge section42 chart: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "live source C2-edge chart solve must remain finite"
    );
}
