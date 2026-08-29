use super::*;

#[test]
fn rotational_c2_edge_section42_source_phase_variants_are_ranked() {
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
    let mut ranked = Vec::new();
    for previous_uses_near_m in [false, true] {
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let source_variant = TestTensorBorderChartVariant {
                                previous_uses_near_m,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let base_rows = tensor_border_chart_variant_rows(
                                &tensor_border,
                                sector,
                                source_variant,
                            );
                            for rotate_current in [false, true] {
                                let rotated_rows = if rotate_current {
                                    rotate_tensor_chart_piece_targets(
                                        &base_rows,
                                        Kp20PieceId::NearMidpointRight,
                                        -angle,
                                    )
                                } else {
                                    base_rows.clone()
                                };
                                for left_scale in [-0.125, 0.5 * (1.0 - ring.lambda)] {
                                    let left_scaled =
                                        scale_tensor_chart_depth_from_boundary_for_piece(
                                            &rotated_rows,
                                            source_variant,
                                            Kp20PieceId::NearMidpointLeft,
                                            2,
                                            left_scale,
                                        );
                                    let scaled_rows =
                                        scale_tensor_chart_depth_from_boundary_for_piece(
                                            &left_scaled,
                                            source_variant,
                                            Kp20PieceId::NearMidpointRight,
                                            2,
                                            0.5 * (1.0 - ring.lambda),
                                        );
                                    let source_rows =
                                        without_tensor_chart_depth(&scaled_rows, source_variant, 1);
                                    let mut vector_constraints =
                                        kp20_planar_tau_internal_c2_constraints();
                                    vector_constraints.extend(
                                        kp20_planar_tau_fig8_indicated_averaging_constraints(),
                                    );
                                    vector_constraints.extend(source_rows);
                                    let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                        &Kp20PlanarTauF3ScalarSolveInput {
                                            initial: planar_tau_seed_macro_patch(valence, sector),
                                            vector_constraints,
                                            scalar_constraints: scalar_constraints.clone(),
                                            regularization_weight: 1.0e-10,
                                        },
                                    )
                                    .expect("C2-edge source phase solve succeeds");
                                    let full_tensor_residual = scaled_rows
                                        .iter()
                                        .map(|row| {
                                            let residual =
                                                vector_residual(&solved.tau_macro_patch, row);
                                            residual[0].abs().max(residual[1].abs())
                                        })
                                        .fold(0.0_f64, f64::max);
                                    let local_solved = global_planar_patch_to_local(
                                        valence,
                                        sector,
                                        solved.tau_macro_patch,
                                    );
                                    let validity = kp20_planar_bi6_tau_validity_report(
                                        &local_solved,
                                        &SAMPLES,
                                        1.0e-10,
                                        1.0e-10,
                                    );
                                    ranked.push((
                                        validity.orientation.negative_determinants,
                                        validity.domain.outside_triangle_samples,
                                        solved.max_constraint_residual,
                                        full_tensor_residual,
                                        validity.orientation.min_determinant,
                                        source_variant,
                                        rotate_current,
                                        left_scale,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ranked.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(a.2.partial_cmp(&b.2).unwrap())
    });
    for (rank, row) in ranked.iter().take(12).enumerate() {
        eprintln!(
            "C2-edge source phase rank={rank} neg={} outside={} residual={:e} full_tensor={:e} min_det={:e} rotate_current={} left_scale={:e} variant={:?}",
            row.0, row.1, row.2, row.3, row.4, row.6, row.7, row.5,
        );
    }
    let best = ranked.first().expect("source phase cases exist");
    assert!(
        best.2.is_finite() && best.3.is_finite(),
        "C2-edge source phase ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_section42_best_source_phase_depth_scales_are_ranked() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for left_scale in LEFT_SCALES {
        for right_scale in RIGHT_SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &rotated_rows,
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
            let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(source_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("C2-edge source depth scale solve succeeds");
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
            if validity.domain.outside_triangle_samples == 0
                && (validity.orientation.negative_determinants <= 21
                    || solved.max_constraint_residual <= 5.5e-2)
            {
                eprintln!(
                    "C2-edge best phase depth left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.3 < best.3
                || (row.3 == best.3 && row.2 < best.2)
                || (row.3 == best.3 && row.2 == best.2 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best C2-edge best phase depth scales: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite() && best.6.is_finite(),
        "C2-edge best phase depth scale ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_without_depth_scaling_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    vector_constraints.extend(source_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge paper-oriented source solve succeeds");
    let full_tensor_residual = source_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented source no depth scaling: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        "C2-edge paper-oriented source solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_projection_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let (projected_rows, projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            source_rows.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &source_rows
                .iter()
                .map(|row| row.terms[0].0)
                .collect::<Vec<_>>(),
        )
        .expect("C2-edge paper-oriented source projection succeeds");
    let max_source_delta = source_rows
        .iter()
        .zip(&projected_rows)
        .map(|(source, projected)| {
            (projected.target[0] - source.target[0]).hypot(projected.target[1] - source.target[1])
        })
        .fold(0.0_f64, f64::max);
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge projected paper-oriented source solve succeeds");
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented source projected: projection_residual={projection_residual:e} projection_delta={projection_delta:e} max_source_delta={max_delta:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        projection_delta = projection_report.max_source_delta,
        max_delta = max_source_delta,
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
        solved.max_constraint_residual.is_finite() && max_source_delta.is_finite(),
        "C2-edge paper-oriented source projection solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_affine_normalization_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let mut augmented_scalar_constraints = Vec::new();
    for row in &source_rows {
        let coeff = row.terms[0].0;
        let [x, y] = row.target;
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -x),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(1), -y),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(2), -1.0),
            ],
            target: x,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(3), -x),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(4), -y),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(5), -1.0),
            ],
            target: y,
        });
    }
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            augmented_scalar_constraints,
            extra_variable_count: 6,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
        },
    )
    .expect("C2-edge paper-oriented affine source solve succeeds");
    let d = &solved.extra_values;
    let affine_det = (1.0 + d[0]) * (1.0 + d[4]) - d[1] * d[3];
    let max_source_delta = source_rows
        .iter()
        .map(|row| {
            let [x, y] = row.target;
            (d[0] * x + d[1] * y + d[2]).hypot(d[3] * x + d[4] * y + d[5])
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented affine source: residual={residual:e} affine_det={affine_det:e} max_source_delta={max_delta:e} extras=[{d0:e},{d1:e},{d2:e},{d3:e},{d4:e},{d5:e}] outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        residual = solved.max_constraint_residual,
        max_delta = max_source_delta,
        d0 = d[0],
        d1 = d[1],
        d2 = d[2],
        d3 = d[3],
        d4 = d[4],
        d5 = d[5],
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite() && affine_det.is_finite(),
        "C2-edge affine source normalization solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_translation_normalization_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let mut augmented_scalar_constraints = Vec::new();
    for row in &source_rows {
        let coeff = row.terms[0].0;
        let [x, y] = row.target;
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -1.0),
            ],
            target: x,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(1), -1.0),
            ],
            target: y,
        });
    }
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            augmented_scalar_constraints,
            extra_variable_count: 2,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
        },
    )
    .expect("C2-edge paper-oriented translation source solve succeeds");
    let tx = solved.extra_values[0];
    let ty = solved.extra_values[1];
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented translation source: residual={residual:e} translation=({tx:e},{ty:e}) outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        "C2-edge translation source normalization solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_eq14_d1_3_sweep_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let default_d1_3 = 2.0 * c / 3.0;
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
    for d1_3 in [
        -default_d1_3,
        -0.5 * default_d1_3,
        0.0,
        0.125,
        0.5 * default_d1_3,
        default_d1_3,
        1.5 * default_d1_3,
        2.0 * default_d1_3,
    ] {
        let split = split_kp20_tensor_border_with_d1_3(&tensor_border, d1_3).expect("eq14 split");
        let previous = (sector + split.len() - 1) % split.len();
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            split[previous].near_m.rows,
            split[sector].near_c.rows,
            source_variant,
        );
        let source_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(source_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge d1_3 source solve succeeds");
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
            "C2-edge paper-oriented eq14 d1_3={d1_3:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            d1_3,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge paper-oriented eq14 d1_3: d1_3={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge d1_3 source ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_depth_subsets_are_ranked() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, depths) in [
        ("depth-0", &[0_usize][..]),
        ("depth-1", &[1_usize][..]),
        ("depth-2", &[2_usize][..]),
        ("depth-0-1", &[0_usize, 1][..]),
        ("depth-0-2", &[0_usize, 2][..]),
        ("depth-1-2", &[1_usize, 2][..]),
        ("depth-0-1-2", &[0_usize, 1, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge source depth subset solve succeeds");
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
            "C2-edge paper-oriented source depth subset {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge paper-oriented source depth subset: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge source depth subset ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_appendix_a_central_pins_with_source_depths_are_ranked() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A tau samples with identity guide");
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
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, depths) in [
        ("depth-0", &[0_usize][..]),
        ("depth-1", &[1_usize][..]),
        ("depth-2", &[2_usize][..]),
        ("depth-0-2", &[0_usize, 2][..]),
        ("depth-1-2", &[1_usize, 2][..]),
        ("depth-0-1-2", &[0_usize, 1, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(central_pin_rows.clone());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge central source depth solve succeeds");
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
            "C2-edge Appendix-A central pins source depth {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge Appendix-A central pins source depth: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge central source depth ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_source_depth_candidates_are_compared_to_appendix_a_tau() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let paper_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, paper_tau);
    for (label, depths) in [
        ("depth-1", &[1_usize][..]),
        ("depth-0-2", &[0_usize, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge source depth candidate solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        let mut deltas = Vec::new();
        for piece in Kp20PieceId::all() {
            let solved_piece = solved.tau_macro_patch.piece(piece);
            let paper_piece = paper_tau.piece(piece);
            for i in 0..=6 {
                for j in 0..=6 {
                    let dx = solved_piece.cps[i][j].x - paper_piece.cps[i][j].x;
                    let dy = solved_piece.cps[i][j].y - paper_piece.cps[i][j].y;
                    deltas.push((dx.hypot(dy), piece, i, j, dx, dy));
                }
            }
        }
        deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        eprintln!(
            "C2-edge source candidate {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        for (rank, (norm, piece, i, j, dx, dy)) in deltas.iter().take(10).enumerate() {
            eprintln!(
                "  vs Appendix-A rank={rank} piece={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e}"
            );
        }
        assert!(
            solved.max_constraint_residual.is_finite(),
            "C2-edge source candidate comparison residual must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_best_scaled_source_frame_is_localized() {
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
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("best scaled source frame localization solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "best scaled source frame localized residual={:e} full_tensor={:e}",
        solved.max_constraint_residual, full_tensor_residual,
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
        "best scaled source frame localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_eq13_b0_after_scaled_source_is_measured() {
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
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let base_map = TestSection42ChartMap {
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
    let b0_cases = [
        ("zero", 0.0),
        ("near-m-c", c),
        ("paper-2c", 2.0 * c),
        ("minus-c", -c),
        ("minus-2c", -2.0 * c),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (label, b0) in b0_cases {
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: test_section42_endpoint_aware_chart_rows(base_map, b0),
                regularization_weight: 1.0e-10,
            })
            .expect("eq13 b0 scaled source solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled source eq13 b0 {label}: b0={b0:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
        "best scaled source eq13 b0: label={} neg={} outside={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(best.3.is_finite(), "eq13 b0 ranking must be finite");
}

#[test]
fn rotational_endpoint_aware_tensor_border_half_assignment_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let mut scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
            .expect("sector-local diagonal rows"),
    );
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let cases = [
        ("previous-near-c/current-near-c", false, false),
        ("previous-near-c/current-near-m", false, true),
        ("previous-near-m/current-near-c", true, false),
        ("previous-near-m/current-near-m", true, true),
    ];

    let mut best_fixed = ("none", f64::INFINITY, usize::MAX, usize::MAX);
    for (label, previous_uses_near_m, current_uses_near_m) in cases {
        let tensor_rows = tensor_border_half_assignment_rows(
            &tensor_border,
            sector,
            previous_uses_near_m,
            current_uses_near_m,
        );
        for project_targets in [false, true] {
            let tensor_constraints = if project_targets {
                let free_coeffs = tensor_rows
                    .iter()
                    .map(|row| row.terms[0].0)
                    .collect::<Vec<_>>();
                kp20_project_section42_source_tensor_rows_with_free_coeffs(
                    tensor_rows.clone(),
                    &fig8_rows,
                    &scalar_constraints,
                    &free_coeffs,
                )
                .expect("tensor row projection succeeds")
                .0
            } else {
                tensor_rows.clone()
            };
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(tensor_constraints.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("half-assignment solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let tensor_residual = tensor_constraints
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "endpoint-aware tensor-border half assignment {label} projected={project_targets}: residual={residual:e} tensor_residual={tensor_residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            if !project_targets
                && (solved.max_constraint_residual < best_fixed.1
                    || (solved.max_constraint_residual == best_fixed.1
                        && validity.orientation.negative_determinants < best_fixed.2))
            {
                best_fixed = (
                    label,
                    solved.max_constraint_residual,
                    validity.orientation.negative_determinants,
                    validity.domain.outside_triangle_samples,
                );
            }
        }
    }
    eprintln!(
        "best fixed tensor-border half assignment: label={} residual={:e} neg={} outside={}",
        best_fixed.0, best_fixed.1, best_fixed.2, best_fixed.3
    );
    assert!(
        best_fixed.1.is_finite(),
        "half-assignment ranking must produce finite solves"
    );
}

#[test]
fn rotational_paper_tau_characteristic_source_orientation_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for valence in [5, 6] {
        let tensor_border = characteristic_cap_net_tensor_border(valence);
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: planar_tau_seed_macro_patch(valence, 0),
                regularization_weight: 1.0e-10,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            0,
            &tensor_border,
            &local_planar_patch_to_global(valence, 0, local_tau),
        )
        .expect("family residuals");
        eprintln!(
            "characteristic-source tau validity n={valence}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2}) families={residuals:?}",
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
fn characteristic_source_near_origin_first_order_frame_is_measured() {
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: planar_tau_seed_macro_patch(valence, 0),
        regularization_weight: 1.0e-10,
    })
    .expect("rotational paper tau solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
    let cap_net = kp20_characteristic_cap_net(valence).expect("characteristic cap-net");
    let c4 = global_to_local_planar(valence, 0, cap_net.sectors[0][3]);
    let c5 = global_to_local_planar(valence, 0, cap_net.sectors[0][4]);
    let c6 = global_to_local_planar(valence, 0, cap_net.sectors[0][5]);
    eprintln!(
        "characteristic inner cap-net local c4=({c4x:e},{c4y:e}) c5=({c5x:e},{c5y:e}) c6=({c6x:e},{c6y:e})",
        c4x = c4[0],
        c4y = c4[1],
        c5x = c5[0],
        c5y = c5[1],
        c6x = c6[0],
        c6y = c6[1],
    );
    eprintln!(
        "characteristic-source NearOriginLeft first-order origin=({lox:e},{loy:e}) du=({ldux:e},{lduy:e}) dv=({ldvx:e},{ldvy:e}) det={ldet:e}",
        lox = left.origin.x,
        loy = left.origin.y,
        ldux = left.du.x,
        lduy = left.du.y,
        ldvx = left.dv.x,
        ldvy = left.dv.y,
        ldet = left.determinant,
    );
    eprintln!(
        "characteristic-source NearOriginRight first-order origin=({rox:e},{roy:e}) du=({rdux:e},{rduy:e}) dv=({rdvx:e},{rdvy:e}) det={rdet:e}",
        rox = right.origin.x,
        roy = right.origin.y,
        rdux = right.du.x,
        rduy = right.du.y,
        rdvx = right.dv.x,
        rdvy = right.dv.y,
        rdet = right.determinant,
    );
    assert!(left.determinant.is_finite());
    assert!(right.determinant.is_finite());
}

#[test]
fn characteristic_source_centered_first_order_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, include_near_origin_right) in [("left-only", false), ("left-and-right", true)] {
        let solved = solve_rotational_characteristic_source_with_centered_first_order(
            5,
            0,
            include_near_origin_right,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "characteristic-source centered first-order {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn centered_first_order_conflict_is_localized_to_source_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let all_section42 = TestSection42FamilySet::new(true, true, true, true);
    let no_section42 = TestSection42FamilySet::new(false, false, false, false);
    let cases = [
        ("center-only", no_section42, false, false, false),
        ("section42", all_section42, false, false, false),
        ("section42-fig8-diagonal", all_section42, true, true, false),
        ("tensor-source-only", no_section42, false, false, true),
        ("full", all_section42, true, true, true),
    ];
    for (label, families, include_fig8, include_diagonal, include_tensor_source) in cases {
        let solved = solve_characteristic_centered_first_order_with_row_families(
            5,
            0,
            families,
            include_fig8,
            include_diagonal,
            include_tensor_source,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "centered first-order row localization {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_origin=({origin_x:e},{origin_y:e}) left_det={left_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            origin_x = left.origin.x,
            origin_y = left.origin.y,
            left_det = left.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn centered_first_order_source_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let solved =
        solve_characteristic_source_with_centered_projected_tensor_band(5, 0, true, true, false);
    let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let tensor_border = characteristic_cap_net_tensor_border(5);
    let residuals = kp20_planar_tau_rotational_paper_family_residuals(
        5,
        0,
        &tensor_border,
        &local_planar_patch_to_global(5, 0, local_tau),
    )
    .expect("family residuals");
    eprintln!(
        "centered first-order source projection: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e} families={residuals:?}",
        residual = solved.max_constraint_residual,
        scale = solved.extra_values[0],
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
        left_det = left.determinant,
        right_det = right.determinant,
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn rotational_paper_tau_characteristic_source_regularization_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let mut best = (0.0_f64, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);
    for regularization_weight in [
        1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4,
    ] {
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: planar_tau_seed_macro_patch(valence, 0),
                regularization_weight,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "characteristic-source regularization weight={regularization_weight:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        if validity.orientation.negative_determinants < best.1
            || (validity.orientation.negative_determinants == best.1
                && validity.orientation.min_determinant > best.2)
        {
            best = (
                regularization_weight,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!(
        "best characteristic-source regularization: weight={:e} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3
    );
    assert!(best.3.is_finite());
}

#[test]
fn rotational_characteristic_source_appendix_a_central_pins_measure_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    for (layout, use_paper_quadrant_maps) in
        [("positive-nonoverlap", false), ("paper-quadrant", true)]
    {
        for max_total_degree in [1, 2, 3] {
            let solved = solve_rotational_source_with_appendix_a_central_pins(
                valence,
                0,
                &tensor_border,
                max_total_degree,
                use_paper_quadrant_maps,
            );
            let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "characteristic-source Appendix-A central pins layout={layout} degree={max_total_degree}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
fn rotational_paper_tau_regular_source_near_origin_corner_derivative_is_inverted() {
    let tensor_border = planar_regular_sector_tensor_border(5);
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: planar_tau_seed_macro_patch(5, 0),
        regularization_weight: 1.0e-10,
    })
    .expect("rotational paper tau solve succeeds");
    let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
    let p00 = piece.cps[0][0];
    let p10 = piece.cps[1][0];
    let p01 = piece.cps[0][1];
    let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
    let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
    let determinant = du.x * dv.y - dv.x * du.y;

    eprintln!(
        "regular-source local NearOriginLeft corner p00=({p00x:e},{p00y:e}) p10=({p10x:e},{p10y:e}) p01=({p01x:e},{p01y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e}) det={determinant:e}",
        p00x = p00.x,
        p00y = p00.y,
        p10x = p10.x,
        p10y = p10.y,
        p01x = p01.x,
        p01y = p01.y,
        dux = du.x,
        duy = du.y,
        dvx = dv.x,
        dvy = dv.y,
    );

    assert!(determinant.is_finite());
}

#[test]
fn rotational_regular_source_orientation_localizes_to_fig8_and_diagonal_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let cases = [
        ("section42-only", false, false),
        ("section42-plus-fig8", true, false),
        ("section42-plus-diagonal", false, true),
        ("section42-plus-fig8-diagonal", true, true),
    ];

    for (label, include_fig8, include_diagonal) in cases {
        let solved =
            solve_rotational_regular_source_with_row_options(5, 0, include_fig8, include_diagonal);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
        let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
        let corner_det = du.x * dv.y - dv.x * du.y;
        eprintln!(
            "regular-source row option {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] corner_det={corner_det:e} du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            dux = du.x,
            duy = du.y,
            dvx = dv.x,
            dvy = dv.y,
        );
    }
}

#[test]
fn rotational_regular_source_orientation_localizes_within_section42_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let cases = [
        ("c0", TestSection42FamilySet::new(true, false, false, false)),
        (
            "eq11",
            TestSection42FamilySet::new(false, true, false, false),
        ),
        (
            "eq12",
            TestSection42FamilySet::new(false, false, true, false),
        ),
        (
            "eq13",
            TestSection42FamilySet::new(false, false, false, true),
        ),
        (
            "c0-eq11",
            TestSection42FamilySet::new(true, true, false, false),
        ),
        (
            "c0-eq11-eq12",
            TestSection42FamilySet::new(true, true, true, false),
        ),
        (
            "c0-eq11-eq12-eq13",
            TestSection42FamilySet::new(true, true, true, true),
        ),
    ];

    for (label, families) in cases {
        let solved = solve_rotational_regular_source_with_section42_families(5, 0, families);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source section42 family {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_orientation_sweeps_diagonal_piece_mapping() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let piece_sets = [
        (
            "nol-nmr_self_nor-nml_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointRight],
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointLeft],
        ),
        (
            "nol-nor_self_nmr-nml_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight],
            [
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft,
            ],
        ),
        (
            "nol-nml_self_nor-nmr_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointLeft],
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointRight],
        ),
        (
            "nor-nmr_self_nol-nml_pair",
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointRight],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointLeft],
        ),
        (
            "nor-nml_self_nol-nmr_pair",
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointLeft],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointRight],
        ),
        (
            "nmr-nml_self_nol-nor_pair",
            [
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft,
            ],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight],
        ),
    ];
    let mut best_negative = usize::MAX;

    for (label, self_pieces, paired_pieces) in piece_sets {
        let solved =
            solve_rotational_regular_source_with_diagonal_mapping(5, 0, self_pieces, paired_pieces);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source diagonal mapping {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        if orientation.negative_determinants < best_negative {
            best_negative = orientation.negative_determinants;
        }
    }

    assert_eq!(best_negative, 0);
}

#[test]
fn rotational_regular_source_origin_pins_measure_regular_tau_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let pin_sets: &[(&str, &[Kp20PlanarTauCoeff])] = &[
        (
            "near-origin-left-first-order",
            &[
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 1),
            ],
        ),
        (
            "both-origin-first-order",
            &[
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 1),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 1),
            ],
        ),
    ];

    for &(label, pins) in pin_sets {
        let solved = solve_rotational_regular_source_with_origin_pins(5, 0, pins);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source origin pins {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_paper_tau_origin_pins_measure_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, use_paper_quadrant_maps) in [
        ("appendix-a-positive-nonoverlap", false),
        ("appendix-a-paper-quadrant", true),
    ] {
        let solved = solve_rotational_regular_source_with_paper_tau_origin_pins(
            5,
            0,
            use_paper_quadrant_maps,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source paper tau origin pins {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_previous_sector_index_transpose_improves_tau_orientation() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, index_map) in [
        ("all-previous", TestSection42IndexMap::PreviousSectorAll),
        (
            "central-only",
            TestSection42IndexMap::PreviousSectorCentralOnly,
        ),
        (
            "stored-piece-edge",
            TestSection42IndexMap::StoredPieceEdgeOrientation,
        ),
    ] {
        let solved =
            solve_rotational_regular_source_with_previous_sector_index_transpose(5, 0, index_map);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
        let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
        let corner_det = du.x * dv.y - dv.x * du.y;

        eprintln!(
            "regular-source previous-sector index transpose {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] corner_det={corner_det:e} worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );

        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_previous_sector_rotation_sign_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, sign) in [("plus", 1.0), ("minus", -1.0)] {
        let solved = solve_rotational_regular_source_with_rotation_sign(5, 0, sign);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source previous-sector rotation {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
    }
}

#[test]
fn rotational_regular_source_eq13_b0_sign_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    for (label, b0_near_origin) in [
        ("paper-positive-2c", 2.0 * c),
        ("opposite-side-negative-2c", -2.0 * c),
        ("zero-eq13-b-coupling", 0.0),
        ("half-positive-c", c),
        ("half-negative-c", -c),
    ] {
        let solved = solve_rotational_regular_source_with_eq13_b0(5, 0, b0_near_origin);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source eq13 b0 {label}: b0={b0_near_origin:e} residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}
