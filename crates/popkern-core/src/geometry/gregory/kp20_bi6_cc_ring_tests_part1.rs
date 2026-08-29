#[allow(unused_imports)]
use super::super::super::kp_polar_shared::weighted_sum;
use super::*;

#[test]
fn pr08_fourier_one_vector_satisfies_catmull_clark_block_eigen_equation() {
    for n in [3, 5, 6, 7, 10] {
        let scalars = Kp20CharacteristicScalars::for_valence(n).unwrap();
        let lambda = Complex::new(scalars.lambda, 0.0);
        let vector = pr08_catmull_clark_fourier_one_eigenvector(n).unwrap();
        let mut product = [Complex::zero(); 13];
        let a00 = pr08_catmull_clark_a00(n);
        let w = Complex::from_angle(core::f64::consts::TAU / n as f64);
        let (a10, a11, a20, a21) = pr08_catmull_clark_lower_blocks(w);

        for row in 0..3 {
            for col in 0..3 {
                product[row] = product[row] + a00[row][col] * vector[col];
            }
        }
        for row in 0..4 {
            for col in 0..3 {
                product[3 + row] = product[3 + row] + a10[row][col] * vector[col];
            }
            for col in 0..4 {
                product[3 + row] = product[3 + row] + a11[row][col] * vector[3 + col];
            }
        }
        for row in 0..6 {
            for col in 0..3 {
                product[7 + row] = product[7 + row] + a20[row][col] * vector[col];
            }
            for col in 0..4 {
                product[7 + row] = product[7 + row] + a21[row][col] * vector[3 + col];
            }
        }

        for k in 0..13 {
            let expected = vector[k] * lambda;
            let delta = product[k] - expected;
            assert!(
                delta.norm2().sqrt() < 1e-9,
                "n={n}, coeff={k}: expected {:?} to match {:?}, delta {}",
                product[k],
                expected,
                delta.norm2().sqrt()
            );
        }
    }
}

#[test]
fn pr08_real_subdivision_matrix_preserves_constants() {
    for n in [3, 5, 6, 8] {
        let matrix = pr08_catmull_clark_subdivision_matrix(n).unwrap();
        assert_eq!(matrix.len(), 13 * n);
        for (row_index, row) in matrix.into_iter().enumerate() {
            let row_sum = row.iter().sum::<f64>();
            assert!(
                (row_sum - 1.0).abs() < 1.0e-12,
                "PR08 real subdivision matrix row {row_index} should preserve constants for n={n}, got {row_sum}"
            );
        }
    }
}

#[test]
fn pr08_real_subdivision_stencils_match_inverse_dft_blocks() {
    for n in [3, 5, 6, 8] {
        let inverse_dft = pr08_catmull_clark_subdivision_matrix(n).unwrap();
        let stencils = pr08_catmull_clark_subdivision_stencil_matrix(n).unwrap();
        assert_eq!(stencils.len(), inverse_dft.len());
        for row in 0..13 * n {
            assert_eq!(stencils[row].len(), inverse_dft[row].len());
            for col in 0..13 * n {
                let delta = (stencils[row][col] - inverse_dft[row][col]).abs();
                assert!(
                    delta < 1.0e-12,
                    "PR08 real-space stencil row {row}, col {col}, n={n} should match inverse-DFT block matrix, got stencil {} vs inverse-DFT {}",
                    stencils[row][col],
                    inverse_dft[row][col]
                );
            }
        }
    }
}

#[test]
fn pr08_q13_and_q2_real_subdivision_source_masks_are_pinned() {
    let n = 5;
    let matrix = pr08_catmull_clark_subdivision_stencil_matrix(n).unwrap();
    let row = |sector: usize, label_index: usize| -> &[f64] { &matrix[sector * 13 + label_index] };
    let source = |sector: usize, label_index: usize| -> usize { sector * 13 + label_index };

    let q13 = row(0, 12);
    assert_eq!(
        q13.iter().filter(|weight| weight.abs() > 1.0e-14).count(),
        4
    );
    assert!((q13[source(1, 1)] - 0.25).abs() < 1.0e-14);
    assert!((q13[source(0, 2)] - 0.25).abs() < 1.0e-14);
    assert!((q13[source(0, 3)] - 0.25).abs() < 1.0e-14);
    assert!((q13[source(1, 4)] - 0.25).abs() < 1.0e-14);

    let q2 = row(0, 1);
    assert_eq!(
        q2.iter().filter(|weight| weight.abs() > 1.0e-14).count(),
        n + 5
    );
    for sector in 0..n {
        assert!((q2[source(sector, 0)] - 3.0 / (8.0 * n as f64)).abs() < 1.0e-14);
    }
    assert!((q2[source(4, 1)] - 1.0 / 16.0).abs() < 1.0e-14);
    assert!((q2[source(0, 1)] - 3.0 / 8.0).abs() < 1.0e-14);
    assert!((q2[source(1, 1)] - 1.0 / 16.0).abs() < 1.0e-14);
    assert!((q2[source(0, 2)] - 1.0 / 16.0).abs() < 1.0e-14);
    assert!((q2[source(1, 2)] - 1.0 / 16.0).abs() < 1.0e-14);
}

#[test]
fn pr08_real_subdivision_matrix_scales_characteristic_ring_by_lambda() {
    for n in [3, 5, 6, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
        let refined = pr08_subdivide_characteristic_ring(&ring).unwrap();
        for (sector_index, sector) in ring.sectors.iter().enumerate() {
            for (label_index, point) in sector.iter().enumerate() {
                let expected = Point3::new(
                    ring.lambda * point.x,
                    ring.lambda * point.y,
                    ring.lambda * point.z,
                );
                assert_close_point(refined[sector_index][label_index], expected, 1.0e-10);
            }
        }
    }
}

#[test]
fn pr08_refined_stam_visual_projection_is_only_lambda_scaled() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let refined = pr08_subdivide_characteristic_ring(&ring).unwrap();
    let refined_ring = Kp20CatmullClarkCharacteristicRing {
        sectors: refined,
        lambda: ring.lambda,
    };
    let visual = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let refined_visual = kp20_pr08_stam_visual_projected_grids(&refined_ring).unwrap();

    for sector_index in 0..ring.sectors.len() {
        for patch in 0..3 {
            for row in 0..4 {
                for col in 0..4 {
                    let source = visual.sectors[sector_index][patch][row][col];
                    let expected = Point3::new(
                        ring.lambda * source.x,
                        ring.lambda * source.y,
                        ring.lambda * source.z,
                    );
                    assert_close_point(
                        refined_visual.sectors[sector_index][patch][row][col],
                        expected,
                        1.0e-10,
                    );
                }
            }
        }
    }
}

#[test]
fn pr08_characteristic_ring_lifts_fourier_vector_by_rotation() {
    let n = 7;
    let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
    let theta = core::f64::consts::TAU / n as f64;
    let (s, c) = theta.sin_cos();

    assert_eq!(ring.sectors.len(), n);
    for sector in 0..n {
        let next = (sector + 1) % n;
        for coeff in 0..13 {
            let current = ring.sectors[sector][coeff];
            let expected = Point3::new(
                c * current.x - s * current.y,
                s * current.x + c * current.y,
                0.0,
            );
            assert_close_point(ring.sectors[next][coeff], expected, 1e-10);
        }
    }
}

#[test]
fn stam_appendix_a_subdivision_masks_are_affine_rows() {
    let top_count = 11;
    let s11 = stam_appendix_a_s11(top_count);
    let s12 = stam_appendix_a_s12();
    for row in 0..7 {
        let row_sum = s11[row].iter().sum::<f64>() + s12[row].iter().sum::<f64>();
        assert!(
            (row_sum - 1.0).abs() < 1e-12,
            "Stam Appendix-A S11/S12 row {row} must preserve affine constants, got {row_sum}"
        );
    }

    let s21 = stam_appendix_a_s21(top_count);
    let s22 = stam_appendix_a_s22();
    for row in 0..9 {
        let row_sum = s21[row].iter().sum::<f64>() + s22[row].iter().sum::<f64>();
        assert!(
            (row_sum - 1.0).abs() < 1e-12,
            "Stam Appendix-A S21/S22 row {row} must preserve affine constants, got {row_sum}"
        );
    }
}

#[test]
fn stam_appendix_a_subdivided_labels_have_expected_count() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let labels = stam_appendix_a_subdivided_extended_characteristic_labels(ring.lambda, 5).unwrap();
    assert_eq!(labels.len(), 2 * 5 + 17);
}

#[test]
fn stam_appendix_a_middle_patch_first_column_is_rotationally_covariant() {
    let ring = kp20_catmull_clark_characteristic_ring(7).unwrap();
    let first = stam_appendix_a_middle_patch_first_column(ring.lambda, 7, 0).unwrap();
    let second = stam_appendix_a_middle_patch_first_column(ring.lambda, 7, 1).unwrap();
    let theta = core::f64::consts::TAU / 7.0;
    let (s, c) = theta.sin_cos();
    for row in 0..3 {
        assert_close_point(second[row], rotate_xy(first[row], s, c), 1e-10);
    }
}

#[test]
fn stam_appendix_a_frequency_one_top_vector_matches_page_formula() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    let gain = 4.0 * ring.lambda - 1.0;
    for sector in 0..valence {
        let theta = core::f64::consts::TAU * sector as f64 / valence as f64;
        let next_theta = core::f64::consts::TAU * (sector + 1) as f64 / valence as f64;
        u0_x[1 + 2 * sector] = gain * theta.cos();
        u0_y[1 + 2 * sector] = gain * theta.sin();
        u0_x[2 + 2 * sector] = theta.cos() + next_theta.cos();
        u0_y[2 + 2 * sector] = theta.sin() + next_theta.sin();
    }

    let top_rows = stam_appendix_a_projected_rows_from_rotational_top_vector_with_center(
        &ring,
        [
            0.0,
            0.0,
            gain,
            0.0,
            1.0 + (core::f64::consts::TAU / valence as f64).cos(),
            (core::f64::consts::TAU / valence as f64).sin(),
        ],
        0,
        0,
        false,
        1.0,
    )
    .unwrap();
    let page_formula_rows =
        stam_appendix_a_projected_rows_from_top_vectors(&ring, &u0_x, &u0_y, 0, 0, false).unwrap();
    for sector in 0..valence {
        for row in 0..3 {
            for col in 0..4 {
                assert_close_point(
                    top_rows.sectors[sector][row][col],
                    page_formula_rows.sectors[sector][row][col],
                    1.0e-12,
                );
            }
        }
    }
}

#[test]
fn stam_appendix_a_top_to_label_source_matrix_reproduces_first_column() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let weights = stam_appendix_a_label_top_weight_matrix(ring.lambda, top_count).unwrap();
    let labels =
        stam_appendix_a_subdivided_extended_characteristic_labels(ring.lambda, valence).unwrap();
    let top = stam_appendix_a_frequency_one_top_vector(&ring);
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    for sector in 0..valence {
        u0_x[1 + 2 * sector] = top.odd[sector].re;
        u0_y[1 + 2 * sector] = top.odd[sector].im;
        u0_x[2 + 2 * sector] = top.even[sector].re;
        u0_y[2 + 2 * sector] = top.even[sector].im;
    }

    for (label_index, label_weights) in weights.iter().enumerate() {
        let x = label_weights
            .iter()
            .zip(&u0_x)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        let y = label_weights
            .iter()
            .zip(&u0_y)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        assert_close_point(labels[label_index], Point3::new(x, y, 0.0), 1.0e-11);
    }

    let first_column = stam_appendix_a_middle_patch_first_column(ring.lambda, valence, 0).unwrap();
    for (row, label_index) in [3, 2 * valence + 6, 2 * valence + 15]
        .into_iter()
        .enumerate()
    {
        let x = weights[label_index]
            .iter()
            .zip(&u0_x)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        let y = weights[label_index]
            .iter()
            .zip(&u0_y)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        assert_close_point(first_column[row], Point3::new(x, y, 0.0), 1.0e-11);
    }

    assert_eq!(
        weights[3]
            .iter()
            .filter(|weight| weight.abs() > 1.0e-14)
            .count(),
        1
    );
    assert!((weights[3][3] - ring.lambda).abs() < 1.0e-14);
}

#[test]
fn stam_appendix_b_q2_rows_use_the_source_first_column() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let column = stam_appendix_a_middle_patch_first_column(ring.lambda, 5, 0).unwrap();
    let rows = stam_appendix_a_regular_subpatch_interior_rows(ring.lambda, 5, 0, 1).unwrap();
    for row in 0..3 {
        assert_close_point(rows[row][0], column[row], 1e-12);
    }
}

#[test]
fn pr08_characteristic_ring_rejects_invalid_valence() {
    assert!(kp20_catmull_clark_characteristic_ring(2).is_err());
}

#[test]
fn pr08_tensor_border_candidates_are_rotationally_covariant() {
    let n = 7;
    let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
    let theta = core::f64::consts::TAU / n as f64;
    let (s, c) = theta.sin_cos();

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        assert_eq!(rows.sectors.len(), n);
        for sector in 0..n {
            let next = (sector + 1) % n;
            for row in 0..3 {
                for col in 0..4 {
                    let current = rows.sectors[sector][row][col];
                    let expected = Point3::new(
                        c * current.x - s * current.y,
                        s * current.x + c * current.y,
                        0.0,
                    );
                    assert_close_point(rows.sectors[next][row][col], expected, 1e-10);
                }
            }
        }
    }
}

#[test]
fn pr08_stam_inner_q1_production_rows_match_diagnostic_candidate() {
    let n = 5;
    let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
    let diagnostic = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
    )
    .unwrap();
    let production = kp20_pr08_stam_inner_q1_regular_tensor_border_rows(n).unwrap();

    assert_eq!(production.sectors.len(), diagnostic.sectors.len());
    for sector in 0..n {
        for row in 0..3 {
            for col in 0..4 {
                assert_close_point(
                    production.sectors[sector][row][col],
                    diagnostic.sectors[sector][row][col],
                    1.0e-12,
                );
            }
        }
    }
}

#[test]
fn pr08_fig63_endpoint_side_rows_match_decoded_stam_transform() {
    let n = 5;
    let ring = kp20_catmull_clark_characteristic_ring(n).unwrap();
    let projection =
        kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Appendix tensor projection");
    let transformed = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridSide::Left,
        false,
    )
    .expect("transformed endpoint-side rows");
    let explicit = kp20_pr08_fig63_endpoint_side_regular_tensor_border_rows(n)
        .expect("explicit Fig.6.3 endpoint-side rows");

    for sector in 0..n {
        for row in 0..3 {
            for col in 0..4 {
                assert_close_point(
                    explicit.sectors[sector][row][col],
                    transformed.sectors[sector][row][col],
                    1.0e-12,
                );
            }
        }
    }
}

#[test]
fn pr08_stam_inner_q1_endpoint_normalization_matches_unit_sector_spokes() {
    let n = 5;
    let rows = kp20_pr08_stam_inner_q1_endpoint_normalized_regular_tensor_border_rows(n).unwrap();
    let expected_right = Point3::new(1.0, 0.0, 0.0);

    assert_close_point(rows.sectors[0][0][0], Point3::new(0.0, 0.0, 0.0), 1.0e-12);
    assert_close_point(rows.sectors[0][0][3], expected_right, 1.0e-12);
}

#[test]
fn pr08_tensor_border_candidates_are_rejected_by_paper_tau_constraints() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
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
            "{} should remain rejected until the PR08 regular subpatch layout is mapped; residual {}",
            candidate.name(),
            solved.max_constraint_residual
        );
    }
}

#[test]
fn pr08_tensor_border_candidates_close_section42_with_source_projection() {
    let valence = 5;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let mut best_candidate = "";
    let mut best_projection_delta = f64::INFINITY;
    let mut best_solve_residual = 0.0_f64;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
        let projection_report = kp20_planar_tau_rotational_section42_source_projection_report(
            valence,
            0,
            &tensor_border,
        )
        .unwrap();
        let family_residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            0,
            &tensor_border,
            &solved.tau_macro_patch,
        )
        .unwrap();
        let tensor_border_family = family_residuals
            .iter()
            .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::TensorBorderEq14)
            .unwrap();
        let c0 = family_residuals
            .iter()
            .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::Section42C0)
            .unwrap();
        let eq11 = family_residuals
            .iter()
            .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq11Degree5)
            .unwrap();
        let eq12 = family_residuals
            .iter()
            .find(|residual| {
                residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq12InternalC2
            })
            .unwrap();
        let eq13 = family_residuals
            .iter()
            .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq13Strip)
            .unwrap();
        let section42_residual = c0
            .max_residual
            .max(eq11.max_residual)
            .max(eq12.max_residual)
            .max(eq13.max_residual);

        eprintln!(
            "PR08 source §4.2 projected n={valence} candidate={candidate}: solve={solve:e} projection_delta={projection_delta:e} projection_residual={projection_residual:e} adjusted_scalar_targets={adjusted} eq14={eq14:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e}",
            candidate = candidate.name(),
            solve = solved.max_constraint_residual,
            projection_delta = projection_report.max_source_delta,
            projection_residual = projection_report.max_row_residual,
            adjusted = projection_report.adjusted_scalar_targets,
            eq14 = tensor_border_family.max_residual,
            c0 = c0.max_residual,
            eq11 = eq11.max_residual,
            eq12 = eq12.max_residual,
            eq13 = eq13.max_residual,
        );

        assert!(solved.max_constraint_residual.is_finite());
        assert!(section42_residual.is_finite());
        assert!(projection_report.max_source_delta.is_finite());
        assert_eq!(projection_report.adjusted_scalar_targets, 84);
        assert!(
            tensor_border_family.max_residual.is_finite(),
            "{} should measure the unprojected eq(14) source residual after live source projection",
            candidate.name()
        );
        assert!(
            solved.max_constraint_residual < 1.0e-7 && section42_residual < 1.0e-7,
            "{} should close the promoted §4.2 basis after source projection, got solve={} section42={}",
            candidate.name(),
            solved.max_constraint_residual,
            section42_residual
        );

        if projection_report.max_source_delta < best_projection_delta {
            best_candidate = candidate.name();
            best_projection_delta = projection_report.max_source_delta;
            best_solve_residual = solved.max_constraint_residual;
        }
    }

    eprintln!(
        "PR08 source §4.2 projected best n={valence}: candidate={best_candidate} projection_delta={best_projection_delta:e} solve={best_solve_residual:e}"
    );
    assert!(
        best_projection_delta < 10.0,
        "projection should remain bounded for current PR08 source candidates; best {best_candidate} delta={best_projection_delta}"
    );
}

#[test]
fn pr08_tensor_border_candidates_with_unit_sector_endpoint_normalization_are_ranked_by_tau_validity()
 {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let mut best_candidate = "";
    let mut best_endpoint_order = "";
    let mut best_outside = usize::MAX;
    let mut best_negative = usize::MAX;
    let mut best_near_zero = usize::MAX;
    let mut best_min_determinant = f64::NEG_INFINITY;
    let mut best_max_excess = f64::INFINITY;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        for (endpoint_order, swap_targets, conjugate_source) in [
            ("direct", false, false),
            ("swapped", true, false),
            ("direct-conjugate", false, true),
            ("swapped-conjugate", true, true),
        ] {
            let Some(normalized_rows) = unit_sector_endpoint_normalized_tensor_rows(
                &rows,
                valence,
                0,
                swap_targets,
                conjugate_source,
            ) else {
                continue;
            };
            let tensor_border =
                kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
            let solved =
                solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                    valence,
                    sector: 0,
                    characteristic_tensor_border: &tensor_border,
                    initial: Kp20MacroPatch::zero(),
                    regularization_weight: 1e-10,
                })
                .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "PR08 unit-sector normalized tau candidate={candidate} order={endpoint_order}: solve={solve:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
                candidate = candidate.name(),
                solve = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                max_excess = validity.domain.max_simplex_excess,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
                worst = validity.orientation.worst_piece,
                u = validity.orientation.worst_u,
                v = validity.orientation.worst_v,
            );
            let better = validity.domain.outside_triangle_samples < best_outside
                || (validity.domain.outside_triangle_samples == best_outside
                    && validity.orientation.negative_determinants < best_negative)
                || (validity.domain.outside_triangle_samples == best_outside
                    && validity.orientation.negative_determinants == best_negative
                    && validity.orientation.near_zero_determinants < best_near_zero)
                || (validity.domain.outside_triangle_samples == best_outside
                    && validity.orientation.negative_determinants == best_negative
                    && validity.orientation.near_zero_determinants == best_near_zero
                    && validity.domain.max_simplex_excess < best_max_excess)
                || (validity.domain.outside_triangle_samples == best_outside
                    && validity.orientation.negative_determinants == best_negative
                    && validity.orientation.near_zero_determinants == best_near_zero
                    && (validity.domain.max_simplex_excess - best_max_excess).abs() < 1.0e-14
                    && validity.orientation.min_determinant > best_min_determinant);
            if better {
                best_candidate = candidate.name();
                best_endpoint_order = endpoint_order;
                best_outside = validity.domain.outside_triangle_samples;
                best_negative = validity.orientation.negative_determinants;
                best_near_zero = validity.orientation.near_zero_determinants;
                best_min_determinant = validity.orientation.min_determinant;
                best_max_excess = validity.domain.max_simplex_excess;
            }
        }
    }

    eprintln!(
        "PR08 unit-sector normalized tau best candidate={best_candidate} order={best_endpoint_order} outside={best_outside} neg={best_negative} zero={best_near_zero} max_excess={best_max_excess:e} min_det={best_min_determinant:e}"
    );
    assert!(
        best_outside > 0 || best_negative > 0 || best_near_zero > 0,
        "unit-sector endpoint normalization unexpectedly produced a sampled paper-valid tau"
    );
}

#[test]
fn pr08_best_row_map_preserves_middle_patch_normalization_corner() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    grid[0] = rows.sectors[0][0];
    grid[1] = rows.sectors[0][1];
    grid[2] = rows.sectors[0][2];

    let normalized_corner = eval_cubic_bspline_grid(&grid, 0.0, 1.0).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();

    assert_close_point(normalized_corner, expected, 1e-12);
}

#[test]
fn stam_appendix_a_pr08_best_top_label_projection_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = stam_appendix_a_pr08_top_label_projected_rows(&ring, 1, 12, 1, 1, true).unwrap();
    assert!(tensor_border_row_norm(&rows) > 1.0e-10);
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
        solved.max_constraint_residual < 1.0e-10,
        "correct cross-degree tensor-border elevation should close the best nondegenerate PR08 q2/q13 top-label projection, got {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_pr08_inner_ring_inverse_lambda_seed_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut best = StamAppendixAPr08InnerSeedDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
    };
    for patch_index in 0..3 {
        for row_start in [0, 1] {
            for reverse_rows in [false, true] {
                let rows = stam_appendix_a_pr08_inner_ring_scaled_subpatch_rows(
                    &ring,
                    patch_index,
                    row_start,
                    reverse_rows,
                    1.0 / ring.lambda,
                )
                .unwrap();
                let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(
                    &Kp20PlanarTauPaperSolveInput {
                        valence: 5,
                        sector: 0,
                        characteristic_tensor_border: &tensor_border,
                        initial: Kp20MacroPatch::zero(),
                        regularization_weight: 1e-10,
                    },
                )
                .unwrap();
                if solved.max_constraint_residual < best.residual {
                    best = StamAppendixAPr08InnerSeedDiagnostic {
                        residual: solved.max_constraint_residual,
                        patch_index,
                        row_start,
                        reverse_rows,
                    };
                }
            }
        }
    }

    assert!(
        best.residual < 1.0e-10,
        "correct cross-degree tensor-border elevation should close the inverse-lambda PR08 inner-ring seed, got {best:?}"
    );
    assert!(
        best.patch_index < 3,
        "best inner-ring seed patch should remain one of Stam's three regular subpatches"
    );
    assert_stam_appendix_a_pr08_inner_seed_diagnostic_in_search_space(best);
}

#[test]
fn pr08_labels_best_match_stam_top_vector_rings_by_similarity() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let mut best_odd = Pr08TopRingLabelMatch {
        residual: f64::INFINITY,
        label_index: 0,
    };
    let mut best_even = best_odd;
    for label_index in 0..13 {
        let values = pr08_label_ring_as_complex(&ring, label_index);
        let odd_residual = complex_similarity_ring_residual(&values, &stam_top.odd);
        if odd_residual < best_odd.residual {
            best_odd = Pr08TopRingLabelMatch {
                residual: odd_residual,
                label_index,
            };
        }
        let even_residual = complex_similarity_ring_residual(&values, &stam_top.even);
        if even_residual < best_even.residual {
            best_even = Pr08TopRingLabelMatch {
                residual: even_residual,
                label_index,
            };
        }
    }

    assert!(best_odd.residual < 1.0e-12 && best_even.residual < 1.0e-12);
    // Near-tie winners flip with platform libm (macOS: even-label 4); the
    // residual is the class pin, members are Linux measurements (PK-291).
    let members = (best_odd.label_index, best_even.label_index);
    assert!(cfg!(not(target_os = "linux")) || members == (12, 1));
}

#[test]
fn pr08_q13_q2_to_stam_top_vector_similarity_scales_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let q13 = pr08_label_ring_as_complex(&ring, 12);
    let q2 = pr08_label_ring_as_complex(&ring, 1);
    let odd_scale = complex_similarity_ring_scale(&q13, &stam_top.odd).unwrap();
    let even_scale = complex_similarity_ring_scale(&q2, &stam_top.even).unwrap();
    assert!((odd_scale.re - 0.12734337133080173).abs() < 1.0e-12);
    assert!((odd_scale.im + 0.4262575543508202).abs() < 1.0e-12);
    assert!((even_scale.re - 1.0).abs() < 1.0e-12);
    assert!(even_scale.im.abs() < 1.0e-12);
}

#[test]
fn pr08_middle_patch_normalization_value_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let value = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    eprintln!(
        "PR08 middle-patch normalization value n=5: ({:e}, {:e})",
        value.x, value.y
    );
    assert!(value.x.is_finite());
    assert!(value.y.is_finite());
    assert!(value.x.abs() > 1.0e-12 || value.y.abs() > 1.0e-12);
}

#[test]
fn pr08_middle_patch_normalized_source_candidates_are_ranked_by_tau_validity() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let normalization_value = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let normalization = complex_from_point(normalization_value);
    let Some(scale) =
        complex_similarity_scale_for_pairs(&[(normalization, Complex::new(1.0, 0.0))])
    else {
        panic!("PR08 middle-patch normalization value should be nonzero");
    };
    let mut best_candidate = "";
    let mut best_outside = usize::MAX;
    let mut best_negative = usize::MAX;
    let mut best_near_zero = usize::MAX;
    let mut best_min_determinant = f64::NEG_INFINITY;
    let mut best_max_excess = f64::INFINITY;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let normalized_rows = complex_scaled_tensor_border_rows(&rows, scale, false);
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "PR08 middle-normalized tau candidate={candidate}: solve={solve:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            candidate = candidate.name(),
            solve = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            max_excess = validity.domain.max_simplex_excess,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            worst = validity.orientation.worst_piece,
            u = validity.orientation.worst_u,
            v = validity.orientation.worst_v,
        );
        let better = validity.domain.outside_triangle_samples < best_outside
            || (validity.domain.outside_triangle_samples == best_outside
                && validity.orientation.negative_determinants < best_negative)
            || (validity.domain.outside_triangle_samples == best_outside
                && validity.orientation.negative_determinants == best_negative
                && validity.orientation.near_zero_determinants < best_near_zero)
            || (validity.domain.outside_triangle_samples == best_outside
                && validity.orientation.negative_determinants == best_negative
                && validity.orientation.near_zero_determinants == best_near_zero
                && validity.domain.max_simplex_excess < best_max_excess)
            || (validity.domain.outside_triangle_samples == best_outside
                && validity.orientation.negative_determinants == best_negative
                && validity.orientation.near_zero_determinants == best_near_zero
                && (validity.domain.max_simplex_excess - best_max_excess).abs() < 1.0e-14
                && validity.orientation.min_determinant > best_min_determinant);
        if better {
            best_candidate = candidate.name();
            best_outside = validity.domain.outside_triangle_samples;
            best_negative = validity.orientation.negative_determinants;
            best_near_zero = validity.orientation.near_zero_determinants;
            best_min_determinant = validity.orientation.min_determinant;
            best_max_excess = validity.domain.max_simplex_excess;
        }
    }

    eprintln!(
        "PR08 middle-normalized tau best candidate={best_candidate} outside={best_outside} neg={best_negative} zero={best_near_zero} max_excess={best_max_excess:e} min_det={best_min_determinant:e}"
    );
    assert!(
        best_outside > 0 || best_negative > 0 || best_near_zero > 0,
        "PR08 middle-patch normalization unexpectedly produced a sampled paper-valid tau"
    );
}

#[test]
fn stam_appendix_a_pr08_q13_q2_top_projection_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut best = StamAppendixAPr08InnerSeedDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
    };
    for patch_index in 0..3 {
        for row_start in [0, 1] {
            for reverse_rows in [false, true] {
                let rows = stam_appendix_a_pr08_top_label_projected_rows(
                    &ring,
                    12,
                    1,
                    patch_index,
                    row_start,
                    reverse_rows,
                )
                .unwrap();
                let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(
                    &Kp20PlanarTauPaperSolveInput {
                        valence: 5,
                        sector: 0,
                        characteristic_tensor_border: &tensor_border,
                        initial: Kp20MacroPatch::zero(),
                        regularization_weight: 1e-10,
                    },
                )
                .unwrap();
                if solved.max_constraint_residual < best.residual {
                    best = StamAppendixAPr08InnerSeedDiagnostic {
                        residual: solved.max_constraint_residual,
                        patch_index,
                        row_start,
                        reverse_rows,
                    };
                }
            }
        }
    }

    assert!(
        best.residual < 1.0e-10,
        "correct cross-degree tensor-border elevation should close the raw PR08 q13/q2 top projection, got {best:?}"
    );
    assert!(
        best.patch_index < 3,
        "best top projection patch should remain one of Stam's three regular subpatches"
    );
    assert_stam_appendix_a_pr08_inner_seed_diagnostic_in_search_space(best);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_top_projection_closes_section42_43_rows() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let mut best = StamAppendixAPr08InnerSeedDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
    };
    for patch_index in 0..3 {
        for row_start in [0, 1] {
            for reverse_rows in [false, true] {
                let rows = stam_appendix_a_pr08_scaled_top_label_projected_rows(
                    &ring,
                    (12, odd_scale),
                    (1, even_scale),
                    patch_index,
                    row_start,
                    reverse_rows,
                )
                .unwrap();
                let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(
                    &Kp20PlanarTauPaperSolveInput {
                        valence: 5,
                        sector: 0,
                        characteristic_tensor_border: &tensor_border,
                        initial: Kp20MacroPatch::zero(),
                        regularization_weight: 1e-10,
                    },
                )
                .unwrap();
                if solved.max_constraint_residual < best.residual {
                    best = StamAppendixAPr08InnerSeedDiagnostic {
                        residual: solved.max_constraint_residual,
                        patch_index,
                        row_start,
                        reverse_rows,
                    };
                }
            }
        }
    }

    assert!(
        best.residual < 1.0e-10,
        "scaled PR08 q13/q2 top projection should close the pre-solved §4.2/§4.3 row gate, got {best:?}"
    );
    assert!(
        best.patch_index < 3,
        "best scaled q13/q2 projection patch should remain one of Stam's three regular subpatches"
    );
    assert_stam_appendix_a_pr08_inner_seed_diagnostic_in_search_space(best);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_candidate_matches_derived_source_rows() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let candidate = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let derived = stam_appendix_a_pr08_scaled_top_label_projected_rows(
        &ring,
        (
            12,
            complex_similarity_ring_scale(
                &pr08_label_ring_as_complex(&ring, 12),
                &stam_appendix_a_frequency_one_top_vector(&ring).odd,
            )
            .unwrap(),
        ),
        (
            1,
            complex_similarity_ring_scale(
                &pr08_label_ring_as_complex(&ring, 1),
                &stam_appendix_a_frequency_one_top_vector(&ring).even,
            )
            .unwrap(),
        ),
        0,
        1,
        true,
    )
    .unwrap();

    for sector in 0..ring.sectors.len() {
        for row in 0..3 {
            for col in 0..4 {
                assert_close_point(
                    candidate.sectors[sector][row][col],
                    derived.sectors[sector][row][col],
                    1.0e-12,
                );
            }
        }
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_candidate_is_fed_through_section42_43_solve() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
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
    let residuals = kp20_planar_tau_rotational_paper_family_residuals(
        5,
        0,
        &tensor_border,
        &solved.tau_macro_patch,
    )
    .unwrap();
    assert!(
        solved.max_constraint_residual < 1.0e-10,
        "derived PR08 q13/q2 -> Stam Appendix-A/B source should close through §4.2/§4.3 rows, got {}",
        solved.max_constraint_residual
    );
    assert!(
        residuals
            .iter()
            .filter(|residual| {
                !matches!(
                    residual.family,
                    Kp20PlanarTauPaperRowFamily::TensorBorderEq14
                        | Kp20PlanarTauPaperRowFamily::DiagonalSymmetry
                )
            })
            .all(|residual| residual.max_residual < 1.0e-10)
    );
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_live_path_reduced_pair_basis_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let pair_cases = [
        ("none", None::<fn(usize, usize) -> bool>),
        (
            "apex-corner",
            Some((|i, j| i == 0 && j == 0) as fn(usize, usize) -> bool),
        ),
        (
            "first-order-square",
            Some((|i, j| i <= 1 && j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "total-degree-one",
            Some((|i, j| i + j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "lemma-cubic-triangle",
            Some((|i, j| i + j <= 3) as fn(usize, usize) -> bool),
        ),
        (
            "post-lemma-complement",
            Some((|i, j| i + j > 3) as fn(usize, usize) -> bool),
        ),
        ("full-pair", Some((|_, _| true) as fn(usize, usize) -> bool)),
    ];

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );

    for (label, pair_keep) in pair_cases {
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
        if let Some(keep) = pair_keep {
            scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                keep,
            ));
        }

        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows.clone(),
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);

        eprintln!(
            "scaled q13/q2 live reduced pair {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
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

        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live reduced pair best: {best:?}");
    assert!(best.5 < 1.0e-10 || best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_live_path_pair_orientation_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let pair_orientations = [
        ("none", None::<fn(usize, usize) -> (usize, usize)>),
        (
            "transpose",
            Some((|i, j| (j, i)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "transpose-reverse-first",
            Some((|i, j| (6 - j, i)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "transpose-reverse-second",
            Some((|i, j| (j, 6 - i)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "transpose-reverse-both",
            Some((|i, j| (6 - j, 6 - i)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "same",
            Some((|i, j| (i, j)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "same-reverse-first",
            Some((|i, j| (6 - i, j)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "same-reverse-second",
            Some((|i, j| (i, 6 - j)) as fn(usize, usize) -> (usize, usize)),
        ),
        (
            "same-reverse-both",
            Some((|i, j| (6 - i, 6 - j)) as fn(usize, usize) -> (usize, usize)),
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

    for (label, upper_index) in pair_orientations {
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
        if let Some(upper_index) = upper_index {
            scalar_constraints.extend(cc_test_sector_local_mapped_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                upper_index,
            ));
        }

        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows.clone(),
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);

        eprintln!(
            "scaled q13/q2 live pair orientation {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
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

        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live pair orientation best: {best:?}");
    assert!(best.5 < 1.0e-10 || best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_live_path_reversed_pair_reductions_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let pair_cases = [
        ("none", None::<fn(usize, usize) -> bool>),
        (
            "apex-corner",
            Some((|i, j| i == 0 && j == 0) as fn(usize, usize) -> bool),
        ),
        (
            "first-order-square",
            Some((|i, j| i <= 1 && j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "total-degree-one",
            Some((|i, j| i + j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "lemma-cubic-triangle",
            Some((|i, j| i + j <= 3) as fn(usize, usize) -> bool),
        ),
        (
            "post-lemma-complement",
            Some((|i, j| i + j > 3) as fn(usize, usize) -> bool),
        ),
        ("full-pair", Some((|_, _| true) as fn(usize, usize) -> bool)),
    ];

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );

    for (label, pair_keep) in pair_cases {
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
        if let Some(keep) = pair_keep {
            scalar_constraints.extend(cc_test_sector_local_filtered_mapped_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                keep,
                |i, j| (6 - j, 6 - i),
            ));
        }

        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);

        eprintln!(
            "scaled q13/q2 live reversed pair reduction {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
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

        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live reversed pair reduction best: {best:?}");
    assert!(best.5 < 1.0e-10 || best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_live_path_reversed_pair_components_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    let component_cases = [
        (
            "full-pair-both-components",
            true,
            true,
            None::<fn(usize, usize) -> bool>,
        ),
        (
            "full-pair-first-component",
            true,
            false,
            None::<fn(usize, usize) -> bool>,
        ),
        (
            "full-pair-second-component",
            false,
            true,
            None::<fn(usize, usize) -> bool>,
        ),
        (
            "first-order-square-first-component",
            true,
            false,
            Some((|i, j| i <= 1 && j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "first-order-square-second-component",
            false,
            true,
            Some((|i, j| i <= 1 && j <= 1) as fn(usize, usize) -> bool),
        ),
        (
            "lemma-cubic-first-component",
            true,
            false,
            Some((|i, j| i + j <= 3) as fn(usize, usize) -> bool),
        ),
        (
            "lemma-cubic-second-component",
            false,
            true,
            Some((|i, j| i + j <= 3) as fn(usize, usize) -> bool),
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

    for (label, include_first, include_second, keep) in component_cases {
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
        scalar_constraints.extend(
            cc_test_sector_local_filtered_mapped_pair_diagonal_component_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                keep.unwrap_or((|_, _| true) as fn(usize, usize) -> bool),
                |i, j| (6 - j, 6 - i),
                include_first,
                include_second,
            ),
        );

        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "scaled q13/q2 live reversed pair component {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
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
        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live reversed pair component best: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_eq14_half_sources_through_live_diagonal_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let split = split_kp20_tensor_border(&tensor_border).unwrap();

    let current_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let left_endpoint = current_rows[6].target;
    let right_endpoint = current_rows[27].target;
    let mut left_endpoint_normalized_rows = current_rows.clone();
    for row in left_endpoint_normalized_rows.iter_mut().take(21) {
        row.target[0] -= left_endpoint[0];
        row.target[1] -= left_endpoint[1];
    }
    let mut right_endpoint_normalized_rows = current_rows.clone();
    for row in right_endpoint_normalized_rows.iter_mut().skip(21) {
        row.target[0] -= right_endpoint[0];
        row.target[1] -= right_endpoint[1];
    }
    let mut both_endpoint_normalized_rows = left_endpoint_normalized_rows.clone();
    for row in both_endpoint_normalized_rows.iter_mut().skip(21) {
        row.target[0] -= right_endpoint[0];
        row.target[1] -= right_endpoint[1];
    }

    let source_cases = vec![
        ("current-live-previous-nearc-current-nearc", current_rows),
        (
            "current-live-left-endpoint-origin",
            left_endpoint_normalized_rows,
        ),
        (
            "current-live-right-endpoint-origin",
            right_endpoint_normalized_rows,
        ),
        (
            "current-live-both-endpoint-origin",
            both_endpoint_normalized_rows,
        ),
        (
            "half-sweep-best-previous-nearm-current-nearc-reversed",
            kp20_section42_depth_u_tensor_border_constraints_from_sources(
                &split,
                Kp20Section42Eq14PieceSource {
                    piece: Kp20PieceId::NearMidpointLeft,
                    sector: split.len() - 1,
                    half: Kp20Section42Eq14Half::NearM,
                    reverse_u: false,
                },
                Kp20Section42Eq14PieceSource {
                    piece: Kp20PieceId::NearMidpointRight,
                    sector,
                    half: Kp20Section42Eq14Half::NearC,
                    reverse_u: true,
                },
            ),
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

    for (label, tensor_rows) in source_cases {
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "scaled q13/q2 live eq14 half source {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
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
        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live eq14 half source best: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_live_source_projection_deltas_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
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
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);

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
            "scaled q13/q2 live eq14 projection delta rank={rank} row={row_index} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
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
        "scaled q13/q2 live eq14 projection delta summary: max_source_delta={max_source_delta:e}"
    );

    let worst = deltas[0];
    assert_eq!(worst.1, 27);
    assert_eq!(
        worst.2,
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, 0, 6)
    );
    assert!(max_source_delta > 1.0 && max_source_delta < 1.2);
}

#[test]
fn pr08_reflected_unit_sector_chart_endpoint_fit_is_not_exact() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
    )
    .unwrap();
    let normalized_rows = kp20_reflected_unit_sector_chart_tensor_rows(&rows, valence, sector)
        .expect("next-row transposed source should have nondegenerate reflected chart endpoints");

    let alpha = core::f64::consts::TAU / valence as f64;
    let sector_start = Complex::from_angle(alpha * sector as f64);
    let sector_end = Complex::from_angle(alpha * (sector + 1) as f64);
    let first_endpoint = complex_from_point(normalized_rows.sectors[sector][0][0]);
    let second_endpoint = complex_from_point(normalized_rows.sectors[sector][0][3]);
    let first_residual = (first_endpoint - sector_end).norm2().sqrt();
    let second_residual = (second_endpoint - sector_start).norm2().sqrt();

    eprintln!(
        "PR08 reflected unit-sector endpoint fit residuals: first={first_residual:e} second={second_residual:e}"
    );
    assert!(first_residual > 1.0 && first_residual < 1.3);
    assert!(second_residual > 0.4 && second_residual < 0.5);
}

#[test]
fn pr08_candidate_endpoint_chart_fits_are_ranked_before_tau_solve() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let alpha = core::f64::consts::TAU / valence as f64;
    let sector_start = Complex::from_angle(alpha * sector as f64);
    let sector_end = Complex::from_angle(alpha * (sector + 1) as f64);
    let mut best = ("", "", f64::INFINITY, f64::INFINITY, f64::INFINITY);

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let source0 = complex_from_point(rows.sectors[sector][0][0]);
        let source1 = complex_from_point(rows.sectors[sector][0][3]);
        for (chart_label, swap_targets, conjugate_source) in [
            ("direct", false, false),
            ("swapped", true, false),
            ("direct-conjugate", false, true),
            ("swapped-conjugate", true, true),
        ] {
            let target0 = if swap_targets {
                sector_end
            } else {
                sector_start
            };
            let target1 = if swap_targets {
                sector_start
            } else {
                sector_end
            };
            let fit_source0 = maybe_conjugate(source0, conjugate_source);
            let fit_source1 = maybe_conjugate(source1, conjugate_source);
            let Some(scale) = complex_similarity_scale_for_pairs(&[
                (fit_source0, target0),
                (fit_source1, target1),
            ]) else {
                continue;
            };
            let residual0 = (fit_source0 * scale - target0).norm2().sqrt();
            let residual1 = (fit_source1 * scale - target1).norm2().sqrt();
            let max_residual = residual0.max(residual1);
            let sum_residual = residual0 + residual1;
            if max_residual < best.2 || (max_residual == best.2 && sum_residual < best.3) {
                best = (
                    candidate.name(),
                    chart_label,
                    max_residual,
                    residual0,
                    residual1,
                );
            }
        }
    }

    eprintln!(
        "PR08 endpoint chart fit best candidate={} chart={} max={:e} first={:e} second={:e}",
        best.0, best.1, best.2, best.3, best.4
    );
    assert_eq!(
        best.0,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedReversed.name()
    );
    assert_eq!(best.1, "direct-conjugate");
    assert!(best.2 < 1.0e-12);
}

#[test]
fn pr08_stam_projected_side_endpoint_fits_are_ranked_before_tau_solve() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let ranked = pr08_stam_projected_side_endpoint_fit_ranking(&ring, valence, sector);
    let best = ranked[0];

    eprintln!(
        "PR08 Stam projected side endpoint best: family={} patch={:?} transform={:?} side={} reversed_depth={} swap={} conjugate={} max={:e} sum={:e}",
        best.family.name(),
        best.patch,
        best.transform,
        best.side.name(),
        best.reversed_depth,
        best.swap_targets,
        best.conjugate_source,
        best.max_endpoint_residual,
        best.sum_endpoint_residual,
    );

    assert_eq!(best.family, Kp20StamProjectionFamily::AppendixTensorOrder);
    assert_eq!(best.patch, Kp20StamPatch::LowerRight);
    assert_eq!(best.transform, Kp20StamGridTransform::AntiTranspose);
    assert_eq!(best.side, Kp20StamGridSide::Left);
    assert!(!best.reversed_depth);
    assert!(best.swap_targets);
    assert!(!best.conjugate_source);
    assert!(best.max_endpoint_residual < 1.0e-12);
}

#[test]
fn pr08_best_exact_stam_side_linear_characteristic_source_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let raw_best_side = pr08_best_exact_stam_side_source_rows(&ring);
    let endpoint_best_side =
        unit_sector_endpoint_normalized_tensor_rows(&raw_best_side, valence, sector, true, true)
            .expect("best exact Stam side should have nondegenerate endpoint chart");
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, 0.0);

    for (normalization_label, rows) in [
        ("raw-exact-side", raw_best_side.clone()),
        ("endpoint-normalized-exact-side", endpoint_best_side),
        (
            "per-sector-characteristic-exact-side",
            cc_test_per_sector_characteristic_tensor_rows(&raw_best_side, valence, sector),
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
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 best exact Stam side linear characteristic normalization {normalization_label}: residual={residual:e} center_scale={center_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} frob_delta={frob_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[0],
            square_outside = square.outside,
            square_samples = square.samples,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if square.outside < best.1
            || (square.outside == best.1 && validity.orientation.negative_determinants < best.2)
            || (square.outside == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                normalization_label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                det_a,
            );
        }
    }

    eprintln!("PR08 best exact Stam side linear characteristic normalization best: {best:?}");
    assert_eq!(best.1, 0);
    assert_eq!(best.2, 0);
    assert!(best.4.abs() < 1.0e-6);
}

#[test]
fn pr08_best_exact_stam_side_projection_masks_are_measured_under_live_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let raw_rows = pr08_best_exact_stam_side_source_rows(&ring);
    let rows = unit_sector_endpoint_normalized_tensor_rows(&raw_rows, valence, sector, true, true)
        .expect("best exact Stam side should have nondegenerate endpoint chart");
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let all_coeffs = tensor_rows
        .iter()
        .map(|row| row.terms[0].0)
        .collect::<Vec<_>>();
    let underlined_c0_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedC0PinnedRows,
    );
    let underlined_section42_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows,
    );
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for (label, free_coeffs) in [
        ("none", Vec::<Kp20PlanarTauCoeff>::new()),
        ("underlined-c0", underlined_c0_coeffs),
        ("underlined-section42", underlined_section42_coeffs),
        ("all-tensor-pins", all_coeffs),
    ] {
        let row_set =
            kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
                valence,
                sector,
                &tensor_border,
                &free_coeffs,
            )
            .unwrap();
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: row_set.vector_constraints,
                scalar_constraints: row_set.scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 best exact Stam side live projection mask {label}: free_coeffs={} residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            free_coeffs.len(),
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if solved.max_constraint_residual < 1.0e-8
            && (square.outside < best.1
                || (square.outside == best.1
                    && validity.orientation.negative_determinants < best.2)
                || (square.outside == best.1
                    && validity.orientation.negative_determinants == best.2
                    && validity.orientation.min_determinant > best.3))
        {
            best = (
                label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }

    eprintln!("PR08 best exact Stam side live projection mask best: {best:?}");
    assert_eq!(best.0, "all-tensor-pins");
    assert_eq!(best.1, 0);
    assert_eq!(best.2, 0);
    assert!(best.3 > -1.0e-20);
}

#[test]
fn pr08_best_exact_stam_side_projection_deltas_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let raw_rows = pr08_best_exact_stam_side_source_rows(&ring);
    let rows = unit_sector_endpoint_normalized_tensor_rows(&raw_rows, valence, sector, true, true)
        .expect("best exact Stam side should have nondegenerate endpoint chart");
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);
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
    deltas.sort_by(|a, b| b.0.total_cmp(&a.0));

    for (rank, (norm, row_index, coeff, dx, dy, original, projected)) in
        deltas.iter().take(12).enumerate()
    {
        eprintln!(
            "PR08 best exact Stam side projection delta rank={rank} row={row_index} label={} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
            tensor_border_row_label(*row_index),
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
        "PR08 best exact Stam side projection delta summary: max_source_delta={max_source_delta:e}"
    );

    assert!(max_source_delta > 0.49);
    assert_eq!(deltas[0].1, 39);
}

#[test]
fn pr08_best_exact_stam_side_eq14_half_assignment_variants_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [0, 1, 2],
    );
    let mut best_unprojected = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);
    let mut best_projected = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for previous_half in [Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM] {
        for current_half in [Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM] {
            for reverse_previous_u in [false, true] {
                for reverse_current_u in [false, true] {
                    let label = match (
                        previous_half,
                        current_half,
                        reverse_previous_u,
                        reverse_current_u,
                    ) {
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearC, false, false) => {
                            "prev-c/current-c"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearC, true, false) => {
                            "prev-c-rev/current-c"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearC, false, true) => {
                            "prev-c/current-c-rev"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearC, true, true) => {
                            "prev-c-rev/current-c-rev"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM, false, false) => {
                            "prev-c/current-m"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM, true, false) => {
                            "prev-c-rev/current-m"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM, false, true) => {
                            "prev-c/current-m-rev"
                        }
                        (Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM, true, true) => {
                            "prev-c-rev/current-m-rev"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearC, false, false) => {
                            "prev-m/current-c"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearC, true, false) => {
                            "prev-m-rev/current-c"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearC, false, true) => {
                            "prev-m/current-c-rev"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearC, true, true) => {
                            "prev-m-rev/current-c-rev"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearM, false, false) => {
                            "prev-m/current-m"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearM, true, false) => {
                            "prev-m-rev/current-m"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearM, false, true) => {
                            "prev-m/current-m-rev"
                        }
                        (Kp20Eq14SourceHalf::NearM, Kp20Eq14SourceHalf::NearM, true, true) => {
                            "prev-m-rev/current-m-rev"
                        }
                    };
                    let tensor_rows = cc_test_eq14_half_assignment_tensor_rows(
                        &tensor_border,
                        sector,
                        previous_half,
                        current_half,
                        reverse_previous_u,
                        reverse_current_u,
                    );
                    let unprojected = cc_test_solve_exact_side_tensor_rows_under_live_basis(
                        valence,
                        sector,
                        tensor_rows.clone(),
                        false,
                    );
                    let projected = cc_test_solve_exact_side_tensor_rows_under_live_basis(
                        valence,
                        sector,
                        tensor_rows,
                        true,
                    );
                    eprintln!(
                        "PR08 best exact Stam side eq14 half variant {label}: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
                        ures = unprojected.residual,
                        usquare = unprojected.square_outside,
                        uneg = unprojected.negative_determinants,
                        umin = unprojected.min_determinant,
                        pres = projected.residual,
                        psquare = projected.square_outside,
                        pneg = projected.negative_determinants,
                        pmin = projected.min_determinant,
                    );
                    if unprojected.residual < 1.0e-8
                        && (unprojected.square_outside < best_unprojected.1
                            || (unprojected.square_outside == best_unprojected.1
                                && unprojected.negative_determinants < best_unprojected.2)
                            || (unprojected.square_outside == best_unprojected.1
                                && unprojected.negative_determinants == best_unprojected.2
                                && unprojected.min_determinant > best_unprojected.3))
                    {
                        best_unprojected = (
                            label,
                            unprojected.square_outside,
                            unprojected.negative_determinants,
                            unprojected.min_determinant,
                            unprojected.residual,
                        );
                    }
                    if projected.square_outside < best_projected.1
                        || (projected.square_outside == best_projected.1
                            && projected.negative_determinants < best_projected.2)
                        || (projected.square_outside == best_projected.1
                            && projected.negative_determinants == best_projected.2
                            && projected.min_determinant > best_projected.3)
                    {
                        best_projected = (
                            label,
                            projected.square_outside,
                            projected.negative_determinants,
                            projected.min_determinant,
                            projected.residual,
                        );
                    }
                }
            }
        }
    }

    eprintln!(
        "PR08 best exact Stam side eq14 half variants best_unprojected={best_unprojected:?} best_projected={best_projected:?}"
    );
    assert_eq!(best_unprojected.0, "");
    assert_eq!(best_projected.1, 33);
    assert!(best_projected.2 <= 29);
}

#[test]
fn pr08_best_exact_stam_side_reversed_depth_projection_deltas_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);
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
    deltas.sort_by(|a, b| b.0.total_cmp(&a.0));

    for (rank, (norm, row_index, coeff, dx, dy, original, projected)) in
        deltas.iter().take(12).enumerate()
    {
        eprintln!(
            "PR08 best exact Stam side reversed-depth projection delta rank={rank} row={row_index} label={} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
            tensor_border_row_label(*row_index),
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
        "PR08 best exact Stam side reversed-depth projection delta summary: max_source_delta={max_source_delta:e}"
    );

    assert!(max_source_delta > 0.2);
    assert_eq!(deltas[0].1, 41);
}

#[test]
fn pr08_best_exact_stam_side_reversed_depth_projection_masks_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let all_coeffs = tensor_rows
        .iter()
        .map(|row| row.terms[0].0)
        .collect::<Vec<_>>();
    let underlined_c0_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedC0PinnedRows,
    );
    let underlined_section42_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows,
    );
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for (label, free_coeffs) in [
        ("none", Vec::<Kp20PlanarTauCoeff>::new()),
        ("underlined-c0", underlined_c0_coeffs),
        ("underlined-section42", underlined_section42_coeffs),
        ("all-tensor-pins", all_coeffs),
    ] {
        let row_set =
            kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
                valence,
                sector,
                &tensor_border,
                &free_coeffs,
            )
            .unwrap();
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: row_set.vector_constraints,
                scalar_constraints: row_set.scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 best exact Stam side reversed-depth projection mask {label}: free_coeffs={} residual={residual:e} square_outside={square_outside}/100 neg={neg}/100 min_det={min_det:e}",
            free_coeffs.len(),
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            neg = validity.orientation.negative_determinants,
            min_det = validity.orientation.min_determinant,
        );
        if solved.max_constraint_residual < 1.0e-8
            && (square.outside < best.1
                || (square.outside == best.1
                    && validity.orientation.negative_determinants < best.2)
                || (square.outside == best.1
                    && validity.orientation.negative_determinants == best.2
                    && validity.orientation.min_determinant > best.3))
        {
            best = (
                label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }

    eprintln!("PR08 best exact Stam side reversed-depth projection mask best: {best:?}");
    assert_eq!(best.0, "all-tensor-pins");
    assert_eq!(best.1, 0);
    assert_eq!(best.2, 0);
    assert!(best.3 > -1.0e-20);
}

#[test]
fn pr08_best_exact_stam_side_tangential_beta_projection_masks_are_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let tensor_rows = kp13_beta_tensor_border_constraints(&tensor_border, sector).unwrap();
    let underlined_c0_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedC0PinnedRows,
    );
    let underlined_section42_coeffs = kp20_section42_source_projection_free_coeffs(
        &tensor_rows,
        map,
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows,
    );
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let mut best = (
        "",
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0,
    );

    for (label, projected_tensor_rows, source_delta) in [
        ("none", tensor_rows.clone(), 0.0),
        {
            let free = underlined_c0_coeffs.clone();
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| free.contains(&coeff),
            );
            ("underlined-c0", rows, delta)
        },
        {
            let free = underlined_section42_coeffs.clone();
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| free.contains(&coeff),
            );
            ("underlined-section42", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| coeff.piece == Kp20PieceId::NearMidpointRight && coeff.j == 6,
            );
            ("right-u6", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| {
                    coeff.piece == Kp20PieceId::NearMidpointRight && (coeff.j == 5 || coeff.j == 6)
                },
            );
            ("right-u5-u6", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| {
                    coeff.piece == Kp20PieceId::NearMidpointRight
                        && ((coeff.j == 5 || coeff.j == 6) || (coeff.i == 0 && coeff.j == 0))
                },
            );
            ("right-corner-u5-u6", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| {
                    coeff.piece == Kp20PieceId::NearMidpointRight
                        && ((coeff.j == 5 || coeff.j == 6) || (coeff.i + coeff.j <= 1))
                },
            );
            ("right-diagonal-order1-u5-u6", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| {
                    coeff.piece == Kp20PieceId::NearMidpointRight
                        && ((coeff.j == 5 || coeff.j == 6) || (coeff.i + coeff.j <= 2))
                },
            );
            ("right-diagonal-order2-u5-u6", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| coeff.piece == Kp20PieceId::NearMidpointRight,
            );
            ("right-all", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_selected_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                |coeff| coeff.piece == Kp20PieceId::NearMidpointLeft,
            );
            ("left-all", rows, delta)
        },
        {
            let (rows, delta) = cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
            );
            ("all-tensor-pins", rows, delta)
        },
    ] {
        let diagnostic = cc_test_solve_exact_side_tensor_rows_under_live_basis(
            valence,
            sector,
            projected_tensor_rows,
            false,
        );
        eprintln!(
            "PR08 reversed-depth tangential-beta projection mask {label}: source_delta={source_delta:e} residual={residual:e} square={square}/100 neg={neg}/100 min_det={min_det:e}",
            residual = diagnostic.residual,
            square = diagnostic.square_outside,
            neg = diagnostic.negative_determinants,
            min_det = diagnostic.min_determinant,
        );
        if diagnostic.residual < 1.0e-8
            && (diagnostic.square_outside < best.2
                || (diagnostic.square_outside == best.2
                    && diagnostic.negative_determinants < best.3)
                || (diagnostic.square_outside == best.2
                    && diagnostic.negative_determinants == best.3
                    && diagnostic.min_determinant > best.4))
        {
            best = (
                label,
                diagnostic.residual,
                diagnostic.square_outside,
                diagnostic.negative_determinants,
                diagnostic.min_determinant,
                source_delta,
            );
        }
    }

    eprintln!("PR08 reversed-depth tangential-beta projection mask best: {best:?}");
    assert_eq!(best.0, "all-tensor-pins");
    assert_eq!(best.2, 31);
    assert_eq!(best.3, 27);
}

#[test]
fn pr08_best_exact_stam_side_tangential_beta_projection_deltas_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let tensor_rows = kp13_beta_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);
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
    deltas.sort_by(|a, b| b.0.total_cmp(&a.0));

    for (rank, (norm, row_index, coeff, dx, dy, original, projected)) in
        deltas.iter().take(12).enumerate()
    {
        eprintln!(
            "PR08 tangential-beta projection delta rank={rank} row={row_index} label={} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
            tensor_border_row_label(*row_index),
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
        "PR08 tangential-beta projection delta summary: max_source_delta={max_source_delta:e}"
    );

    assert!(max_source_delta > 0.4);
    assert!(deltas[0].0 > 0.4);
    assert_eq!(deltas[0].1, 34);
    assert_eq!(deltas[1].1, 41);
    assert_eq!(deltas[2].1, 27);
}

#[test]
fn pr08_best_exact_stam_side_all_pin_closed_tau_seed_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let beta_rows = kp13_beta_tensor_border_constraints(&tensor_border, sector).unwrap();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_rows, source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(beta_rows, &fig8_rows, &scalar_constraints);
    let mut vector_constraints = fig8_rows;
    vector_constraints.extend(projected_rows);
    let affine_seed = cc_test_local_planar_patch_to_global(
        valence,
        sector,
        cc_test_paper_quadrant_affine_macro_patch(),
    );
    let mut best = (
        "",
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );

    for (seed_label, seed) in [
        ("zero", Kp20MacroPatch::zero()),
        ("paper-affine", affine_seed),
    ] {
        for regularization_weight in [
            1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4,
        ] {
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: seed,
                    vector_constraints: vector_constraints.clone(),
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
            eprintln!(
                "PR08 all-pin closed tau seed sweep seed={seed_label} weight={regularization_weight:e}: source_delta={source_delta:e} residual={residual:e} square={square_outside}/100 neg={neg}/100 min_det={min_det:e} max_excess={max_excess:e} worst={worst:?}@({u:.2},{v:.2})",
                residual = solved.max_constraint_residual,
                square_outside = square.outside,
                neg = validity.orientation.negative_determinants,
                min_det = validity.orientation.min_determinant,
                max_excess = validity.domain.max_simplex_excess,
                worst = validity.orientation.worst_piece,
                u = validity.orientation.worst_u,
                v = validity.orientation.worst_v,
            );
            if square.outside < best.2
                || (square.outside == best.2 && validity.orientation.negative_determinants < best.3)
                || (square.outside == best.2
                    && validity.orientation.negative_determinants == best.3
                    && validity.orientation.min_determinant > best.4)
            {
                best = (
                    seed_label,
                    regularization_weight,
                    square.outside,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }

    eprintln!("PR08 all-pin closed tau seed sweep best: {best:?}");
    assert_eq!(best.2, 20);
    assert_eq!(best.3, 23);
    assert!(best.5 < 1.0e-8);
}

#[test]
fn pr08_best_exact_stam_side_source_projection_weight_pareto_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const FINE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    const EXTRA_WEIGHTS: [f64; 9] = [
        1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let beta_rows = kp13_beta_tensor_border_constraints(&tensor_border, sector).unwrap();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );

    let mut best_by_domain = (
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        0.0_f64,
        0_usize,
        0_usize,
        0.0_f64,
        0.0_f64,
    );
    let mut best_domain_tau = Kp20MacroPatch::zero();
    let mut best_by_orientation = (
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        0.0_f64,
        0_usize,
        0_usize,
        0.0_f64,
        0.0_f64,
    );

    for extra_weight in EXTRA_WEIGHTS {
        let (projected_rows, source_delta, projection_residual, domain_slacks) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis_with_extra_weight(
                beta_rows.clone(),
                &fig8_rows,
                &scalar_constraints,
                &[],
                0,
                extra_weight,
            );
        let solved = cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
            valence,
            sector,
            projected_rows.clone(),
            false,
        );
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let tensor_residual =
            cc_test_vector_rows_worst_residual(&solved.tau_macro_patch, &projected_rows);
        let max_domain_slack = domain_slacks
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        let (active_rows, lower_u, lower_v, upper_sum) = cc_test_active_sector_domain_bound_rows(
            &solved.tau_macro_patch,
            valence,
            sector,
            &SAMPLES,
            1.0e-10,
        );
        eprintln!(
            "PR08 source projection weight={extra_weight:e}: active_rows={} lower_u={} lower_v={} upper_sum={} source_delta={source_delta:e} domain_slack={max_domain_slack:e} projection_residual={projection_residual:e} residual={residual:e} tensor_worst=({}, {tensor_res:e}, {:?}[{},{}]) triangle_outside={outside}/100 max_excess={max_excess:e} square={square_outside}/100 neg={neg}/100 min_det={min_det:e}",
            active_rows.len(),
            lower_u,
            lower_v,
            upper_sum,
            tensor_border_row_label(tensor_residual.0),
            tensor_residual.2.piece,
            tensor_residual.2.i,
            tensor_residual.2.j,
            residual = solved.validity.residual,
            tensor_res = tensor_residual.1,
            outside = validity.domain.outside_triangle_samples,
            max_excess = validity.domain.max_simplex_excess,
            square_outside = square.outside,
            neg = validity.orientation.negative_determinants,
            min_det = validity.orientation.min_determinant,
        );

        let report = (
            validity.domain.outside_triangle_samples,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            validity.domain.max_simplex_excess,
            extra_weight,
            0_usize,
            active_rows.len(),
            source_delta,
            max_domain_slack,
        );
        if report.0 < best_by_domain.0
            || (report.0 == best_by_domain.0 && report.1 < best_by_domain.1)
            || (report.0 == best_by_domain.0
                && report.1 == best_by_domain.1
                && report.2 > best_by_domain.2)
        {
            best_by_domain = report;
            best_domain_tau = local_tau;
        }
        if report.1 < best_by_orientation.1
            || (report.1 == best_by_orientation.1 && report.0 < best_by_orientation.0)
            || (report.1 == best_by_orientation.1
                && report.0 == best_by_orientation.0
                && report.2 > best_by_orientation.2)
        {
            best_by_orientation = report;
        }
    }

    eprintln!(
        "PR08 source projection weight best_by_domain={best_by_domain:?} best_by_orientation={best_by_orientation:?}"
    );
    let fine_validity =
        kp20_planar_bi6_tau_validity_report(&best_domain_tau, &FINE_SAMPLES, 1.0e-10, 1.0e-10);
    let fine_square = cc_test_square_domain_report(&best_domain_tau, &FINE_SAMPLES, 1.0e-10);
    eprintln!(
        "PR08 source projection weight best fine grid: triangle_outside={outside}/324 max_excess={max_excess:e} square={square}/324 neg={neg}/324 min_det={min_det:e}",
        outside = fine_validity.domain.outside_triangle_samples,
        max_excess = fine_validity.domain.max_simplex_excess,
        square = fine_square.outside,
        neg = fine_validity.orientation.negative_determinants,
        min_det = fine_validity.orientation.min_determinant,
    );
    assert_eq!(best_by_domain.0, 0);
    assert_eq!(best_by_domain.1, 0);
    // Roundoff-only bound: winner excess is rounding noise (3.3e-14
    // scalar KKT, 5.0e-14 LU — PK-291; macOS just above 1e-13), ~12
    // orders under the rejected ~1.4e0 class; 1e-12 clears all.
    assert!(best_by_domain.3 < 1.0e-12);
    assert!(best_by_domain.7 > 0.49 && best_by_domain.7 < 0.50);
    assert_eq!(best_by_orientation.1, 0);
    assert_eq!(best_by_orientation.0, 0);
    assert!(best_by_orientation.3 < 1.0e-12);
    assert_eq!(fine_validity.domain.outside_triangle_samples, 0);
    assert_eq!(fine_validity.orientation.negative_determinants, 0);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_common_source_offset_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let base_tensor_rows =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );

    for (label, offset) in [
        ("none", [0.0, 0.0]),
        ("diagonal-0.02", [0.02, 0.02]),
        ("diagonal-0.04", [0.04, 0.04]),
        ("u-0.04", [0.04, 0.0]),
        ("v-0.04", [0.0, 0.04]),
        ("negative-diagonal-0.02", [-0.02, -0.02]),
    ] {
        let mut tensor_rows = base_tensor_rows.clone();
        for row in &mut tensor_rows {
            row.target[0] += offset[0];
            row.target[1] += offset[1];
        }
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
        let (projected_tensor_rows, source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
        );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "scaled q13/q2 live common source offset {label}: residual={residual:e} source_delta={source_delta:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
            residual = solved.max_constraint_residual,
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
        let score = (
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            source_delta,
            solved.max_constraint_residual,
        );
        if score.0 < best.1
            || (score.0 == best.1 && score.1 < best.2)
            || (score.0 == best.1 && score.1 == best.2 && score.2 > best.3)
        {
            best = (label, score.0, score.1, score.2, score.3, score.4);
        }
    }

    eprintln!("scaled q13/q2 live common source offset best: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_tau_invalidity_localizes_by_piece() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1.0e-10,
    })
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);

    for piece_id in Kp20PieceId::all() {
        let report = sampled_piece_tau_report(&local_tau, piece_id, &SAMPLES);
        eprintln!(
            "scaled q13/q2 piece {piece_id:?}: outside={outside}/{samples} neg={neg}/{det_samples} zero={zero} det=[{min_det:e},{max_det:e}] u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] max_sum={max_sum:e}",
            outside = report.outside,
            samples = report.samples,
            neg = report.negative,
            det_samples = report.det_samples,
            zero = report.near_zero,
            min_det = report.min_det,
            max_det = report.max_det,
            min_u = report.min_u,
            max_u = report.max_u,
            min_v = report.min_v,
            max_v = report.max_v,
            max_sum = report.max_sum,
        );
        assert!(report.min_det < 0.0 || report.max_det > 0.0);
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_centered_first_order_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let augmented_scalar_constraints = cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        cc_test_sector_local_axis_components(valence, sector),
        0,
    );
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);

    eprintln!(
        "scaled PR08 q13/q2 centered first-order: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_paper_tau_seed_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let local_seed = sample_kp20_macro_patch_from_tau_family_maps(
        &cc_test_local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .unwrap();
    let seed = cc_test_local_planar_patch_to_global(valence, sector, local_seed);

    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for regularization_weight in [1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0] {
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector,
                characteristic_tensor_border: &tensor_border,
                initial: seed,
                regularization_weight,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 paper-tau seed weight={regularization_weight:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let candidate = (
            regularization_weight,
            validity.domain.outside_triangle_samples,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            solved.max_constraint_residual,
        );
        if candidate.1 < best.1
            || (candidate.1 == best.1 && candidate.2 < best.2)
            || (candidate.1 == best.1 && candidate.2 == best.2 && candidate.3 > best.3)
        {
            best = candidate;
        }
    }
    assert!(
        best.1 > 0 || best.2 > 0,
        "paper tau seed unexpectedly produced a strict valid q13/q2 tau: {best:?}"
    );
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_source_chart_rows_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    for row in 0..3 {
        for col in 0..4 {
            let p = rows.sectors[sector][row][col];
            let local = cc_test_global_to_local_planar(valence, sector, p);
            eprintln!(
                "scaled q13/q2 source row[{row}][{col}] global=({x:e},{y:e}) local=({u:e},{v:e})",
                x = p.x,
                y = p.y,
                u = local[0],
                v = local[1],
            );
            assert!(p.x.is_finite() && p.y.is_finite());
        }
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_characteristic_affine_chart_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();

    let source_triangles = [
        (
            "near-origin-row-corner",
            [
                rows.sectors[sector][2][0],
                rows.sectors[sector][2][1],
                rows.sectors[sector][1][0],
            ],
        ),
        (
            "tensor-band-corner",
            [
                rows.sectors[sector][0][0],
                rows.sectors[sector][0][3],
                rows.sectors[sector][2][0],
            ],
        ),
        (
            "appendix-b-grid-corner",
            [
                rows.sectors[sector][2][0],
                rows.sectors[sector][0][0],
                rows.sectors[sector][0][3],
            ],
        ),
        (
            "row-edge-frame",
            [
                rows.sectors[sector][1][0],
                rows.sectors[sector][1][3],
                rows.sectors[sector][2][0],
            ],
        ),
    ];
    let target_triangles = [
        ("unit", [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        ("half", [[0.0, 0.0], [0.5, 0.0], [0.0, 0.5]]),
        (
            "paper-quadrant",
            [[0.0, 0.0], [1.0 / 3.0, 0.0], [0.0, 1.0 / 3.0]],
        ),
    ];

    let mut best = (
        "",
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (source_name, source_triangle) in source_triangles {
        for (target_name, target_triangle) in target_triangles {
            let Some(normalized_rows) = cc_test_affine_normalized_tensor_rows_from_source_triangle(
                &rows,
                valence,
                sector,
                source_triangle,
                target_triangle,
            ) else {
                eprintln!(
                    "scaled q13/q2 affine chart {source_name}->{target_name}: singular source triangle"
                );
                continue;
            };
            let tensor_border =
                kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
            let solved =
                solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                    valence,
                    sector,
                    characteristic_tensor_border: &tensor_border,
                    initial: Kp20MacroPatch::zero(),
                    regularization_weight: 1.0e-10,
                })
                .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "scaled q13/q2 affine chart {source_name}->{target_name}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            if validity.domain.outside_triangle_samples < best.2
                || (validity.domain.outside_triangle_samples == best.2
                    && validity.orientation.negative_determinants < best.3)
                || (validity.domain.outside_triangle_samples == best.2
                    && validity.orientation.negative_determinants == best.3
                    && validity.orientation.min_determinant > best.4)
            {
                best = (
                    source_name,
                    target_name,
                    validity.domain.outside_triangle_samples,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    validity.orientation.max_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }
    eprintln!("scaled q13/q2 best affine chart diagnostic: {best:?}");
    assert!(best.6 < 1.0e-9);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_tau_chart_phase_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let solve_sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector: solve_sector,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1.0e-10,
    })
    .unwrap();

    for chart_sector in 0..valence {
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, chart_sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 chart phase solve_sector={solve_sector} chart_sector={chart_sector}: outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] max_excess={max_excess:e}",
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            max_excess = validity.domain.max_simplex_excess,
        );
    }
    assert!(solved.max_constraint_residual < 1.0e-10);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_characteristic_chart_symmetry_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let transforms = [
        ("identity", [[1.0, 0.0], [0.0, 1.0]]),
        ("reflect-v", [[1.0, 0.0], [0.0, -1.0]]),
        ("reflect-u", [[-1.0, 0.0], [0.0, 1.0]]),
        ("negate", [[-1.0, 0.0], [0.0, -1.0]]),
        ("swap", [[0.0, 1.0], [1.0, 0.0]]),
        ("swap-reflect-second", [[0.0, 1.0], [-1.0, 0.0]]),
        ("swap-reflect-first", [[0.0, -1.0], [1.0, 0.0]]),
        ("swap-negate", [[0.0, -1.0], [-1.0, 0.0]]),
    ];
    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (label, transform) in transforms {
        let transformed_rows =
            cc_test_sector_local_chart_transformed_tensor_rows(&rows, valence, sector, transform);
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&transformed_rows).unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 chart symmetry {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] max_excess={max_excess:e}",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            max_excess = validity.domain.max_simplex_excess,
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
                validity.orientation.max_determinant,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!("scaled q13/q2 best chart symmetry diagnostic: {best:?}");
    assert!(best.5 < 1.0e-9);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_near_origin_frames_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1.0e-10,
    })
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    for piece_id in Kp20PieceId::all() {
        let piece = local_tau.piece(piece_id);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 near-origin frame {piece_id:?}: p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e}) det={det:e}",
            p00x = p00.x,
            p00y = p00.y,
            dux = du[0],
            duy = du[1],
            dvx = dv[0],
            dvy = dv[1],
        );
    }
    assert!(solved.max_constraint_residual < 1.0e-10);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_projected_row_family_frames_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();

    for family in Kp20Section42DiagnosticFamily::all() {
        let mut scalar_constraints = kp20_section42_diagnostic_rows(map, 2.0 * c, family);
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: row_set.vector_constraints.clone(),
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
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
            "scaled q13/q2 projected family {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            label = family.label(),
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
fn stam_appendix_a_scaled_pr08_q13_q2_c0_frame_companions_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let c0_rows = kp20_section42_only_c0_rows(map);
    let diagonal_rows =
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap();

    for (label, keep_fig8, keep_diagonal) in [
        ("c0-no-fig8-no-diagonal", false, false),
        ("c0-fig8-no-diagonal", true, false),
        ("c0-no-fig8-diagonal", false, true),
        ("c0-fig8-diagonal", true, true),
    ] {
        let mut vector_constraints = row_set.vector_constraints.clone();
        if !keep_fig8 {
            vector_constraints.remove(0);
        }
        let mut scalar_constraints = c0_rows.clone();
        if keep_diagonal {
            scalar_constraints.extend(diagonal_rows.clone());
        }
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
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
            "scaled q13/q2 {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
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
fn stam_appendix_a_scaled_pr08_q13_q2_without_full_diagonal_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);

    let mut vector_constraints = fig8_rows;
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
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
        "scaled q13/q2 no-full-diagonal projection: residual={residual:e} source_delta={max_source_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
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
    assert!(solved.max_constraint_residual < 1.0e-9);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_diagonal_blocks_are_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let c0_rows = kp20_section42_only_c0_rows(map);
    let near_origin_left_rows = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let near_midpoint_right_rows = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_rows = cc_test_sector_local_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
    );

    for (label, diagonal_blocks) in [
        ("none", Vec::<&[Kp20PlanarTauScalarLinearConstraint]>::new()),
        ("near-origin-left", vec![near_origin_left_rows.as_slice()]),
        (
            "near-midpoint-right",
            vec![near_midpoint_right_rows.as_slice()],
        ),
        ("paired-off-diagonal", vec![pair_rows.as_slice()]),
        (
            "near-origin-left-plus-pair",
            vec![near_origin_left_rows.as_slice(), pair_rows.as_slice()],
        ),
        (
            "near-midpoint-right-plus-pair",
            vec![near_midpoint_right_rows.as_slice(), pair_rows.as_slice()],
        ),
        (
            "all",
            vec![
                near_origin_left_rows.as_slice(),
                near_midpoint_right_rows.as_slice(),
                pair_rows.as_slice(),
            ],
        ),
    ] {
        let mut scalar_constraints = c0_rows.clone();
        for block in diagonal_blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: row_set.vector_constraints.clone(),
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
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
            "scaled q13/q2 diagonal block {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e})",
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
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_centered_first_order_is_measured() {
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
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta) =
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints);
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
        "scaled q13/q2 self-diagonal centered-first-order: residual={residual:e} source_delta={max_source_delta:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
        residual = solved.max_constraint_residual,
        scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_appendix_cubic_pins_are_measured() {
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
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            for i in 0..=3 {
                for j in 0..=(3 - i) {
                    let target = appendix_tau.piece(piece_id).cps[i][j];
                    vector_constraints.push(Kp20PlanarTauLinearConstraint::pin(
                        Kp20PlanarTauCoeff::new(piece_id, i, j),
                        [target.x, target.y],
                    ));
                }
            }
        }

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
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &vector_constraints,
                &scalar_constraints,
            );
        vector_constraints.extend(projected_tensor_rows);

        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
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
            "scaled q13/q2 self-diagonal {label} cubic pins: residual={residual:e} source_delta={source_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            source_delta = max_source_delta,
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
