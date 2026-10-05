#[allow(unused_imports)]
use super::super::super::kp_polar_shared::weighted_sum;
use super::*;

#[test]
fn stam_appendix_a_rotational_top_vector_corrected_eq12_endpoint_frame_is_measured() {
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
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let corrected_plain_endpoint = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::UnderlinedP,
        Kp20Section42PaperSymbol::P,
        Kp20Section42Eq12Split::PlainEndpoint,
    );
    let corrected_plain_plus = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::UnderlinedP,
        Kp20Section42PaperSymbol::P,
        Kp20Section42Eq12Split::PlainPlusRow,
    );
    let corrected_plain_minus = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::UnderlinedP,
        Kp20Section42PaperSymbol::P,
        Kp20Section42Eq12Split::PlainMinusRow,
    );
    let corrected_arrow_endpoint = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::ArrowUnderlinedP,
        Kp20Section42PaperSymbol::ArrowP,
        Kp20Section42Eq12Split::ArrowEndpoint,
    );
    let corrected_arrow_plus = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::ArrowUnderlinedP,
        Kp20Section42PaperSymbol::ArrowP,
        Kp20Section42Eq12Split::ArrowPlusRow,
    );
    let corrected_arrow_minus = kp20_section42_eq12_split_rows_with_symbols(
        map,
        Kp20Section42PaperSymbol::ArrowUnderlinedP,
        Kp20Section42PaperSymbol::ArrowP,
        Kp20Section42Eq12Split::ArrowMinusRow,
    );

    let empty_rows: &[Kp20PlanarTauScalarLinearConstraint] = &[];
    let cases: [(&str, [&[Kp20PlanarTauScalarLinearConstraint]; 9]); 4] = [
        (
            "full-corrected-section42",
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                corrected_plain_endpoint.as_slice(),
                corrected_plain_plus.as_slice(),
                corrected_plain_minus.as_slice(),
                corrected_arrow_endpoint.as_slice(),
                corrected_arrow_plus.as_slice(),
                corrected_arrow_minus.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "without-plain-endpoint",
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                empty_rows,
                corrected_plain_plus.as_slice(),
                corrected_plain_minus.as_slice(),
                corrected_arrow_endpoint.as_slice(),
                corrected_arrow_plus.as_slice(),
                corrected_arrow_minus.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "without-arrow-endpoint",
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                corrected_plain_endpoint.as_slice(),
                corrected_plain_plus.as_slice(),
                corrected_plain_minus.as_slice(),
                empty_rows,
                corrected_arrow_plus.as_slice(),
                corrected_arrow_minus.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "without-both-endpoints",
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                empty_rows,
                corrected_plain_plus.as_slice(),
                corrected_plain_minus.as_slice(),
                empty_rows,
                corrected_arrow_plus.as_slice(),
                corrected_arrow_minus.as_slice(),
                eq13_rows.as_slice(),
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
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let origin_left = cc_test_endpoint_frame(&local_tau, Kp20PieceId::NearOriginLeft, 0, 0);
        let midpoint_right =
            cc_test_endpoint_frame(&local_tau, Kp20PieceId::NearMidpointRight, 6, 0);
        let endpoint_delta = origin_left.origin.distance_to(&midpoint_right.origin);
        eprintln!(
            "rotational top corrected eq12 endpoint frame {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] endpoint_delta={endpoint_delta:e} nol_origin=({nox:e},{noy:e}) nol_du=({nodux:e},{noduy:e}) nol_dv=({nodvx:e},{nodvy:e}) nol_det={nodet:e} nmr_origin=({nmx:e},{nmy:e}) nmr_du=({nmdux:e},{nmduy:e}) nmr_dv=({nmdvx:e},{nmdvy:e}) nmr_det={nmdet:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
            residual = solved.max_constraint_residual,
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            nox = origin_left.origin.x,
            noy = origin_left.origin.y,
            nodux = origin_left.du.x,
            noduy = origin_left.du.y,
            nodvx = origin_left.dv.x,
            nodvy = origin_left.dv.y,
            nodet = origin_left.determinant,
            nmx = midpoint_right.origin.x,
            nmy = midpoint_right.origin.y,
            nmdux = midpoint_right.du.x,
            nmduy = midpoint_right.du.y,
            nmdvx = midpoint_right.dv.x,
            nmdvy = midpoint_right.dv.y,
            nmdet = midpoint_right.determinant,
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
        );
        assert!(solved.max_constraint_residual < 1.0e-10);
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_endpoint_j_orientation_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut base_rows = Vec::new();
    base_rows.extend(kp20_section42_only_c0_rows(map));
    base_rows.extend(kp20_section42_only_eq11_rows(map));
    base_rows.extend(kp20_section42_only_eq13_rows(map, 2.0 * c));
    for split in [
        Kp20Section42Eq12Split::PlainPlusRow,
        Kp20Section42Eq12Split::PlainMinusRow,
    ] {
        base_rows.extend(kp20_section42_eq12_split_rows_with_symbols(
            map,
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            split,
        ));
    }
    for split in [
        Kp20Section42Eq12Split::ArrowPlusRow,
        Kp20Section42Eq12Split::ArrowMinusRow,
    ] {
        base_rows.extend(kp20_section42_eq12_split_rows_with_symbols(
            map,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            split,
        ));
    }

    let endpoint_rows = |upper: Kp20Section42PaperSymbol,
                         lower: Kp20Section42PaperSymbol,
                         reverse_upper_j: bool|
     -> Vec<Kp20PlanarTauScalarLinearConstraint> {
        let mut rows = Vec::new();
        for j in 0..=2 {
            let upper_j = if reverse_upper_j { 2 - j } else { j };
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, upper_j, -1.0),
                ],
            );
        }
        rows
    };
    let mut best = ("", usize::MAX, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for plain_reverse in [false, true] {
        for arrow_reverse in [false, true] {
            let mut scalar_constraints = base_rows.clone();
            scalar_constraints.extend(endpoint_rows(
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
                plain_reverse,
            ));
            scalar_constraints.extend(endpoint_rows(
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
                arrow_reverse,
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
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let orientation =
                kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
            let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
            let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
            let label = match (plain_reverse, arrow_reverse) {
                (false, false) => "plain-same/arrow-same",
                (true, false) => "plain-reversed/arrow-same",
                (false, true) => "plain-same/arrow-reversed",
                (true, true) => "plain-reversed/arrow-reversed",
            };
            eprintln!(
                "rotational top corrected eq12 endpoint j-orientation {label}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e} nmr_det={nmr_det:e}",
                residual = solved.max_constraint_residual,
                neg = orientation.negative_determinants,
                samples = orientation.sample_count,
                worst = orientation.worst_piece,
                u = orientation.worst_u,
                v = orientation.worst_v,
                nml_det = nml.determinant,
                nmr_det = nmr.determinant,
            );
            assert!(solved.max_constraint_residual < 1.0e-10);
            if orientation.negative_determinants < best.1
                || (orientation.negative_determinants == best.1 && nml.determinant > best.2)
            {
                best = (
                    label,
                    orientation.negative_determinants,
                    nml.determinant,
                    nmr.determinant,
                );
            }
        }
    }

    assert_eq!(
        best.0, "plain-reversed/arrow-reversed",
        "reversing both endpoint j pairings should be the best endpoint-only orientation diagnostic, got {best:?}"
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_coherent_j_orientation_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut base_rows = Vec::new();
    base_rows.extend(kp20_section42_only_c0_rows(map));
    base_rows.extend(kp20_section42_only_eq11_rows(map));
    base_rows.extend(kp20_section42_only_eq13_rows(map, 2.0 * c));

    let coherent_eq12_rows = |upper: Kp20Section42PaperSymbol,
                              lower: Kp20Section42PaperSymbol,
                              reverse_upper_j: bool|
     -> Vec<Kp20PlanarTauScalarLinearConstraint> {
        let mut rows = Vec::new();
        for j in 0..=2 {
            let upper_j = if reverse_upper_j { 2 - j } else { j };
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, upper_j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            let upper_j = if reverse_upper_j { 2 - j } else { j };
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, upper_j, -1.0),
                    section42_term(upper, 4, upper_j, -0.25),
                    section42_term(upper, 2, upper_j, 0.25),
                ],
            );
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, upper_j, -1.0),
                    section42_term(upper, 4, upper_j, 0.25),
                    section42_term(upper, 2, upper_j, -0.25),
                ],
            );
        }
        rows
    };

    let mut best = ("", usize::MAX, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for plain_reverse in [false, true] {
        for arrow_reverse in [false, true] {
            let mut scalar_constraints = base_rows.clone();
            scalar_constraints.extend(coherent_eq12_rows(
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
                plain_reverse,
            ));
            scalar_constraints.extend(coherent_eq12_rows(
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
                arrow_reverse,
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
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let orientation =
                kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
            let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
            let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
            let label = match (plain_reverse, arrow_reverse) {
                (false, false) => "plain-same/arrow-same",
                (true, false) => "plain-reversed/arrow-same",
                (false, true) => "plain-same/arrow-reversed",
                (true, true) => "plain-reversed/arrow-reversed",
            };
            eprintln!(
                "rotational top corrected eq12 coherent j-orientation {label}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e} nmr_det={nmr_det:e}",
                residual = solved.max_constraint_residual,
                neg = orientation.negative_determinants,
                samples = orientation.sample_count,
                worst = orientation.worst_piece,
                u = orientation.worst_u,
                v = orientation.worst_v,
                nml_det = nml.determinant,
                nmr_det = nmr.determinant,
            );
            assert!(solved.max_constraint_residual < 1.0e-10);
            if orientation.negative_determinants < best.1
                || (orientation.negative_determinants == best.1 && nml.determinant > best.2)
            {
                best = (
                    label,
                    orientation.negative_determinants,
                    nml.determinant,
                    nmr.determinant,
                );
            }
        }
    }

    assert_eq!(
        best.0, "plain-reversed/arrow-reversed",
        "coherent eq12 j-orientation reversal should be the best measured endpoint orientation diagnostic, got {best:?}"
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_outer_corner_constraints_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let vector_constraints = cc_test_fig8_indicated_average_variant_rows(false, false);
    let augmented_scalar_constraints =
        stam_appendix_a_rotational_top_vector_tensor_border_constraints(&ring, 0, 0, false, -1.0);
    let corner_coeffs = [
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 6, 6),
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 5, 6),
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 6, 5),
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, 6, 6),
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, 5, 6),
        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, 6, 5),
    ];

    for coeff in corner_coeffs {
        let vector_count = vector_constraints
            .iter()
            .filter(|row| row.terms.iter().any(|(row_coeff, _)| *row_coeff == coeff))
            .count();
        let scalar_count = scalar_constraints
            .iter()
            .filter(|row| {
                row.terms
                    .iter()
                    .any(|(row_coeff, _)| row_coeff.coeff == coeff)
            })
            .count();
        let augmented_count = augmented_scalar_constraints
            .iter()
            .filter(|row| {
                row.terms.iter().any(|(row_coeff, _)| {
                    matches!(
                        row_coeff,
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(tau_coeff)
                            if tau_coeff.coeff == coeff
                    )
                })
            })
            .count();
        eprintln!(
            "rotational top corrected §4.2 outer-corner constraint ownership {coeff:?}: vector={vector_count} scalar={scalar_count} augmented={augmented_count}"
        );
        assert_eq!(
            vector_count + scalar_count + augmented_count,
            0,
            "near-midpoint outer-corner coefficient should be an unconstrained F3 completion DOF, got vector={vector_count} scalar={scalar_count} augmented={augmented_count} for {coeff:?}"
        );
    }
}

#[test]
fn stam_appendix_a_source_span_under_corrected_section42_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let anchored_scale = ring.lambda.powi(4);

    let rotational_extra_count = 4;
    let mut rotational_source_rows =
        stam_appendix_a_rotational_top_vector_tensor_border_constraints(&ring, 0, 0, false, -1.0);
    rotational_source_rows.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        rotational_extra_count,
    ));
    rotational_source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(rotational_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    rotational_source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: 1.0,
    });
    let rotational = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: section42_rows.clone(),
            augmented_scalar_constraints: rotational_source_rows,
            extra_variable_count: rotational_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();

    let top_count = 2 * valence + 1;
    let unknown_extra_count = 2 * top_count;
    let mut unknown_source_rows =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);
    unknown_source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(1), 1.0)],
        target: ring.sectors[0][1].x,
    });
    let unknown = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: section42_rows.clone(),
            augmented_scalar_constraints: unknown_source_rows,
            extra_variable_count: unknown_extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();

    let mut unknown_first_order_source_rows =
        stam_appendix_a_unknown_top_vector_tensor_border_constraints(&ring, 0, 0, false);
    unknown_first_order_source_rows.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        unknown_extra_count,
    ));
    unknown_first_order_source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(unknown_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    let unknown_first_order = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: section42_rows,
            augmented_scalar_constraints: unknown_first_order_source_rows,
            extra_variable_count: unknown_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();

    for (label, solved) in [
        ("rotational-frequency-one", &rotational),
        ("unknown-top-anchor", &unknown),
        ("unknown-top-first-order", &unknown_first_order),
    ] {
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let coarse =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let dense = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let square =
            cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
        eprintln!(
            "corrected §4.2 source span {label}: residual={residual:e} coarse_neg={coarse_neg}/{coarse_samples} dense_neg={dense_neg}/{dense_samples} triangle_outside={triangle_outside}/{domain_samples} square_outside={square_outside}/{square_samples} square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) det=[{min_det:e},{max_det:e}] nml66={nml66:e} nml_det={nml_det:e} nmr_det={nmr_det:e}",
            residual = solved.max_constraint_residual,
            coarse_neg = coarse.orientation.negative_determinants,
            coarse_samples = coarse.orientation.sample_count,
            dense_neg = dense.orientation.negative_determinants,
            dense_samples = dense.orientation.sample_count,
            triangle_outside = dense.domain.outside_triangle_samples,
            domain_samples = dense.domain.sample_count,
            square_outside = square.outside,
            square_samples = square.samples,
            min_u = square.min_u,
            max_u = square.max_u,
            min_v = square.min_v,
            max_v = square.max_v,
            min_det = dense.orientation.min_determinant,
            max_det = dense.orientation.max_determinant,
            nml66 = nml.p66.x.hypot(nml.p66.y),
            nml_det = nml.determinant,
            nmr_det = nmr.determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-8,
            "{label} should close the corrected §4.2 algebra to be comparable"
        );
    }

    let unknown_first_order_odd = (0..valence)
        .map(|sector| {
            Complex::new(
                unknown_first_order.extra_values[1 + 2 * sector],
                unknown_first_order.extra_values[unknown_extra_count / 2 + 1 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let unknown_first_order_even = (0..valence)
        .map(|sector| {
            Complex::new(
                unknown_first_order.extra_values[2 + 2 * sector],
                unknown_first_order.extra_values[unknown_extra_count / 2 + 2 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let odd_magnitudes = fourier_mode_magnitudes(&unknown_first_order_odd);
    let even_magnitudes = fourier_mode_magnitudes(&unknown_first_order_even);
    eprintln!(
        "corrected §4.2 unknown first-order top-vector Fourier magnitudes: odd={odd_magnitudes:?} even={even_magnitudes:?}"
    );

    let unknown_first_order_local =
        cc_test_global_planar_patch_to_local(valence, 0, unknown_first_order.tau_macro_patch);
    let unknown_first_order_validity = kp20_planar_bi6_tau_validity_report(
        &unknown_first_order_local,
        &[0.0, 0.5, 1.0],
        1.0e-10,
        1.0e-10,
    );
    assert!(
        !unknown_first_order_validity.is_paper_valid_reparameterization(),
        "unknown source span with Lemma-4 first-order normalization unexpectedly became promotable; inspect the source construction before accepting"
    );
}

#[test]
fn stam_appendix_a_rotational_source_section42_map_square_domain_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        let mut source_rows = stam_appendix_a_rotational_top_vector_tensor_border_constraints(
            &ring, 0, 0, false, -1.0,
        );
        source_rows.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: section42_rows,
                augmented_scalar_constraints: source_rows,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let square =
            cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
        eprintln!(
            "rotational source §4.2 square-domain map {map_label}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
            triangle_outside = validity.domain.outside_triangle_samples,
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
        assert!(solved.max_constraint_residual.is_finite());
        if square.outside < best.1
            || (square.outside == best.1 && validity.orientation.negative_determinants < best.2)
            || (square.outside == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                map_label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }

    eprintln!(
        "rotational source §4.2 square-domain best map {} square_outside={} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4
    );
    assert!(
        best.1 > 0 || best.2 > 0 || best.3 <= 1.0e-10,
        "a §4.2 map unexpectedly produced a square-bounded, orientation-preserving rotational source tau"
    );
}

#[test]
fn stam_appendix_a_rotational_source_kp15_beta_half_strip_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut source_rows = stam_appendix_a_rotational_top_vector_kp15_beta_tensor_border_constraints(
        &ring, 0, 0, false, -1.0,
    );
    source_rows.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        source_extra_count,
    ));
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: 1.0,
    });
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
            augmented_scalar_constraints: source_rows,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let dense = kp20_planar_bi6_tau_validity_report(
        &local_tau,
        &[0.0, 0.25, 0.5, 0.75, 1.0],
        1.0e-10,
        1.0e-10,
    );
    let square = cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
    let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
    let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
    eprintln!(
        "rotational source KP15 β half-strip: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) nml_det={nml_det:e} nmr_det={nmr_det:e}",
        residual = solved.max_constraint_residual,
        square_outside = square.outside,
        square_samples = square.samples,
        triangle_outside = dense.domain.outside_triangle_samples,
        domain_samples = dense.domain.sample_count,
        neg = dense.orientation.negative_determinants,
        orientation_samples = dense.orientation.sample_count,
        zero = dense.orientation.near_zero_determinants,
        min_det = dense.orientation.min_determinant,
        max_det = dense.orientation.max_determinant,
        min_u = square.min_u,
        max_u = square.max_u,
        min_v = square.min_v,
        max_v = square.max_v,
        nml_det = nml.determinant,
        nmr_det = nmr.determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite(),
        "KP15 β half-strip diagnostic should produce finite algebra"
    );
}

#[test]
fn stam_appendix_a_rotational_source_kp15_beta_section42_map_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        let mut source_rows =
            stam_appendix_a_rotational_top_vector_kp15_beta_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        source_rows.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
                augmented_scalar_constraints: source_rows,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let square =
            cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
        eprintln!(
            "rotational source KP15 β §4.2 map {map_label}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
            residual = solved.max_constraint_residual,
            square_outside = square.outside,
            square_samples = square.samples,
            triangle_outside = validity.domain.outside_triangle_samples,
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
        assert!(solved.max_constraint_residual.is_finite());
        if square.outside < best.1
            || (square.outside == best.1 && validity.orientation.negative_determinants < best.2)
            || (square.outside == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                map_label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }

    eprintln!(
        "rotational source KP15 β §4.2 best map {} square_outside={} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4
    );
    assert!(
        best.1 > 0 || best.2 > 0 || best.3 <= 1.0e-10,
        "KP15 β source unexpectedly produced a square-bounded, orientation-preserving rotational tau"
    );
}

#[test]
fn stam_appendix_a_rotational_source_kp15_beta_full_then_split_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut source_rows =
        stam_appendix_a_rotational_top_vector_kp15_beta_full_then_split_tensor_border_constraints(
            &ring, 0, 0, false, -1.0,
        );
    source_rows.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        source_extra_count,
    ));
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: 1.0,
    });
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
            augmented_scalar_constraints: source_rows,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let dense = kp20_planar_bi6_tau_validity_report(
        &local_tau,
        &[0.0, 0.25, 0.5, 0.75, 1.0],
        1.0e-10,
        1.0e-10,
    );
    let square = cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
    let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
    let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
    eprintln!(
        "rotational source KP15 β full-then-split: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) nml_det={nml_det:e} nmr_det={nmr_det:e}",
        residual = solved.max_constraint_residual,
        square_outside = square.outside,
        square_samples = square.samples,
        triangle_outside = dense.domain.outside_triangle_samples,
        domain_samples = dense.domain.sample_count,
        neg = dense.orientation.negative_determinants,
        orientation_samples = dense.orientation.sample_count,
        zero = dense.orientation.near_zero_determinants,
        min_det = dense.orientation.min_determinant,
        max_det = dense.orientation.max_determinant,
        min_u = square.min_u,
        max_u = square.max_u,
        min_v = square.min_v,
        max_v = square.max_v,
        nml_det = nml.determinant,
        nmr_det = nmr.determinant,
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_rotational_source_kp13_beta_two_piece_row2_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut source_rows = stam_appendix_a_rotational_top_vector_kp13_beta_tensor_border_constraints(
        &ring, 0, 0, false, -1.0,
    );
    source_rows.extend(cc_test_first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes,
        source_extra_count,
    ));
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    });
    source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: 1.0,
    });
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
            scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
            augmented_scalar_constraints: source_rows,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let dense = kp20_planar_bi6_tau_validity_report(
        &local_tau,
        &[0.0, 0.25, 0.5, 0.75, 1.0],
        1.0e-10,
        1.0e-10,
    );
    let square = cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
    let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
    let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
    eprintln!(
        "rotational source KP13 β two-piece row2: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) nml_det={nml_det:e} nmr_det={nmr_det:e}",
        residual = solved.max_constraint_residual,
        square_outside = square.outside,
        square_samples = square.samples,
        triangle_outside = dense.domain.outside_triangle_samples,
        domain_samples = dense.domain.sample_count,
        neg = dense.orientation.negative_determinants,
        orientation_samples = dense.orientation.sample_count,
        zero = dense.orientation.near_zero_determinants,
        min_det = dense.orientation.min_determinant,
        max_det = dense.orientation.max_determinant,
        min_u = square.min_u,
        max_u = square.max_u,
        min_v = square.min_v,
        max_v = square.max_v,
        nml_det = nml.determinant,
        nmr_det = nmr.determinant,
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_rotational_source_kp13_beta_map_and_center_frame_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = (
        "",
        Kp20PieceId::NearOriginLeft,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for center_piece in Kp20PieceId::all() {
            let mut source_rows =
                stam_appendix_a_rotational_top_vector_kp13_beta_tensor_border_constraints(
                    &ring, 0, 0, false, -1.0,
                );
            source_rows.extend(cc_test_first_order_centered_source_rows(
                center_piece,
                local_axes,
                source_extra_count,
            ));
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                    1.0,
                )],
                target: anchored_scale,
            });
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: 1.0,
            });
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                    scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
                    augmented_scalar_constraints: source_rows,
                    extra_variable_count: source_extra_count + 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-12,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity = kp20_planar_bi6_tau_validity_report(
                &local_tau,
                &[0.0, 0.25, 0.5, 0.75, 1.0],
                1.0e-10,
                1.0e-10,
            );
            let square =
                cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
            eprintln!(
                "rotational source KP13 β map/frame map={map_label} center={center_piece:?}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
                residual = solved.max_constraint_residual,
                square_outside = square.outside,
                square_samples = square.samples,
                triangle_outside = validity.domain.outside_triangle_samples,
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
            assert!(solved.max_constraint_residual.is_finite());
            let algebra_closes = solved.max_constraint_residual < 1.0e-8;
            let nonsingular_samples = validity.orientation.near_zero_determinants == 0;
            if algebra_closes
                && nonsingular_samples
                && (square.outside < best.2
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants < best.3)
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants == best.3
                        && validity.orientation.min_determinant > best.4))
            {
                best = (
                    map_label,
                    center_piece,
                    square.outside,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }

    eprintln!(
        "rotational source KP13 β map/frame best map={} center={:?} square_outside={} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5
    );
    assert!(
        best.2 > 0 || best.3 > 0 || best.4 <= 1.0e-10,
        "KP13 β row2 map/frame sweep unexpectedly produced a square-bounded, orientation-preserving rotational tau"
    );
}

#[test]
fn stam_appendix_a_rotational_source_kp13_beta_exact_pullback_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = (
        "",
        Kp20PieceId::NearOriginLeft,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for center_piece in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            let mut source_rows =
                stam_appendix_a_rotational_top_vector_kp13_beta_exact_tensor_border_constraints(
                    &ring, 0, 0, false, -1.0,
                );
            source_rows.extend(cc_test_first_order_centered_source_rows(
                center_piece,
                local_axes,
                source_extra_count,
            ));
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                    1.0,
                )],
                target: anchored_scale,
            });
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: 1.0,
            });
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                    scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
                    augmented_scalar_constraints: source_rows,
                    extra_variable_count: source_extra_count + 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-12,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity = kp20_planar_bi6_tau_validity_report(
                &local_tau,
                &[0.0, 0.25, 0.5, 0.75, 1.0],
                1.0e-10,
                1.0e-10,
            );
            let square =
                cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
            eprintln!(
                "rotational source KP13 β exact pullback map={map_label} center={center_piece:?}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={triangle_outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
                residual = solved.max_constraint_residual,
                square_outside = square.outside,
                square_samples = square.samples,
                triangle_outside = validity.domain.outside_triangle_samples,
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
            assert!(solved.max_constraint_residual.is_finite());
            let algebra_closes = solved.max_constraint_residual < 1.0e-8;
            let nonsingular_samples = validity.orientation.near_zero_determinants == 0;
            if algebra_closes
                && nonsingular_samples
                && (square.outside < best.2
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants < best.3)
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants == best.3
                        && validity.orientation.min_determinant > best.4))
            {
                best = (
                    map_label,
                    center_piece,
                    square.outside,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }

    eprintln!(
        "rotational source KP13 β exact pullback best map={} center={:?} square_outside={} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5
    );
    assert!(
        best.2 > 0 || best.3 > 0 || best.4 <= 1.0e-10,
        "KP13 β exact pullback unexpectedly produced a square-bounded, orientation-preserving rotational tau"
    );
}

#[test]
fn kp13_beta_equation15_boundary_row_is_derived_and_c2() {
    let c = (core::f64::consts::TAU / 5.0).cos();
    let row = kp13_beta_equation15_boundary_row(c);
    let tol = 1.0e-14;

    assert_close_xy(
        row.middle[2],
        [1.0 / 3.0 + c / 216.0, 1.0 / 3.0 + c / 216.0],
        tol,
    );
    assert_close_xy(
        row.middle[3],
        [0.5 + 7.0 * c / 1296.0, 1.0 / 3.0 + c / 72.0],
        tol,
    );
    assert_close_xy(row.right[3], [1.0, 1.0 / 3.0 + c / 18.0], tol);

    assert_close_xy(row.right[0], row.middle[3], tol);
    assert_close_xy(
        sub_xy(row.right[1], row.right[0]),
        sub_xy(row.middle[3], row.middle[2]),
        tol,
    );
    assert_close_xy(
        add_xy(
            sub_xy(row.right[2], scale_xy(row.right[1], 2.0)),
            row.right[0],
        ),
        add_xy(
            sub_xy(row.middle[3], scale_xy(row.middle[2], 2.0)),
            row.middle[1],
        ),
        tol,
    );
    assert_close_xy(row.middle[2], [row.middle[2][1], row.middle[2][0]], tol);

    let (middle_e, middle_d, right_e, right_d) = kp13_beta_two_piece_row2_controls(c);
    for i in 0..4 {
        assert_close_scalar(middle_e[i], 6.0 * (row.middle[i][0] - i as f64 / 6.0), tol);
        assert_close_scalar(middle_d[i], 6.0 * (row.middle[i][1] - 1.0 / 3.0), tol);
        assert_close_scalar(
            right_e[i],
            6.0 * (row.right[i][0] - (0.5 + i as f64 / 6.0)),
            tol,
        );
        assert_close_scalar(right_d[i], 6.0 * (row.right[i][1] - 1.0 / 3.0), tol);
    }
}

#[test]
fn stam_appendix_a_rotational_source_kp13_beta_exact_phase_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let selectors = [
        Kp13BetaBoundaryHalfSelector::PreviousMiddle,
        Kp13BetaBoundaryHalfSelector::PreviousRight,
        Kp13BetaBoundaryHalfSelector::CurrentMiddle,
        Kp13BetaBoundaryHalfSelector::CurrentRight,
    ];
    let mut best = (
        Kp13BetaBoundaryHalfSelector::PreviousMiddle,
        Kp13BetaBoundaryHalfSelector::CurrentMiddle,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );

    for &left_selector in &selectors {
        for &right_selector in &selectors {
            let mut source_rows =
                stam_appendix_a_rotational_top_vector_kp13_beta_exact_phase_constraints(
                    &ring,
                    0,
                    0,
                    false,
                    -1.0,
                    left_selector,
                    right_selector,
                );
            source_rows.extend(cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                source_extra_count,
            ));
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                    1.0,
                )],
                target: anchored_scale,
            });
            source_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: 1.0,
            });
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                    scalar_constraints: kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
                    augmented_scalar_constraints: source_rows,
                    extra_variable_count: source_extra_count + 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-12,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity = kp20_planar_bi6_tau_validity_report(
                &local_tau,
                &[0.0, 0.25, 0.5, 0.75, 1.0],
                1.0e-10,
                1.0e-10,
            );
            let square =
                cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
            eprintln!(
                "rotational source KP13 β exact phase left={left:?} right={right:?}: residual={residual:e} square_outside={square_outside}/{square_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                left = left_selector,
                right = right_selector,
                residual = solved.max_constraint_residual,
                square_outside = square.outside,
                square_samples = square.samples,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(solved.max_constraint_residual.is_finite());
            let algebra_closes = solved.max_constraint_residual < 1.0e-8;
            let nonsingular_samples = validity.orientation.near_zero_determinants == 0;
            if algebra_closes
                && nonsingular_samples
                && (square.outside < best.2
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants < best.3)
                    || (square.outside == best.2
                        && validity.orientation.negative_determinants == best.3
                        && validity.orientation.min_determinant > best.4))
            {
                best = (
                    left_selector,
                    right_selector,
                    square.outside,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }

    eprintln!(
        "rotational source KP13 β exact phase best left={:?} right={:?} square_outside={} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5
    );
    assert_eq!(
        (best.0, best.1),
        (
            Kp13BetaBoundaryHalfSelector::PreviousMiddle,
            Kp13BetaBoundaryHalfSelector::CurrentMiddle
        ),
        "recalibrate if the paper beta boundary phase improves on the current previous/current middle assignment"
    );
}

#[test]
fn stam_appendix_a_kp13_beta_exact_source_localizes_by_section42_family() {
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
        ("source-only", vec![]),
        ("c0", vec![c0_rows.as_slice()]),
        ("c0-eq11", vec![c0_rows.as_slice(), eq11_rows.as_slice()]),
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
            "section42",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
    ];
    let mut full_invalid = None;

    for (label, blocks) in cases {
        let mut scalar_constraints = Vec::new();
        for block in blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_kp13_beta_exact_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
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
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
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
        let validity = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let square =
            cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
        eprintln!(
            "KP13 β exact source §4.2 family {label}: residual={residual:e} scale={scale:e} square_outside={outside}/{square_samples} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[source_extra_count],
            outside = square.outside,
            square_samples = square.samples,
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
        if label == "section42" {
            full_invalid = Some((
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
            ));
        }
    }

    let (outside, negative, min_det) = full_invalid.expect("section42 case should be measured");
    assert!(
        outside > 0 || negative > 0 || min_det <= 1.0e-10,
        "corrected KP13 beta exact source unexpectedly made full §4.2 tau valid"
    );
}

#[test]
fn stam_appendix_a_kp13_beta_exact_source_eq12_split_rows_are_localized() {
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
            stam_appendix_a_rotational_top_vector_kp13_beta_exact_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
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
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
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
        let validity = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let square =
            cc_test_square_domain_report(&local_tau, &[0.0, 0.25, 0.5, 0.75, 1.0], 1.0e-10);
        eprintln!(
            "KP13 β exact source eq12 split {label}: residual={residual:e} square_outside={outside}/{square_samples} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = square.outside,
            square_samples = square.samples,
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected KP13 beta eq12 split {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
    }
}

#[test]
fn section42_eq12_split_rows_reconstruct_full_row_builder() {
    let angle = core::f64::consts::TAU / 5.0;
    let map = kp20_section42_map_candidates(angle)[1].1;
    let mut split_rows = Vec::new();
    for split in Kp20Section42Eq12Split::ALL {
        split_rows.extend(kp20_section42_eq12_split_rows(map, split));
    }
    let mut split_rows = split_rows
        .iter()
        .map(|row| format!("{row:?}"))
        .collect::<Vec<_>>();
    let mut full_rows = kp20_section42_only_eq12_rows(map)
        .iter()
        .map(|row| format!("{row:?}"))
        .collect::<Vec<_>>();
    split_rows.sort();
    full_rows.sort();

    assert_eq!(
        split_rows, full_rows,
        "eq12 split diagnostics must use the same upper/lower symbol order as the full §4.2 row builder"
    );
}

#[test]
fn stam_appendix_a_source_span_central_frame_piece_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let anchored_scale = ring.lambda.powi(4);
    let frame_cases: [(&str, &[Kp20PieceId]); 3] = [
        ("near-origin-left", &[Kp20PieceId::NearOriginLeft]),
        ("near-midpoint-right", &[Kp20PieceId::NearMidpointRight]),
        (
            "both-self-diagonal",
            &[Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointRight],
        ),
    ];

    for (label, frame_pieces) in frame_cases {
        let source_extra_count = 4;
        let mut augmented_scalar_constraints =
            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                &ring, 0, 0, false, -1.0,
            );
        for &piece_id in frame_pieces {
            augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
                piece_id,
                local_axes,
                source_extra_count,
            ));
        }
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(
                Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                1.0,
            )],
            target: anchored_scale,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
            target: 1.0,
        });
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: section42_rows.clone(),
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let coarse =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let dense = kp20_planar_bi6_tau_validity_report(
            &local_tau,
            &[0.0, 0.25, 0.5, 0.75, 1.0],
            1.0e-10,
            1.0e-10,
        );
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
        eprintln!(
            "rotational source central-frame piece {label}: residual={residual:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) scale={scale:e} coarse_neg={coarse_neg}/{coarse_samples} dense_neg={dense_neg}/{dense_samples} outside={outside}/{domain_samples} det=[{min_det:e},{max_det:e}] nml_det={nml_det:e} nmr_det={nmr_det:e}",
            residual = solved.max_constraint_residual,
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            scale = solved.extra_values[source_extra_count],
            coarse_neg = coarse.orientation.negative_determinants,
            coarse_samples = coarse.orientation.sample_count,
            dense_neg = dense.orientation.negative_determinants,
            dense_samples = dense.orientation.sample_count,
            outside = dense.domain.outside_triangle_samples,
            domain_samples = dense.domain.sample_count,
            min_det = dense.orientation.min_determinant,
            max_det = dense.orientation.max_determinant,
            nml_det = nml.determinant,
            nmr_det = nmr.determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "central-frame piece sweep should produce finite diagnostics"
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_corrected_section42_chart_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let mut augmented_scalar_constraints =
        stam_appendix_a_rotational_top_vector_tensor_border_constraints(&ring, 0, 0, false, -1.0);
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
        usize::MAX,
        "",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );

    for chart_sector in 0..valence {
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, chart_sector, solved.tau_macro_patch);
        for (label, transform) in transforms {
            let transformed_tau = cc_test_transform_tau_chart(local_tau, transform);
            let validity =
                kp20_planar_bi6_tau_validity_report(&transformed_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "rotational corrected §4.2 chart sector={chart_sector} transform={label}: outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                max_excess = validity.domain.max_simplex_excess,
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
                    && validity.domain.max_simplex_excess < best.4)
                || (validity.domain.outside_triangle_samples == best.2
                    && validity.orientation.negative_determinants == best.3
                    && (validity.domain.max_simplex_excess - best.4).abs() < 1.0e-14
                    && validity.orientation.min_determinant > best.5)
            {
                best = (
                    chart_sector,
                    label,
                    validity.domain.outside_triangle_samples,
                    validity.orientation.negative_determinants,
                    validity.domain.max_simplex_excess,
                    validity.orientation.min_determinant,
                );
            }
        }
    }

    eprintln!(
        "rotational corrected §4.2 best chart diagnostic: sector={} transform={} outside={} neg={} max_excess={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5
    );
    assert!(
        best.2 > 0 || best.3 > 0,
        "a chart transform unexpectedly made the corrected rotational tau paper-valid"
    );
}

#[test]
fn stam_appendix_a_paper_quadrant_affine_skeleton_section42_residual_is_measured() {
    let valence = 5;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[1].1;
    let tau = cc_test_paper_quadrant_affine_macro_patch();
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let breakdown =
        kp20_section42_residual_breakdown(&tau, &c0_rows, &eq11_rows, &eq12_rows, &eq13_rows);
    let fig8 = max_kp20_vector_rows_residual(&tau, &fig8_rows);
    let triangle_validity =
        kp20_planar_bi6_tau_validity_report(&tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
    let square = cc_test_square_domain_report(&tau, &[0.0, 0.5, 1.0], 1.0e-10);

    eprintln!(
        "paper-quadrant affine skeleton: fig8={fig8:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} triangle_outside={triangle_outside}/{triangle_samples} square_outside={square_outside}/{square_samples} square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        c0 = breakdown.c0.0,
        eq11 = breakdown.eq11.0,
        eq12 = breakdown.eq12.0,
        eq13 = breakdown.eq13.0,
        triangle_outside = triangle_validity.domain.outside_triangle_samples,
        triangle_samples = triangle_validity.domain.sample_count,
        square_outside = square.outside,
        square_samples = square.samples,
        min_u = square.min_u,
        max_u = square.max_u,
        min_v = square.min_v,
        max_v = square.max_v,
        neg = triangle_validity.orientation.negative_determinants,
        orientation_samples = triangle_validity.orientation.sample_count,
        zero = triangle_validity.orientation.near_zero_determinants,
        min_det = triangle_validity.orientation.min_determinant,
        max_det = triangle_validity.orientation.max_determinant,
    );
    assert_eq!(square.outside, 0);
    assert_eq!(triangle_validity.orientation.negative_determinants, 0);
    assert_eq!(triangle_validity.orientation.near_zero_determinants, 0);
    assert_eq!(
        triangle_validity.domain.outside_triangle_samples, 14,
        "the Fig.8 paper-quadrant skeleton is a unit-square/quadrant layout, not a unit-simplex layout"
    );
    assert!(
        breakdown.eq12.0 > 0.1 || breakdown.eq13.0 > 0.1,
        "if the paper-quadrant affine skeleton closes §4.2, promote it as a domain anchor candidate"
    );
}

#[test]
fn stam_appendix_a_paper_quadrant_affine_skeleton_sweeps_section42_symbol_maps() {
    let valence = 5;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tau = cc_test_paper_quadrant_affine_macro_patch();
    let fig8 = max_kp20_vector_rows_residual(
        &tau,
        &kp20_planar_tau_fig8_indicated_averaging_constraints(),
    );
    let square = cc_test_square_domain_report(&tau, &[0.0, 0.5, 1.0], 1.0e-10);
    let mut best = (
        "",
        f64::INFINITY,
        Kp20Section42ResidualBreakdown {
            c0: (f64::INFINITY, 0),
            eq11: (f64::INFINITY, 0),
            eq12: (f64::INFINITY, 0),
            eq13: (f64::INFINITY, 0),
        },
    );

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        let c0_rows = kp20_section42_only_c0_rows(map);
        let eq11_rows = kp20_section42_only_eq11_rows(map);
        let eq12_rows = kp20_section42_only_eq12_rows(map);
        let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
        let breakdown =
            kp20_section42_residual_breakdown(&tau, &c0_rows, &eq11_rows, &eq12_rows, &eq13_rows);
        let worst_family = breakdown
            .c0
            .0
            .max(breakdown.eq11.0)
            .max(breakdown.eq12.0)
            .max(breakdown.eq13.0);
        eprintln!(
            "paper-quadrant affine skeleton §4.2 map {map_label}: fig8={fig8:e} square_outside={square_outside}/{square_samples} c0={c0:e}@{c0_row} eq11={eq11:e}@{eq11_row} eq12={eq12:e}@{eq12_row} eq13={eq13:e}@{eq13_row} worst={worst:e}",
            square_outside = square.outside,
            square_samples = square.samples,
            c0 = breakdown.c0.0,
            c0_row = breakdown.c0.1,
            eq11 = breakdown.eq11.0,
            eq11_row = breakdown.eq11.1,
            eq12 = breakdown.eq12.0,
            eq12_row = breakdown.eq12.1,
            eq13 = breakdown.eq13.0,
            eq13_row = breakdown.eq13.1,
            worst = worst_family,
        );
        if worst_family < best.1 {
            best = (map_label, worst_family, breakdown);
        }
    }

    eprintln!(
        "paper-quadrant affine skeleton §4.2 best map {} worst={:e} breakdown={:?}",
        best.0, best.1, best.2
    );
    assert_eq!(
        square.outside, 0,
        "Fig. 8 paper-quadrant affine skeleton must remain square-bounded"
    );
    assert!(
        fig8 < 1.0e-12,
        "Fig. 8 affine skeleton should satisfy the indicated averaging rows before §4.2 is interpreted"
    );
    assert!(
        best.1 > 0.1,
        "a candidate §4.2 symbol map unexpectedly closed the paper-quadrant affine skeleton; promote that map as the row-symbol anchor"
    );
}

#[test]
fn stam_appendix_a_unknown_top_vector_pair_first_order_rows_are_localized() {
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
    let first_order_terms = [
        ("corner-00", 0, 0, true, true),
        ("u-10-both-components", 1, 0, true, true),
        ("v-01-both-components", 0, 1, true, true),
        ("mixed-11-both-components", 1, 1, true, true),
        ("u-10-local0", 1, 0, true, false),
        ("u-10-local1", 1, 0, false, true),
        ("v-01-local0", 0, 1, true, false),
        ("v-01-local1", 0, 1, false, true),
        ("mixed-11-local0", 1, 1, true, false),
        ("mixed-11-local1", 1, 1, false, true),
    ];

    for (label, coeff_i, coeff_j, keep_local0, keep_local1) in first_order_terms {
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
        if keep_local0 {
            scalar_constraints.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, coeff_i, coeff_j),
                    local_axes[0],
                    Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, coeff_j, coeff_i),
                    local_axes[1],
                ),
            ));
        }
        if keep_local1 {
            scalar_constraints.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, coeff_i, coeff_j),
                    local_axes[1],
                    Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, coeff_j, coeff_i),
                    local_axes[0],
                ),
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
            "unknown top pair first-order row {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_pair_first_difference_rows_are_measured() {
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
    let derivative_rows = [
        ("none", false, false, false, false),
        ("corner-only", true, false, false, false),
        ("u10-delta-local0", true, true, false, false),
        ("u10-delta-local1", true, false, true, false),
        ("u10-delta-both", true, true, true, false),
        ("first-order-delta", true, true, true, true),
    ];

    for (label, keep_corner, keep_u_local0, keep_u_local1, keep_v_both) in derivative_rows {
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
        if keep_corner {
            scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 0 && j == 0,
            ));
        }
        if keep_u_local0 {
            scalar_constraints.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                        0,
                    ),
                    local_axes[0][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                        1,
                    ),
                    local_axes[0][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                        0,
                    ),
                    -local_axes[0][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                        1,
                    ),
                    -local_axes[0][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 1),
                        0,
                    ),
                    -local_axes[1][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 1),
                        1,
                    ),
                    -local_axes[1][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 0),
                        0,
                    ),
                    local_axes[1][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 0),
                        1,
                    ),
                    local_axes[1][1],
                ),
            ]));
        }
        if keep_u_local1 {
            scalar_constraints.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                        0,
                    ),
                    local_axes[1][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                        1,
                    ),
                    local_axes[1][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                        0,
                    ),
                    -local_axes[1][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                        1,
                    ),
                    -local_axes[1][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 1),
                        0,
                    ),
                    -local_axes[0][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 1),
                        1,
                    ),
                    -local_axes[0][1],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 0),
                        0,
                    ),
                    local_axes[0][0],
                ),
                (
                    Kp20PlanarTauScalarCoeff::new(
                        Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 0),
                        1,
                    ),
                    local_axes[0][1],
                ),
            ]));
        }
        if keep_v_both {
            scalar_constraints.extend(cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 0 && j == 1,
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
            "unknown top pair first-difference {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_pair_bad_derivative_sign_is_measured() {
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
    let lower_origin = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0);
    let lower_step = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0);
    let upper_origin = Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 0);
    let upper_step = Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 0, 1);

    for (label, upper_sign) in [("same-sign", 1.0), ("opposite-sign", -1.0)] {
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
            |i, j| i == 0 && j == 0,
        ));
        scalar_constraints.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
            cc_test_sector_local_pair_first_difference_terms(
                lower_step,
                lower_origin,
                local_axes[1],
                upper_step,
                upper_origin,
                local_axes[0],
                upper_sign,
            ),
        ));

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
            "unknown top pair bad derivative sign {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_composed_pr08_subdivision_source_is_one_lambda_level() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let composed = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let raw = stam_appendix_a_pr08_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let expected = scale_tensor_border_rows(&raw, ring.lambda);
    for sector in 0..ring.sectors.len() {
        for row in 0..3 {
            for col in 0..4 {
                assert_close_point(
                    composed.sectors[sector][row][col],
                    expected.sectors[sector][row][col],
                    1.0e-10,
                );
            }
        }
    }
}

#[test]
fn stam_appendix_a_composed_pr08_subdivision_source_residual_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
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
        solved.max_constraint_residual < 1.0e-10,
        "correct cross-degree tensor-border elevation should close the composed PR08 subdivision source rows, got {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_composed_pr08_subdivision_source_sweep_is_rejected() {
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
                let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
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
        "correct cross-degree tensor-border elevation should close the best composed PR08 subdivision source sweep, got {best:?}"
    );
    assert!(
        best.patch_index < 3,
        "best composed subdivision source patch should remain one of Stam's three regular subpatches"
    );
    assert_stam_appendix_a_pr08_inner_seed_diagnostic_in_search_space(best);
}

#[test]
fn stam_appendix_a_composed_source_full_degree8_eq9_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints
        .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(
        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8(
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
            kp20_b_coeffs_near_m(c),
            angle,
        ),
    );
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1e-10,
        })
        .unwrap();

    assert!(
        solved.max_constraint_residual > 10.6 && solved.max_constraint_residual < 10.8,
        "full-degree-8 near-midpoint eq9 should drop the composed source to the near-origin frontier, got {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_composed_source_full_degree8_eq9_at_both_ends_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints
        .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(
        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8(
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginLeft,
            [c, 2.0 * c, 2.0 * c],
            kp20_b_coeffs_near_o(c),
            angle,
        ),
    );
    scalar_constraints.extend(
        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8(
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
            kp20_b_coeffs_near_m(c),
            angle,
        ),
    );
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1e-10,
        })
        .unwrap();

    assert!(
        solved.max_constraint_residual > 10.6 && solved.max_constraint_residual < 10.8,
        "full-degree-8 eq9 at both ends should remain on the near-origin frontier, got {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_eq13_near_m_strip_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints
        .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_eq13_rotated_near_m_strip_assignment_rows(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c)[0],
        angle,
    ));
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1e-10,
        })
        .unwrap();

    assert!(
        solved.max_constraint_residual > 1.13 && solved.max_constraint_residual < 1.15,
        "KP20 eq(13) near-m strip diagnostic residual changed: {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_map_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();

    let candidates = kp20_section42_map_candidates(angle);

    let mut best = ("", "", "", f64::INFINITY);
    for (label, map) in candidates {
        for include_diagonal in [false, true] {
            let diagonal_label = if include_diagonal {
                "with-diagonal"
            } else {
                "without-diagonal"
            };
            for family in Kp20Section42DiagnosticFamily::all() {
                let mut vector_constraints = Vec::new();
                vector_constraints
                    .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
                let mut scalar_constraints = if include_diagonal {
                    kp20_planar_tau_diagonal_symmetry_constraints()
                } else {
                    Vec::new()
                };
                if family == Kp20Section42DiagnosticFamily::Full {
                    scalar_constraints.extend(kp20_section42_presolved_strip_basis_rows(
                        map,
                        kp20_b_coeffs_near_o(c)[0],
                    ));
                } else {
                    scalar_constraints.extend(kp20_section42_diagnostic_rows(
                        map,
                        kp20_b_coeffs_near_o(c)[0],
                        family,
                    ));
                }
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: Kp20MacroPatch::zero(),
                        vector_constraints,
                        scalar_constraints,
                        regularization_weight: 1e-10,
                    },
                )
                .unwrap();
                eprintln!(
                    "KP20 §4.2 map sweep {label}/{diagonal_label}/{}: residual {}",
                    family.label(),
                    solved.max_constraint_residual
                );
                if family == Kp20Section42DiagnosticFamily::Full
                    && solved.max_constraint_residual < best.3
                {
                    best = (
                        label,
                        diagonal_label,
                        family.label(),
                        solved.max_constraint_residual,
                    );
                }
            }
        }
    }

    eprintln!(
        "KP20 §4.2 map sweep best full basis {}/{}/{}: residual {}",
        best.0, best.1, best.2, best.3
    );
    assert!(
        best.3.is_finite(),
        "full KP20 §4.2 pre-solved strip basis must produce a finite solve"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_eq11_residual_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let with_diagonal =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], true);
    let without_diagonal =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], false);

    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, kp20_b_coeffs_near_o(c)[0]);
    let with_diag_breakdown = kp20_section42_residual_breakdown(
        &with_diagonal.tau_macro_patch,
        &c0_rows,
        &eq11_rows,
        &eq12_rows,
        &eq13_rows,
    );
    let without_diag_breakdown = kp20_section42_residual_breakdown(
        &without_diagonal.tau_macro_patch,
        &c0_rows,
        &eq11_rows,
        &eq12_rows,
        &eq13_rows,
    );

    eprintln!("KP20 §4.2 with diagonal: {with_diag_breakdown:?}");
    eprintln!("KP20 §4.2 without diagonal: {without_diag_breakdown:?}");
    assert!(
        without_diag_breakdown.eq11.0 < 1.0e-10,
        "eq(11) should be exactly satisfiable in the native §4.2 subsystem, got {:?}",
        without_diag_breakdown
    );
    assert!(
        without_diag_breakdown.c0.0 < 1.0e-10
            && without_diag_breakdown.eq12.0 < 1.0e-10
            && without_diag_breakdown.eq13.0 < 1.0e-10,
        "native §4.2 rows should be feasible after the eq(14) adjacent-half source assignment, got {:?}",
        without_diag_breakdown
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_c0_projection_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let maps = kp20_section42_map_candidates(angle);
    let mut best = Kp20Section42ProjectionSweepDiagnostic {
        residual: f64::INFINITY,
        map_label: "",
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
    };

    for patch_index in 0..3 {
        for row_start in [0, 1] {
            for reverse_rows in [false, true] {
                let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
                    &ring,
                    (12, odd_scale),
                    (1, even_scale),
                    patch_index,
                    row_start,
                    reverse_rows,
                )
                .unwrap();
                let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                for &(map_label, map) in &maps {
                    let residual = solve_kp20_section42_family_residual(
                        &tensor_border,
                        map,
                        0.0,
                        Kp20Section42DiagnosticFamily::C0Only,
                        false,
                    );
                    if residual < best.residual {
                        best = Kp20Section42ProjectionSweepDiagnostic {
                            residual,
                            map_label,
                            patch_index,
                            row_start,
                            reverse_rows,
                        };
                    }
                }
            }
        }
    }

    eprintln!("KP20 §4.2 C0 source projection sweep best: {best:?}");
    assert!(
        best.residual < 1.0e-10,
        "corrected eq(14) adjacent-half source assignment should make the composed Stam/PR08 C0 source feasible, got {best:?}"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_endpoint_projection_exposes_eq13_frontier() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let mut original_tensor_rows =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap();
    let projected_delta =
        project_kp20_section42_underlined_c0_endpoint_tensor_pins(&mut original_tensor_rows, map);
    let mut scalar_constraints =
        kp20_section42_presolved_strip_basis_rows(map, kp20_b_coeffs_near_o(c)[0]);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: original_tensor_rows,
            scalar_constraints: std::mem::take(&mut scalar_constraints),
            regularization_weight: 1e-10,
        })
        .unwrap();
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, kp20_b_coeffs_near_o(c)[0]);
    let breakdown = kp20_section42_residual_breakdown(
        &solved.tau_macro_patch,
        &c0_rows,
        &eq11_rows,
        &eq12_rows,
        &eq13_rows,
    );

    eprintln!(
        "KP20 §4.2 endpoint-projected solve: max adjusted endpoint by {projected_delta:e}, residual {}, breakdown {breakdown:?}",
        solved.max_constraint_residual
    );
    assert!(
        projected_delta < 1.0e-12,
        "the former endpoint projection should be a no-op after the eq(14) adjacent-half source assignment"
    );
    assert!(
        solved.max_constraint_residual < 1.0e-10
            && breakdown.c0.0 < 1.0e-10
            && breakdown.eq11.0 < 1.0e-10
            && breakdown.eq12.0 < 1.0e-10
            && breakdown.eq13.0 < 1.0e-10,
        "corrected source rows should close the full native §4.2 basis, got {breakdown:?}"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_source_projection_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let underlined_c0_projection = solve_kp20_section42_source_projection(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        Kp20Section42SourceProjectionScope::UnderlinedC0PinnedRows,
    );
    let full_tensor_projection = solve_kp20_section42_source_projection(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        Kp20Section42SourceProjectionScope::AllTensorPins,
    );
    let underlined_section42_projection = solve_kp20_section42_source_projection(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows,
    );

    eprintln!("KP20 §4.2 source projection underlined C0: {underlined_c0_projection:?}");
    eprintln!(
        "KP20 §4.2 source projection underlined §4.2 pins: {underlined_section42_projection:?}"
    );
    eprintln!("KP20 §4.2 source projection all tensor pins: {full_tensor_projection:?}");
    assert!(underlined_c0_projection.max_row_residual < 1.0e-8);
    assert!(underlined_c0_projection.max_source_delta < 0.03);
    assert!(underlined_section42_projection.max_row_residual < 1.0e-8);
    assert!(full_tensor_projection.max_row_residual < 1.0e-8);
    assert!(
        underlined_section42_projection.max_source_delta < full_tensor_projection.max_source_delta,
        "full tensor projection should remain the broadest source move after cross-degree row elevation"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_projected_tensor_rows_are_feasible() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let (projected_tensor_rows, projection_report) = project_kp20_section42_source_tensor_rows(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows,
    );
    let scalar_constraints =
        kp20_section42_presolved_strip_basis_rows(map, kp20_b_coeffs_near_o(c)[0]);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: projected_tensor_rows,
            scalar_constraints,
            regularization_weight: 1e-10,
        })
        .unwrap();

    eprintln!(
        "KP20 §4.2 projected tensor rows report {projection_report:?}, residual {}",
        solved.max_constraint_residual
    );
    assert!(projection_report.max_row_residual < 1.0e-8);
    assert!(solved.max_constraint_residual < 1.0e-8);
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_eq14_uv_orientation_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let split = split_kp20_tensor_border(&tensor_border).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let current_depth_u =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], false)
            .max_constraint_residual;
    let transposed_u_depth = solve_kp20_section42_rows_with_tensor_constraints(
        kp20_section42_uv_tensor_border_constraints_from_split(&split, 0).unwrap(),
        map,
        kp20_b_coeffs_near_o(c)[0],
        false,
    )
    .max_constraint_residual;

    eprintln!(
        "KP20 §4.2 eq(14) orientation depth-u residual {current_depth_u}, u-depth residual {transposed_u_depth}"
    );
    assert!(
        current_depth_u < transposed_u_depth,
        "simple transposition was expected to remain rejected by the §4.2 diagnostic; depth-u={current_depth_u}, u-depth={transposed_u_depth}"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_eq14_half_source_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let split = split_kp20_tensor_border(&tensor_border).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let mut best = Kp20Section42Eq14HalfSourceSweep {
        residual: f64::INFINITY,
        left_sector_offset: 0,
        left_half: Kp20Section42Eq14Half::NearC,
        left_reverse: false,
        right_sector_offset: 0,
        right_half: Kp20Section42Eq14Half::NearM,
        right_reverse: false,
    };

    for left_sector_offset in [0, split.len() - 1] {
        for right_sector_offset in [0, split.len() - 1] {
            for left_half in [Kp20Section42Eq14Half::NearC, Kp20Section42Eq14Half::NearM] {
                for right_half in [Kp20Section42Eq14Half::NearC, Kp20Section42Eq14Half::NearM] {
                    for left_reverse in [false, true] {
                        for right_reverse in [false, true] {
                            let tensor_constraints =
                                kp20_section42_depth_u_tensor_border_constraints_from_sources(
                                    &split,
                                    Kp20Section42Eq14PieceSource {
                                        piece: Kp20PieceId::NearMidpointLeft,
                                        sector: left_sector_offset,
                                        half: left_half,
                                        reverse_u: left_reverse,
                                    },
                                    Kp20Section42Eq14PieceSource {
                                        piece: Kp20PieceId::NearMidpointRight,
                                        sector: right_sector_offset,
                                        half: right_half,
                                        reverse_u: right_reverse,
                                    },
                                );
                            let residual = solve_kp20_section42_rows_with_tensor_constraints(
                                tensor_constraints,
                                map,
                                kp20_b_coeffs_near_o(c)[0],
                                false,
                            )
                            .max_constraint_residual;
                            if residual < best.residual {
                                best = Kp20Section42Eq14HalfSourceSweep {
                                    residual,
                                    left_sector_offset,
                                    left_half,
                                    left_reverse,
                                    right_sector_offset,
                                    right_half,
                                    right_reverse,
                                };
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("KP20 §4.2 eq(14) half-source sweep best: {best:?}");
    assert!(
        best.residual < 0.61,
        "half-source sweep should at least recover the current eq(14) assignment, got {best:?}"
    );
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_diagonal_conflict_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };

    let solved =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], false);
    let diagonal_rows = kp20_planar_tau_diagonal_symmetry_constraints();
    let worst = worst_kp20_scalar_row_residual(&solved.tau_macro_patch, &diagonal_rows);

    eprintln!(
        "KP20 §4.2 diagonal conflict worst row {} residual {} terms {:?}",
        worst.row_index, worst.residual, worst.terms
    );
    assert!(worst.residual > 0.98 && worst.residual < 0.99);
}

#[test]
fn stam_appendix_a_composed_source_kp20_section42_diagonal_subset_residuals_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let diagonal_rows = kp20_planar_tau_diagonal_symmetry_constraints();
    let diagonal_near_origin_left_rows = diagonal_rows[0..49].to_vec();
    let diagonal_near_midpoint_right_rows = diagonal_rows[49..98].to_vec();
    let diagonal_self_rows = diagonal_rows[0..98].to_vec();
    let diagonal_pair_rows = diagonal_rows[98..].to_vec();

    let no_diagonal =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], false);
    let self_diagonal = solve_kp20_section42_rows_with_extra_scalar_constraints(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        diagonal_self_rows,
    );
    let near_origin_left_diagonal = solve_kp20_section42_rows_with_extra_scalar_constraints(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        diagonal_near_origin_left_rows,
    );
    let near_midpoint_right_diagonal = solve_kp20_section42_rows_with_extra_scalar_constraints(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        diagonal_near_midpoint_right_rows,
    );
    let pair_diagonal = solve_kp20_section42_rows_with_extra_scalar_constraints(
        &tensor_border,
        map,
        kp20_b_coeffs_near_o(c)[0],
        diagonal_pair_rows,
    );
    let full_diagonal =
        solve_kp20_section42_native_rows(&tensor_border, map, kp20_b_coeffs_near_o(c)[0], true);

    eprintln!(
        "KP20 §4.2 diagonal subset residuals: none {}, self {}, near-origin-left {}, near-midpoint-right {}, pair {}, full {}",
        no_diagonal.max_constraint_residual,
        self_diagonal.max_constraint_residual,
        near_origin_left_diagonal.max_constraint_residual,
        near_midpoint_right_diagonal.max_constraint_residual,
        pair_diagonal.max_constraint_residual,
        full_diagonal.max_constraint_residual
    );
    assert!(no_diagonal.max_constraint_residual < 1.0e-10);
    assert!(self_diagonal.max_constraint_residual > pair_diagonal.max_constraint_residual);
}

#[test]
fn stam_appendix_a_composed_source_eq13_near_m_with_full_degree_origin_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints
        .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(
        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8(
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginLeft,
            [c, 2.0 * c, 2.0 * c],
            kp20_b_coeffs_near_o(c),
            angle,
        ),
    );
    scalar_constraints.extend(kp20_eq13_rotated_near_m_strip_assignment_rows(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c)[0],
        angle,
    ));
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1e-10,
        })
        .unwrap();

    assert!(
        solved.max_constraint_residual > 1.13 && solved.max_constraint_residual < 1.15,
        "KP20 eq(13) near-m strip with full-degree origin diagnostic changed: {}",
        solved.max_constraint_residual
    );
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_dyadic_level_sweep_closes_finite_level() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let base_rows = stam_appendix_a_pr08_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    let mut finite_level_best = StamAppendixADyadicLevelDiagnostic {
        residual: f64::INFINITY,
        level_exponent: 0,
    };
    for level_exponent in -1..=4 {
        let rows = scale_tensor_border_rows(&base_rows, ring.lambda.powi(level_exponent));
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: 5,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
        assert!(solved.max_constraint_residual < 1.0e-12, "{level_exponent}");
        if solved.max_constraint_residual < finite_level_best.residual {
            finite_level_best = StamAppendixADyadicLevelDiagnostic {
                residual: solved.max_constraint_residual,
                level_exponent,
            };
        }
    }
    assert!(
        finite_level_best.residual < 1.0e-10,
        "§4.2/§4.3 closure: {:e} at level {}",
        finite_level_best.residual,
        finite_level_best.level_exponent
    ); // Every level closes (1.9e-15..3.9e-14); the argmin is a last-bits tie (PK-676).

    let collapsed_rows = scale_tensor_border_rows(&base_rows, ring.lambda.powi(12));
    let collapsed_tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&collapsed_rows).unwrap();
    let collapsed =
        solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
            valence: 5,
            sector: 0,
            characteristic_tensor_border: &collapsed_tensor_border,
            initial: Kp20MacroPatch::zero(),
            regularization_weight: 1e-10,
        })
        .unwrap();
    assert!(
        collapsed.max_constraint_residual < 0.12,
        "very deep dyadic scaling should expose only collapse-to-zero feasibility under the accepted paper row set, got {}",
        collapsed.max_constraint_residual
    );
    assert!(
        tensor_border_row_norm(&collapsed_rows) < 0.02,
        "deep dyadic source should be visibly collapsed, got norm {}",
        tensor_border_row_norm(&collapsed_rows)
    );
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_sector_coupling_sweep_is_rejected() {
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_pairs = [
        (
            "origin-right-left",
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginLeft,
        ),
        (
            "origin-left-right",
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
        ),
    ];
    let midpoint_pairs = [
        (
            "midpoint-right-left",
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "midpoint-left-right",
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        ),
    ];
    let origin_b = [
        ("origin-paper", kp20_b_coeffs_near_o(c)),
        ("origin-reversed", [c, 2.0 * c, 2.0 * c]),
    ];
    let midpoint_b = [
        ("midpoint-paper", kp20_b_coeffs_near_m(c)),
        ("midpoint-reversed", [0.0, 0.0, c]),
        ("midpoint-middle", [0.0, c, 0.0]),
    ];
    let mut best = Kp20SectorCouplingDiagnostic {
        residual: f64::INFINITY,
        origin_pair: "",
        midpoint_pair: "",
        origin_b: "",
        midpoint_b: "",
        rotation: 0.0,
    };

    for &(origin_pair, origin_previous, origin_current) in &origin_pairs {
        for &(midpoint_pair, midpoint_previous, midpoint_current) in &midpoint_pairs {
            for &(origin_b_label, origin_b_coeffs) in &origin_b {
                for &(midpoint_b_label, midpoint_b_coeffs) in &midpoint_b {
                    for rotation in [angle, -angle] {
                        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                        vector_constraints.extend(
                            kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap(),
                        );
                        let mut scalar_constraints =
                            kp20_planar_tau_diagonal_symmetry_constraints();
                        scalar_constraints.extend(
                            kp20_planar_tau_rotated_sector_g1_g2_constraints(
                                origin_previous,
                                origin_current,
                                origin_b_coeffs,
                                rotation,
                            ),
                        );
                        scalar_constraints.extend(
                            kp20_planar_tau_rotated_sector_g1_g2_constraints(
                                midpoint_previous,
                                midpoint_current,
                                midpoint_b_coeffs,
                                rotation,
                            ),
                        );
                        let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                            &Kp20PlanarTauF3ScalarSolveInput {
                                initial: Kp20MacroPatch::zero(),
                                vector_constraints,
                                scalar_constraints,
                                regularization_weight: 1e-10,
                            },
                        )
                        .unwrap();
                        if solved.max_constraint_residual < best.residual {
                            best = Kp20SectorCouplingDiagnostic {
                                residual: solved.max_constraint_residual,
                                origin_pair,
                                midpoint_pair,
                                origin_b: origin_b_label,
                                midpoint_b: midpoint_b_label,
                                rotation,
                            };
                        }
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 11.8 && best.residual < 11.9,
        "coherent sector-coupling sweep should remain rejected near 22.14 after the eq(14) adjacent-half source assignment, got {best:?}"
    );
    assert_eq!(best.origin_pair, "origin-right-left");
    assert_eq!(best.midpoint_pair, "midpoint-left-right");
    assert_eq!(best.origin_b, "origin-paper");
    assert_eq!(best.midpoint_b, "midpoint-reversed");
    assert!((best.rotation - angle).abs() < 1.0e-12);
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_eq9_source_sweep_remains_rejected() {
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_b = kp20_b_coeffs_near_o(c);
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let mut best = Kp20Eq9TermSourceDiagnostic {
        residual: f64::INFINITY,
        origin_terms: Kp20Eq9DiagnosticTerms::paper(),
        midpoint_terms: Kp20Eq9DiagnosticTerms::paper(),
    };
    let mut least_bad_class: Vec<Kp20Eq9DiagnosticTerms> = Vec::new();

    for midpoint_signs in kp20_eq9_term_sign_variants() {
        for midpoint_terms in kp20_eq9_term_source_variants(midpoint_signs) {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    origin_b,
                    origin_b,
                    Kp20Eq9DiagnosticTerms::paper(),
                    angle,
                ),
            );
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
                    Kp20PieceId::NearMidpointLeft,
                    Kp20PieceId::NearMidpointRight,
                    midpoint_b,
                    midpoint_b,
                    midpoint_terms,
                    angle,
                ),
            );
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints,
                    regularization_weight: 1e-10,
                },
            )
            .unwrap();
            eprintln!(
                "eq9 source sweep variant duv={:?} duu={:?} du={:?} signs=({},{},{}) residual={:e}",
                midpoint_terms.duv_source,
                midpoint_terms.duu_source,
                midpoint_terms.du_source,
                midpoint_terms.signs.two_b_duv,
                midpoint_terms.signs.b_sq_duu,
                midpoint_terms.signs.b_bprime_du,
                solved.max_constraint_residual
            );
            if solved.max_constraint_residual > 19.4 && solved.max_constraint_residual < 19.5 {
                least_bad_class.push(midpoint_terms);
            }
            if solved.max_constraint_residual < best.residual {
                best = Kp20Eq9TermSourceDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_terms: Kp20Eq9DiagnosticTerms::paper(),
                    midpoint_terms,
                };
            }
        }
    }

    assert!(
        best.residual > 19.4 && best.residual < 19.5,
        "corrected scaled q13/q2 source should remain rejected after midpoint eq9 source/sign sweep, got {best:?}"
    );
    // The least-bad class is an 8-way tie (identical residuals to
    // ~11 significant digits — measured 2026-08-16, PK-291): the
    // argmin's identity within it is solver rounding noise, so the
    // pinned facts are the class's structure, not which member the
    // sweep happens to visit at the lowest noise. Every member
    // shares duu/du sourced from the previous sector with the
    // b²·d_uu and b·b′·d_u sign flips anti-correlated; the d_uv
    // source and the 2b·d_uv sign are fully degenerate. The
    // remaining sign combinations of the same source class are
    // rejected far away (≈41.2), as are the other source classes
    // (≈26.3, ≈34.0).
    assert_eq!(least_bad_class.len(), 8, "{least_bad_class:?}");
    for terms in &least_bad_class {
        assert_eq!(terms.duu_source, Kp20Eq9TermSource::Previous);
        assert_eq!(terms.du_source, Kp20Eq9TermSource::Previous);
        assert_eq!(terms.signs.b_sq_duu * terms.signs.b_bprime_du, -1.0);
    }
    for duv in [Kp20Eq9TermSource::Previous, Kp20Eq9TermSource::Current] {
        assert!(
            least_bad_class.iter().any(|terms| terms.duv_source == duv),
            "d_uv source {duv:?} missing from the tied class"
        );
    }
    for sign in [1.0, -1.0] {
        assert!(
            least_bad_class
                .iter()
                .any(|terms| terms.signs.two_b_duv == sign),
            "2b·d_uv sign {sign} missing from the tied class"
        );
    }
}

#[test]
fn stam_appendix_a_unknown_top_vector_source_span_orientation_sweep_is_near_feasible() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let top_count = 2 * ring.sectors.len() + 1;
    let mut best = StamAppendixAUnknownTopDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
    };
    for patch_index in 0..3 {
        for row_start in [0, 1] {
            for reverse_rows in [false, true] {
                let tensor_border_constraints =
                    stam_appendix_a_unknown_top_vector_tensor_border_constraints(
                        &ring,
                        patch_index,
                        row_start,
                        reverse_rows,
                    );
                let vector_constraints = kp20_planar_tau_internal_c2_constraints();
                let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                scalar_constraints.extend(
                    kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
                        Kp20PieceId::NearOriginRight,
                        Kp20PieceId::NearOriginLeft,
                        [c, 2.0 * c, 2.0 * c],
                        kp20_b_coeffs_near_o(c),
                        angle,
                    ),
                );
                scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                    Kp20PieceId::NearMidpointLeft,
                    Kp20PieceId::NearMidpointRight,
                    kp20_b_coeffs_near_m(c),
                    angle,
                ));
                let mut augmented_scalar_constraints = tensor_border_constraints;
                augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                    terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(1), 1.0)],
                    target: ring.sectors[0][1].x,
                });
                let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                    &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                        initial: Kp20MacroPatch::zero(),
                        vector_constraints,
                        scalar_constraints,
                        augmented_scalar_constraints,
                        extra_variable_count: 2 * top_count,
                        regularization_weight: 1e-10,
                        extra_regularization_weight: 1e-12,
                    },
                )
                .unwrap();
                if solved.max_constraint_residual < best.residual {
                    best = StamAppendixAUnknownTopDiagnostic {
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
        best.residual < 1.0e-6,
        "unknown Appendix-A top-vector source span should satisfy the paper gate, got {best:?}"
    );
    assert_eq!(best.patch_index, 0);
    assert_eq!(best.row_start, 0);
    assert!(!best.reverse_rows);
}

#[test]
fn stam_appendix_a_feasible_unknown_top_vector_is_not_pure_frequency_one() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let solved = solve_stam_appendix_a_unknown_top_vector_source_span(&ring, 0, 0, false);
    assert!(solved.max_constraint_residual < 1.0e-6);

    let top_count = 2 * ring.sectors.len() + 1;
    let odd = (0..ring.sectors.len())
        .map(|sector| {
            Complex::new(
                solved.extra_values[1 + 2 * sector],
                solved.extra_values[top_count + 1 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let even = (0..ring.sectors.len())
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
    assert!(
        odd_residual > 1.04 && odd_residual < 1.05,
        "feasible unknown Appendix-A top vector should expose non-frequency-one odd ring, got {odd_residual}"
    );
    assert!(
        even_residual < 1.0e-12,
        "feasible unknown Appendix-A top vector even ring should remain frequency-one, got {even_residual}"
    );
    for magnitude in odd_magnitudes {
        assert!(
            magnitude > 0.26 && magnitude < 0.27,
            "feasible unknown Appendix-A odd top vector should expose sector-local freedom, got magnitude {magnitude}"
        );
    }
    for magnitude in even_magnitudes {
        assert!(
            magnitude < 1.0e-12,
            "feasible unknown Appendix-A even top vector should be zero in all DFT modes, got {magnitude}"
        );
    }
}

#[test]
fn stam_appendix_a_feasible_unknown_top_vector_tau_orientation_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let solved = solve_stam_appendix_a_unknown_top_vector_source_span(&ring, 0, 0, false);
    assert!(solved.max_constraint_residual < 1.0e-6);

    let local_tau = cc_test_global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "feasible unknown Appendix-A top-vector tau validity: residual={residual:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
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

    assert!(
        !validity.is_paper_valid_reparameterization(),
        "unknown top-vector row closure unexpectedly produced a sampled paper-valid tau"
    );
    assert_eq!(
        validity.orientation.near_zero_determinants,
        validity.orientation.sample_count
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_source_span_orientation_sweep_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut best = StamAppendixARotationalTopDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
        rotation_direction: 1.0,
    };
    for rotation_direction in [-1.0, 1.0] {
        for patch_index in 0..3 {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    let tensor_border_constraints =
                        stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                            &ring,
                            patch_index,
                            row_start,
                            reverse_rows,
                            rotation_direction,
                        );
                    let vector_constraints = kp20_planar_tau_internal_c2_constraints();
                    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
                            Kp20PieceId::NearOriginRight,
                            Kp20PieceId::NearOriginLeft,
                            [c, 2.0 * c, 2.0 * c],
                            kp20_b_coeffs_near_o(c),
                            angle,
                        ),
                    );
                    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                        Kp20PieceId::NearMidpointLeft,
                        Kp20PieceId::NearMidpointRight,
                        kp20_b_coeffs_near_m(c),
                        angle,
                    ));
                    let mut augmented_scalar_constraints = tensor_border_constraints;
                    augmented_scalar_constraints.push(
                        Kp20PlanarTauAugmentedScalarLinearConstraint {
                            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                            target: ring.sectors[0][1].x,
                        },
                    );
                    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                            initial: Kp20MacroPatch::zero(),
                            vector_constraints,
                            scalar_constraints,
                            augmented_scalar_constraints,
                            extra_variable_count: 4,
                            regularization_weight: 1e-10,
                            extra_regularization_weight: 1e-12,
                        },
                    )
                    .unwrap();
                    if solved.max_constraint_residual < best.residual {
                        best = StamAppendixARotationalTopDiagnostic {
                            residual: solved.max_constraint_residual,
                            patch_index,
                            row_start,
                            reverse_rows,
                            rotation_direction,
                        };
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 1.3 && best.residual < 1.32,
        "current rotational Appendix-A top-vector source span frontier changed: {best:?}"
    );
    assert_eq!(best.patch_index, 0);
    assert_eq!(best.row_start, 0);
    assert!(!best.reverse_rows);
    assert_eq!(best.rotation_direction, -1.0);
}

#[test]
fn stam_appendix_a_rotational_top_vector_with_center_source_span_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut best = StamAppendixARotationalTopDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        row_start: 0,
        reverse_rows: false,
        rotation_direction: 1.0,
    };
    for rotation_direction in [-1.0, 1.0] {
        for patch_index in 0..3 {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    let tensor_border_constraints =
                        stam_appendix_a_rotational_top_vector_with_center_tensor_border_constraints(
                            &ring,
                            patch_index,
                            row_start,
                            reverse_rows,
                            rotation_direction,
                        );
                    let vector_constraints = kp20_planar_tau_internal_c2_constraints();
                    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
                            Kp20PieceId::NearOriginRight,
                            Kp20PieceId::NearOriginLeft,
                            [c, 2.0 * c, 2.0 * c],
                            kp20_b_coeffs_near_o(c),
                            angle,
                        ),
                    );
                    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                        Kp20PieceId::NearMidpointLeft,
                        Kp20PieceId::NearMidpointRight,
                        kp20_b_coeffs_near_m(c),
                        angle,
                    ));
                    let mut augmented_scalar_constraints = tensor_border_constraints;
                    augmented_scalar_constraints.push(
                        Kp20PlanarTauAugmentedScalarLinearConstraint {
                            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(2), 1.0)],
                            target: ring.sectors[0][1].x,
                        },
                    );
                    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                            initial: Kp20MacroPatch::zero(),
                            vector_constraints,
                            scalar_constraints,
                            augmented_scalar_constraints,
                            extra_variable_count: 6,
                            regularization_weight: 1e-10,
                            extra_regularization_weight: 1e-12,
                        },
                    )
                    .unwrap();
                    if solved.max_constraint_residual < best.residual {
                        best = StamAppendixARotationalTopDiagnostic {
                            residual: solved.max_constraint_residual,
                            patch_index,
                            row_start,
                            reverse_rows,
                            rotation_direction,
                        };
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 1.3 && best.residual < 1.32,
        "center-augmented rotational Appendix-A top-vector source span frontier changed: {best:?}"
    );
    assert_eq!(best.patch_index, 0);
    assert_eq!(best.row_start, 0);
    assert!(!best.reverse_rows);
    assert_eq!(best.rotation_direction, -1.0);
}

#[test]
fn pr08_tensor_border_candidate_rejects_invalid_ring() {
    let ring = Kp20CatmullClarkCharacteristicRing {
        sectors: vec![[Point3::new(0.0, 0.0, 0.0); 13]; 2],
        lambda: 0.5,
    };
    assert!(
        kp20_pr08_characteristic_tensor_border_candidate_rows(
            &ring,
            Kp20Pr08TensorBorderCandidate::VisibleStrip
        )
        .is_err()
    );
}
