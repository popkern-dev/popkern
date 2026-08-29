use super::*;

#[test]
fn paper_row_set_assembles_all_kp20_tau_constraint_families() {
    let tensor_border = constant_tensor_border(5, p(0.25, 0.25));
    let row_set = kp20_planar_tau_paper_row_set(5, 2, &tensor_border).expect("row set");

    assert_eq!(row_set.vector_constraints.len(), 84 + 4 + 42 + 15 + 15);
    assert_eq!(
        row_set.scalar_constraints.len(),
        0,
        "KP20 Fig. 8 indicated averaging is not the same as the full-grid diagonal diagnostic"
    );
}

#[test]
fn fig8_indicated_average_maps_to_four_quadrant_bottom_edges() {
    let rows = kp20_planar_tau_fig8_indicated_averaging_constraints();

    assert_eq!(rows.len(), 4);
    for (row, piece) in rows.iter().zip(Kp20PieceId::all()) {
        assert_eq!(
            row.terms,
            vec![
                (Kp20PlanarTauCoeff::new(piece, 1, 0), 1.0),
                (Kp20PlanarTauCoeff::new(piece, 0, 0), -0.5),
                (Kp20PlanarTauCoeff::new(piece, 2, 0), -0.5),
            ]
        );
        assert_eq!(row.target, [0.0, 0.0]);
    }
}

#[test]
fn sector_local_diagonal_rows_separate_self_and_reflected_pair_blocks() {
    let self_rows = kp20_planar_tau_sector_local_self_diagonal_symmetry_constraints(5, 0).unwrap();
    let full_diagonal_rows =
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(5, 0).unwrap();

    assert_eq!(self_rows.len(), 98);
    assert_eq!(full_diagonal_rows.len(), 196);
    assert_eq!(&full_diagonal_rows[..self_rows.len()], self_rows.as_slice());

    for row in &self_rows {
        for &(scalar_coeff, _) in &row.terms {
            assert!(
                matches!(
                    scalar_coeff.coeff.piece,
                    Kp20PieceId::NearOriginLeft | Kp20PieceId::NearMidpointRight
                ),
                "self-diagonal row touched {:?}",
                scalar_coeff.coeff.piece
            );
        }
    }

    let pair_rows = &full_diagonal_rows[self_rows.len()..];
    assert!(
        pair_rows[0]
            .terms
            .iter()
            .any(|(scalar_coeff, coefficient)| {
                scalar_coeff.coeff == Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0)
                    && *coefficient != 0.0
            })
    );
    assert!(
        pair_rows[0]
            .terms
            .iter()
            .any(|(scalar_coeff, coefficient)| {
                scalar_coeff.coeff == Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 6, 6)
                    && *coefficient != 0.0
            })
    );
    assert!(pair_rows.iter().any(|row| {
        row.terms
            .iter()
            .any(|(scalar_coeff, _)| scalar_coeff.coeff.piece == Kp20PieceId::NearOriginRight)
    }));
    assert!(pair_rows.iter().any(|row| {
        row.terms
            .iter()
            .any(|(scalar_coeff, _)| scalar_coeff.coeff.piece == Kp20PieceId::NearMidpointLeft)
    }));
}

#[test]
fn paper_row_set_rejects_invalid_valence_and_sector() {
    let tensor_border = constant_tensor_border(5, p(0.0, 0.0));

    assert_eq!(
        kp20_planar_tau_paper_row_set(2, 0, &tensor_border).unwrap_err(),
        Kp20PlanarTauPaperSolveError::InvalidValence(2)
    );
    assert!(matches!(
        kp20_planar_tau_paper_row_set(5, 5, &tensor_border),
        Err(Kp20PlanarTauPaperSolveError::TensorBorder(_))
    ));
}

#[test]
fn paper_tau_f3_solve_preserves_constant_characteristic_border() {
    let value = p(0.25, 0.25);
    let tensor_border = constant_tensor_border(5, value);
    let initial = constant_tau(value);
    let solved = solve_kp20_planar_tau_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial,
        regularization_weight: 1.0e-12,
    })
    .expect("constant paper tau solve succeeds");
    let row_set = kp20_planar_tau_paper_row_set(5, 0, &tensor_border).expect("row set");

    assert!(solved.max_constraint_residual <= 1.0e-8);
    assert!(kp20_planar_bi6_tau_f3_energy(&solved.tau_macro_patch) <= 1.0e-16);
    for (row_index, constraint) in row_set.vector_constraints.iter().enumerate() {
        let residual = vector_residual(&solved.tau_macro_patch, constraint);
        assert!(
            residual[0].abs() <= 1.0e-8 && residual[1].abs() <= 1.0e-8,
            "vector row {row_index} residual {residual:?}"
        );
    }
    for (row_index, constraint) in row_set.scalar_constraints.iter().enumerate() {
        let residual = scalar_residual(&solved.tau_macro_patch, constraint);
        assert!(
            residual.abs() <= 1.0e-8,
            "scalar row {row_index} residual {residual}"
        );
    }
}

#[test]
fn explicit_appendix_a_tau_row_family_residuals_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let tensor_border = planar_regular_sector_tensor_border(valence);

    for (label, maps) in [
        (
            "positive-nonoverlap",
            kp20_positive_nonoverlap_c2_sector_piece_maps(),
        ),
        ("paper-quadrant", kp20_paper_quadrant_sector_piece_maps()),
    ] {
        let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
            &local_identity_tri6_guide(),
            valence,
            sector,
            Kp20TauFamily::Paper,
            maps,
        )
        .expect("paper tau samples with identity guide");
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            sector,
            &tensor_border,
            &global_tau,
        )
        .expect("family residuals");

        eprintln!(
            "explicit Appendix-A tau family residuals {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] {:?}",
            residuals,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
        );
        assert!(orientation.is_orientation_preserving());
        assert!(
            residuals
                .iter()
                .all(|residual| residual.max_residual.is_finite())
        );
    }
}

#[test]
fn explicit_appendix_a_tau_implied_eq14_tensor_band_delta_is_measured() {
    let valence = 5;
    let sector = 0;
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let eq14_rows =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).expect("eq14 rows");
    let mut max_delta = 0.0_f64;
    let mut worst = (Kp20PieceId::NearOriginLeft, 0_usize, 0_usize);

    for row in eq14_rows {
        let coeff = row.terms[0].0;
        let actual = global_tau.piece(coeff.piece).cps[coeff.i][coeff.j];
        let delta =
            ((actual.x - row.target[0]).powi(2) + (actual.y - row.target[1]).powi(2)).sqrt();
        if delta > max_delta {
            max_delta = delta;
            worst = (coeff.piece, coeff.i, coeff.j);
        }
    }

    eprintln!(
        "explicit Appendix-A tau implied eq14 tensor-band delta: max_delta={max_delta:e} worst={piece:?}[{i},{j}]",
        piece = worst.0,
        i = worst.1,
        j = worst.2,
    );
    assert!(max_delta.is_finite());
    assert!(max_delta > 1.0);
}

#[test]
fn rotational_regular_source_with_appendix_a_tau_seed_measures_orientation() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let local_seed = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let initial = local_planar_patch_to_global(valence, sector, local_seed);

    for regularization_weight in [1.0e-10, 1.0e-6, 1.0e-2, 1.0, 1.0e2, 1.0e4] {
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector,
                characteristic_tensor_border: &tensor_border,
                initial,
                regularization_weight,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "regular-source Appendix-A tau seed weight={regularization_weight:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
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
fn appendix_a_implied_eq14_source_through_section42_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);

    let mut appendix_a_source_rows = Vec::new();
    for piece_id in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let piece = global_tau.piece(piece_id);
        for depth in 0..3 {
            for u_index in 0..=6 {
                let cp = piece.cps[depth][u_index];
                appendix_a_source_rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece_id, depth, u_index),
                    [cp.x, cp.y],
                ));
            }
        }
    }
    let map = kp20_section42_standard_map(angle);
    let section42_family_cases = [
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
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("diagonal rows");

    let mut cases: Vec<(String, bool, TestSection42FamilySet, bool)> = vec![
        (
            "source-only".to_string(),
            false,
            TestSection42FamilySet::new(false, false, false, false),
            false,
        ),
        (
            "source-fig8".to_string(),
            true,
            TestSection42FamilySet::new(false, false, false, false),
            false,
        ),
    ];
    for (family_label, families) in section42_family_cases {
        cases.push((format!("source-{family_label}"), false, families, false));
        cases.push((format!("source-{family_label}-fig8"), true, families, false));
    }
    cases.push((
        "source-section42-diagonal".to_string(),
        false,
        TestSection42FamilySet::new(true, true, true, true),
        true,
    ));
    cases.push((
        "source-section42-fig8-diagonal".to_string(),
        true,
        TestSection42FamilySet::new(true, true, true, true),
        true,
    ));

    for (label, include_fig8, families, include_diagonal) in cases {
        let mut vector_constraints = Vec::new();
        if include_fig8 {
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        }
        vector_constraints.extend(appendix_a_source_rows.clone());
        let mut scalar_constraints = Vec::new();
        if families.c0 {
            scalar_constraints.extend(kp20_section42_only_c0_rows(map));
        }
        if families.eq11 {
            scalar_constraints.extend(kp20_section42_only_eq11_rows(map));
        }
        if families.eq12 {
            scalar_constraints.extend(kp20_section42_only_eq12_rows(map));
        }
        if families.eq13 {
            scalar_constraints.extend(kp20_section42_only_eq13_rows(map, 2.0 * c));
        }
        if include_diagonal {
            scalar_constraints.extend(diagonal_rows.clone());
        }

        let explicit_vector_residual = vector_constraints
            .iter()
            .map(|row| {
                let r = vector_residual(&global_tau, row);
                r[0].abs().max(r[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let explicit_scalar_residual = scalar_constraints
            .iter()
            .map(|row| scalar_residual(&global_tau, row).abs())
            .fold(0.0_f64, f64::max);

        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: global_tau,
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .expect("Appendix-A source solve succeeds");
        let solved_local = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&solved_local, &SAMPLES, 1.0e-10, 1.0e-10);
        let source_validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let source_left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let solved_left = near_origin_first_order_frame(&solved_local, Kp20PieceId::NearOriginLeft);

        eprintln!(
            "Appendix-A implied eq14 source row split {label}: explicit_vector_residual={explicit_vector_residual:e} explicit_scalar_residual={explicit_scalar_residual:e} solve_residual={residual:e} source_neg={source_neg}/{samples} source_zero={source_zero} source_det=[{source_min:e},{source_max:e}] solved_neg={neg}/{samples} solved_zero={zero} solved_det=[{min_det:e},{max_det:e}] source_apex_det={source_apex_det:e} solved_apex_det={solved_apex_det:e}",
            residual = solved.max_constraint_residual,
            source_neg = source_validity.orientation.negative_determinants,
            source_zero = source_validity.orientation.near_zero_determinants,
            source_min = source_validity.orientation.min_determinant,
            source_max = source_validity.orientation.max_determinant,
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            source_apex_det = source_left.determinant,
            solved_apex_det = solved_left.determinant,
        );

        assert!(source_validity.is_paper_valid_reparameterization());
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn appendix_a_section42_c0_symbol_map_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mapped_pieces = Kp20PieceId::all()
        .into_iter()
        .flat_map(|piece| {
            [
                Kp20Section42MappedPiece::current(piece),
                Kp20Section42MappedPiece::previous(piece),
            ]
        })
        .collect::<Vec<_>>();
    let mut best = (
        f64::INFINITY,
        "unmeasured",
        1.0_f64,
        None,
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
    );

    for &p in &mapped_pieces {
        for &arrow_p in &mapped_pieces {
            if p.piece == arrow_p.piece || p.is_previous_sector == arrow_p.is_previous_sector {
                continue;
            }
            for &underlined_p in &mapped_pieces {
                for &arrow_underlined_p in &mapped_pieces {
                    if underlined_p.piece == arrow_underlined_p.piece
                        || underlined_p.is_previous_sector == arrow_underlined_p.is_previous_sector
                    {
                        continue;
                    }
                    for &(orientation_label, index_map) in &[
                        ("standard", None),
                        (
                            "previous-sector-all-transposed",
                            Some(TestSection42IndexMap::PreviousSectorAll),
                        ),
                        (
                            "previous-sector-central-transposed",
                            Some(TestSection42IndexMap::PreviousSectorCentralOnly),
                        ),
                        (
                            "stored-piece-edge-orientation",
                            Some(TestSection42IndexMap::StoredPieceEdgeOrientation),
                        ),
                    ] {
                        for rotation_sign in [1.0_f64, -1.0] {
                            let map = Kp20Section42PaperSymbolMap {
                                p,
                                arrow_p,
                                underlined_p,
                                arrow_underlined_p,
                                previous_sector_rotation_radians: rotation_sign * angle,
                            };
                            let rows = if let Some(index_map) = index_map {
                                let mut rows = Vec::new();
                                test_add_section42_c0_rows_previous_transposed(
                                    &mut rows, map, index_map,
                                );
                                rows
                            } else {
                                kp20_section42_only_c0_rows(map)
                            };
                            let residual = rows
                                .iter()
                                .map(|row| scalar_residual(&global_tau, row).abs())
                                .fold(0.0_f64, f64::max);
                            if residual < best.0 {
                                best = (
                                    residual,
                                    orientation_label,
                                    rotation_sign,
                                    index_map,
                                    p,
                                    arrow_p,
                                    underlined_p,
                                    arrow_underlined_p,
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "Appendix-A §4.2 C0 nondegenerate symbol map best residual={residual:e} orientation={orientation} rotation_sign={rotation_sign:+} index_map={index_map:?} p={p:?} arrow_p={arrow_p:?} underlined_p={underlined_p:?} arrow_underlined_p={arrow_underlined_p:?}",
        residual = best.0,
        orientation = best.1,
        rotation_sign = best.2,
        index_map = best.3,
        p = best.4,
        arrow_p = best.5,
        underlined_p = best.6,
        arrow_underlined_p = best.7,
    );

    assert!(
        best.0.is_finite(),
        "Appendix-A tau C0 symbol correspondence diagnostic must produce a finite residual"
    );
}

#[test]
fn appendix_a_macro_edge_correspondence_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mapped_pieces = Kp20PieceId::all()
        .into_iter()
        .flat_map(|piece| {
            [
                Kp20Section42MappedPiece::current(piece),
                Kp20Section42MappedPiece::previous(piece),
            ]
        })
        .collect::<Vec<_>>();
    let mut best = (
        f64::INFINITY,
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        TestMacroEdge::I0,
        Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        TestMacroEdge::I0,
        false,
        1.0_f64,
    );
    let mut exact_correspondences = Vec::new();

    for &left in &mapped_pieces {
        for &right in &mapped_pieces {
            if left.piece == right.piece
                || (left.is_previous_sector == right.is_previous_sector && left.is_previous_sector)
            {
                continue;
            }
            for left_edge in TestMacroEdge::ALL {
                for right_edge in TestMacroEdge::ALL {
                    for reversed in [false, true] {
                        for rotation_sign in [1.0_f64, -1.0] {
                            let residual = macro_edge_pair_residual(
                                &global_tau,
                                left,
                                left_edge,
                                right,
                                right_edge,
                                reversed,
                                rotation_sign * angle,
                            );
                            if residual < best.0 {
                                best = (
                                    residual,
                                    left,
                                    left_edge,
                                    right,
                                    right_edge,
                                    reversed,
                                    rotation_sign,
                                );
                            }
                            if residual < 1.0e-12 {
                                exact_correspondences.push((
                                    left,
                                    left_edge,
                                    right,
                                    right_edge,
                                    reversed,
                                    rotation_sign,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "Appendix-A macro edge correspondence best residual={residual:e} left={left:?}:{left_edge:?} right={right:?}:{right_edge:?} reversed={reversed} rotation_sign={rotation_sign:+}",
        residual = best.0,
        left = best.1,
        left_edge = best.2,
        right = best.3,
        right_edge = best.4,
        reversed = best.5,
        rotation_sign = best.6,
    );
    for (index, correspondence) in exact_correspondences.iter().take(24).enumerate() {
        eprintln!(
            "  exact[{index}] left={left:?}:{left_edge:?} right={right:?}:{right_edge:?} reversed={reversed} rotation_sign={rotation_sign:+}",
            left = correspondence.0,
            left_edge = correspondence.1,
            right = correspondence.2,
            right_edge = correspondence.3,
            reversed = correspondence.4,
            rotation_sign = correspondence.5,
        );
    }
    eprintln!(
        "Appendix-A macro edge exact correspondence count={}",
        exact_correspondences.len()
    );

    assert!(
        best.0.is_finite(),
        "Appendix-A macro edge correspondence diagnostic must produce a finite residual"
    );
}

#[test]
fn appendix_a_section42_strip_orientation_map_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mapped_pieces = Kp20PieceId::all()
        .into_iter()
        .flat_map(|piece| {
            [
                Kp20Section42MappedPiece::current(piece),
                Kp20Section42MappedPiece::previous(piece),
            ]
        })
        .collect::<Vec<_>>();
    let strip_maps = mapped_pieces
        .iter()
        .copied()
        .flat_map(|mapped| {
            TestMacroEdge::ALL.into_iter().flat_map(move |edge| {
                [false, true]
                    .into_iter()
                    .map(move |reversed| TestSection42StripMap {
                        mapped,
                        edge,
                        reversed,
                        sector_rotation_steps: if mapped.is_previous_sector { 1 } else { 0 },
                    })
            })
        })
        .collect::<Vec<_>>();
    let mut c0_pairs = Vec::new();
    for &left in &strip_maps {
        for &right in &strip_maps {
            if left.mapped.piece == right.mapped.piece
                || (left.mapped.is_previous_sector && right.mapped.is_previous_sector)
            {
                continue;
            }
            let residual = strip_boundary_pair_residual(&global_tau, left, right, angle);
            if residual < 1.0e-12 {
                c0_pairs.push((left, right));
            }
        }
    }

    let mut best = (
        f64::INFINITY,
        f64::INFINITY,
        0_usize,
        0_usize,
        c0_pairs[0].0,
        c0_pairs[0].1,
        c0_pairs[0].0,
        c0_pairs[0].1,
    );
    for (p_index, &(p, arrow_p)) in c0_pairs.iter().enumerate() {
        for (underlined_index, &(underlined_p, arrow_underlined_p)) in c0_pairs.iter().enumerate() {
            let symbol_map = TestSection42StripSymbolMap {
                p,
                arrow_p,
                underlined_p,
                arrow_underlined_p,
                previous_sector_rotation_radians: angle,
            };
            let rows = test_section42_strip_basis_rows(symbol_map, 2.0 * c);
            let residual = rows
                .iter()
                .map(|row| scalar_residual(&global_tau, row).abs())
                .fold(0.0_f64, f64::max);
            let c0_rows = test_section42_strip_c0_rows(symbol_map);
            let c0_residual = c0_rows
                .iter()
                .map(|row| scalar_residual(&global_tau, row).abs())
                .fold(0.0_f64, f64::max);
            if residual < best.0 {
                best = (
                    residual,
                    c0_residual,
                    p_index,
                    underlined_index,
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                );
            }
        }
    }

    eprintln!(
        "Appendix-A §4.2 strip orientation best residual={residual:e} c0={c0:e} c0_pairs={c0_pairs} p_pair={p_pair} underlined_pair={underlined_pair} p={p:?} arrow_p={arrow_p:?} underlined_p={underlined_p:?} arrow_underlined_p={arrow_underlined_p:?}",
        residual = best.0,
        c0 = best.1,
        c0_pairs = c0_pairs.len(),
        p_pair = best.2,
        underlined_pair = best.3,
        p = best.4,
        arrow_p = best.5,
        underlined_p = best.6,
        arrow_underlined_p = best.7,
    );
    let best_map = TestSection42StripSymbolMap {
        p: best.4,
        arrow_p: best.5,
        underlined_p: best.6,
        arrow_underlined_p: best.7,
        previous_sector_rotation_radians: angle,
    };
    let mut eq11_rows = Vec::new();
    test_add_section42_strip_eq11_degree5_rows(&mut eq11_rows, best_map);
    let mut eq12_rows = Vec::new();
    test_add_section42_strip_eq12_internal_c2_rows(&mut eq12_rows, best_map);
    let mut eq13_rows = Vec::new();
    test_add_section42_strip_eq13_rows(&mut eq13_rows, best_map, 2.0 * c);
    eprintln!(
        "  family residuals c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e}",
        c0 = best.1,
        eq11 = scalar_row_set_residual(&global_tau, &eq11_rows),
        eq12 = scalar_row_set_residual(&global_tau, &eq12_rows),
        eq13 = scalar_row_set_residual(&global_tau, &eq13_rows),
    );
    let mut best_eq12_variant = (f64::INFINITY, false, 6_usize, 0_usize);
    for source_is_underlined in [false, true] {
        for source_endpoint in [0_usize, 6] {
            for receiving_endpoint in [0_usize, 6] {
                let mut rows = test_section42_strip_c0_rows(best_map);
                test_add_section42_strip_eq11_degree5_rows(&mut rows, best_map);
                test_add_section42_strip_eq12_variant_rows(
                    &mut rows,
                    best_map,
                    source_is_underlined,
                    source_endpoint,
                    receiving_endpoint,
                );
                test_add_section42_strip_eq13_rows(&mut rows, best_map, 2.0 * c);
                let residual = scalar_row_set_residual(&global_tau, &rows);
                if residual < best_eq12_variant.0 {
                    best_eq12_variant = (
                        residual,
                        source_is_underlined,
                        source_endpoint,
                        receiving_endpoint,
                    );
                }
            }
        }
    }
    eprintln!(
        "  best eq12 variant full residual={residual:e} source_is_underlined={source_is_underlined} source_endpoint={source_endpoint} receiving_endpoint={receiving_endpoint}",
        residual = best_eq12_variant.0,
        source_is_underlined = best_eq12_variant.1,
        source_endpoint = best_eq12_variant.2,
        receiving_endpoint = best_eq12_variant.3,
    );

    assert!(
        best.0.is_finite(),
        "Appendix-A §4.2 strip orientation diagnostic must produce a finite residual"
    );
}

#[test]
fn appendix_a_fig6c_semantic_strip_graph_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let piece_maps = kp20_paper_quadrant_sector_piece_maps();
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        piece_maps,
    )
    .expect("paper tau samples with identity guide");
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let paper_square = test_square_domain_report(&local_tau, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A Fig.6c paper-quadrant paper tau square-domain: outside={outside}/{samples} excess={excess:e} worst={worst} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
        outside = paper_square.outside,
        samples = paper_square.samples,
        excess = paper_square.max_excess,
        worst = paper_square.worst_sample,
        min_u = paper_square.min_u,
        max_u = paper_square.max_u,
        min_v = paper_square.min_v,
        max_v = paper_square.max_v,
    );
    let global_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let refined_local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        piece_maps,
    )
    .expect("refined paper tau samples with identity guide");
    let refined_square =
        test_square_domain_report(&refined_local_tau, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let refined_validity = kp20_planar_bi6_tau_validity_report(
        &refined_local_tau,
        &TAU_VALIDITY_SAMPLES,
        1.0e-10,
        1.0e-10,
    );
    eprintln!(
        "Appendix-A Fig.6c paper-quadrant refined-apex tau validity: outside={outside}/{samples} max_excess={max_excess:e} worst_domain={worst_domain} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        outside = refined_validity.domain.outside_triangle_samples,
        samples = refined_validity.domain.sample_count,
        max_excess = refined_validity.domain.max_simplex_excess,
        worst_domain = refined_validity.domain.worst_sample,
        neg = refined_validity.orientation.negative_determinants,
        orientation_samples = refined_validity.orientation.sample_count,
        zero = refined_validity.orientation.near_zero_determinants,
        min_det = refined_validity.orientation.min_determinant,
        max_det = refined_validity.orientation.max_determinant,
    );
    eprintln!(
        "Appendix-A Fig.6c paper-quadrant refined-apex square-domain: outside={outside}/{samples} excess={excess:e} worst={worst} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
        outside = refined_square.outside,
        samples = refined_square.samples,
        excess = refined_square.max_excess,
        worst = refined_square.worst_sample,
        min_u = refined_square.min_u,
        max_u = refined_square.max_u,
        min_v = refined_square.min_v,
        max_v = refined_square.max_v,
    );
    let refined_global_tau = local_planar_patch_to_global(valence, sector, refined_local_tau);
    let strips = appendix_a_candidate_section42_strips();
    let sector_join_pairs = exact_strip_boundary_pairs(&global_tau, &strips, angle);
    let classified_sector_pairs = sector_join_pairs
        .iter()
        .filter(|&&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::SectorJoin)
            )
        })
        .count();
    let classified_two_piece_pairs = sector_join_pairs
        .iter()
        .filter(|&&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::TwoPieceJoin)
            )
        })
        .count();
    eprintln!(
        "Appendix-A exact strip pairs: total={} sector_kind={} two_piece_kind={}",
        sector_join_pairs.len(),
        classified_sector_pairs,
        classified_two_piece_pairs,
    );
    for (index, &(left, right)) in sector_join_pairs.iter().enumerate() {
        eprintln!(
            "  exact_strip_pair[{index}] kind={kind:?} left={left:?} right={right:?}",
            kind = fig6c_exact_pair_kind(left, right),
        );
    }
    let mut candidates = Vec::new();

    for &(p, arrow_p) in &sector_join_pairs {
        if !matches!(
            fig6c_exact_pair_kind(p, arrow_p),
            Some(TestFig6cExactPairKind::SectorJoin)
        ) {
            continue;
        }
        for &(underlined_p, arrow_underlined_p) in &sector_join_pairs {
            if !matches!(
                fig6c_exact_pair_kind(underlined_p, arrow_underlined_p),
                Some(TestFig6cExactPairKind::SectorJoin)
            ) {
                continue;
            }
            if !fig6c_semantic_symbols_are_distinct([p, arrow_p, underlined_p, arrow_underlined_p])
            {
                continue;
            }
            let map = TestSection42StripSymbolMap {
                p,
                arrow_p,
                underlined_p,
                arrow_underlined_p,
                previous_sector_rotation_radians: angle,
            };
            let plain_join = best_strip_eq12_join_residual(&global_tau, p, underlined_p, angle);
            let arrow_join =
                best_strip_eq12_join_residual(&global_tau, arrow_p, arrow_underlined_p, angle);
            candidates.push((map, plain_join, arrow_join));
        }
    }

    let mut best = (
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        None,
    );
    for &(map, plain_join, arrow_join) in &candidates {
        let rows = test_section42_strip_basis_rows(map, 2.0 * c);
        let full_residual = scalar_row_set_residual(&global_tau, &rows);
        let c0_residual = scalar_row_set_residual(&global_tau, &test_section42_strip_c0_rows(map));
        let mut eq11_rows = Vec::new();
        test_add_section42_strip_eq11_degree5_rows(&mut eq11_rows, map);
        let mut eq12_rows = Vec::new();
        test_add_section42_strip_eq12_internal_c2_rows(&mut eq12_rows, map);
        let mut eq13_rows = Vec::new();
        test_add_section42_strip_eq13_rows(&mut eq13_rows, map, 2.0 * c);
        let eq11_residual = scalar_row_set_residual(&global_tau, &eq11_rows);
        let eq12_residual = scalar_row_set_residual(&global_tau, &eq12_rows);
        let eq13_residual = scalar_row_set_residual(&global_tau, &eq13_rows);
        if full_residual < best.0 {
            best = (
                full_residual,
                c0_residual.max(eq11_residual),
                eq12_residual,
                eq13_residual,
                Some((map, plain_join, arrow_join)),
            );
        }
    }

    if let Some((map, plain_join, arrow_join)) = best.4 {
        eprintln!(
            "Appendix-A Fig.6c semantic strip graph best full={full:e} c0_eq11={c0_eq11:e} eq12={eq12:e} eq13={eq13:e} candidates={candidates} exact_internal_candidates={exact_internal_candidates} map={map:?} plain_join={plain_join:?} plain_source_endpoint={plain_source_endpoint} plain_receiving_endpoint={plain_receiving_endpoint} arrow_join={arrow_join:?} arrow_source_endpoint={arrow_source_endpoint} arrow_receiving_endpoint={arrow_receiving_endpoint}",
            full = best.0,
            c0_eq11 = best.1,
            eq12 = best.2,
            eq13 = best.3,
            candidates = candidates.len(),
            exact_internal_candidates = candidates
                .iter()
                .filter(|(_, plain, arrow)| {
                    plain.best_endpoint_residual < 1.0e-12 && arrow.best_endpoint_residual < 1.0e-12
                })
                .count(),
            plain_source_endpoint = plain_join.source_endpoint,
            plain_receiving_endpoint = plain_join.receiving_endpoint,
            arrow_source_endpoint = arrow_join.source_endpoint,
            arrow_receiving_endpoint = arrow_join.receiving_endpoint,
        );
    } else {
        eprintln!(
            "Appendix-A Fig.6c semantic strip graph found no candidate with exact sector joins; sector_join_pairs={}",
            sector_join_pairs.len()
        );
    }
    if let Some((sector_residual, p, arrow_p, underlined_p, arrow_underlined_p)) =
        best_fig6c_internal_join_first_graph(&global_tau, &sector_join_pairs, angle)
    {
        eprintln!(
            "Appendix-A Fig.6c internal-first graph best sector_residual={sector_residual:e} p={p:?} arrow_p={arrow_p:?} underlined_p={underlined_p:?} arrow_underlined_p={arrow_underlined_p:?}",
        );
    }
    if let Some((
        report,
        best_c0_report,
        best_eq13_report,
        best_c0_closed_report,
        c0_closed_count,
    )) = best_fig6c_eq12_normalized_source_graph(&global_tau, &strips, &sector_join_pairs, angle)
    {
        eprintln!(
            "Appendix-A Fig.6c eq12-eq13-normalized source graph best full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} candidates={candidates} map={map:?}",
            full = report.full_residual,
            c0 = report.c0_residual,
            eq11 = report.eq11_residual,
            eq12 = report.eq12_residual,
            eq13 = report.eq13_residual,
            candidates = report.candidate_count,
            map = report.map,
        );
        eprintln!(
            "  eq12-eq13-normalized source graph best_c0={best_c0:e} best_c0_full={best_c0_full:e} best_c0_eq13={best_c0_eq13:e} best_eq13={best_eq13:e} best_eq13_full={best_eq13_full:e} best_c0_closed_full={best_c0_closed_full:e} best_c0_closed_eq13={best_c0_closed_eq13:e} c0_closed_count={c0_closed_count}",
            best_c0 = best_c0_report.c0_residual,
            best_c0_full = best_c0_report.full_residual,
            best_c0_eq13 = best_c0_report.eq13_residual,
            best_eq13 = best_eq13_report.eq13_residual,
            best_eq13_full = best_eq13_report.full_residual,
            best_c0_closed_full = best_c0_closed_report.full_residual,
            best_c0_closed_eq13 = best_c0_closed_report.eq13_residual,
        );
        let eq12_patch = fig6c_eq12_normalized_source_patch(&global_tau, report.map);
        eprintln!(
            "  eq12/eq13 overlap after eq12-only normalization={overlap:e}",
            overlap = fig6c_eq12_eq13_overlap_residual(&eq12_patch, report.map),
        );
    }
    if let Some((
        report,
        source_delta,
        overlap,
        best_overlap_report,
        best_overlap_delta,
        best_overlap,
        overlap_closed_count,
        independent_count,
        independent_best_report,
        independent_best_overlap,
        c0_eq11_eq12_closed_count,
        best_c0_eq11_eq12_closed_report,
    )) = best_fig6c_shared_row_projected_source_graph(
        &global_tau,
        &strips,
        &sector_join_pairs,
        angle,
    ) {
        eprintln!(
            "Appendix-A Fig.6c shared-row projected source graph best full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} candidates={candidates} max_source_delta={source_delta:e} overlap={overlap:e} map={map:?}",
            full = report.full_residual,
            c0 = report.c0_residual,
            eq11 = report.eq11_residual,
            eq12 = report.eq12_residual,
            eq13 = report.eq13_residual,
            candidates = report.candidate_count,
            map = report.map,
        );
        eprintln!(
            "  shared-row projection best_overlap={best_overlap:e} best_overlap_full={best_overlap_full:e} best_overlap_delta={best_overlap_delta:e} overlap_closed_count={overlap_closed_count} independent_count={independent_count} independent_best_full={independent_best_full:e} independent_best_overlap={independent_best_overlap:e} c0_eq11_eq12_closed_count={c0_eq11_eq12_closed_count} best_c0_eq11_eq12_closed_full={closed_full:e} best_c0_eq11_eq12_closed_eq13={closed_eq13:e}",
            best_overlap_full = best_overlap_report.full_residual,
            independent_best_full = independent_best_report.full_residual,
            closed_full = best_c0_eq11_eq12_closed_report.full_residual,
            closed_eq13 = best_c0_eq11_eq12_closed_report.eq13_residual,
        );
    }
    if let Some((report, max_source_delta, row_lsq_residual)) =
        best_fig6c_source_row_basis_projected_graph(&global_tau, &strips, &sector_join_pairs, angle)
    {
        eprintln!(
            "Appendix-A Fig.6c source-row basis projected graph best full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} candidates={candidates} max_source_delta={max_source_delta:e} row_lsq_residual={row_lsq_residual:e} map={map:?}",
            full = report.full_residual,
            c0 = report.c0_residual,
            eq11 = report.eq11_residual,
            eq12 = report.eq12_residual,
            eq13 = report.eq13_residual,
            candidates = report.candidate_count,
            map = report.map,
        );
        let (boundary_projected, boundary_delta, boundary_lsq_residual) =
            fig6c_boundary_row_basis_projected_source_patch(&global_tau, report.map);
        let boundary_report = measure_fig6c_paper_direction_rows(&boundary_projected, report.map);
        let boundary_local = global_planar_patch_to_local(valence, sector, boundary_projected);
        let boundary_validity = kp20_planar_bi6_tau_validity_report(
            &boundary_local,
            &TAU_VALIDITY_SAMPLES,
            1.0e-10,
            1.0e-10,
        );
        let boundary_square =
            test_square_domain_report(&boundary_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  boundary-row basis projection on best source-row map full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} max_source_delta={boundary_delta:e} row_lsq_residual={boundary_lsq_residual:e} triangle_outside={outside}/{domain_samples} max_excess={max_excess:e} worst_domain={worst_domain} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            full = boundary_report.full_residual,
            c0 = boundary_report.c0_residual,
            eq11 = boundary_report.eq11_residual,
            eq12 = boundary_report.eq12_residual,
            eq13 = boundary_report.eq13_residual,
            outside = boundary_validity.domain.outside_triangle_samples,
            domain_samples = boundary_validity.domain.sample_count,
            max_excess = boundary_validity.domain.max_simplex_excess,
            worst_domain = boundary_validity.domain.worst_sample,
            square_outside = boundary_square.outside,
            square_samples = boundary_square.samples,
            square_excess = boundary_square.max_excess,
            square_worst = boundary_square.worst_sample,
            neg = boundary_validity.orientation.negative_determinants,
            orientation_samples = boundary_validity.orientation.sample_count,
            zero = boundary_validity.orientation.near_zero_determinants,
            min_det = boundary_validity.orientation.min_determinant,
            max_det = boundary_validity.orientation.max_determinant,
        );
        for boundary_weight in [0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let blended =
                blend_macro_patch_points(&global_tau, &boundary_projected, boundary_weight);
            let blended_report = measure_fig6c_paper_direction_rows(&blended, report.map);
            let blended_local = global_planar_patch_to_local(valence, sector, blended);
            let blended_validity = kp20_planar_bi6_tau_validity_report(
                &blended_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            let blended_square =
                test_square_domain_report(&blended_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            eprintln!(
                "    boundary blend weight={boundary_weight:e}: full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} triangle_outside={outside}/{domain_samples} max_excess={max_excess:e} worst_domain={worst_domain} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                full = blended_report.full_residual,
                c0 = blended_report.c0_residual,
                eq11 = blended_report.eq11_residual,
                eq12 = blended_report.eq12_residual,
                eq13 = blended_report.eq13_residual,
                outside = blended_validity.domain.outside_triangle_samples,
                domain_samples = blended_validity.domain.sample_count,
                max_excess = blended_validity.domain.max_simplex_excess,
                worst_domain = blended_validity.domain.worst_sample,
                square_outside = blended_square.outside,
                square_samples = blended_square.samples,
                square_excess = blended_square.max_excess,
                square_worst = blended_square.worst_sample,
                neg = blended_validity.orientation.negative_determinants,
                orientation_samples = blended_validity.orientation.sample_count,
                zero = blended_validity.orientation.near_zero_determinants,
                min_det = blended_validity.orientation.min_determinant,
                max_det = blended_validity.orientation.max_determinant,
            );
        }
    }
    let refined_exact_pairs = exact_strip_boundary_pairs(&refined_global_tau, &strips, angle);
    let refined_classified_sector_pairs = refined_exact_pairs
        .iter()
        .filter(|&&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::SectorJoin)
            )
        })
        .count();
    let refined_classified_two_piece_pairs = refined_exact_pairs
        .iter()
        .filter(|&&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::TwoPieceJoin)
            )
        })
        .count();
    eprintln!(
        "Appendix-A Fig.6c refined-apex exact strip pairs: total={} sector_kind={} two_piece_kind={}",
        refined_exact_pairs.len(),
        refined_classified_sector_pairs,
        refined_classified_two_piece_pairs,
    );
    if let Some((report, max_source_delta, row_lsq_residual)) =
        best_fig6c_source_row_basis_projected_graph(
            &refined_global_tau,
            &strips,
            &refined_exact_pairs,
            angle,
        )
    {
        eprintln!(
            "Appendix-A Fig.6c refined-apex source-row basis projected graph best full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} candidates={candidates} max_source_delta={max_source_delta:e} row_lsq_residual={row_lsq_residual:e} map={map:?}",
            full = report.full_residual,
            c0 = report.c0_residual,
            eq11 = report.eq11_residual,
            eq12 = report.eq12_residual,
            eq13 = report.eq13_residual,
            candidates = report.candidate_count,
            map = report.map,
        );
        let (boundary_projected, boundary_delta, boundary_lsq_residual) =
            fig6c_boundary_row_basis_projected_source_patch(&refined_global_tau, report.map);
        let boundary_report = measure_fig6c_paper_direction_rows(&boundary_projected, report.map);
        let boundary_local = global_planar_patch_to_local(valence, sector, boundary_projected);
        let boundary_validity = kp20_planar_bi6_tau_validity_report(
            &boundary_local,
            &TAU_VALIDITY_SAMPLES,
            1.0e-10,
            1.0e-10,
        );
        let boundary_square =
            test_square_domain_report(&boundary_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  refined-apex boundary-row basis projection full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} max_source_delta={boundary_delta:e} row_lsq_residual={boundary_lsq_residual:e} triangle_outside={outside}/{domain_samples} max_excess={max_excess:e} worst_domain={worst_domain} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            full = boundary_report.full_residual,
            c0 = boundary_report.c0_residual,
            eq11 = boundary_report.eq11_residual,
            eq12 = boundary_report.eq12_residual,
            eq13 = boundary_report.eq13_residual,
            outside = boundary_validity.domain.outside_triangle_samples,
            domain_samples = boundary_validity.domain.sample_count,
            max_excess = boundary_validity.domain.max_simplex_excess,
            worst_domain = boundary_validity.domain.worst_sample,
            square_outside = boundary_square.outside,
            square_samples = boundary_square.samples,
            square_excess = boundary_square.max_excess,
            square_worst = boundary_square.worst_sample,
            neg = boundary_validity.orientation.negative_determinants,
            orientation_samples = boundary_validity.orientation.sample_count,
            zero = boundary_validity.orientation.near_zero_determinants,
            min_det = boundary_validity.orientation.min_determinant,
            max_det = boundary_validity.orientation.max_determinant,
        );
        let mut corner_pinned = boundary_projected;
        corner_pinned.piece_mut(Kp20PieceId::NearMidpointRight).cps[6][6] =
            local_to_global_planar(valence, sector, [1.0, 1.0]);
        let corner_pinned_report = measure_fig6c_paper_direction_rows(&corner_pinned, report.map);
        let corner_pinned_local = global_planar_patch_to_local(valence, sector, corner_pinned);
        let corner_pinned_validity = kp20_planar_bi6_tau_validity_report(
            &corner_pinned_local,
            &TAU_VALIDITY_SAMPLES,
            1.0e-10,
            1.0e-10,
        );
        let corner_pinned_square =
            test_square_domain_report(&corner_pinned_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  refined-apex boundary-row plus Fig.8 NMR[6,6]=(1,1) pin full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            full = corner_pinned_report.full_residual,
            c0 = corner_pinned_report.c0_residual,
            eq11 = corner_pinned_report.eq11_residual,
            eq12 = corner_pinned_report.eq12_residual,
            eq13 = corner_pinned_report.eq13_residual,
            square_outside = corner_pinned_square.outside,
            square_samples = corner_pinned_square.samples,
            square_excess = corner_pinned_square.max_excess,
            square_worst = corner_pinned_square.worst_sample,
            neg = corner_pinned_validity.orientation.negative_determinants,
            orientation_samples = corner_pinned_validity.orientation.sample_count,
            zero = corner_pinned_validity.orientation.near_zero_determinants,
            min_det = corner_pinned_validity.orientation.min_determinant,
            max_det = corner_pinned_validity.orientation.max_determinant,
        );
        for (label, outer_boundary_only) in [("outer-boundary", true), ("full-quadrant", false)] {
            let mut normalized = boundary_projected;
            set_fig8_affine_square_piece(
                &mut normalized,
                valence,
                sector,
                Kp20PieceId::NearMidpointRight,
                outer_boundary_only,
            );
            let normalized_report = measure_fig6c_paper_direction_rows(&normalized, report.map);
            let normalized_local = global_planar_patch_to_local(valence, sector, normalized);
            let normalized_validity = kp20_planar_bi6_tau_validity_report(
                &normalized_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            let normalized_square =
                test_square_domain_report(&normalized_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            eprintln!(
                "  refined-apex boundary-row plus Fig.8 NMR {label} square normalization full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                full = normalized_report.full_residual,
                c0 = normalized_report.c0_residual,
                eq11 = normalized_report.eq11_residual,
                eq12 = normalized_report.eq12_residual,
                eq13 = normalized_report.eq13_residual,
                square_outside = normalized_square.outside,
                square_samples = normalized_square.samples,
                square_excess = normalized_square.max_excess,
                square_worst = normalized_square.worst_sample,
                neg = normalized_validity.orientation.negative_determinants,
                orientation_samples = normalized_validity.orientation.sample_count,
                zero = normalized_validity.orientation.near_zero_determinants,
                min_det = normalized_validity.orientation.min_determinant,
                max_det = normalized_validity.orientation.max_determinant,
            );
        }
        let midpoint_pieces = [
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ];
        let all_pieces = Kp20PieceId::all();
        for (label, pieces) in [
            ("midpoint-quadrants", &midpoint_pieces[..]),
            ("all-quadrants", &all_pieces[..]),
        ] {
            let mut normalized = boundary_projected;
            for &piece_id in pieces {
                set_fig8_affine_square_piece(&mut normalized, valence, sector, piece_id, false);
            }
            let normalized_report = measure_fig6c_paper_direction_rows(&normalized, report.map);
            let normalized_local = global_planar_patch_to_local(valence, sector, normalized);
            let normalized_validity = kp20_planar_bi6_tau_validity_report(
                &normalized_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            let normalized_square =
                test_square_domain_report(&normalized_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            eprintln!(
                "  refined-apex boundary-row plus Fig.8 {label} square normalization full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                full = normalized_report.full_residual,
                c0 = normalized_report.c0_residual,
                eq11 = normalized_report.eq11_residual,
                eq12 = normalized_report.eq12_residual,
                eq13 = normalized_report.eq13_residual,
                square_outside = normalized_square.outside,
                square_samples = normalized_square.samples,
                square_excess = normalized_square.max_excess,
                square_worst = normalized_square.worst_sample,
                neg = normalized_validity.orientation.negative_determinants,
                orientation_samples = normalized_validity.orientation.sample_count,
                zero = normalized_validity.orientation.near_zero_determinants,
                min_det = normalized_validity.orientation.min_determinant,
                max_det = normalized_validity.orientation.max_determinant,
            );
        }
        let affine_square = fig8_affine_square_patch(valence, sector);
        for affine_weight in [0.0_f64, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let blended =
                blend_macro_patch_points(&boundary_projected, &affine_square, affine_weight);
            let blended_report = measure_fig6c_paper_direction_rows(&blended, report.map);
            let blended_local = global_planar_patch_to_local(valence, sector, blended);
            let blended_validity = kp20_planar_bi6_tau_validity_report(
                &blended_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            let blended_square =
                test_square_domain_report(&blended_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            eprintln!(
                "  refined-apex boundary-row to Fig.8 affine-square blend weight={affine_weight:e} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                full = blended_report.full_residual,
                c0 = blended_report.c0_residual,
                eq11 = blended_report.eq11_residual,
                eq12 = blended_report.eq12_residual,
                eq13 = blended_report.eq13_residual,
                square_outside = blended_square.outside,
                square_samples = blended_square.samples,
                square_excess = blended_square.max_excess,
                square_worst = blended_square.worst_sample,
                neg = blended_validity.orientation.negative_determinants,
                orientation_samples = blended_validity.orientation.sample_count,
                zero = blended_validity.orientation.near_zero_determinants,
                min_det = blended_validity.orientation.min_determinant,
                max_det = blended_validity.orientation.max_determinant,
            );
        }
        let active_set = fig6c_square_domain_active_set_projection(
            &boundary_projected,
            report.map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
        );
        let active_set_report = measure_fig6c_paper_direction_rows(&active_set.patch, report.map);
        let active_set_local = global_planar_patch_to_local(valence, sector, active_set.patch);
        let active_set_validity = kp20_planar_bi6_tau_validity_report(
            &active_set_local,
            &TAU_VALIDITY_SAMPLES,
            1.0e-10,
            1.0e-10,
        );
        let active_set_orientation =
            test_orientation_detail_report(&active_set_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let active_set_square =
            test_square_domain_report(&active_set_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  refined-apex active-set square projection iterations={iterations} active={active} max_delta={max_delta:e} linear_residual={linear_residual:e} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            iterations = active_set.iterations,
            active = active_set.active_constraints,
            max_delta = active_set.max_delta,
            linear_residual = active_set.linear_residual,
            full = active_set_report.full_residual,
            c0 = active_set_report.c0_residual,
            eq11 = active_set_report.eq11_residual,
            eq12 = active_set_report.eq12_residual,
            eq13 = active_set_report.eq13_residual,
            square_outside = active_set_square.outside,
            square_samples = active_set_square.samples,
            square_excess = active_set_square.max_excess,
            square_worst = active_set_square.worst_sample,
            neg = active_set_validity.orientation.negative_determinants,
            orientation_samples = active_set_validity.orientation.sample_count,
            zero = active_set_validity.orientation.near_zero_determinants,
            min_det = active_set_validity.orientation.min_determinant,
            max_det = active_set_validity.orientation.max_determinant,
            det_worst = active_set_orientation.worst_sample,
        );
        for active_weight in [0.0_f64, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let blended =
                blend_macro_patch_points(&boundary_projected, &active_set.patch, active_weight);
            let blended_report = measure_fig6c_paper_direction_rows(&blended, report.map);
            let blended_local = global_planar_patch_to_local(valence, sector, blended);
            let blended_validity = kp20_planar_bi6_tau_validity_report(
                &blended_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            let blended_orientation =
                test_orientation_detail_report(&blended_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            let blended_square =
                test_square_domain_report(&blended_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            eprintln!(
                "  refined-apex active-set blend weight={active_weight:e} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
                full = blended_report.full_residual,
                c0 = blended_report.c0_residual,
                eq11 = blended_report.eq11_residual,
                eq12 = blended_report.eq12_residual,
                eq13 = blended_report.eq13_residual,
                square_outside = blended_square.outside,
                square_samples = blended_square.samples,
                square_excess = blended_square.max_excess,
                square_worst = blended_square.worst_sample,
                neg = blended_validity.orientation.negative_determinants,
                orientation_samples = blended_validity.orientation.sample_count,
                zero = blended_validity.orientation.near_zero_determinants,
                min_det = blended_validity.orientation.min_determinant,
                max_det = blended_validity.orientation.max_determinant,
                det_worst = blended_orientation.worst_sample,
            );
        }
    }

    if !best.0.is_finite() {
        eprintln!(
            "Appendix-A Fig.6c semantic strip graph has no finite sector-first candidate under the current exact-pair classifier"
        );
    }
}

#[test]
fn appendix_a_fig6c_active_set_source_is_measured_through_sampled_cap_assembly() {
    let valence = 5;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point: Point3::new(0.0, 0.0, 0.0),
        tensor_border: tensor_border.clone(),
    })
    .expect("guide builds");
    let mut worst_iterations = 0_usize;
    let mut worst_active = 0_usize;
    let mut worst_orientation_active = 0_usize;
    let mut worst_max_delta = 0.0_f64;
    let mut worst_linear_residual = 0.0_f64;
    let mut worst_square_outside = 0_usize;
    let mut worst_square_samples = 0_usize;
    let mut worst_square_excess = 0.0_f64;
    let mut worst_square_sample = String::new();
    let mut worst_negative = 0_usize;
    let mut worst_orientation_samples = 0_usize;
    let mut worst_zero = 0_usize;
    let mut min_det = f64::INFINITY;
    let mut max_det = f64::NEG_INFINITY;
    let mut worst_det_sample = String::new();
    let mut tau_geometry_patches = Vec::with_capacity(valence);
    let sampled_macro_patches = (0..valence)
        .map(|sector| {
            let active_set =
                fig6c_active_set_source_projection(valence, sector, &TAU_VALIDITY_SAMPLES)
                    .expect("active-set source projection");
            tau_geometry_patches.push(active_set.patch);
            let active_local = global_planar_patch_to_local(valence, sector, active_set.patch);
            let active_square =
                test_square_domain_report(&active_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            let active_orientation =
                test_orientation_detail_report(&active_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            let active_validity = kp20_planar_bi6_tau_validity_report(
                &active_local,
                &TAU_VALIDITY_SAMPLES,
                1.0e-10,
                1.0e-10,
            );
            worst_iterations = worst_iterations.max(active_set.iterations);
            worst_active = worst_active.max(active_set.active_constraints);
            worst_orientation_active =
                worst_orientation_active.max(active_set.orientation_constraints);
            worst_max_delta = worst_max_delta.max(active_set.max_delta);
            worst_linear_residual = worst_linear_residual.max(active_set.linear_residual);
            if active_square.max_excess > worst_square_excess {
                worst_square_excess = active_square.max_excess;
                worst_square_sample = active_square.worst_sample;
            }
            worst_square_outside = worst_square_outside.max(active_square.outside);
            worst_square_samples = worst_square_samples.max(active_square.samples);
            worst_negative = worst_negative.max(active_validity.orientation.negative_determinants);
            worst_orientation_samples =
                worst_orientation_samples.max(active_validity.orientation.sample_count);
            worst_zero = worst_zero.max(active_validity.orientation.near_zero_determinants);
            if active_validity.orientation.min_determinant < min_det {
                min_det = active_validity.orientation.min_determinant;
                worst_det_sample = active_orientation.worst_sample;
            }
            max_det = max_det.max(active_validity.orientation.max_determinant);
            sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
                guide.sector(sector),
                &active_local,
            )
        })
        .collect::<Vec<_>>();
    let tau_geometry_cap = test_macro_cap_from_sampled(valence, tau_geometry_patches);
    let tau_geometry_report = measure_kp20_assembly_constraints(&tau_geometry_cap);
    let tau_geometry_folds = test_summarize_cap_folds(&tau_geometry_cap);
    let tau_geometry_c0 = test_cap_v0_pair_c0_summary(&tau_geometry_cap);
    let pre_cap = test_macro_cap_from_sampled(valence, sampled_macro_patches.clone());
    let pre_report = measure_kp20_assembly_constraints(&pre_cap);
    let pre_folds = test_summarize_cap_folds(&pre_cap);
    let pre_c0 = test_cap_v0_pair_c0_summary(&pre_cap);
    let assembled =
        assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source(&Kp20Bi6MacroCapInput {
            tensor_border,
            sampled_macro_patches,
        })
        .expect("sampled-source assembly succeeds");
    let assembled_report = measure_kp20_assembly_constraints(&assembled);
    let assembled_folds = test_summarize_cap_folds(&assembled);
    let assembled_c0 = test_cap_v0_pair_c0_summary(&assembled);

    eprintln!(
        "Appendix-A Fig.6c active-set cap assembly n={valence}: tau_iterations={iterations} tau_active={active} tau_orientation_active={orientation_active} tau_max_delta={max_delta:e} tau_linear_residual={linear_residual:e} square_outside={square_outside}/{square_samples} square_excess={square_excess:e} square_worst={square_worst} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] det_worst={det_worst} tau_geom_c0_o={tau_geom_c0_o:e} tau_geom_c0_m={tau_geom_c0_m:e} tau_geom_eq8={tau_geom_eq8:e} tau_geom_eq8_o={tau_geom_eq8_o:e} tau_geom_eq8_m={tau_geom_eq8_m:e} tau_geom_eq9={tau_geom_eq9:e} tau_geom_eq9_o={tau_geom_eq9_o:e} tau_geom_eq9_m={tau_geom_eq9_m:e} tau_geom_folds={tau_geom_folds} tau_geom_min_fold_cos={tau_geom_fold_cos:e} tau_geom_piece_folds={tau_geom_piece_folds:?} pre_c0_o={pre_c0_o:e} pre_c0_m={pre_c0_m:e} pre_eq8={pre_eq8:e} pre_eq8_o={pre_eq8_o:e} pre_eq8_m={pre_eq8_m:e} pre_eq9={pre_eq9:e} pre_eq9_o={pre_eq9_o:e} pre_eq9_m={pre_eq9_m:e} pre_folds={pre_folds} pre_min_fold_cos={pre_fold_cos:e} pre_piece_folds={pre_piece_folds:?} post_c0_o={post_c0_o:e} post_c0_m={post_c0_m:e} post_eq8={post_eq8:e} post_eq8_o={post_eq8_o:e} post_eq8_m={post_eq8_m:e} post_eq9={post_eq9:e} post_eq9_o={post_eq9_o:e} post_eq9_m={post_eq9_m:e} post_folds={post_folds} post_min_fold_cos={post_fold_cos:e} post_piece_folds={post_piece_folds:?}",
        iterations = worst_iterations,
        active = worst_active,
        orientation_active = worst_orientation_active,
        max_delta = worst_max_delta,
        linear_residual = worst_linear_residual,
        square_outside = worst_square_outside,
        square_samples = worst_square_samples,
        square_excess = worst_square_excess,
        square_worst = worst_square_sample,
        neg = worst_negative,
        orientation_samples = worst_orientation_samples,
        zero = worst_zero,
        min_det = min_det,
        max_det = max_det,
        det_worst = worst_det_sample,
        tau_geom_c0_o = tau_geometry_c0.near_origin,
        tau_geom_c0_m = tau_geometry_c0.near_midpoint,
        tau_geom_eq8 = tau_geometry_report.worst_eq8(),
        tau_geom_eq8_o = max_slice(&tau_geometry_report.worst_eq8_near_o_per_sector),
        tau_geom_eq8_m = max_slice(&tau_geometry_report.worst_eq8_near_m_per_sector),
        tau_geom_eq9 = tau_geometry_report.worst_eq9(),
        tau_geom_eq9_o = max_slice(&tau_geometry_report.worst_eq9_near_o_per_sector),
        tau_geom_eq9_m = max_slice(&tau_geometry_report.worst_eq9_near_m_per_sector),
        tau_geom_folds = tau_geometry_folds.fold_count,
        tau_geom_fold_cos = tau_geometry_folds.min_jacobian_normal_cos,
        tau_geom_piece_folds = tau_geometry_folds.piece_fold_counts,
        pre_eq8 = pre_report.worst_eq8(),
        pre_c0_o = pre_c0.near_origin,
        pre_c0_m = pre_c0.near_midpoint,
        pre_eq8_o = max_slice(&pre_report.worst_eq8_near_o_per_sector),
        pre_eq8_m = max_slice(&pre_report.worst_eq8_near_m_per_sector),
        pre_eq9 = pre_report.worst_eq9(),
        pre_eq9_o = max_slice(&pre_report.worst_eq9_near_o_per_sector),
        pre_eq9_m = max_slice(&pre_report.worst_eq9_near_m_per_sector),
        pre_folds = pre_folds.fold_count,
        pre_fold_cos = pre_folds.min_jacobian_normal_cos,
        pre_piece_folds = pre_folds.piece_fold_counts,
        post_eq8 = assembled_report.worst_eq8(),
        post_c0_o = assembled_c0.near_origin,
        post_c0_m = assembled_c0.near_midpoint,
        post_eq8_o = max_slice(&assembled_report.worst_eq8_near_o_per_sector),
        post_eq8_m = max_slice(&assembled_report.worst_eq8_near_m_per_sector),
        post_eq9 = assembled_report.worst_eq9(),
        post_eq9_o = max_slice(&assembled_report.worst_eq9_near_o_per_sector),
        post_eq9_m = max_slice(&assembled_report.worst_eq9_near_m_per_sector),
        post_folds = assembled_folds.fold_count,
        post_fold_cos = assembled_folds.min_jacobian_normal_cos,
        post_piece_folds = assembled_folds.piece_fold_counts,
    );

    assert_eq!(worst_negative, 0);
    assert!(tau_geometry_report.all_finite());
    assert!(pre_report.all_finite());
    assert!(assembled_report.all_finite());
}

#[test]
fn appendix_a_fig6c_active_set_composed_cap_projection_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point: Point3::new(0.0, 0.0, 0.0),
        tensor_border,
    })
    .expect("guide builds");
    let composed = fig6c_composed_cap_residual_projection(
        &active_set.patch,
        map,
        guide.sector(sector),
        valence,
        sector,
    );
    let before_tau = measure_fig6c_paper_direction_rows(&active_set.patch, map);
    let after_tau = measure_fig6c_paper_direction_rows(&composed.patch, map);
    let before_cap = fig6c_composed_cap_residual_vector(
        &active_set.patch,
        guide.sector(sector),
        valence,
        sector,
    );
    let after_cap =
        fig6c_composed_cap_residual_vector(&composed.patch, guide.sector(sector), valence, sector);
    let before_local = global_planar_patch_to_local(valence, sector, active_set.patch);
    let after_local = global_planar_patch_to_local(valence, sector, composed.patch);
    let before_square = test_square_domain_report(&before_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let before_orientation =
        test_orientation_detail_report(&before_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);

    eprintln!(
        "Appendix-A Fig.6c composed cap projection n={valence} sector={sector}: linear_residual={linear_residual:e} alpha={alpha:e} max_delta={max_delta:e} cap_before={cap_before:e} cap_after={cap_after:e} tau_before_full={tau_before:e} tau_after_full={tau_after:e} tau_after_c0={tau_after_c0:e} tau_after_eq11={tau_after_eq11:e} tau_after_eq12={tau_after_eq12:e} tau_after_eq13={tau_after_eq13:e} square_before={square_before}/{square_samples} square_after={square_after}/{after_square_samples} square_after_excess={square_after_excess:e} orient_before_neg={before_neg}/{before_orientation_samples} orient_after_neg={after_neg}/{after_orientation_samples} det_before=[{before_min_det:e},{before_max_det:e}] det_after=[{after_min_det:e},{after_max_det:e}] after_det_worst={after_det_worst}",
        linear_residual = composed.linear_residual,
        alpha = composed.alpha,
        max_delta = composed.max_delta,
        cap_before = max_abs(&before_cap),
        cap_after = max_abs(&after_cap),
        tau_before = before_tau.full_residual,
        tau_after = after_tau.full_residual,
        tau_after_c0 = after_tau.c0_residual,
        tau_after_eq11 = after_tau.eq11_residual,
        tau_after_eq12 = after_tau.eq12_residual,
        tau_after_eq13 = after_tau.eq13_residual,
        square_before = before_square.outside,
        square_samples = before_square.samples,
        square_after = after_square.outside,
        after_square_samples = after_square.samples,
        square_after_excess = after_square.max_excess,
        before_neg = before_orientation.negative,
        before_orientation_samples = before_orientation.samples,
        after_neg = after_orientation.negative,
        after_orientation_samples = after_orientation.samples,
        before_min_det = before_orientation.min_det,
        before_max_det = before_orientation.max_det,
        after_min_det = after_orientation.min_det,
        after_max_det = after_orientation.max_det,
        after_det_worst = after_orientation.worst_sample,
    );

    assert!(composed.linear_residual.is_finite());
    assert!(max_abs(&after_cap).is_finite());
}

#[test]
fn appendix_a_fig6c_coupled_composed_projection_weight_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point: Point3::new(0.0, 0.0, 0.0),
        tensor_border,
    })
    .expect("guide builds");
    let base_cap = max_abs(&fig6c_composed_cap_residual_vector(
        &active_set.patch,
        guide.sector(sector),
        valence,
        sector,
    ));
    let base_tau = measure_fig6c_paper_direction_rows(&active_set.patch, map);
    let base_local = global_planar_patch_to_local(valence, sector, active_set.patch);
    let base_square = test_square_domain_report(&base_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let base_orientation =
        test_orientation_detail_report(&base_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A Fig.6c coupled composed projection baseline n={valence} sector={sector}: cap={base_cap:e} tau_full={tau_full:e} square={square}/{samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        tau_full = base_tau.full_residual,
        square = base_square.outside,
        samples = base_square.samples,
        neg = base_orientation.negative,
        orientation_samples = base_orientation.samples,
        min_det = base_orientation.min_det,
        max_det = base_orientation.max_det,
    );

    for cap_weight in [1.0e-4_f64, 1.0e-3, 1.0e-2, 1.0e-1, 1.0] {
        let projected = fig6c_composed_cap_residual_projection_with_weights(
            &active_set.patch,
            map,
            guide.sector(sector),
            valence,
            sector,
            100.0,
            cap_weight,
            false,
        );
        let after_cap = max_abs(&fig6c_composed_cap_residual_vector(
            &projected.patch,
            guide.sector(sector),
            valence,
            sector,
        ));
        let after_tau = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  coupled cap_weight={cap_weight:e}: alpha={alpha:e} linear_residual={linear_residual:e} max_delta={max_delta:e} cap={after_cap:e} cap_ratio={cap_ratio:e} tau_full={tau_full:e} tau_c0={tau_c0:e} tau_eq11={tau_eq11:e} tau_eq12={tau_eq12:e} tau_eq13={tau_eq13:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} valid_alpha={valid_alpha:e} valid_cap={valid_cap:e} valid_cap_ratio={valid_cap_ratio:e} valid_tau_full={valid_tau_full:e}",
            alpha = projected.alpha,
            linear_residual = projected.linear_residual,
            max_delta = projected.max_delta,
            cap_ratio = after_cap / base_cap,
            tau_full = after_tau.full_residual,
            tau_c0 = after_tau.c0_residual,
            tau_eq11 = after_tau.eq11_residual,
            tau_eq12 = after_tau.eq12_residual,
            tau_eq13 = after_tau.eq13_residual,
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
            valid_alpha = projected.valid_alpha,
            valid_cap = projected.valid_cap_residual,
            valid_cap_ratio = projected.valid_cap_residual / base_cap,
            valid_tau_full = projected.valid_tau_full,
        );
        assert!(after_cap.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_guarded_composed_projection_weight_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point: Point3::new(0.0, 0.0, 0.0),
        tensor_border,
    })
    .expect("guide builds");
    let base_cap = max_abs(&fig6c_composed_cap_residual_vector(
        &active_set.patch,
        guide.sector(sector),
        valence,
        sector,
    ));
    for cap_weight in [1.0e-4_f64, 1.0e-3, 1.0e-2, 1.0e-1, 1.0] {
        let projected = fig6c_composed_cap_residual_projection_with_weights(
            &active_set.patch,
            map,
            guide.sector(sector),
            valence,
            sector,
            100.0,
            cap_weight,
            true,
        );
        let after_cap = max_abs(&fig6c_composed_cap_residual_vector(
            &projected.patch,
            guide.sector(sector),
            valence,
            sector,
        ));
        let after_tau = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "  guarded coupled cap_weight={cap_weight:e}: alpha={alpha:e} valid_alpha={valid_alpha:e} linear_residual={linear_residual:e} max_delta={max_delta:e} cap={after_cap:e} cap_ratio={cap_ratio:e} tau_full={tau_full:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} valid_cap={valid_cap:e} valid_cap_ratio={valid_cap_ratio:e} valid_tau_full={valid_tau_full:e}",
            alpha = projected.alpha,
            valid_alpha = projected.valid_alpha,
            linear_residual = projected.linear_residual,
            max_delta = projected.max_delta,
            cap_ratio = after_cap / base_cap,
            tau_full = after_tau.full_residual,
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
            valid_cap = projected.valid_cap_residual,
            valid_cap_ratio = projected.valid_cap_residual / base_cap,
            valid_tau_full = projected.valid_tau_full,
        );
        assert!(after_cap.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_c0_source_projection_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut cap_c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut cap_c0_rows);
    let projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_c0_rows,
    );
    let before_fig6c = measure_fig6c_paper_direction_rows(&active_set.patch, map);
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let before_cap = test_macro_cap_from_sampled(valence, vec![active_set.patch; valence]);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let before_c0 = test_cap_v0_pair_c0_summary(&before_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);

    eprintln!(
        "Appendix-A Fig.6c cap-pair C0 source projection n={valence} sector={sector}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} before_fig6c={before_fig6c:e} after_fig6c={after_fig6c:e} before_c0_o={before_c0_o:e} before_c0_m={before_c0_m:e} after_c0_o={after_c0_o:e} after_c0_m={after_c0_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        before_fig6c = before_fig6c.full_residual,
        after_fig6c = after_fig6c.full_residual,
        before_c0_o = before_c0.near_origin,
        before_c0_m = before_c0.near_midpoint,
        after_c0_o = after_c0.near_origin,
        after_c0_m = after_c0.near_midpoint,
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );

    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_eq9_source_projection_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut cap_rows);
    add_cap_pair_eq8_eq9_rows(&mut cap_rows, valence);
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_rows,
        Some(0.5),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);

    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8/eq9 source projection n={valence} sector={sector}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        c0_o = after_c0.near_origin,
        c0_m = after_c0.near_midpoint,
        eq8 = after_report.worst_eq8(),
        eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
        eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
        eq9 = after_report.worst_eq9(),
        eq9_o = max_slice(&after_report.worst_eq9_near_o_per_sector),
        eq9_m = max_slice(&after_report.worst_eq9_near_m_per_sector),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );

    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_eq9_trust_region_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut cap_rows);
    add_cap_pair_eq8_eq9_rows(&mut cap_rows, valence);
    for step_limit in [0.25_f64, 0.1, 0.05] {
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &active_set.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &cap_rows,
            Some(step_limit),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
        let after_report = measure_kp20_assembly_constraints(&after_cap);
        let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "Appendix-A Fig.6c cap-pair eq8/eq9 trust step={step_limit:e}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            fig6c = after_fig6c.full_residual,
            c0_o = after_c0.near_origin,
            c0_m = after_c0.near_midpoint,
            eq8 = after_report.worst_eq8(),
            eq9 = after_report.worst_eq9(),
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
        );
        assert!(projected.linear_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_eq9_trust_region_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut cap_rows);
    add_cap_pair_residual_rows(&mut cap_rows, valence, TestCapResidualFamily::Eq9);
    for step_limit in [0.25_f64, 0.1, 0.05] {
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &active_set.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &cap_rows,
            Some(step_limit),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
        let after_report = measure_kp20_assembly_constraints(&after_cap);
        let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        eprintln!(
            "Appendix-A Fig.6c cap-pair eq9 trust step={step_limit:e}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            fig6c = after_fig6c.full_residual,
            c0_o = after_c0.near_origin,
            c0_m = after_c0.near_midpoint,
            eq8 = after_report.worst_eq8(),
            eq9 = after_report.worst_eq9(),
            eq9_o = max_slice(&after_report.worst_eq9_near_o_per_sector),
            eq9_m = max_slice(&after_report.worst_eq9_near_m_per_sector),
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
        );
        assert!(projected.linear_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_eq9_from_c0_bridge_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut eq9_rows);
    add_cap_pair_residual_rows(&mut eq9_rows, valence, TestCapResidualFamily::Eq9);
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &eq9_rows,
        Some(0.5),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq9 from C0 bridge: c0_iterations={c0_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        c0_iterations = c0_projected.iterations,
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        c0_o = after_c0.near_origin,
        c0_m = after_c0.near_midpoint,
        eq8 = after_report.worst_eq8(),
        eq9 = after_report.worst_eq9(),
        eq9_o = max_slice(&after_report.worst_eq9_near_o_per_sector),
        eq9_m = max_slice(&after_report.worst_eq9_near_m_per_sector),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_from_c0_eq9_bridge_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let mut full_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut full_rows);
    add_cap_pair_eq8_eq9_rows(&mut full_rows, valence);
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq9_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &full_rows,
        Some(0.25),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 from C0+eq9 bridge: c0_iterations={c0_iterations} eq9_iterations={eq9_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        c0_iterations = c0_projected.iterations,
        eq9_iterations = c0_eq9_projected.iterations,
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        c0_o = after_c0.near_origin,
        c0_m = after_c0.near_midpoint,
        eq8 = after_report.worst_eq8(),
        eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
        eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
        eq9 = after_report.worst_eq9(),
        eq9_o = max_slice(&after_report.worst_eq9_near_o_per_sector),
        eq9_m = max_slice(&after_report.worst_eq9_near_m_per_sector),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_orientation_variants_are_ranked() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let mut best = (f64::INFINITY, 1.0_f64, 1.0_f64, false, false);
    for origin_b_sign in [-1.0_f64, 1.0] {
        for midpoint_b_sign in [-1.0_f64, 1.0] {
            for origin_swap in [false, true] {
                for midpoint_swap in [false, true] {
                    let mut rows = Vec::new();
                    add_cap_pair_eq8_rows_with_orientation(
                        &mut rows,
                        valence,
                        origin_b_sign,
                        midpoint_b_sign,
                        origin_swap,
                        midpoint_swap,
                    );
                    let residual = scalar_row_set_residual(&c0_eq9_projected.patch, &rows);
                    if residual < best.0 {
                        best = (
                            residual,
                            origin_b_sign,
                            midpoint_b_sign,
                            origin_swap,
                            midpoint_swap,
                        );
                    }
                }
            }
        }
    }
    let mut paper_rows = Vec::new();
    add_cap_pair_residual_rows(&mut paper_rows, valence, TestCapResidualFamily::Eq8);
    let paper_residual = scalar_row_set_residual(&c0_eq9_projected.patch, &paper_rows);
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 orientation variants on C0+eq9 bridge: paper={paper_residual:e} best={best_residual:e} origin_b_sign={origin_b_sign:+e} midpoint_b_sign={midpoint_b_sign:+e} origin_swap={origin_swap} midpoint_swap={midpoint_swap}",
        best_residual = best.0,
        origin_b_sign = best.1,
        midpoint_b_sign = best.2,
        origin_swap = best.3,
        midpoint_swap = best.4,
    );
    assert!(best.0.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_row_groups_are_localized() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let bridge_cap = test_macro_cap_from_sampled(valence, vec![c0_eq9_projected.patch; valence]);
    let bridge_report = measure_kp20_assembly_constraints(&bridge_cap);
    let mut best_valid_direction = (f64::INFINITY, TestCapEq8Join::NearOrigin, 0_usize);
    let mut worst_domain_damage = (0_usize, 0_usize, TestCapEq8Join::NearOrigin, 0_usize);

    for join in [TestCapEq8Join::NearOrigin, TestCapEq8Join::NearMidpoint] {
        for row_index in 0..7 {
            let mut subset_rows = Vec::new();
            add_cap_pair_eq8_index_rows(&mut subset_rows, valence, join, row_index);
            let bridge_subset_residual =
                scalar_row_set_residual(&c0_eq9_projected.patch, &subset_rows);
            let mut rows = c0_eq9_rows.clone();
            rows.extend(subset_rows.iter().cloned());
            let (one_step_patch, max_delta, linear_residual) = solve_square_domain_active_step(
                &c0_eq9_projected.patch,
                map,
                valence,
                sector,
                &[],
                &[],
                &rows,
                Some(0.25),
                true,
            );
            let after_subset_residual = scalar_row_set_residual(&one_step_patch, &subset_rows);
            let after_cap = test_macro_cap_from_sampled(valence, vec![one_step_patch; valence]);
            let after_report = measure_kp20_assembly_constraints(&after_cap);
            let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
            let after_local = global_planar_patch_to_local(valence, sector, one_step_patch);
            let after_square =
                test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            let after_orientation =
                test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
            if after_square.outside == 0
                && after_orientation.negative == 0
                && after_subset_residual < best_valid_direction.0
            {
                best_valid_direction = (after_subset_residual, join, row_index);
            }
            let damage = after_square.outside + after_orientation.negative;
            if damage > worst_domain_damage.0 + worst_domain_damage.1 {
                worst_domain_damage = (
                    after_square.outside,
                    after_orientation.negative,
                    join,
                    row_index,
                );
            }
            eprintln!(
                "Appendix-A Fig.6c cap-pair eq8 row group {join:?}[{row_index}]: bridge_subset={bridge_subset_residual:e} after_subset={after_subset_residual:e} linear_residual={linear_residual:e} max_delta={max_delta:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
                c0_o = after_c0.near_origin,
                c0_m = after_c0.near_midpoint,
                eq8 = after_report.worst_eq8(),
                eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
                eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
                eq9 = after_report.worst_eq9(),
                square = after_square.outside,
                samples = after_square.samples,
                square_excess = after_square.max_excess,
                neg = after_orientation.negative,
                orientation_samples = after_orientation.samples,
                min_det = after_orientation.min_det,
                max_det = after_orientation.max_det,
                det_worst = after_orientation.worst_sample,
            );
            assert!(bridge_subset_residual.is_finite());
        }
    }
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 row-group summary: bridge_eq8={bridge_eq8:e} bridge_eq8_o={bridge_eq8_o:e} bridge_eq8_m={bridge_eq8_m:e} best_valid={best_valid_residual:e}@{best_valid_join:?}[{best_valid_row}] worst_damage_square={worst_square} worst_damage_neg={worst_neg}@{worst_join:?}[{worst_row}]",
        bridge_eq8 = bridge_report.worst_eq8(),
        bridge_eq8_o = max_slice(&bridge_report.worst_eq8_near_o_per_sector),
        bridge_eq8_m = max_slice(&bridge_report.worst_eq8_near_m_per_sector),
        best_valid_residual = best_valid_direction.0,
        best_valid_join = best_valid_direction.1,
        best_valid_row = best_valid_direction.2,
        worst_square = worst_domain_damage.0,
        worst_neg = worst_domain_damage.1,
        worst_join = worst_domain_damage.2,
        worst_row = worst_domain_damage.3,
    );
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_selected_row_groups_are_projected() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    for (join, row_index) in [
        (TestCapEq8Join::NearOrigin, 6_usize),
        (TestCapEq8Join::NearMidpoint, 5_usize),
        (TestCapEq8Join::NearOrigin, 4_usize),
    ] {
        let mut subset_rows = Vec::new();
        add_cap_pair_eq8_index_rows(&mut subset_rows, valence, join, row_index);
        let mut rows = c0_eq9_rows.clone();
        rows.extend(subset_rows.iter().cloned());
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_eq9_projected.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &rows,
            Some(0.25),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
        let after_report = measure_kp20_assembly_constraints(&after_cap);
        let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_subset_residual = scalar_row_set_residual(&projected.patch, &subset_rows);
        eprintln!(
            "Appendix-A Fig.6c cap-pair eq8 selected row group {join:?}[{row_index}]: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} subset={subset:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            subset = after_subset_residual,
            fig6c = after_fig6c.full_residual,
            c0_o = after_c0.near_origin,
            c0_m = after_c0.near_midpoint,
            eq8 = after_report.worst_eq8(),
            eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
            eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
            eq9 = after_report.worst_eq9(),
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
        );
        assert!(after_subset_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_join_groups_are_projected() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    for join in [TestCapEq8Join::NearOrigin, TestCapEq8Join::NearMidpoint] {
        let mut join_rows = Vec::new();
        add_cap_pair_eq8_join_rows(&mut join_rows, valence, join);
        let mut rows = c0_eq9_rows.clone();
        rows.extend(join_rows.iter().cloned());
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_eq9_projected.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &rows,
            Some(0.25),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
        let after_report = measure_kp20_assembly_constraints(&after_cap);
        let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_join_residual = scalar_row_set_residual(&projected.patch, &join_rows);
        eprintln!(
            "Appendix-A Fig.6c cap-pair eq8 join group {join:?}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} join_residual={join_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            join_residual = after_join_residual,
            fig6c = after_fig6c.full_residual,
            c0_o = after_c0.near_origin,
            c0_m = after_c0.near_midpoint,
            eq8 = after_report.worst_eq8(),
            eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
            eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
            eq9 = after_report.worst_eq9(),
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
        );
        assert!(after_join_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_near_origin_frames_are_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    log_near_origin_frame_summary("active-set", valence, sector, &active_set.patch);
    log_near_origin_frame_summary("c0-bridge", valence, sector, &c0_projected.patch);
    log_near_origin_frame_summary("c0-eq9-bridge", valence, sector, &c0_eq9_projected.patch);
    for (label, join, row_index) in [
        ("near-origin-row6", TestCapEq8Join::NearOrigin, 6_usize),
        ("near-origin-row4", TestCapEq8Join::NearOrigin, 4_usize),
        ("near-midpoint-row5", TestCapEq8Join::NearMidpoint, 5_usize),
    ] {
        let mut subset_rows = c0_eq9_rows.clone();
        add_cap_pair_eq8_index_rows(&mut subset_rows, valence, join, row_index);
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_eq9_projected.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &subset_rows,
            Some(0.25),
            true,
        );
        log_near_origin_frame_summary(label, valence, sector, &projected.patch);
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_with_near_origin_frame_rows_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let mut frame_rows = Vec::new();
    add_near_origin_right_transverse_matches_left_rows(&mut frame_rows);
    let mut rows = c0_eq9_rows.clone();
    rows.extend(frame_rows.iter().cloned());
    add_cap_pair_eq8_eq9_rows(&mut rows, valence);
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq9_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &rows,
        Some(0.25),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let frame_residual = scalar_row_set_residual(&projected.patch, &frame_rows);
    let left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 with near-origin frame rows: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} frame={frame:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        frame = frame_residual,
        fig6c = after_fig6c.full_residual,
        c0_o = after_c0.near_origin,
        c0_m = after_c0.near_midpoint,
        eq8 = after_report.worst_eq8(),
        eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
        eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
        eq9 = after_report.worst_eq9(),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
        left_det = left.determinant,
        right_det = right.determinant,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_with_fixed_near_origin_frame_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let mut frame_rows = Vec::new();
    add_fixed_near_origin_first_order_frame_rows(&mut frame_rows, &c0_eq9_projected.patch);
    let mut rows = c0_eq9_rows.clone();
    rows.extend(frame_rows.iter().cloned());
    add_cap_pair_eq8_eq9_rows(&mut rows, valence);
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq9_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &rows,
        Some(0.25),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let frame_residual = scalar_row_set_residual(&projected.patch, &frame_rows);
    let left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 with fixed near-origin frame: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} frame={frame:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        frame = frame_residual,
        fig6c = after_fig6c.full_residual,
        c0_o = after_c0.near_origin,
        c0_m = after_c0.near_midpoint,
        eq8 = after_report.worst_eq8(),
        eq8_o = max_slice(&after_report.worst_eq8_near_o_per_sector),
        eq8_m = max_slice(&after_report.worst_eq8_near_m_per_sector),
        eq9 = after_report.worst_eq9(),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
        left_det = left.determinant,
        right_det = right.determinant,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_endpoint_frame_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let local = global_planar_patch_to_local(valence, sector, c0_eq9_projected.patch);
    let previous = near_origin_first_order_frame(&local, Kp20PieceId::NearOriginLeft);
    let current = near_origin_first_order_frame(&local, Kp20PieceId::NearOriginRight);
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let b0 = kp20_b_coeffs_near_o(c)[0];
    let implied_current_dv = Point3::new(
        b0 * previous.du.x - previous.dv.x,
        b0 * previous.du.y - previous.dv.y,
        0.0,
    );
    let implied_det = current.du.x * implied_current_dv.y - current.du.y * implied_current_dv.x;
    let current_det = current.determinant;
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 endpoint frame: b0={b0:e} prev_du=({pdux:e},{pduy:e}) prev_dv=({pdvx:e},{pdvy:e}) current_du=({cdux:e},{cduy:e}) current_dv=({cdvx:e},{cdvy:e}) current_det={current_det:e} implied_current_dv=({idvx:e},{idvy:e}) implied_det={implied_det:e}",
        pdux = previous.du.x,
        pduy = previous.du.y,
        pdvx = previous.dv.x,
        pdvy = previous.dv.y,
        cdux = current.du.x,
        cduy = current.du.y,
        cdvx = current.dv.x,
        cdvy = current.dv.y,
        idvx = implied_current_dv.x,
        idvy = implied_current_dv.y,
    );
    assert!(implied_det.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_endpoint_orientation_variants_are_ranked() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_rows);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut c0_eq9_rows);
    add_cap_pair_residual_rows(&mut c0_eq9_rows, valence, TestCapResidualFamily::Eq9);
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let local = global_planar_patch_to_local(valence, sector, c0_eq9_projected.patch);
    let left = near_origin_first_order_frame(&local, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&local, Kp20PieceId::NearOriginRight);
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let paper_b0 = kp20_b_coeffs_near_o(c)[0];
    let mut best = (f64::NEG_INFINITY, false, 1.0_f64);
    for swap in [false, true] {
        for b_sign in [-1.0_f64, 1.0] {
            let (previous, current) = if swap { (right, left) } else { (left, right) };
            let b0 = b_sign * paper_b0;
            let implied_current_dv = Point3::new(
                b0 * previous.du.x - previous.dv.x,
                b0 * previous.du.y - previous.dv.y,
                0.0,
            );
            let implied_det =
                current.du.x * implied_current_dv.y - current.du.y * implied_current_dv.x;
            if implied_det > best.0 {
                best = (implied_det, swap, b_sign);
            }
            eprintln!(
                "Appendix-A Fig.6c cap-pair eq8 endpoint variant swap={swap} b_sign={b_sign:+e}: implied_dv=({dvx:e},{dvy:e}) implied_det={implied_det:e}",
                dvx = implied_current_dv.x,
                dvy = implied_current_dv.y,
            );
        }
    }
    eprintln!(
        "Appendix-A Fig.6c cap-pair eq8 endpoint best orientation variant: det={det:e} swap={swap} b_sign={b_sign:+e}",
        det = best.0,
        swap = best.1,
        b_sign = best.2,
    );
    assert!(best.0.is_finite());
}

#[test]
fn appendix_a_fig6c_reversed_current_cap_pair_source_projection_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
        Some(0.5),
        true,
    );
    let mut c0_eq9_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_eq9_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut c0_eq9_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let c0_eq9_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq9_rows,
        Some(0.5),
        true,
    );
    let mut full_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut full_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut full_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    add_cap_pair_residual_rows_with_chart(
        &mut full_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq9_projected.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &full_rows,
        Some(0.25),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![projected.patch; valence]);
    let same_direction_report = measure_kp20_assembly_constraints(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let reversed_c0 = scalar_row_set_residual(&projected.patch, &c0_rows);
    let reversed_eq8 = {
        let mut rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut rows,
            valence,
            TestCapResidualFamily::Eq8,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        scalar_row_set_residual(&projected.patch, &rows)
    };
    let reversed_eq9 = {
        let mut rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut rows,
            valence,
            TestCapResidualFamily::Eq9,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        scalar_row_set_residual(&projected.patch, &rows)
    };
    let local_left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
    let local_right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);

    eprintln!(
        "Appendix-A Fig.6c reversed-current cap-pair projection: c0_iterations={c0_iterations} eq9_iterations={eq9_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} reversed_c0={reversed_c0:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} same_eq8={same_eq8:e} same_eq9={same_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
        c0_iterations = c0_projected.iterations,
        eq9_iterations = c0_eq9_projected.iterations,
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        same_eq8 = same_direction_report.worst_eq8(),
        same_eq9 = same_direction_report.worst_eq9(),
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
        left_det = local_left.determinant,
        right_det = local_right.determinant,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_reversed_current_cap_pair_stages_are_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let stage_specs = [
        ("reversed-c0", TestCapResidualFamilySet::new(false, false)),
        (
            "reversed-c0-eq9",
            TestCapResidualFamilySet::new(false, true),
        ),
        (
            "reversed-c0-eq8",
            TestCapResidualFamilySet::new(true, false),
        ),
        (
            "reversed-c0-eq8-eq9",
            TestCapResidualFamilySet::new(true, true),
        ),
    ];
    for (label, families) in stage_specs {
        let mut rows = Vec::new();
        add_cap_pair_v0_c0_rows_with_chart(&mut rows, TestCapPairBoundaryChart::ReverseCurrentU);
        if families.eq8 {
            add_cap_pair_residual_rows_with_chart(
                &mut rows,
                valence,
                TestCapResidualFamily::Eq8,
                TestCapPairBoundaryChart::ReverseCurrentU,
            );
        }
        if families.eq9 {
            add_cap_pair_residual_rows_with_chart(
                &mut rows,
                valence,
                TestCapResidualFamily::Eq9,
                TestCapPairBoundaryChart::ReverseCurrentU,
            );
        }
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &active_set.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &rows,
            Some(0.25),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let mut reversed_c0_rows = Vec::new();
        add_cap_pair_v0_c0_rows_with_chart(
            &mut reversed_c0_rows,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let mut reversed_eq8_rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut reversed_eq8_rows,
            valence,
            TestCapResidualFamily::Eq8,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let mut reversed_eq9_rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut reversed_eq9_rows,
            valence,
            TestCapResidualFamily::Eq9,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);
        eprintln!(
            "Appendix-A Fig.6c reversed-current stage {label}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} reversed_c0={reversed_c0:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            fig6c = after_fig6c.full_residual,
            reversed_c0 = scalar_row_set_residual(&projected.patch, &reversed_c0_rows),
            reversed_eq8 = scalar_row_set_residual(&projected.patch, &reversed_eq8_rows),
            reversed_eq9 = scalar_row_set_residual(&projected.patch, &reversed_eq9_rows),
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(projected.linear_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_reversed_current_cap_pair_staged_orders_are_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
        Some(0.25),
        true,
    );
    for first_family in [TestCapResidualFamily::Eq8, TestCapResidualFamily::Eq9] {
        let second_family = match first_family {
            TestCapResidualFamily::Eq8 => TestCapResidualFamily::Eq9,
            TestCapResidualFamily::Eq9 => TestCapResidualFamily::Eq8,
        };
        let mut first_rows = c0_rows.clone();
        add_cap_pair_residual_rows_with_chart(
            &mut first_rows,
            valence,
            first_family,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let first_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_projected.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &first_rows,
            Some(0.25),
            true,
        );
        let mut second_rows = first_rows.clone();
        add_cap_pair_residual_rows_with_chart(
            &mut second_rows,
            valence,
            second_family,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &first_projected.patch,
            map,
            valence,
            sector,
            &TAU_VALIDITY_SAMPLES,
            &second_rows,
            Some(0.25),
            true,
        );
        let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
        let oriented_cap = test_reversed_current_cap_residual_summary(
            &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
            valence,
        );
        let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
        let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let after_orientation =
            test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
        let mut reversed_eq8_rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut reversed_eq8_rows,
            valence,
            TestCapResidualFamily::Eq8,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let mut reversed_eq9_rows = Vec::new();
        add_cap_pair_residual_rows_with_chart(
            &mut reversed_eq9_rows,
            valence,
            TestCapResidualFamily::Eq9,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);
        eprintln!(
            "Appendix-A Fig.6c reversed-current staged order first={first_family:?}: c0_iterations={c0_iterations} first_iterations={first_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
            c0_iterations = c0_projected.iterations,
            first_iterations = first_projected.iterations,
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            fig6c = after_fig6c.full_residual,
            reversed_eq8 = scalar_row_set_residual(&projected.patch, &reversed_eq8_rows),
            reversed_eq9 = scalar_row_set_residual(&projected.patch, &reversed_eq9_rows),
            oriented_cap_c0 = oriented_cap.c0,
            oriented_cap_eq8 = oriented_cap.eq8,
            oriented_cap_eq9 = oriented_cap.eq9,
            square = after_square.outside,
            samples = after_square.samples,
            square_excess = after_square.max_excess,
            neg = after_orientation.negative,
            orientation_samples = after_orientation.samples,
            min_det = after_orientation.min_det,
            max_det = after_orientation.max_det,
            det_worst = after_orientation.worst_sample,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(projected.linear_residual.is_finite());
    }
}

#[test]
fn appendix_a_reversed_current_without_fig6c_base_rows_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let empty_base_rows: Vec<Kp20PlanarTauScalarLinearConstraint> = Vec::new();
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    let c0_projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
        Some(0.25),
        true,
    );
    let mut c0_eq8_rows = c0_rows.clone();
    add_cap_pair_residual_rows_with_chart(
        &mut c0_eq8_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let c0_eq8_projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &c0_projected.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq8_rows,
        Some(0.25),
        true,
    );
    let mut full_rows = c0_eq8_rows.clone();
    add_cap_pair_residual_rows_with_chart(
        &mut full_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &c0_eq8_projected.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &full_rows,
        Some(0.25),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&projected.patch, map);
    let oriented_cap = test_reversed_current_cap_residual_summary(
        &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
        valence,
    );
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let mut reversed_eq8_rows = Vec::new();
    add_cap_pair_residual_rows_with_chart(
        &mut reversed_eq8_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let mut reversed_eq9_rows = Vec::new();
    add_cap_pair_residual_rows_with_chart(
        &mut reversed_eq9_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let left = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&after_local, Kp20PieceId::NearOriginRight);
    eprintln!(
        "Appendix-A reversed-current without Fig.6c base rows: c0_iterations={c0_iterations} eq8_iterations={eq8_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst} left_det={left_det:e} right_det={right_det:e}",
        c0_iterations = c0_projected.iterations,
        eq8_iterations = c0_eq8_projected.iterations,
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        reversed_eq8 = scalar_row_set_residual(&projected.patch, &reversed_eq8_rows),
        reversed_eq9 = scalar_row_set_residual(&projected.patch, &reversed_eq9_rows),
        oriented_cap_c0 = oriented_cap.c0,
        oriented_cap_eq8 = oriented_cap.eq8,
        oriented_cap_eq9 = oriented_cap.eq9,
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
        left_det = left.determinant,
        right_det = right.determinant,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_empty_base_tau_fig6c_strip_graph_is_ranked() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, _) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let empty_base_rows: Vec<Kp20PlanarTauScalarLinearConstraint> = Vec::new();
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    let c0_projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_rows,
        Some(0.25),
        true,
    );
    let mut c0_eq8_rows = c0_rows.clone();
    add_cap_pair_residual_rows_with_chart(
        &mut c0_eq8_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let c0_eq8_projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &c0_projected.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &c0_eq8_rows,
        Some(0.25),
        true,
    );
    let mut full_rows = c0_eq8_rows.clone();
    add_cap_pair_residual_rows_with_chart(
        &mut full_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &c0_eq8_projected.patch,
        &empty_base_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &full_rows,
        Some(0.25),
        true,
    );

    let strips = appendix_a_candidate_section42_strips();
    let exact_pairs = exact_strip_boundary_pairs(&projected.patch, &strips, angle);
    let sector_pairs = exact_pairs
        .iter()
        .copied()
        .filter(|&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::SectorJoin)
            )
        })
        .collect::<Vec<_>>();
    let mut best = (
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        None,
    );
    let mut candidate_count = 0_usize;
    let mut exact_internal_count = 0_usize;
    for &(p, arrow_p) in &sector_pairs {
        for &(underlined_p, arrow_underlined_p) in &sector_pairs {
            if !fig6c_semantic_symbols_are_distinct([p, arrow_p, underlined_p, arrow_underlined_p])
            {
                continue;
            }
            candidate_count += 1;
            let map = TestSection42StripSymbolMap {
                p,
                arrow_p,
                underlined_p,
                arrow_underlined_p,
                previous_sector_rotation_radians: angle,
            };
            let plain_join =
                best_strip_eq12_join_residual(&projected.patch, p, underlined_p, angle);
            let arrow_join =
                best_strip_eq12_join_residual(&projected.patch, arrow_p, arrow_underlined_p, angle);
            if plain_join.best_endpoint_residual < 1.0e-12
                && arrow_join.best_endpoint_residual < 1.0e-12
            {
                exact_internal_count += 1;
            }
            let rows = test_section42_strip_basis_rows(map, 2.0 * c);
            let full_residual = scalar_row_set_residual(&projected.patch, &rows);
            let c0_residual =
                scalar_row_set_residual(&projected.patch, &test_section42_strip_c0_rows(map));
            let mut eq11_rows = Vec::new();
            test_add_section42_strip_eq11_degree5_rows(&mut eq11_rows, map);
            let mut eq12_rows = Vec::new();
            test_add_section42_strip_paper_direction_eq12_rows(&mut eq12_rows, map);
            let mut eq13_rows = Vec::new();
            test_add_section42_strip_eq13_rows(&mut eq13_rows, map, 2.0 * c);
            let eq11_residual = scalar_row_set_residual(&projected.patch, &eq11_rows);
            let eq12_residual = scalar_row_set_residual(&projected.patch, &eq12_rows);
            let eq13_residual = scalar_row_set_residual(&projected.patch, &eq13_rows);
            if full_residual < best.0 {
                best = (
                    full_residual,
                    c0_residual.max(eq11_residual),
                    eq12_residual,
                    eq13_residual,
                    Some((map, plain_join, arrow_join)),
                );
            }
        }
    }
    if let Some((map, plain_join, arrow_join)) = best.4 {
        eprintln!(
            "Appendix-A empty-base tau Fig.6c strip graph best full={full:e} c0_eq11={c0_eq11:e} eq12={eq12:e} eq13={eq13:e} exact_pairs={exact_pairs} sector_pairs={sector_pairs} candidates={candidates} exact_internal={exact_internal} map={map:?} plain_join={plain_join:?} arrow_join={arrow_join:?}",
            full = best.0,
            c0_eq11 = best.1,
            eq12 = best.2,
            eq13 = best.3,
            exact_pairs = exact_pairs.len(),
            sector_pairs = sector_pairs.len(),
            candidates = candidate_count,
            exact_internal = exact_internal_count,
        );
    } else {
        eprintln!(
            "Appendix-A empty-base tau Fig.6c strip graph found no sector-pair candidate: exact_pairs={} sector_pairs={}",
            exact_pairs.len(),
            sector_pairs.len(),
        );
    }
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_oriented_boundary_section42_basis_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, paper_map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let oriented_rows = test_section42_strip_basis_rows(oriented_map, 2.0 * c);
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &oriented_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &[],
        Some(0.25),
        true,
    );
    let report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
    let paper_report = measure_fig6c_paper_direction_rows(&projected.patch, paper_map);
    let oriented_cap = test_reversed_current_cap_residual_summary(
        &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
        valence,
    );
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A oriented-boundary Section 4.2 basis: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} paper_full={paper_full:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        full = report.full_residual,
        c0 = report.c0_residual,
        eq11 = report.eq11_residual,
        eq12 = report.eq12_residual,
        eq13 = report.eq13_residual,
        paper_full = paper_report.full_residual,
        oriented_cap_c0 = oriented_cap.c0,
        oriented_cap_eq8 = oriented_cap.eq8,
        oriented_cap_eq9 = oriented_cap.eq9,
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_oriented_boundary_eq12_variants_are_ranked() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, _) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let mut best = (
        f64::INFINITY,
        false,
        0_usize,
        0_usize,
        0_usize,
        0_usize,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for source_is_underlined in [false, true] {
        for source_endpoint in [0_usize, 6] {
            for receiving_endpoint in [0_usize, 6] {
                let mut rows = test_section42_strip_c0_rows(oriented_map);
                test_add_section42_strip_eq11_degree5_rows(&mut rows, oriented_map);
                test_add_section42_strip_eq12_variant_rows(
                    &mut rows,
                    oriented_map,
                    source_is_underlined,
                    source_endpoint,
                    receiving_endpoint,
                );
                test_add_section42_strip_eq13_rows(&mut rows, oriented_map, 2.0 * c);
                let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
                    &active_set.patch,
                    &rows,
                    valence,
                    sector,
                    &TAU_VALIDITY_SAMPLES,
                    &[],
                    Some(0.25),
                    true,
                );
                let report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
                let oriented_cap = test_reversed_current_cap_residual_summary(
                    &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
                    valence,
                );
                let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
                let after_square =
                    test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
                let after_orientation =
                    test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
                eprintln!(
                    "Appendix-A oriented-boundary eq12 variant source_is_underlined={source_is_underlined} source_endpoint={source_endpoint} receiving_endpoint={receiving_endpoint}: iterations={iterations} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} cap_eq8={cap_eq8:e} cap_eq9={cap_eq9:e} square={square}/{samples} neg={neg}/{orientation_samples} min_det={min_det:e}",
                    iterations = projected.iterations,
                    full = report.full_residual,
                    c0 = report.c0_residual,
                    eq11 = report.eq11_residual,
                    eq12 = report.eq12_residual,
                    eq13 = report.eq13_residual,
                    cap_eq8 = oriented_cap.eq8,
                    cap_eq9 = oriented_cap.eq9,
                    square = after_square.outside,
                    samples = after_square.samples,
                    neg = after_orientation.negative,
                    orientation_samples = after_orientation.samples,
                    min_det = after_orientation.min_det,
                );
                let validity_penalty = (after_square.outside + after_orientation.negative) as f64;
                let score = validity_penalty + report.full_residual;
                if score < best.0 {
                    best = (
                        score,
                        source_is_underlined,
                        source_endpoint,
                        receiving_endpoint,
                        after_square.outside,
                        after_orientation.negative,
                        report.full_residual,
                        oriented_cap.eq8,
                        oriented_cap.eq9,
                    );
                }
                assert!(projected.linear_residual.is_finite());
            }
        }
    }
    eprintln!(
        "Appendix-A oriented-boundary eq12 variant best score={score:e} source_is_underlined={source_is_underlined} source_endpoint={source_endpoint} receiving_endpoint={receiving_endpoint} square={square} neg={neg} full={full:e} cap_eq8={cap_eq8:e} cap_eq9={cap_eq9:e}",
        score = best.0,
        source_is_underlined = best.1,
        source_endpoint = best.2,
        receiving_endpoint = best.3,
        square = best.4,
        neg = best.5,
        full = best.6,
        cap_eq8 = best.7,
        cap_eq9 = best.8,
    );
}

#[test]
fn appendix_a_oriented_boundary_section42_with_reversed_cap_rows_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, _) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let mut section42_rows = test_section42_strip_c0_rows(oriented_map);
    test_add_section42_strip_eq11_degree5_rows(&mut section42_rows, oriented_map);
    test_add_section42_strip_eq12_variant_rows(&mut section42_rows, oriented_map, false, 6, 0);
    test_add_section42_strip_eq13_rows(&mut section42_rows, oriented_map, 2.0 * c);

    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut cap_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &section42_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_rows,
        Some(0.25),
        true,
    );
    let report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
    let oriented_cap = test_reversed_current_cap_residual_summary(
        &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
        valence,
    );
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A oriented-boundary Section 4.2 plus reversed cap rows: iterations={iterations} active={active} orientation_active={orientation_active} linear_residual={linear_residual:e} full={full:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} cap_rows={cap_rows:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        linear_residual = projected.linear_residual,
        full = report.full_residual,
        c0 = report.c0_residual,
        eq11 = report.eq11_residual,
        eq12 = report.eq12_residual,
        eq13 = report.eq13_residual,
        cap_rows = scalar_row_set_residual(&projected.patch, &cap_rows),
        oriented_cap_c0 = oriented_cap.c0,
        oriented_cap_eq8 = oriented_cap.eq8,
        oriented_cap_eq9 = oriented_cap.eq9,
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_oriented_boundary_eq12_variants_with_reversed_cap_rows_are_ranked() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, _) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut cap_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let mut best = (
        f64::INFINITY,
        false,
        0_usize,
        0_usize,
        0_usize,
        0_usize,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for source_is_underlined in [false, true] {
        for source_endpoint in [0_usize, 6] {
            for receiving_endpoint in [0_usize, 6] {
                let mut section42_rows = test_section42_strip_c0_rows(oriented_map);
                test_add_section42_strip_eq11_degree5_rows(&mut section42_rows, oriented_map);
                test_add_section42_strip_eq12_variant_rows(
                    &mut section42_rows,
                    oriented_map,
                    source_is_underlined,
                    source_endpoint,
                    receiving_endpoint,
                );
                test_add_section42_strip_eq13_rows(&mut section42_rows, oriented_map, 2.0 * c);
                let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
                    &active_set.patch,
                    &section42_rows,
                    valence,
                    sector,
                    &TAU_VALIDITY_SAMPLES,
                    &cap_rows,
                    Some(0.25),
                    true,
                );
                let report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
                let oriented_cap = test_reversed_current_cap_residual_summary(
                    &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
                    valence,
                );
                let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
                let after_square =
                    test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
                let after_orientation =
                    test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
                let cap_row_residual = scalar_row_set_residual(&projected.patch, &cap_rows);
                eprintln!(
                    "Appendix-A oriented-boundary eq12+cap variant source_is_underlined={source_is_underlined} source_endpoint={source_endpoint} receiving_endpoint={receiving_endpoint}: iterations={iterations} full={full:e} eq12={eq12:e} cap_rows={cap_rows:e} cap_eq8={cap_eq8:e} cap_eq9={cap_eq9:e} square={square}/{samples} neg={neg}/{orientation_samples} min_det={min_det:e}",
                    iterations = projected.iterations,
                    full = report.full_residual,
                    eq12 = report.eq12_residual,
                    cap_rows = cap_row_residual,
                    cap_eq8 = oriented_cap.eq8,
                    cap_eq9 = oriented_cap.eq9,
                    square = after_square.outside,
                    samples = after_square.samples,
                    neg = after_orientation.negative,
                    orientation_samples = after_orientation.samples,
                    min_det = after_orientation.min_det,
                );
                let validity_penalty = (after_square.outside + after_orientation.negative) as f64;
                let score = validity_penalty + report.full_residual + cap_row_residual;
                if score < best.0 {
                    best = (
                        score,
                        source_is_underlined,
                        source_endpoint,
                        receiving_endpoint,
                        after_square.outside,
                        after_orientation.negative,
                        report.full_residual,
                        cap_row_residual,
                        after_orientation.min_det,
                    );
                }
                assert!(projected.linear_residual.is_finite());
            }
        }
    }
    eprintln!(
        "Appendix-A oriented-boundary eq12+cap variant best score={score:e} source_is_underlined={source_is_underlined} source_endpoint={source_endpoint} receiving_endpoint={receiving_endpoint} square={square} neg={neg} full={full:e} cap_rows={cap_rows:e} min_det={min_det:e}",
        score = best.0,
        source_is_underlined = best.1,
        source_endpoint = best.2,
        receiving_endpoint = best.3,
        square = best.4,
        neg = best.5,
        full = best.6,
        cap_rows = best.7,
        min_det = best.8,
    );
}

#[test]
fn appendix_a_oriented_boundary_cap_compatible_eq12_rows_are_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, _) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let mut section42_rows = test_section42_strip_c0_rows(oriented_map);
    test_add_section42_strip_eq11_degree5_rows(&mut section42_rows, oriented_map);
    test_add_section42_strip_eq12_variant_rows(&mut section42_rows, oriented_map, false, 6, 0);
    test_add_section42_strip_eq13_rows(&mut section42_rows, oriented_map, 2.0 * c);
    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut cap_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &section42_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_rows,
        Some(0.25),
        true,
    );
    let mut eq12_rows = Vec::new();
    let mut labels = Vec::new();
    for (source, receiving, family) in [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
            "plain",
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            "arrow",
        ),
    ] {
        let start = eq12_rows.len();
        for j in 0..=2 {
            test_add_section42_strip_vector_rows(
                &mut eq12_rows,
                oriented_map,
                &[
                    section42_term(receiving, 0, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                ],
            );
            labels.push(format!("{family} endpoint j={j} x"));
            labels.push(format!("{family} endpoint j={j} y"));
        }
        for j in 1..=2 {
            test_add_section42_strip_vector_rows(
                &mut eq12_rows,
                oriented_map,
                &[
                    section42_term(receiving, 5, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, -0.25),
                    section42_term(source, 2, j, 0.25),
                ],
            );
            labels.push(format!("{family} row5 j={j} x"));
            labels.push(format!("{family} row5 j={j} y"));
            test_add_section42_strip_vector_rows(
                &mut eq12_rows,
                oriented_map,
                &[
                    section42_term(receiving, 1, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, 0.25),
                    section42_term(source, 2, j, -0.25),
                ],
            );
            labels.push(format!("{family} row1 j={j} x"));
            labels.push(format!("{family} row1 j={j} y"));
        }
        assert_eq!(eq12_rows.len() - start, 14);
    }
    let mut ranked = eq12_rows
        .iter()
        .enumerate()
        .map(|(index, row)| {
            (
                scalar_residual(&projected.patch, row).abs(),
                labels[index].clone(),
            )
        })
        .collect::<Vec<_>>();
    ranked.sort_by(|left, right| {
        right
            .0
            .partial_cmp(&left.0)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    for (rank, (residual, label)) in ranked.iter().take(12).enumerate() {
        eprintln!(
            "Appendix-A oriented-boundary cap-compatible eq12 row rank={rank} residual={residual:e} label={label}"
        );
    }
    let report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
    eprintln!(
        "Appendix-A oriented-boundary cap-compatible eq12 localized full={full:e} eq12={eq12:e} cap_rows={cap_rows:e}",
        full = report.full_residual,
        eq12 = report.eq12_residual,
        cap_rows = scalar_row_set_residual(&projected.patch, &cap_rows),
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_endpoint_aware_section42_with_reversed_cap_rows_closes() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, paper_map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let oriented_map = TestSection42StripSymbolMap {
        p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
            edge: TestMacroEdge::J0,
            reversed: false,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: TestSection42StripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
            edge: TestMacroEdge::J0,
            reversed: true,
            sector_rotation_steps: 0,
        },
        previous_sector_rotation_radians: angle,
    };
    let mut section42_rows = test_section42_strip_c0_rows(oriented_map);
    test_add_section42_strip_eq11_degree5_rows(&mut section42_rows, oriented_map);
    test_add_section42_strip_endpoint_aware_eq12_rows(&mut section42_rows, oriented_map);
    test_add_section42_strip_eq13_rows(&mut section42_rows, oriented_map, 2.0 * c);

    let mut cap_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut cap_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    add_cap_pair_residual_rows_with_chart(
        &mut cap_rows,
        valence,
        TestCapResidualFamily::Eq9,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let projected = square_domain_active_set_projection_with_base_and_extra_rows_solver(
        &active_set.patch,
        &section42_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_rows,
        Some(0.25),
        true,
    );
    let endpoint_report =
        measure_fig6c_endpoint_aware_section42_rows(&projected.patch, oriented_map);
    let paper_report = measure_fig6c_paper_direction_rows(&projected.patch, paper_map);
    let fixed_oriented_report = measure_fig6c_paper_direction_rows(&projected.patch, oriented_map);
    let oriented_cap = test_reversed_current_cap_residual_summary(
        &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
        valence,
    );
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A endpoint-aware Section 4.2 plus reversed cap rows: iterations={iterations} active={active} orientation_active={orientation_active} linear_residual={linear_residual:e} endpoint_full={endpoint_full:e} endpoint_c0={endpoint_c0:e} endpoint_eq11={endpoint_eq11:e} endpoint_eq12={endpoint_eq12:e} endpoint_eq13={endpoint_eq13:e} fixed_oriented_full={fixed_oriented_full:e} paper_full={paper_full:e} cap_rows={cap_rows:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        linear_residual = projected.linear_residual,
        endpoint_full = endpoint_report.full_residual,
        endpoint_c0 = endpoint_report.c0_residual,
        endpoint_eq11 = endpoint_report.eq11_residual,
        endpoint_eq12 = endpoint_report.eq12_residual,
        endpoint_eq13 = endpoint_report.eq13_residual,
        fixed_oriented_full = fixed_oriented_report.full_residual,
        paper_full = paper_report.full_residual,
        cap_rows = scalar_row_set_residual(&projected.patch, &cap_rows),
        oriented_cap_c0 = oriented_cap.c0,
        oriented_cap_eq8 = oriented_cap.eq8,
        oriented_cap_eq9 = oriented_cap.eq9,
        square = after_square.outside,
        samples = after_square.samples,
        square_excess = after_square.max_excess,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
        det_worst = after_orientation.worst_sample,
    );
    assert!(projected.linear_residual.is_finite());
}
