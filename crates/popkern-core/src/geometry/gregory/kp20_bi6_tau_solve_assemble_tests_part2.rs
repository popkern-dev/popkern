use super::*;

#[test]
fn appendix_a_promoted_oriented_section42_rows_match_endpoint_diagnostic() {
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
    let promoted_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let mut diagnostic_rows = test_section42_strip_c0_rows(oriented_map);
    test_add_section42_strip_eq11_degree5_rows(&mut diagnostic_rows, oriented_map);
    test_add_section42_strip_endpoint_aware_eq12_rows(&mut diagnostic_rows, oriented_map);
    test_add_section42_strip_eq13_rows(&mut diagnostic_rows, oriented_map, 2.0 * c);
    assert_eq!(promoted_rows.len(), diagnostic_rows.len());

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
        &promoted_rows,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &cap_rows,
        Some(0.25),
        true,
    );
    let endpoint_report =
        measure_fig6c_endpoint_aware_section42_rows(&projected.patch, oriented_map);
    let promoted_residual = scalar_row_set_residual(&projected.patch, &promoted_rows);
    let diagnostic_residual = scalar_row_set_residual(&projected.patch, &diagnostic_rows);
    let oriented_cap = test_reversed_current_cap_residual_summary(
        &test_macro_cap_from_sampled(valence, vec![projected.patch; valence]),
        valence,
    );
    let after_local = global_planar_patch_to_local(valence, sector, projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    eprintln!(
        "Appendix-A promoted oriented Section 4.2 rows: iterations={iterations} promoted={promoted:e} diagnostic={diagnostic:e} endpoint_full={endpoint_full:e} endpoint_eq12={endpoint_eq12:e} cap_rows={cap_rows:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        iterations = projected.iterations,
        promoted = promoted_residual,
        diagnostic = diagnostic_residual,
        endpoint_full = endpoint_report.full_residual,
        endpoint_eq12 = endpoint_report.eq12_residual,
        cap_rows = scalar_row_set_residual(&projected.patch, &cap_rows),
        oriented_cap_eq8 = oriented_cap.eq8,
        oriented_cap_eq9 = oriented_cap.eq9,
        square = after_square.outside,
        samples = after_square.samples,
        neg = after_orientation.negative,
        orientation_samples = after_orientation.samples,
        min_det = after_orientation.min_det,
        max_det = after_orientation.max_det,
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_reversed_current_paper_quadrant_seed_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (_, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let seed = test_paper_quadrant_affine_macro_patch();
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    let c0_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &seed,
        map,
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
    let c0_eq8_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
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
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq8_projected.patch,
        map,
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
    eprintln!(
        "Appendix-A Fig.6c reversed-current paper-quadrant seed: c0_iterations={c0_iterations} eq8_iterations={eq8_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
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
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_reversed_current_with_quadrant_edges_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    let mut c0_rows = Vec::new();
    add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
    add_paper_quadrant_internal_c0_rows(&mut c0_rows);
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
    let mut c0_eq8_rows = c0_rows.clone();
    add_cap_pair_residual_rows_with_chart(
        &mut c0_eq8_rows,
        valence,
        TestCapResidualFamily::Eq8,
        TestCapPairBoundaryChart::ReverseCurrentU,
    );
    let c0_eq8_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_projected.patch,
        map,
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
    let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &c0_eq8_projected.patch,
        map,
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
    let mut quadrant_rows = Vec::new();
    add_paper_quadrant_internal_c0_rows(&mut quadrant_rows);
    eprintln!(
        "Appendix-A Fig.6c reversed-current with quadrant edges: c0_iterations={c0_iterations} eq8_iterations={eq8_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} quadrant_c0={quadrant_c0:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        c0_iterations = c0_projected.iterations,
        eq8_iterations = c0_eq8_projected.iterations,
        iterations = projected.iterations,
        active = projected.active_constraints,
        orientation_active = projected.orientation_constraints,
        max_delta = projected.max_delta,
        linear_residual = projected.linear_residual,
        fig6c = after_fig6c.full_residual,
        quadrant_c0 = scalar_row_set_residual(&projected.patch, &quadrant_rows),
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
    );
    assert!(projected.linear_residual.is_finite());
}

#[test]
fn appendix_a_fig6c_reversed_current_quadrant_edge_families_are_ranked() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");

    for edge in TestPaperQuadrantInternalEdge::ALL {
        let mut c0_rows = Vec::new();
        add_cap_pair_v0_c0_rows_with_chart(&mut c0_rows, TestCapPairBoundaryChart::ReverseCurrentU);
        add_paper_quadrant_internal_c0_edge_rows(&mut c0_rows, edge);
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
        let mut c0_eq8_rows = c0_rows.clone();
        add_cap_pair_residual_rows_with_chart(
            &mut c0_eq8_rows,
            valence,
            TestCapResidualFamily::Eq8,
            TestCapPairBoundaryChart::ReverseCurrentU,
        );
        let c0_eq8_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_projected.patch,
            map,
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
        let projected = fig6c_square_domain_active_set_projection_with_rows_solver(
            &c0_eq8_projected.patch,
            map,
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
        let mut edge_rows = Vec::new();
        add_paper_quadrant_internal_c0_edge_rows(&mut edge_rows, edge);
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
        eprintln!(
            "Appendix-A Fig.6c reversed-current quadrant edge {edge:?}: c0_iterations={c0_iterations} eq8_iterations={eq8_iterations} iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} edge_c0={edge_c0:e} reversed_eq8={reversed_eq8:e} reversed_eq9={reversed_eq9:e} oriented_cap_c0={oriented_cap_c0:e} oriented_cap_eq8={oriented_cap_eq8:e} oriented_cap_eq9={oriented_cap_eq9:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
            c0_iterations = c0_projected.iterations,
            eq8_iterations = c0_eq8_projected.iterations,
            iterations = projected.iterations,
            active = projected.active_constraints,
            orientation_active = projected.orientation_constraints,
            max_delta = projected.max_delta,
            linear_residual = projected.linear_residual,
            fig6c = after_fig6c.full_residual,
            edge_c0 = scalar_row_set_residual(&projected.patch, &edge_rows),
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
        );
        assert!(projected.linear_residual.is_finite());
    }
}

#[test]
fn appendix_a_fig6c_cap_pair_eq8_eq9_source_projection_split_is_measured() {
    let valence = 5;
    let sector = 0;
    const TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let (active_set, map) =
        fig6c_active_set_source_projection_with_map(valence, sector, &TAU_VALIDITY_SAMPLES)
            .expect("active-set source projection");
    for family in [TestCapResidualFamily::Eq8, TestCapResidualFamily::Eq9] {
        let mut cap_rows = Vec::new();
        add_cap_pair_v0_c0_rows(&mut cap_rows);
        add_cap_pair_residual_rows(&mut cap_rows, valence, family);
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
            "Appendix-A Fig.6c cap-pair {family:?} source projection n={valence} sector={sector}: iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq8_o={eq8_o:e} eq8_m={eq8_m:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
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
}

#[test]
fn appendix_a_fig6c_cap_pair_eq9_diagnostic_terms_are_ranked() {
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
    let variants = test_eq9_term_source_variants();
    let mut best = (
        f64::INFINITY,
        Kp20Eq9DiagnosticTerms::paper(),
        Kp20Eq9DiagnosticTerms::paper(),
    );
    let mut best_origin = (f64::INFINITY, Kp20Eq9DiagnosticTerms::paper());
    let mut best_midpoint = (f64::INFINITY, Kp20Eq9DiagnosticTerms::paper());
    for &origin_terms in &variants {
        let origin_residual =
            cap_pair_eq9_terms_residual(&c0_projected.patch, valence, origin_terms, None);
        if origin_residual < best_origin.0 {
            best_origin = (origin_residual, origin_terms);
        }
        for &midpoint_terms in &variants {
            if origin_terms == variants[0] {
                let midpoint_residual = cap_pair_eq9_terms_residual(
                    &c0_projected.patch,
                    valence,
                    Kp20Eq9DiagnosticTerms::paper(),
                    Some(midpoint_terms),
                );
                if midpoint_residual < best_midpoint.0 {
                    best_midpoint = (midpoint_residual, midpoint_terms);
                }
            }
            let mut rows = Vec::new();
            add_cap_pair_eq9_rows_with_terms(&mut rows, valence, origin_terms, midpoint_terms);
            let residual = scalar_row_set_residual(&c0_projected.patch, &rows);
            if residual < best.0 {
                best = (residual, origin_terms, midpoint_terms);
            }
        }
    }
    let mut best_rows = Vec::new();
    add_cap_pair_v0_c0_rows(&mut best_rows);
    add_cap_pair_eq9_rows_with_terms(&mut best_rows, valence, best.1, best.2);
    let best_projected = fig6c_square_domain_active_set_projection_with_rows_solver(
        &active_set.patch,
        map,
        valence,
        sector,
        &TAU_VALIDITY_SAMPLES,
        &best_rows,
        Some(0.5),
        true,
    );
    let after_fig6c = measure_fig6c_paper_direction_rows(&best_projected.patch, map);
    let after_cap = test_macro_cap_from_sampled(valence, vec![best_projected.patch; valence]);
    let after_report = measure_kp20_assembly_constraints(&after_cap);
    let after_c0 = test_cap_v0_pair_c0_summary(&after_cap);
    let after_local = global_planar_patch_to_local(valence, sector, best_projected.patch);
    let after_square = test_square_domain_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);
    let after_orientation =
        test_orientation_detail_report(&after_local, &TAU_VALIDITY_SAMPLES, 1.0e-10);

    eprintln!(
        "Appendix-A Fig.6c cap-pair eq9 diagnostic terms n={valence} sector={sector}: c0_iterations={c0_iterations} best_residual={best_residual:e} best_origin={best_origin:?} best_midpoint={best_midpoint:?} best_origin_only={best_origin_only:e} best_midpoint_only={best_midpoint_only:e} projection_iterations={iterations} active={active} orientation_active={orientation_active} max_delta={max_delta:e} linear_residual={linear_residual:e} fig6c={fig6c:e} c0_o={c0_o:e} c0_m={c0_m:e} eq8={eq8:e} eq9={eq9:e} eq9_o={eq9_o:e} eq9_m={eq9_m:e} square={square}/{samples} square_excess={square_excess:e} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] det_worst={det_worst}",
        c0_iterations = c0_projected.iterations,
        best_residual = best.0,
        best_origin = best.1,
        best_midpoint = best.2,
        best_origin_only = best_origin.0,
        best_midpoint_only = best_midpoint.0,
        iterations = best_projected.iterations,
        active = best_projected.active_constraints,
        orientation_active = best_projected.orientation_constraints,
        max_delta = best_projected.max_delta,
        linear_residual = best_projected.linear_residual,
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

    assert!(best.0.is_finite());
}

#[test]
fn appendix_a_fig6c_cap_pair_eq9_row_space_is_audited() {
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
    add_cap_pair_residual_rows(&mut eq9_rows, valence, TestCapResidualFamily::Eq9);
    let active_audit = test_scalar_row_space_audit(&active_set.patch, &eq9_rows, 100.0);
    let c0_audit = test_scalar_row_space_audit(&c0_projected.patch, &eq9_rows, 100.0);
    let c0_eq9_residual = scalar_row_set_residual(&c0_projected.patch, &eq9_rows);
    let c0_summary = test_cap_v0_pair_c0_summary(&test_macro_cap_from_sampled(
        valence,
        vec![c0_projected.patch; valence],
    ));

    eprintln!(
        "Appendix-A Fig.6c cap-pair eq9 row-space audit n={valence} sector={sector}: active={active:?} c0_bridge={c0:?} c0_eq9_residual={c0_eq9_residual:e} c0_o={c0_o:e} c0_m={c0_m:e}",
        active = active_audit,
        c0 = c0_audit,
        c0_o = c0_summary.near_origin,
        c0_m = c0_summary.near_midpoint,
    );

    assert!(active_audit.max_rhs.is_finite());
    assert!(c0_audit.max_rhs.is_finite());
}

#[test]
fn rotational_endpoint_aware_tensor_projection_collapse_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let tensor_border = planar_regular_sector_tensor_border(valence);

    let projected =
        solve_kp20_planar_tau_rotational_paper_sector_f3_diagnosed(&Kp20PlanarTauPaperSolveInput {
            valence,
            sector,
            characteristic_tensor_border: &tensor_border,
            initial: planar_tau_seed_macro_patch(valence, sector),
            regularization_weight: 1.0e-10,
        })
        .expect("projected tensor-band solve succeeds");
    let fixed_row_set =
        kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
            valence,
            sector,
            &tensor_border,
            &[],
        )
        .expect("fixed tensor-target row set builds");
    let fixed = solve_kp20_planar_tau_f3_with_scalar_constraints_diagnosed(
        &Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints: fixed_row_set.vector_constraints,
            scalar_constraints: fixed_row_set.scalar_constraints,
            regularization_weight: 1.0e-10,
        },
    )
    .expect("fixed tensor-target solve succeeds");
    let projection_report = kp20_planar_tau_rotational_section42_source_projection_report(
        valence,
        sector,
        &tensor_border,
    )
    .expect("source projection report");

    for (label, solved) in [("projected", &projected), ("fixed", &fixed)] {
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            sector,
            &tensor_border,
            &local_planar_patch_to_global(valence, sector, local_tau),
        )
        .expect("family residuals");
        let max_family_residual = residuals
            .iter()
            .map(|residual| residual.max_residual)
            .fold(0.0_f64, f64::max);
        eprintln!(
            "endpoint-aware tensor-target policy {label}: residual={residual:e} kkt={kkt:e} max_family={max_family:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2}) families={residuals:?}",
            residual = solved.max_constraint_residual,
            kkt = solved.kkt_condition_number.unwrap_or(f64::NAN),
            max_family = max_family_residual,
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
    }
    eprintln!(
        "endpoint-aware source projection report: adjusted_targets={} max_row_residual={:e} max_source_delta={:e}",
        projection_report.adjusted_scalar_targets,
        projection_report.max_row_residual,
        projection_report.max_source_delta,
    );

    let projected_local = global_planar_patch_to_local(valence, sector, projected.tau_macro_patch);
    let fixed_local = global_planar_patch_to_local(valence, sector, fixed.tau_macro_patch);
    let projected_validity =
        kp20_planar_bi6_tau_validity_report(&projected_local, &SAMPLES, 1.0e-10, 1.0e-10);
    let fixed_validity =
        kp20_planar_bi6_tau_validity_report(&fixed_local, &SAMPLES, 1.0e-10, 1.0e-10);

    assert!(
        projected.max_constraint_residual.is_finite() && fixed.max_constraint_residual.is_finite(),
        "both tensor-target policies must produce finite solves"
    );
    assert_eq!(
        projected_validity.orientation.near_zero_determinants, 100,
        "baseline projected tensor-target path documents the current collapsed live solve"
    );
    assert!(
        fixed_validity.orientation.max_determinant > projected_validity.orientation.max_determinant,
        "fixed tensor targets should reveal whether source projection is suppressing scale"
    );
}

#[test]
fn rotational_endpoint_aware_tensor_projection_groups_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector)
        .expect("tensor-border rows build");
    let all_free = tensor_rows
        .iter()
        .map(|row| row.terms[0].0)
        .collect::<Vec<_>>();

    let mut masks: Vec<(&'static str, Vec<Kp20PlanarTauCoeff>)> =
        vec![("none", Vec::new()), ("all", all_free.clone())];
    for piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let label = match piece {
            Kp20PieceId::NearMidpointLeft => "piece-near-midpoint-left",
            Kp20PieceId::NearMidpointRight => "piece-near-midpoint-right",
            _ => unreachable!("eq14 tensor-border rows only pin near-midpoint pieces"),
        };
        masks.push((
            label,
            all_free
                .iter()
                .copied()
                .filter(|coeff| coeff.piece == piece)
                .collect(),
        ));
    }
    for depth in 0..=2 {
        let label = match depth {
            0 => "depth-0",
            1 => "depth-1",
            2 => "depth-2",
            _ => unreachable!(),
        };
        masks.push((
            label,
            all_free
                .iter()
                .copied()
                .filter(|coeff| coeff.j == depth)
                .collect(),
        ));
    }
    for piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        for depth in 0..=2 {
            let label = match (piece, depth) {
                (Kp20PieceId::NearMidpointLeft, 0) => "left-depth-0",
                (Kp20PieceId::NearMidpointLeft, 1) => "left-depth-1",
                (Kp20PieceId::NearMidpointLeft, 2) => "left-depth-2",
                (Kp20PieceId::NearMidpointRight, 0) => "right-depth-0",
                (Kp20PieceId::NearMidpointRight, 1) => "right-depth-1",
                (Kp20PieceId::NearMidpointRight, 2) => "right-depth-2",
                _ => unreachable!("eq14 tensor-border rows only pin near-midpoint pieces"),
            };
            masks.push((
                label,
                all_free
                    .iter()
                    .copied()
                    .filter(|coeff| coeff.piece == piece && coeff.j == depth)
                    .collect(),
            ));
        }
    }

    let mut best_valid = ("none", f64::NEG_INFINITY, 0.0_f64, 0_usize);
    let mut best_row_closure = ("none", f64::INFINITY, usize::MAX, usize::MAX);

    for (label, free_coeffs) in masks {
        let row_set =
            kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
                valence,
                sector,
                &tensor_border,
                &free_coeffs,
            )
            .expect("tensor projection mask row set builds");
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints: row_set.vector_constraints,
                scalar_constraints: row_set.scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .expect("masked tensor-target solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            sector,
            &tensor_border,
            &local_planar_patch_to_global(valence, sector, local_tau),
        )
        .expect("family residuals");
        let max_family_residual = residuals
            .iter()
            .map(|residual| residual.max_residual)
            .fold(0.0_f64, f64::max);
        eprintln!(
            "endpoint-aware tensor projection mask {label}: free={} residual={residual:e} max_family={max_family:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] families={residuals:?}",
            free_coeffs.len(),
            residual = solved.max_constraint_residual,
            max_family = max_family_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if validity.domain.outside_triangle_samples == 0
            && validity.orientation.negative_determinants == 0
            && validity.orientation.near_zero_determinants < 100
            && validity.orientation.min_determinant > best_valid.1
        {
            best_valid = (
                label,
                validity.orientation.min_determinant,
                max_family_residual,
                free_coeffs.len(),
            );
        }
        if max_family_residual < best_row_closure.1
            || (max_family_residual == best_row_closure.1
                && validity.orientation.near_zero_determinants < best_row_closure.2)
        {
            best_row_closure = (
                label,
                max_family_residual,
                validity.orientation.near_zero_determinants,
                free_coeffs.len(),
            );
        }
    }

    eprintln!(
        "best valid noncollapsed tensor projection mask: label={} min_det={:e} max_family={:e} free={}",
        best_valid.0, best_valid.1, best_valid.2, best_valid.3
    );
    eprintln!(
        "best row-closing tensor projection mask: label={} max_family={:e} near_zero={} free={}",
        best_row_closure.0, best_row_closure.1, best_row_closure.2, best_row_closure.3
    );
    assert!(
        best_row_closure.1.is_finite(),
        "masked tensor-target sweep must produce finite row residuals"
    );
}

#[test]
fn rotational_endpoint_aware_centered_source_scale_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let mut projected_scalar_constraints = scalar_constraints.clone();
    projected_scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
            .expect("sector-local diagonal rows"),
    );
    let fixed_vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector)
        .expect("tensor-border rows");
    let axes = sector_local_axis_components(valence, sector);

    let scale_cases = [
        ("free-scale", None),
        ("lambda4-scale", Some(0.0914985_f64)),
        ("unit-scale", Some(1.0_f64)),
    ];
    for (label, fixed_scale) in scale_cases {
        let mut centered_rows =
            first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, axes, 0);
        if let Some(scale) = fixed_scale {
            centered_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: scale,
            });
        }
        let projected_tensor_constraints = project_tensor_band_with_centered_first_order(
            tensor_constraints.clone(),
            &fixed_vector_constraints,
            &projected_scalar_constraints,
            &centered_rows,
        );
        let mut vector_constraints = fixed_vector_constraints.clone();
        vector_constraints.extend(projected_tensor_constraints);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: projected_scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("endpoint-aware centered source solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            sector,
            &tensor_border,
            &local_planar_patch_to_global(valence, sector, local_tau),
        )
        .expect("family residuals");
        eprintln!(
            "endpoint-aware centered source scale {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] left_origin=({origin_x:e},{origin_y:e}) left_det={left_det:e} families={residuals:?}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            origin_x = left.origin.x,
            origin_y = left.origin.y,
            left_det = left.determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "centered source scale case {label} must produce a finite solve"
        );
    }
}

#[test]
fn rotational_endpoint_aware_tensor_border_index_frame_is_checked() {
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
    let standard_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector)
        .expect("standard tensor rows");
    let transposed_rows = transposed_tensor_border_rows(&tensor_border, sector);

    for (label, tensor_rows, project_targets) in [
        ("standard-fixed", standard_rows.clone(), false),
        ("standard-projected", standard_rows, true),
        ("transposed-fixed", transposed_rows.clone(), false),
        ("transposed-projected", transposed_rows, true),
    ] {
        let tensor_constraints = if project_targets {
            let free_coeffs = tensor_rows
                .iter()
                .map(|row| row.terms[0].0)
                .collect::<Vec<_>>();
            kp20_project_section42_source_tensor_rows_with_free_coeffs(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &free_coeffs,
            )
            .expect("tensor row projection succeeds")
            .0
        } else {
            tensor_rows
        };
        let mut vector_constraints = fig8_rows.clone();
        vector_constraints.extend(tensor_constraints.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("index-frame solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let tensor_residual = tensor_constraints
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "endpoint-aware tensor-border index frame {label}: residual={residual:e} tensor_residual={tensor_residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            "index-frame case {label} must produce a finite solve"
        );
    }
}

#[test]
fn rotational_endpoint_aware_tensor_border_chart_variants_are_ranked() {
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

    let mut best_residual = (
        TestTensorBorderChartVariant {
            previous_uses_near_m: false,
            current_uses_near_m: false,
            previous_reverse_u: false,
            current_reverse_u: false,
            transpose_indices: false,
            swap_pieces: false,
        },
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
    );
    let mut best_validity = best_residual;

    for previous_uses_near_m in [false, true] {
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let mut vector_constraints = fig8_rows.clone();
                            let tensor_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            vector_constraints.extend(tensor_rows.clone());
                            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                &Kp20PlanarTauF3ScalarSolveInput {
                                    initial: planar_tau_seed_macro_patch(valence, sector),
                                    vector_constraints,
                                    scalar_constraints: scalar_constraints.clone(),
                                    regularization_weight: 1.0e-10,
                                },
                            )
                            .expect("chart variant solve succeeds");
                            let local_tau = global_planar_patch_to_local(
                                valence,
                                sector,
                                solved.tau_macro_patch,
                            );
                            let validity = kp20_planar_bi6_tau_validity_report(
                                &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                            );
                            if solved.max_constraint_residual < best_residual.1
                                || (solved.max_constraint_residual == best_residual.1
                                    && validity.orientation.negative_determinants < best_residual.2)
                            {
                                best_residual = (
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                );
                            }
                            if validity.orientation.negative_determinants < best_validity.2
                                || (validity.orientation.negative_determinants == best_validity.2
                                    && validity.domain.outside_triangle_samples < best_validity.3)
                                || (validity.orientation.negative_determinants == best_validity.2
                                    && validity.domain.outside_triangle_samples == best_validity.3
                                    && solved.max_constraint_residual < best_validity.1)
                            {
                                best_validity = (
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                );
                            }
                            let should_report = solved.max_constraint_residual < 0.75
                                || (validity.orientation.negative_determinants <= 20
                                    && validity.domain.outside_triangle_samples <= 50);
                            if should_report {
                                eprintln!(
                                    "endpoint-aware tensor-border chart variant {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                                    label = variant.label(),
                                    residual = solved.max_constraint_residual,
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
                }
            }
        }
    }

    eprintln!(
        "best tensor-border chart residual: label={} residual={:e} neg={} outside={} min_det={:e}",
        best_residual.0.label(),
        best_residual.1,
        best_residual.2,
        best_residual.3,
        best_residual.4,
    );
    eprintln!(
        "best tensor-border chart validity: label={} residual={:e} neg={} outside={} min_det={:e}",
        best_validity.0.label(),
        best_validity.1,
        best_validity.2,
        best_validity.3,
        best_validity.4,
    );
    assert!(
        best_residual.1.is_finite() && best_validity.1.is_finite(),
        "chart variant ranking must produce finite solves"
    );
}

#[test]
fn rotational_endpoint_aware_best_tensor_chart_family_conflict_is_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut section42_c0_rows = Vec::new();
    add_kp20_section42_oriented_c0_rows(&mut section42_c0_rows, map);
    let mut section42_eq11_rows = Vec::new();
    add_kp20_section42_oriented_eq11_degree5_rows(&mut section42_eq11_rows, map);
    let mut section42_eq12_rows = Vec::new();
    add_kp20_section42_oriented_endpoint_aware_eq12_rows(&mut section42_eq12_rows, map);
    let mut section42_eq13_rows = Vec::new();
    add_kp20_section42_oriented_eq13_rows(&mut section42_eq13_rows, map, 2.0 * c);
    let tensor_rows = tensor_border_chart_variant_rows(
        &tensor_border,
        sector,
        TestTensorBorderChartVariant {
            previous_uses_near_m: true,
            current_uses_near_m: true,
            previous_reverse_u: true,
            current_reverse_u: false,
            transpose_indices: false,
            swap_pieces: false,
        },
    );

    for (label, include_fig8, include_tensor, include_diagonal) in [
        ("section42-only", false, false, false),
        ("section42-fig8", true, false, false),
        ("section42-diagonal", false, false, true),
        ("section42-tensor", false, true, false),
        ("section42-fig8-tensor", true, true, false),
        ("section42-tensor-diagonal", false, true, true),
        ("section42-fig8-diagonal", true, false, true),
        ("section42-fig8-tensor-diagonal", true, true, true),
    ] {
        let mut vector_constraints = Vec::new();
        if include_fig8 {
            vector_constraints.extend(fig8_rows.clone());
        }
        if include_tensor {
            vector_constraints.extend(tensor_rows.clone());
        }
        let mut scalar_constraints = section42_rows.clone();
        if include_diagonal {
            scalar_constraints.extend(diagonal_rows.clone());
        }
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .expect("best-chart family solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let fig8_residual = fig8_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let tensor_residual = tensor_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let section42_residual = section42_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        let section42_c0_residual = section42_c0_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        let section42_eq11_residual = section42_eq11_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        let section42_eq12_residual = section42_eq12_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        let section42_eq13_residual = section42_eq13_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        let mut eq13_row_reports = section42_eq13_rows
            .iter()
            .enumerate()
            .map(|(row_index, row)| {
                (
                    row_index,
                    scalar_residual(&solved.tau_macro_patch, row).abs(),
                )
            })
            .collect::<Vec<_>>();
        eq13_row_reports.sort_by(|a, b| b.1.total_cmp(&a.1));
        let diagonal_residual = diagonal_rows
            .iter()
            .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
            .fold(0.0_f64, f64::max);
        eprintln!(
            "endpoint-aware best tensor chart family {label}: residual={residual:e} fig8={fig8:e} tensor={tensor:e} section42={section42:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} diagonal={diagonal:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            fig8 = fig8_residual,
            tensor = tensor_residual,
            section42 = section42_residual,
            c0 = section42_c0_residual,
            eq11 = section42_eq11_residual,
            eq12 = section42_eq12_residual,
            eq13 = section42_eq13_residual,
            diagonal = diagonal_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if include_tensor && section42_eq13_residual > 1.0e-8 {
            for (row_index, row_residual) in eq13_row_reports.iter().take(4) {
                eprintln!(
                    "  best-chart solved eq13 row {row_index} semantic=\"{semantic}\" residual={row_residual:e}",
                    semantic = eq13_scalar_row_label(*row_index),
                );
            }
        }
        assert!(
            solved.max_constraint_residual.is_finite(),
            "best-chart family case {label} must produce a finite solve"
        );
    }
}

#[test]
fn rotational_endpoint_aware_tensor_targets_eq13_consistency_is_ranked() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut eq13_rows = Vec::new();
    add_kp20_section42_oriented_eq13_rows(&mut eq13_rows, map, 2.0 * c);

    let mut best = (
        TestTensorBorderChartVariant {
            previous_uses_near_m: false,
            current_uses_near_m: false,
            previous_reverse_u: false,
            current_reverse_u: false,
            transpose_indices: false,
            swap_pieces: false,
        },
        f64::INFINITY,
        usize::MAX,
    );

    for previous_uses_near_m in [false, true] {
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let mut target_patch = Kp20MacroPatch::zero();
                            let mut pinned_coeffs = Vec::new();
                            for row in
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant)
                            {
                                let coeff = row.terms[0].0;
                                pinned_coeffs.push(coeff);
                                target_patch.piece_mut(coeff.piece).cps[coeff.i][coeff.j] =
                                    Point3::new(row.target[0], row.target[1], 0.0);
                            }
                            let mut worst_row = 0_usize;
                            let mut residual = 0.0_f64;
                            for (row_index, row) in eq13_rows.iter().enumerate() {
                                let row_residual = scalar_residual(&target_patch, row).abs();
                                if row_residual > residual {
                                    residual = row_residual;
                                    worst_row = row_index;
                                }
                            }
                            if residual < best.1 {
                                best = (variant, residual, worst_row);
                            }
                            if residual < 0.1 {
                                let (covered, total, missing_weight) =
                                    scalar_row_coverage(&eq13_rows[worst_row], &pinned_coeffs);
                                eprintln!(
                                    "endpoint-aware tensor target eq13 direct {label}: residual={residual:e} worst_row={worst_row} semantic=\"{semantic}\" coverage={covered}/{total} missing_weight={missing_weight:e}",
                                    label = variant.label(),
                                    semantic = eq13_scalar_row_label(worst_row),
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "best tensor target eq13 direct consistency: label={} residual={:e} worst_row={}",
        best.0.label(),
        best.1,
        best.2,
    );
    let mut best_target_patch = Kp20MacroPatch::zero();
    let mut best_pinned_coeffs = Vec::new();
    for row in tensor_border_chart_variant_rows(&tensor_border, sector, best.0) {
        let coeff = row.terms[0].0;
        best_pinned_coeffs.push(coeff);
        best_target_patch.piece_mut(coeff.piece).cps[coeff.i][coeff.j] =
            Point3::new(row.target[0], row.target[1], 0.0);
    }
    let mut row_reports = eq13_rows
        .iter()
        .enumerate()
        .map(|(row_index, row)| {
            let residual = scalar_residual(&best_target_patch, row).abs();
            let (covered, total, missing_weight) = scalar_row_coverage(row, &best_pinned_coeffs);
            (row_index, residual, covered, total, missing_weight)
        })
        .collect::<Vec<_>>();
    row_reports.sort_by(|a, b| b.1.total_cmp(&a.1));
    for (row_index, residual, covered, total, missing_weight) in row_reports.iter().take(8) {
        eprintln!(
            "  best direct eq13 row {row_index} semantic=\"{semantic}\" residual={residual:e} coverage={covered}/{total} missing_weight={missing_weight:e}",
            semantic = eq13_scalar_row_label(*row_index),
        );
    }
    assert!(
        best.1.is_finite(),
        "direct tensor-target eq13 consistency ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_eq14_second_order_factor_is_ranked() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut eq13_rows = Vec::new();
    add_kp20_section42_oriented_eq13_rows(&mut eq13_rows, map, 2.0 * c);
    let factors = [-0.2, -0.1, 0.0, 0.05, 0.1, 0.2, 0.25, 0.5];

    let mut best = (
        0.0_f64,
        TestTensorBorderChartVariant {
            previous_uses_near_m: false,
            current_uses_near_m: false,
            previous_reverse_u: false,
            current_reverse_u: false,
            transpose_indices: false,
            swap_pieces: false,
        },
        f64::INFINITY,
        usize::MAX,
    );

    for second_order_factor in factors {
        for previous_uses_near_m in [false, true] {
            for current_uses_near_m in [false, true] {
                let (previous_rows, current_rows) =
                    split_regular_sector_rows_with_second_order_factor(
                        &tensor_border,
                        sector,
                        previous_uses_near_m,
                        current_uses_near_m,
                        second_order_factor,
                    );
                for previous_reverse_u in [false, true] {
                    for current_reverse_u in [false, true] {
                        for transpose_indices in [false, true] {
                            for swap_pieces in [false, true] {
                                let variant = TestTensorBorderChartVariant {
                                    previous_uses_near_m,
                                    current_uses_near_m,
                                    previous_reverse_u,
                                    current_reverse_u,
                                    transpose_indices,
                                    swap_pieces,
                                };
                                let mut target_patch = Kp20MacroPatch::zero();
                                for row in tensor_border_chart_variant_rows_from_split_rows(
                                    previous_rows,
                                    current_rows,
                                    variant,
                                ) {
                                    let coeff = row.terms[0].0;
                                    target_patch.piece_mut(coeff.piece).cps[coeff.i][coeff.j] =
                                        Point3::new(row.target[0], row.target[1], 0.0);
                                }
                                let mut worst_row = 0_usize;
                                let mut residual = 0.0_f64;
                                for (row_index, row) in eq13_rows.iter().enumerate() {
                                    let row_residual = scalar_residual(&target_patch, row).abs();
                                    if row_residual > residual {
                                        residual = row_residual;
                                        worst_row = row_index;
                                    }
                                }
                                if residual < best.2 {
                                    best = (second_order_factor, variant, residual, worst_row);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "best eq14 second-order factor for direct eq13 consistency: factor={:e} label={} residual={:e} worst_row={}",
        best.0,
        best.1.label(),
        best.2,
        best.3,
    );
    assert!(
        best.2.is_finite(),
        "second-order factor sweep must produce finite residuals"
    );
}

#[test]
fn rotational_endpoint_aware_eq13_source_row_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut section42_eq13_rows = Vec::new();
    add_kp20_section42_oriented_eq13_rows(&mut section42_eq13_rows, map, 2.0 * c);

    let variants = [
        (
            "best-fixed-row-closing",
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: true,
                previous_reverse_u: true,
                current_reverse_u: false,
                transpose_indices: false,
                swap_pieces: false,
            },
        ),
        (
            "best-direct-eq13",
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: false,
                transpose_indices: true,
                swap_pieces: false,
            },
        ),
    ];

    for (case_label, variant) in variants {
        let tensor_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
        for include_diagonal in [false, true] {
            let projected_tensor_rows = without_tensor_chart_depth(&tensor_rows, variant, 1);
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(projected_tensor_rows.clone());
            let mut scalar_constraints = section42_rows.clone();
            if include_diagonal {
                scalar_constraints.extend(diagonal_rows.clone());
            }
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("source-row projected tensor solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let pinned_tensor_residual = projected_tensor_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let full_tensor_residual = tensor_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let eq13_residual = section42_eq13_rows
                .iter()
                .map(|row| scalar_residual(&solved.tau_macro_patch, row).abs())
                .fold(0.0_f64, f64::max);
            eprintln!(
                "endpoint-aware eq13 source-row projection {case_label} chart={} diagonal={include_diagonal}: residual={residual:e} pinned_tensor={pinned_tensor:e} full_tensor={full_tensor:e} eq13={eq13:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                variant.label(),
                residual = solved.max_constraint_residual,
                pinned_tensor = pinned_tensor_residual,
                full_tensor = full_tensor_residual,
                eq13 = eq13_residual,
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
                "source-row projected case {case_label} must produce finite residual"
            );
        }
    }
}

#[test]
fn rotational_endpoint_aware_eq13_source_row_projection_chart_space_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");

    let mut best_row_closing = (
        TestTensorBorderChartVariant {
            previous_uses_near_m: false,
            current_uses_near_m: false,
            previous_reverse_u: false,
            current_reverse_u: false,
            transpose_indices: false,
            swap_pieces: false,
        },
        false,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_validity = best_row_closing;

    for previous_uses_near_m in [false, true] {
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let tensor_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            let projected_tensor_rows =
                                without_tensor_chart_depth(&tensor_rows, variant, 1);
                            for include_diagonal in [false, true] {
                                let mut vector_constraints = fig8_rows.clone();
                                vector_constraints.extend(projected_tensor_rows.clone());
                                let mut scalar_constraints = section42_rows.clone();
                                if include_diagonal {
                                    scalar_constraints.extend(diagonal_rows.clone());
                                }
                                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                    &Kp20PlanarTauF3ScalarSolveInput {
                                        initial: planar_tau_seed_macro_patch(valence, sector),
                                        vector_constraints,
                                        scalar_constraints,
                                        regularization_weight: 1.0e-10,
                                    },
                                )
                                .expect("source-row projected chart solve succeeds");
                                let local_tau = global_planar_patch_to_local(
                                    valence,
                                    sector,
                                    solved.tau_macro_patch,
                                );
                                let validity = kp20_planar_bi6_tau_validity_report(
                                    &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                                );
                                let full_tensor_residual = tensor_rows
                                    .iter()
                                    .map(|row| {
                                        let residual =
                                            vector_residual(&solved.tau_macro_patch, row);
                                        residual[0].abs().max(residual[1].abs())
                                    })
                                    .fold(0.0_f64, f64::max);
                                let candidate = (
                                    variant,
                                    include_diagonal,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                    full_tensor_residual,
                                );
                                if candidate.2 < best_row_closing.2
                                    || (candidate.2 == best_row_closing.2
                                        && candidate.3 < best_row_closing.3)
                                {
                                    best_row_closing = candidate;
                                }
                                if candidate.3 < best_validity.3
                                    || (candidate.3 == best_validity.3
                                        && candidate.4 < best_validity.4)
                                    || (candidate.3 == best_validity.3
                                        && candidate.4 == best_validity.4
                                        && candidate.2 < best_validity.2)
                                {
                                    best_validity = candidate;
                                }
                                if validity.orientation.negative_determinants <= 20
                                    || solved.max_constraint_residual < 1.0e-8
                                {
                                    eprintln!(
                                        "endpoint-aware eq13 source-row projected chart {} diagonal={include_diagonal}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                                        variant.label(),
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
                    }
                }
            }
        }
    }

    eprintln!(
        "best eq13 source-row projected row closure: label={} diagonal={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_row_closing.0.label(),
        best_row_closing.1,
        best_row_closing.2,
        best_row_closing.3,
        best_row_closing.4,
        best_row_closing.5,
        best_row_closing.6,
    );
    eprintln!(
        "best eq13 source-row projected validity: label={} diagonal={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.label(),
        best_validity.1,
        best_validity.2,
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
    );
    assert!(
        best_row_closing.2.is_finite() && best_validity.2.is_finite(),
        "source-row projected chart ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_eq13_source_row_projection_active_set_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: true,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let tensor_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
    let projected_tensor_rows = without_tensor_chart_depth(&tensor_rows, variant, 1);

    let mut vector_constraints = fig8_rows.clone();
    vector_constraints.extend(projected_tensor_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: section42_rows.clone(),
            regularization_weight: 1.0e-10,
        })
        .expect("source-row projected F3 solve succeeds");

    let mut scalar_rows = section42_rows.clone();
    scalar_rows.extend(vector_constraints_as_scalar_rows(&projected_tensor_rows));
    for (label, normalize_rows, max_step_delta) in [
        ("raw", false, None),
        ("row-normalized", true, None),
        ("row-normalized-step-limited", true, Some(0.25)),
    ] {
        let active_set = square_domain_active_set_projection_with_base_and_extra_rows_solver(
            &solved.tau_macro_patch,
            &scalar_rows,
            valence,
            sector,
            &SAMPLES,
            &[],
            max_step_delta,
            normalize_rows,
        );
        let local_tau = global_planar_patch_to_local(valence, sector, active_set.patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let row_residual = scalar_row_set_residual(&active_set.patch, &scalar_rows);
        let full_tensor_residual = tensor_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&active_set.patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "endpoint-aware eq13 source-row active-set {label}: chart={} iterations={} active={} orientation_active={} max_delta={:e} linear_residual={:e} row_residual={row_residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            variant.label(),
            active_set.iterations,
            active_set.active_constraints,
            active_set.orientation_constraints,
            active_set.max_delta,
            active_set.linear_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            row_residual.is_finite(),
            "active-set source-row projection must produce finite residual"
        );
    }
}

#[test]
fn rotational_endpoint_aware_q1_eq13_source_row_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border =
        kp20_pr08_stam_q1_characteristic_tensor_border(valence).expect("Q1 tensor-border");
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");
    let variants = [
        (
            "regular-best-row-closing",
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: true,
                current_reverse_u: false,
                transpose_indices: false,
                swap_pieces: false,
            },
        ),
        (
            "regular-best-validity",
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: true,
                transpose_indices: true,
                swap_pieces: true,
            },
        ),
    ];

    for (case_label, variant) in variants {
        let tensor_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
        let projected_tensor_rows = without_tensor_chart_depth(&tensor_rows, variant, 1);
        for include_diagonal in [false, true] {
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(projected_tensor_rows.clone());
            let mut scalar_constraints = section42_rows.clone();
            if include_diagonal {
                scalar_constraints.extend(diagonal_rows.clone());
            }
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("Q1 source-row projected solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = tensor_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "endpoint-aware Q1 eq13 source-row projection {case_label} chart={} diagonal={include_diagonal}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                variant.label(),
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
            assert!(
                solved.max_constraint_residual.is_finite(),
                "Q1 source-row projected case must produce finite residual"
            );
        }
    }
}

#[test]
fn rotational_endpoint_aware_pr08_candidates_eq13_source_frame_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();

    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best_row_closing = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_validity = best_row_closing;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let cnet_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&cnet_rows).expect("candidate tensor rows");
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m: true,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let tensor_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            let projected_tensor_rows =
                                without_tensor_chart_depth(&tensor_rows, variant, 1);
                            let mut vector_constraints = fig8_rows.clone();
                            vector_constraints.extend(projected_tensor_rows);
                            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                &Kp20PlanarTauF3ScalarSolveInput {
                                    initial: planar_tau_seed_macro_patch(valence, sector),
                                    vector_constraints,
                                    scalar_constraints: section42_rows.clone(),
                                    regularization_weight: 1.0e-10,
                                },
                            )
                            .expect("PR08 candidate source-frame solve succeeds");
                            let local_tau = global_planar_patch_to_local(
                                valence,
                                sector,
                                solved.tau_macro_patch,
                            );
                            let validity = kp20_planar_bi6_tau_validity_report(
                                &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                            );
                            let full_tensor_residual = tensor_rows
                                .iter()
                                .map(|row| {
                                    let residual = vector_residual(&solved.tau_macro_patch, row);
                                    residual[0].abs().max(residual[1].abs())
                                })
                                .fold(0.0_f64, f64::max);
                            let row = (
                                candidate,
                                variant,
                                solved.max_constraint_residual,
                                validity.orientation.negative_determinants,
                                validity.domain.outside_triangle_samples,
                                validity.orientation.min_determinant,
                                full_tensor_residual,
                            );
                            if row.2 < best_row_closing.2
                                || (row.2 == best_row_closing.2 && row.3 < best_row_closing.3)
                            {
                                best_row_closing = row;
                            }
                            if row.3 < best_validity.3
                                || (row.3 == best_validity.3 && row.4 < best_validity.4)
                                || (row.3 == best_validity.3
                                    && row.4 == best_validity.4
                                    && row.2 < best_validity.2)
                            {
                                best_validity = row;
                            }
                            if validity.orientation.negative_determinants <= 20
                                || (solved.max_constraint_residual < 1.0e-8
                                    && full_tensor_residual < 0.1)
                            {
                                eprintln!(
                                    "endpoint-aware PR08 candidate eq13 source frame candidate={} chart={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                                    candidate.name(),
                                    variant.label(),
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
                }
            }
        }
    }

    eprintln!(
        "best PR08 eq13 source-frame row closure: candidate={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_row_closing.0.name(),
        best_row_closing.1.label(),
        best_row_closing.2,
        best_row_closing.3,
        best_row_closing.4,
        best_row_closing.5,
        best_row_closing.6,
    );
    eprintln!(
        "best PR08 eq13 source-frame validity: candidate={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.name(),
        best_validity.1.label(),
        best_validity.2,
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
    );
    assert!(
        best_row_closing.2.is_finite() && best_validity.2.is_finite(),
        "PR08 source-frame ranking must produce finite candidates"
    );
}

#[test]
fn rotational_endpoint_aware_endpoint_normalized_pr08_eq13_source_frame_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();

    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best_row_closing = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_validity = best_row_closing;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let Some(cnet_rows) = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        ) else {
            continue;
        };
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&cnet_rows).expect("candidate tensor rows");
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m: true,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let tensor_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            let projected_tensor_rows =
                                without_tensor_chart_depth(&tensor_rows, variant, 1);
                            let mut vector_constraints = fig8_rows.clone();
                            vector_constraints.extend(projected_tensor_rows);
                            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                &Kp20PlanarTauF3ScalarSolveInput {
                                    initial: planar_tau_seed_macro_patch(valence, sector),
                                    vector_constraints,
                                    scalar_constraints: section42_rows.clone(),
                                    regularization_weight: 1.0e-10,
                                },
                            )
                            .expect("endpoint-normalized PR08 source-frame solve succeeds");
                            let local_tau = global_planar_patch_to_local(
                                valence,
                                sector,
                                solved.tau_macro_patch,
                            );
                            let validity = kp20_planar_bi6_tau_validity_report(
                                &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                            );
                            let full_tensor_residual = tensor_rows
                                .iter()
                                .map(|row| {
                                    let residual = vector_residual(&solved.tau_macro_patch, row);
                                    residual[0].abs().max(residual[1].abs())
                                })
                                .fold(0.0_f64, f64::max);
                            let row = (
                                candidate,
                                variant,
                                solved.max_constraint_residual,
                                validity.orientation.negative_determinants,
                                validity.domain.outside_triangle_samples,
                                validity.orientation.min_determinant,
                                full_tensor_residual,
                            );
                            if row.2 < best_row_closing.2
                                || (row.2 == best_row_closing.2 && row.3 < best_row_closing.3)
                            {
                                best_row_closing = row;
                            }
                            if row.3 < best_validity.3
                                || (row.3 == best_validity.3 && row.4 < best_validity.4)
                                || (row.3 == best_validity.3
                                    && row.4 == best_validity.4
                                    && row.2 < best_validity.2)
                            {
                                best_validity = row;
                            }
                            if validity.orientation.negative_determinants <= 20
                                || (solved.max_constraint_residual < 1.0e-8
                                    && full_tensor_residual < 0.1)
                            {
                                eprintln!(
                                    "endpoint-aware normalized PR08 eq13 source frame candidate={} chart={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                                    candidate.name(),
                                    variant.label(),
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
                }
            }
        }
    }

    eprintln!(
        "best normalized PR08 eq13 source-frame row closure: candidate={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_row_closing.0.name(),
        best_row_closing.1.label(),
        best_row_closing.2,
        best_row_closing.3,
        best_row_closing.4,
        best_row_closing.5,
        best_row_closing.6,
    );
    eprintln!(
        "best normalized PR08 eq13 source-frame validity: candidate={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.name(),
        best_validity.1.label(),
        best_validity.2,
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
    );
    assert!(
        best_row_closing.2.is_finite() && best_validity.2.is_finite(),
        "endpoint-normalized PR08 source-frame ranking must produce finite candidates"
    );
}

#[test]
fn rotational_endpoint_aware_normalized_pr08_depth_order_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const DEPTH_ORDERS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let candidate_set = [
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    ];

    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        [0_usize, 1, 2],
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );

    for candidate in candidate_set {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        for order in DEPTH_ORDERS {
            let ordered_rows =
                reorder_regular_tensor_border_depth_rows_for_test(&normalized_rows, order);
            let tensor_border =
                kp20_tensor_border_from_regular_cnet_rows(&ordered_rows).expect("tensor rows");
            for current_uses_near_m in [false, true] {
                for previous_reverse_u in [false, true] {
                    for current_reverse_u in [false, true] {
                        for transpose_indices in [false, true] {
                            for swap_pieces in [false, true] {
                                let variant = TestTensorBorderChartVariant {
                                    previous_uses_near_m: true,
                                    current_uses_near_m,
                                    previous_reverse_u,
                                    current_reverse_u,
                                    transpose_indices,
                                    swap_pieces,
                                };
                                let tensor_rows = tensor_border_chart_variant_rows(
                                    &tensor_border,
                                    sector,
                                    variant,
                                );
                                let projected_tensor_rows =
                                    without_tensor_chart_depth(&tensor_rows, variant, 1);
                                let mut vector_constraints = fig8_rows.clone();
                                vector_constraints.extend(projected_tensor_rows);
                                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                    &Kp20PlanarTauF3ScalarSolveInput {
                                        initial: planar_tau_seed_macro_patch(valence, sector),
                                        vector_constraints,
                                        scalar_constraints: section42_rows.clone(),
                                        regularization_weight: 1.0e-10,
                                    },
                                )
                                .expect("normalized PR08 depth-order solve succeeds");
                                let local_tau = global_planar_patch_to_local(
                                    valence,
                                    sector,
                                    solved.tau_macro_patch,
                                );
                                let validity = kp20_planar_bi6_tau_validity_report(
                                    &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                                );
                                let full_tensor_residual = tensor_rows
                                    .iter()
                                    .map(|row| {
                                        let residual =
                                            vector_residual(&solved.tau_macro_patch, row);
                                        residual[0].abs().max(residual[1].abs())
                                    })
                                    .fold(0.0_f64, f64::max);
                                let row = (
                                    candidate,
                                    order,
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                    full_tensor_residual,
                                );
                                if row.4 < best.4
                                    || (row.4 == best.4 && row.5 < best.5)
                                    || (row.4 == best.4 && row.5 == best.5 && row.3 < best.3)
                                {
                                    best = row;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!(
        "best normalized PR08 depth order: candidate={} order={:?} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best.0.name(),
        best.1,
        best.2.label(),
        best.3,
        best.4,
        best.5,
        best.6,
        best.7,
    );
    assert!(best.3.is_finite(), "depth-order ranking must be finite");
}

#[test]
fn rotational_endpoint_aware_best_normalized_pr08_source_frame_is_localized() {
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
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: true,
        previous_reverse_u: true,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let tensor_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
    let projected_tensor_rows = without_tensor_chart_depth(&tensor_rows, variant, 1);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("best normalized PR08 localization solve succeeds");
    let full_tensor_residual = tensor_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "best normalized PR08 localized source frame candidate={} chart={} residual={:e} full_tensor={:e}",
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows.name(),
        variant.label(),
        solved.max_constraint_residual,
        full_tensor_residual,
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
    for piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let patch = local_tau.piece(piece);
        let row0_left = patch.cps[0][0];
        let row0_right = patch.cps[0][6];
        let row2_left = patch.cps[2][0];
        let row2_right = patch.cps[2][6];
        let tangent = sub_point(row0_right, row0_left);
        let left_cross = sub_point(row2_left, row0_left);
        let right_cross = sub_point(row2_right, row0_right);
        let left_det = tangent.x * left_cross.y - tangent.y * left_cross.x;
        let right_det = tangent.x * right_cross.y - tangent.y * right_cross.x;
        eprintln!(
            "piece {:?} pinned depth frame: row0_left=({:.6},{:.6}) row0_right=({:.6},{:.6}) row2_left=({:.6},{:.6}) row2_right=({:.6},{:.6}) tangent=({:.6},{:.6}) left_cross=({:.6},{:.6}) right_cross=({:.6},{:.6}) endpoint_det=[{:.6},{:.6}]",
            piece,
            row0_left.x,
            row0_left.y,
            row0_right.x,
            row0_right.y,
            row2_left.x,
            row2_left.y,
            row2_right.x,
            row2_right.y,
            tangent.x,
            tangent.y,
            left_cross.x,
            left_cross.y,
            right_cross.x,
            right_cross.y,
            left_det,
            right_det,
        );
    }
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "localized source-frame solve must remain finite"
    );
}

#[test]
fn rotational_endpoint_aware_best_normalized_pr08_cross_depth_orientation_is_measured() {
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
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: true,
        previous_reverse_u: true,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let tensor_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
    let cases = [
        ("baseline", tensor_rows.clone()),
        (
            "reflect-current-depth2",
            reflect_tensor_chart_depth_across_boundary_for_piece(
                &tensor_rows,
                variant,
                Kp20PieceId::NearMidpointRight,
                2,
            ),
        ),
        (
            "reflect-both-depth2",
            reflect_tensor_chart_depth_across_boundary_for_piece(
                &reflect_tensor_chart_depth_across_boundary_for_piece(
                    &tensor_rows,
                    variant,
                    Kp20PieceId::NearMidpointRight,
                    2,
                ),
                variant,
                Kp20PieceId::NearMidpointLeft,
                2,
            ),
        ),
    ];
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, f64::INFINITY);
    for (label, source_rows) in cases {
        let projected_tensor_rows = without_tensor_chart_depth(&source_rows, variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: kp20_section42_oriented_endpoint_aware_strip_basis_rows(
                    2.0 * c,
                ),
                regularization_weight: 1.0e-10,
            })
            .expect("cross-depth orientation solve succeeds");
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
            "normalized PR08 cross-depth orientation {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        "best normalized PR08 cross-depth orientation: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "cross-depth orientation ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_beta_source_reparameterizations_are_measured() {
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
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: true,
        previous_reverse_u: true,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let kp20_split = split_kp20_tensor_border(&tensor_border).expect("KP20 split");
    let kp15_split_then_apply =
        kp15_beta_split_then_apply_for_test(&tensor_border).expect("KP15 split-then-apply");
    let kp15_full_then_split =
        kp15_beta_full_then_split_for_test(&tensor_border).expect("KP15 full-then-split");
    let cases = [
        ("kp20-eq14-current", kp20_split.as_slice()),
        (
            "kp15-beta-split-then-apply",
            kp15_split_then_apply.as_slice(),
        ),
        ("kp15-beta-full-then-split", kp15_full_then_split.as_slice()),
    ];
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, f64::INFINITY);
    for (label, split) in cases {
        let source_rows = tensor_border_split_chart_variant_rows(split, sector, variant);
        let projected_tensor_rows = without_tensor_chart_depth(&source_rows, variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: kp20_section42_oriented_endpoint_aware_strip_basis_rows(
                    2.0 * c,
                ),
                regularization_weight: 1.0e-10,
            })
            .expect("beta source reparameterization solve succeeds");
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
            "endpoint-aware beta source {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        "best endpoint-aware beta source: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "beta source reparameterization ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_previous_source_sector_rotation_is_measured() {
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
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: true,
        previous_reverse_u: true,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let baseline_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
    let cases = [
        ("unrotated", baseline_rows.clone()),
        (
            "previous-plus-alpha",
            rotate_tensor_chart_piece_targets(&baseline_rows, Kp20PieceId::NearMidpointLeft, angle),
        ),
        (
            "previous-minus-alpha",
            rotate_tensor_chart_piece_targets(
                &baseline_rows,
                Kp20PieceId::NearMidpointLeft,
                -angle,
            ),
        ),
    ];
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, f64::INFINITY);
    for (label, source_rows) in cases {
        let projected_tensor_rows = without_tensor_chart_depth(&source_rows, variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: kp20_section42_oriented_endpoint_aware_strip_basis_rows(
                    2.0 * c,
                ),
                regularization_weight: 1.0e-10,
            })
            .expect("previous-sector source rotation solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let nml = local_tau.piece(Kp20PieceId::NearMidpointLeft);
        eprintln!(
            "endpoint-aware previous-source rotation {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} square={square_outside}/{square_samples} square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] nml_row0=({x0:e},{y0:e})->({x1:e},{y1:e})",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            square_outside = square.outside,
            square_samples = square.samples,
            min_u = square.min_u,
            max_u = square.max_u,
            min_v = square.min_v,
            max_v = square.max_v,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            x0 = nml.cps[0][0].x,
            y0 = nml.cps[0][0].y,
            x1 = nml.cps[0][6].x,
            y1 = nml.cps[0][6].y,
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
        "best endpoint-aware previous-source rotation: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "previous-source rotation ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_sector_rotated_normalized_source_charts_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let section42_rows = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let candidate_set = [
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    ];
    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best_validity = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        "none",
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_domain = best_validity;

    for candidate in candidate_set {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m: true,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let previous_piece = if variant.swap_pieces {
                                Kp20PieceId::NearMidpointRight
                            } else {
                                Kp20PieceId::NearMidpointLeft
                            };
                            let base_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            for (rotation_label, source_rows) in [
                                ("unrotated", base_rows.clone()),
                                (
                                    "previous-minus-alpha",
                                    rotate_tensor_chart_piece_targets(
                                        &base_rows,
                                        previous_piece,
                                        -angle,
                                    ),
                                ),
                            ] {
                                let projected_tensor_rows =
                                    without_tensor_chart_depth(&source_rows, variant, 1);
                                let mut vector_constraints = fig8_rows.clone();
                                vector_constraints.extend(projected_tensor_rows);
                                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                    &Kp20PlanarTauF3ScalarSolveInput {
                                        initial: planar_tau_seed_macro_patch(valence, sector),
                                        vector_constraints,
                                        scalar_constraints: section42_rows.clone(),
                                        regularization_weight: 1.0e-10,
                                    },
                                )
                                .expect("sector-rotated normalized source chart solve succeeds");
                                let local_tau = global_planar_patch_to_local(
                                    valence,
                                    sector,
                                    solved.tau_macro_patch,
                                );
                                let validity = kp20_planar_bi6_tau_validity_report(
                                    &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                                );
                                let full_tensor_residual = source_rows
                                    .iter()
                                    .map(|row| {
                                        let residual =
                                            vector_residual(&solved.tau_macro_patch, row);
                                        residual[0].abs().max(residual[1].abs())
                                    })
                                    .fold(0.0_f64, f64::max);
                                let row = (
                                    candidate,
                                    rotation_label,
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                    full_tensor_residual,
                                );
                                if row.4 < best_validity.4
                                    || (row.4 == best_validity.4 && row.5 < best_validity.5)
                                    || (row.4 == best_validity.4
                                        && row.5 == best_validity.5
                                        && row.3 < best_validity.3)
                                {
                                    best_validity = row;
                                }
                                if row.5 < best_domain.5
                                    || (row.5 == best_domain.5 && row.4 < best_domain.4)
                                    || (row.5 == best_domain.5
                                        && row.4 == best_domain.4
                                        && row.3 < best_domain.3)
                                {
                                    best_domain = row;
                                }
                                if validity.domain.outside_triangle_samples <= 30
                                    || validity.orientation.negative_determinants <= 20
                                {
                                    eprintln!(
                                        "sector-rotated normalized source candidate={} rotation={} chart={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                                        candidate.name(),
                                        rotation_label,
                                        variant.label(),
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
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "best sector-rotated normalized source validity: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.name(),
        best_validity.1,
        best_validity.2.label(),
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
        best_validity.7,
    );
    eprintln!(
        "best sector-rotated normalized source domain: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_domain.0.name(),
        best_domain.1,
        best_domain.2.label(),
        best_domain.3,
        best_domain.4,
        best_domain.5,
        best_domain.6,
        best_domain.7,
    );
    assert!(
        best_validity.3.is_finite() && best_domain.3.is_finite(),
        "sector-rotated normalized source ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_best_domain_sector_rotated_source_is_localized() {
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
    let variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, variant, 1);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("best domain sector-rotated source solve succeeds");
    let full_tensor_residual = source_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "best domain sector-rotated source candidate={} chart={} residual={:e} full_tensor={:e}",
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows.name(),
        variant.label(),
        solved.max_constraint_residual,
        full_tensor_residual,
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
    for piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let patch = local_tau.piece(piece);
        let row0_left = patch.cps[0][0];
        let row0_right = patch.cps[0][6];
        let row2_left = patch.cps[2][0];
        let row2_right = patch.cps[2][6];
        let tangent = sub_point(row0_right, row0_left);
        let left_cross = sub_point(row2_left, row0_left);
        let right_cross = sub_point(row2_right, row0_right);
        let left_det = tangent.x * left_cross.y - tangent.y * left_cross.x;
        let right_det = tangent.x * right_cross.y - tangent.y * right_cross.x;
        eprintln!(
            "piece {:?} pinned depth frame: row0_left=({:.6},{:.6}) row0_right=({:.6},{:.6}) row2_left=({:.6},{:.6}) row2_right=({:.6},{:.6}) endpoint_det=[{:.6},{:.6}]",
            piece,
            row0_left.x,
            row0_left.y,
            row0_right.x,
            row0_right.y,
            row2_left.x,
            row2_left.y,
            row2_right.x,
            row2_right.y,
            left_det,
            right_det,
        );
    }
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "best domain sector-rotated source localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_near_origin_chart_handedness_is_ranked() {
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
    let base_source_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows = rotate_tensor_chart_piece_targets(
        &base_source_rows,
        Kp20PieceId::NearMidpointRight,
        -angle,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
    let underlined_p = TestSection42ChartStrip {
        piece: Kp20PieceId::NearMidpointLeft,
        reverse_i: false,
        reverse_j: false,
        transpose: false,
    };
    let arrow_underlined_p = TestSection42ChartStrip {
        piece: Kp20PieceId::NearMidpointRight,
        reverse_i: true,
        reverse_j: false,
        transpose: false,
    };
    let seed_strip = TestSection42ChartStrip {
        piece: Kp20PieceId::NearOriginLeft,
        reverse_i: false,
        reverse_j: false,
        transpose: false,
    };
    let mut best_validity = (
        seed_strip,
        seed_strip,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_domain = best_validity;

    for p_piece in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        let arrow_piece = match p_piece {
            Kp20PieceId::NearOriginLeft => Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginRight => Kp20PieceId::NearOriginLeft,
            _ => unreachable!("near-origin piece sweep"),
        };
        for transpose in [false, true] {
            for reverse_j in [false, true] {
                for p_reverse_i in [false, true] {
                    for arrow_reverse_i in [false, true] {
                        let p = TestSection42ChartStrip {
                            piece: p_piece,
                            reverse_i: p_reverse_i,
                            reverse_j,
                            transpose,
                        };
                        let arrow_p = TestSection42ChartStrip {
                            piece: arrow_piece,
                            reverse_i: arrow_reverse_i,
                            reverse_j,
                            transpose,
                        };
                        let map = TestSection42ChartMap {
                            p,
                            arrow_p,
                            underlined_p,
                            arrow_underlined_p,
                        };
                        let mut vector_constraints =
                            kp20_planar_tau_fig8_indicated_averaging_constraints();
                        vector_constraints.extend(projected_tensor_rows.clone());
                        let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                            &Kp20PlanarTauF3ScalarSolveInput {
                                initial: planar_tau_seed_macro_patch(valence, sector),
                                vector_constraints,
                                scalar_constraints: test_section42_endpoint_aware_chart_rows(
                                    map,
                                    2.0 * c,
                                ),
                                regularization_weight: 1.0e-10,
                            },
                        )
                        .expect("near-origin chart handedness solve succeeds");
                        let local_tau =
                            global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                        );
                        let full_tensor_residual = source_rows
                            .iter()
                            .map(|row| {
                                let residual = vector_residual(&solved.tau_macro_patch, row);
                                residual[0].abs().max(residual[1].abs())
                            })
                            .fold(0.0_f64, f64::max);
                        let row = (
                            p,
                            arrow_p,
                            solved.max_constraint_residual,
                            validity.orientation.negative_determinants,
                            validity.domain.outside_triangle_samples,
                            validity.orientation.min_determinant,
                            full_tensor_residual,
                        );
                        if row.3 < best_validity.3
                            || (row.3 == best_validity.3 && row.4 < best_validity.4)
                            || (row.3 == best_validity.3
                                && row.4 == best_validity.4
                                && row.2 < best_validity.2)
                        {
                            best_validity = row;
                        }
                        if row.4 < best_domain.4
                            || (row.4 == best_domain.4 && row.3 < best_domain.3)
                            || (row.4 == best_domain.4
                                && row.3 == best_domain.3
                                && row.2 < best_domain.2)
                        {
                            best_domain = row;
                        }
                        if validity.orientation.negative_determinants <= 20
                            || validity.domain.outside_triangle_samples <= 10
                        {
                            eprintln!(
                                "near-origin chart p={} arrow={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                                p.label(),
                                arrow_p.label(),
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
                    }
                }
            }
        }
    }
    eprintln!(
        "best near-origin chart validity: p={} arrow={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.label(),
        best_validity.1.label(),
        best_validity.2,
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
    );
    eprintln!(
        "best near-origin chart domain: p={} arrow={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_domain.0.label(),
        best_domain.1.label(),
        best_domain.2,
        best_domain.3,
        best_domain.4,
        best_domain.5,
        best_domain.6,
    );
    assert!(
        best_validity.2.is_finite() && best_domain.2.is_finite(),
        "near-origin chart handedness ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_best_near_origin_chart_is_localized() {
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
    let base_source_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows = rotate_tensor_chart_piece_targets(
        &base_source_rows,
        Kp20PieceId::NearMidpointRight,
        -angle,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
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
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("best near-origin chart localization solve succeeds");
    let full_tensor_residual = source_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "best near-origin chart localized residual={:e} full_tensor={:e}",
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
        "best near-origin chart localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_underlined_chart_handedness_after_near_origin_fix_is_ranked() {
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
    let base_source_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows = rotate_tensor_chart_piece_targets(
        &base_source_rows,
        Kp20PieceId::NearMidpointRight,
        -angle,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
    let p = TestSection42ChartStrip {
        piece: Kp20PieceId::NearOriginLeft,
        reverse_i: true,
        reverse_j: false,
        transpose: false,
    };
    let arrow_p = TestSection42ChartStrip {
        piece: Kp20PieceId::NearOriginRight,
        reverse_i: false,
        reverse_j: false,
        transpose: false,
    };
    let seed_strip = TestSection42ChartStrip {
        piece: Kp20PieceId::NearMidpointLeft,
        reverse_i: false,
        reverse_j: false,
        transpose: false,
    };
    let mut best = (
        seed_strip,
        seed_strip,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );

    for underlined_piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let arrow_underlined_piece = match underlined_piece {
            Kp20PieceId::NearMidpointLeft => Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointRight => Kp20PieceId::NearMidpointLeft,
            _ => unreachable!("near-midpoint piece sweep"),
        };
        for transpose in [false, true] {
            for reverse_j in [false, true] {
                for underlined_reverse_i in [false, true] {
                    for arrow_underlined_reverse_i in [false, true] {
                        let underlined_p = TestSection42ChartStrip {
                            piece: underlined_piece,
                            reverse_i: underlined_reverse_i,
                            reverse_j,
                            transpose,
                        };
                        let arrow_underlined_p = TestSection42ChartStrip {
                            piece: arrow_underlined_piece,
                            reverse_i: arrow_underlined_reverse_i,
                            reverse_j,
                            transpose,
                        };
                        let map = TestSection42ChartMap {
                            p,
                            arrow_p,
                            underlined_p,
                            arrow_underlined_p,
                        };
                        let mut vector_constraints =
                            kp20_planar_tau_fig8_indicated_averaging_constraints();
                        vector_constraints.extend(projected_tensor_rows.clone());
                        let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                            &Kp20PlanarTauF3ScalarSolveInput {
                                initial: planar_tau_seed_macro_patch(valence, sector),
                                vector_constraints,
                                scalar_constraints: test_section42_endpoint_aware_chart_rows(
                                    map,
                                    2.0 * c,
                                ),
                                regularization_weight: 1.0e-10,
                            },
                        )
                        .expect("underlined chart handedness solve succeeds");
                        let local_tau =
                            global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                        );
                        let full_tensor_residual = source_rows
                            .iter()
                            .map(|row| {
                                let residual = vector_residual(&solved.tau_macro_patch, row);
                                residual[0].abs().max(residual[1].abs())
                            })
                            .fold(0.0_f64, f64::max);
                        let row = (
                            underlined_p,
                            arrow_underlined_p,
                            solved.max_constraint_residual,
                            validity.orientation.negative_determinants,
                            validity.domain.outside_triangle_samples,
                            validity.orientation.min_determinant,
                            full_tensor_residual,
                        );
                        if row.3 < best.3
                            || (row.3 == best.3 && row.4 < best.4)
                            || (row.3 == best.3 && row.4 == best.4 && row.2 < best.2)
                        {
                            best = row;
                        }
                        if validity.orientation.negative_determinants <= 18
                            || validity.domain.outside_triangle_samples <= 9
                        {
                            eprintln!(
                                "underlined chart underlined={} arrow_underlined={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                                underlined_p.label(),
                                arrow_underlined_p.label(),
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
                    }
                }
            }
        }
    }
    eprintln!(
        "best underlined chart after near-origin fix: underlined={} arrow_underlined={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best.0.label(),
        best.1.label(),
        best.2,
        best.3,
        best.4,
        best.5,
        best.6,
    );
    assert!(
        best.2.is_finite(),
        "underlined chart ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_charts_after_near_origin_fix_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
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
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let candidate_set = [
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    ];
    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best_validity = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        "none",
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_domain = best_validity;

    for candidate in candidate_set {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m: true,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let previous_piece = if variant.swap_pieces {
                                Kp20PieceId::NearMidpointRight
                            } else {
                                Kp20PieceId::NearMidpointLeft
                            };
                            let base_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            for (rotation_label, source_rows) in [
                                ("unrotated", base_rows.clone()),
                                (
                                    "previous-minus-alpha",
                                    rotate_tensor_chart_piece_targets(
                                        &base_rows,
                                        previous_piece,
                                        -angle,
                                    ),
                                ),
                            ] {
                                let projected_tensor_rows =
                                    without_tensor_chart_depth(&source_rows, variant, 1);
                                let mut vector_constraints = fig8_rows.clone();
                                vector_constraints.extend(projected_tensor_rows);
                                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                    &Kp20PlanarTauF3ScalarSolveInput {
                                        initial: planar_tau_seed_macro_patch(valence, sector),
                                        vector_constraints,
                                        scalar_constraints: scalar_constraints.clone(),
                                        regularization_weight: 1.0e-10,
                                    },
                                )
                                .expect("source chart after near-origin fix solve succeeds");
                                let local_tau = global_planar_patch_to_local(
                                    valence,
                                    sector,
                                    solved.tau_macro_patch,
                                );
                                let validity = kp20_planar_bi6_tau_validity_report(
                                    &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                                );
                                let full_tensor_residual = source_rows
                                    .iter()
                                    .map(|row| {
                                        let residual =
                                            vector_residual(&solved.tau_macro_patch, row);
                                        residual[0].abs().max(residual[1].abs())
                                    })
                                    .fold(0.0_f64, f64::max);
                                let row = (
                                    candidate,
                                    rotation_label,
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                    full_tensor_residual,
                                );
                                if row.4 < best_validity.4
                                    || (row.4 == best_validity.4 && row.5 < best_validity.5)
                                    || (row.4 == best_validity.4
                                        && row.5 == best_validity.5
                                        && row.3 < best_validity.3)
                                {
                                    best_validity = row;
                                }
                                if row.5 < best_domain.5
                                    || (row.5 == best_domain.5 && row.4 < best_domain.4)
                                    || (row.5 == best_domain.5
                                        && row.4 == best_domain.4
                                        && row.3 < best_domain.3)
                                {
                                    best_domain = row;
                                }
                                if validity.orientation.negative_determinants <= 18
                                    || validity.domain.outside_triangle_samples == 0
                                {
                                    eprintln!(
                                        "source after near-origin fix candidate={} rotation={} chart={} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                                        candidate.name(),
                                        rotation_label,
                                        variant.label(),
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
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "best source after near-origin fix validity: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.name(),
        best_validity.1,
        best_validity.2.label(),
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
        best_validity.7,
    );
    eprintln!(
        "best source after near-origin fix domain: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_domain.0.name(),
        best_domain.1,
        best_domain.2.label(),
        best_domain.3,
        best_domain.4,
        best_domain.5,
        best_domain.6,
        best_domain.7,
    );
    assert!(
        best_validity.3.is_finite() && best_domain.3.is_finite(),
        "source chart ranking after near-origin fix must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_charts_with_spectral_depth_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let candidate_set = [
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    ];
    let seed_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: false,
        swap_pieces: false,
    };
    let mut best_validity = (
        Kp20Pr08TensorBorderCandidate::VisibleStrip,
        "none",
        seed_variant,
        f64::INFINITY,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    let mut best_domain = best_validity;

    for candidate in candidate_set {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let variant = TestTensorBorderChartVariant {
                                previous_uses_near_m: true,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let previous_piece = if variant.swap_pieces {
                                Kp20PieceId::NearMidpointRight
                            } else {
                                Kp20PieceId::NearMidpointLeft
                            };
                            let current_piece = if variant.swap_pieces {
                                Kp20PieceId::NearMidpointLeft
                            } else {
                                Kp20PieceId::NearMidpointRight
                            };
                            let base_rows =
                                tensor_border_chart_variant_rows(&tensor_border, sector, variant);
                            for (rotation_label, source_rows) in [
                                ("unrotated", base_rows.clone()),
                                (
                                    "previous-minus-alpha",
                                    rotate_tensor_chart_piece_targets(
                                        &base_rows,
                                        previous_piece,
                                        -angle,
                                    ),
                                ),
                            ] {
                                let previous_scaled =
                                    scale_tensor_chart_depth_from_boundary_for_piece(
                                        &source_rows,
                                        variant,
                                        previous_piece,
                                        2,
                                        spectral_cross_depth_scale,
                                    );
                                let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                                    &previous_scaled,
                                    variant,
                                    current_piece,
                                    2,
                                    spectral_cross_depth_scale,
                                );
                                let projected_tensor_rows =
                                    without_tensor_chart_depth(&scaled_rows, variant, 1);
                                let mut vector_constraints = fig8_rows.clone();
                                vector_constraints.extend(projected_tensor_rows);
                                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                    &Kp20PlanarTauF3ScalarSolveInput {
                                        initial: planar_tau_seed_macro_patch(valence, sector),
                                        vector_constraints,
                                        scalar_constraints: scalar_constraints.clone(),
                                        regularization_weight: 1.0e-10,
                                    },
                                )
                                .expect("spectral source chart solve succeeds");
                                let local_tau = global_planar_patch_to_local(
                                    valence,
                                    sector,
                                    solved.tau_macro_patch,
                                );
                                let validity = kp20_planar_bi6_tau_validity_report(
                                    &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                                );
                                let full_tensor_residual = scaled_rows
                                    .iter()
                                    .map(|row| {
                                        let residual =
                                            vector_residual(&solved.tau_macro_patch, row);
                                        residual[0].abs().max(residual[1].abs())
                                    })
                                    .fold(0.0_f64, f64::max);
                                let row = (
                                    candidate,
                                    rotation_label,
                                    variant,
                                    solved.max_constraint_residual,
                                    validity.orientation.negative_determinants,
                                    validity.domain.outside_triangle_samples,
                                    validity.orientation.min_determinant,
                                    full_tensor_residual,
                                );
                                if row.4 < best_validity.4
                                    || (row.4 == best_validity.4 && row.5 < best_validity.5)
                                    || (row.4 == best_validity.4
                                        && row.5 == best_validity.5
                                        && row.3 < best_validity.3)
                                {
                                    best_validity = row;
                                }
                                if row.5 < best_domain.5
                                    || (row.5 == best_domain.5 && row.4 < best_domain.4)
                                    || (row.5 == best_domain.5
                                        && row.4 == best_domain.4
                                        && row.3 < best_domain.3)
                                {
                                    best_domain = row;
                                }
                                if validity.orientation.negative_determinants <= 17
                                    || validity.domain.outside_triangle_samples == 0
                                {
                                    eprintln!(
                                        "spectral source chart candidate={} rotation={} chart={} scale={scale:e} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                                        candidate.name(),
                                        rotation_label,
                                        variant.label(),
                                        scale = spectral_cross_depth_scale,
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
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "best spectral source chart validity: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_validity.0.name(),
        best_validity.1,
        best_validity.2.label(),
        best_validity.3,
        best_validity.4,
        best_validity.5,
        best_validity.6,
        best_validity.7,
    );
    eprintln!(
        "best spectral source chart domain: candidate={} rotation={} chart={} residual={:e} neg={} outside={} min_det={:e} full_tensor={:e}",
        best_domain.0.name(),
        best_domain.1,
        best_domain.2.label(),
        best_domain.3,
        best_domain.4,
        best_domain.5,
        best_domain.6,
        best_domain.7,
    );
    assert!(
        best_validity.3.is_finite() && best_domain.3.is_finite(),
        "spectral source chart ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_cross_depth_orientation_after_near_origin_fix_is_measured() {
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
    let cases = [
        ("baseline", sector_rotated_rows.clone()),
        (
            "reflect-near-midpoint-right-depth2",
            reflect_tensor_chart_depth_across_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
            ),
        ),
        (
            "reflect-near-midpoint-left-depth2",
            reflect_tensor_chart_depth_across_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
            ),
        ),
        (
            "reflect-both-depth2",
            reflect_tensor_chart_depth_across_boundary_for_piece(
                &reflect_tensor_chart_depth_across_boundary_for_piece(
                    &sector_rotated_rows,
                    source_variant,
                    Kp20PieceId::NearMidpointRight,
                    2,
                ),
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
            ),
        ),
    ];
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
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, f64::INFINITY);
    for (label, source_rows) in cases {
        let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("cross-depth orientation after near-origin fix solve succeeds");
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
            "cross-depth after near-origin fix {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
        "best cross-depth after near-origin fix: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "cross-depth orientation after near-origin fix must be finite"
    );
}
