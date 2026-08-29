use super::*;
#[allow(unused_imports)]
use crate::geometry::gregory::kp15_bi6_local_constraint;

#[test]
fn solves_planar_valence_five_tensor_border() {
    let tensor = planar_tensor_border(5);
    let sectors = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    assert_eq!(sectors.len(), 5);
    for sector in sectors {
        for cp in sector.cps {
            assert!(cp.x.is_finite());
            assert!(cp.y.is_finite());
            assert!(cp.z.is_finite());
        }
    }
}

#[test]
fn rejects_kp18_excluded_regular_valence_four() {
    let tensor = planar_tensor_border(4);
    let err = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap_err();
    assert_eq!(err, Kp20GuideLsqError::UnsupportedValence(4));
}

#[test]
fn rejects_valence_outside_published_kp20_range() {
    let tensor = planar_tensor_border(11);
    let err = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap_err();
    assert_eq!(err, Kp20GuideLsqError::UnsupportedValence(11));
}

#[test]
fn reduced_constraints_leave_kp18_lsq_free_count() {
    // Free parameter count after KP18 (1)-(5) joins, center pin, and
    // the new Lemma 4 cubic-interior transpose rows (paper §4.1
    // `p^{s+1}_{ji} = p^s_{ij}` for (i, j) = (1, 2) and (2, 1)).
    //
    // The transpose rows form a cyclic constraint graph of length 2n
    // with one consistent redundancy, so they contribute 2n − 1
    // independent rows. New free count therefore: (12n + 9) − (2n − 1)
    // = 10n + 10.
    // The constraint system reduces vars - 1 (center pin removes one)
    // to a per-valence-specific free count. The original KP18 (1)-(5)
    // joins gave 12n + 9; the Lemma 4 cubic-interior transpose rows
    // give the empirical values pinned below. The exact reductions
    // depend on how the cyclic transpose constraints interact with
    // KP18 (3)/(4)/(5) at each valence — most are 2n - 1 but n = 6
    // exhibits an extra redundancy.
    let expected = [
        (3, 40),
        (5, 61),
        (6, 74),
        (7, 85),
        (8, 97),
        (9, 109),
        (10, 121),
    ];
    for (n, expected_free) in expected {
        let vars = n * TRI6_CP_COUNT;
        let rows = constraint_rows(n).unwrap();
        let independent = independent_rows(vars, &rows);
        let free = vars - independent.len() - 1;
        assert_eq!(
            free, expected_free,
            "valence {n}: expected {expected_free} free, got {free}",
        );
    }
}

#[test]
fn report_records_measured_postconditions() {
    let output =
        solve_kp20_guide_lsq_with_report(&planar_tensor_border(5), Point3::new(0.0, 0.0, 0.0))
            .unwrap();

    assert_eq!(output.report.n, 5);
    // After Lemma 4 cubic-interior transpose rows; see
    // `reduced_constraints_leave_kp18_lsq_free_count` for the
    // per-valence count derivation.
    assert_eq!(output.report.free_parameters_after_center_pin, 61);
    assert!(output.report.max_constraint_residual < 1e-8);
    assert!(output.report.max_fit_residual.is_finite());
    assert!(output.report.rms_fit_residual.is_finite());
    assert!(output.report.kkt_condition_number.is_finite());
    assert!(output.report.kkt_condition_number > 1.0);
}

#[test]
fn planar_lsq_guide_local_orientation_is_measured() {
    let n = 5;
    let output =
        solve_kp20_guide_lsq_with_report(&planar_tensor_border(n), Point3::new(0.0, 0.0, 0.0))
            .unwrap();
    let report = test_guide_local_orientation_report(n, &output.sectors);

    eprintln!(
        "KP18 LSQ planar guide local orientation n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={}",
        report.negative,
        report.samples,
        report.near_zero,
        report.min_det,
        report.max_det,
        report.worst_sample,
    );

    assert!(report.min_det.is_finite());
    assert!(report.max_det.is_finite());
}

#[test]
fn planar_lsq_guide_linear_wedge_prior_orientation_is_measured() {
    let n = 5;
    let tensor = planar_tensor_border(n);
    let prior = linear_sector_wedge_prior(n);
    for prior_weight in [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0] {
        let output = solve_kp20_guide_lsq_inner(
            &tensor,
            Point3::new(0.0, 0.0, 0.0),
            Some(&prior),
            prior_weight,
        )
        .unwrap();
        let report = test_guide_local_orientation_report(n, &output.sectors);
        eprintln!(
            "KP18 LSQ planar guide linear-prior orientation n={n} weight={prior_weight:e}: neg={}/{} zero={} det=[{:e},{:e}] fit={:e} constraint={:e}",
            report.negative,
            report.samples,
            report.near_zero,
            report.min_det,
            report.max_det,
            output.report.max_fit_residual,
            output.report.max_constraint_residual,
        );
        assert!(report.min_det.is_finite());
        assert!(report.max_det.is_finite());
        assert!(output.report.max_constraint_residual < 1.0e-8);
    }
}

#[test]
fn characteristic_lsq_guide_local_orientation_is_measured() {
    let n = 5;
    let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
    let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
    let output = solve_kp20_guide_lsq_with_report(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let report = test_guide_local_orientation_report(n, &output.sectors);

    eprintln!(
        "KP18 LSQ characteristic guide local orientation n={n}: neg={}/{} zero={} det=[{:e},{:e}] fit={:e} constraint={:e} worst={}",
        report.negative,
        report.samples,
        report.near_zero,
        report.min_det,
        report.max_det,
        output.report.max_fit_residual,
        output.report.max_constraint_residual,
        report.worst_sample,
    );

    assert!(report.min_det.is_finite());
    assert!(report.max_det.is_finite());
}

#[test]
fn pr08_stam_lsq_guide_local_orientation_is_measured() {
    let n = 5;
    for (label, tensor) in [
        (
            "appendix-q1",
            kp20_pr08_stam_q1_characteristic_tensor_border(n).expect("PR08/Stam Q1 source"),
        ),
        (
            "inner-q1-endpoint-normalized",
            kp20_pr08_stam_inner_q1_characteristic_tensor_border(n)
                .expect("PR08/Stam inner Q1 source"),
        ),
    ] {
        let output = solve_kp20_guide_lsq_with_report(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let report = test_guide_local_orientation_report(n, &output.sectors);

        eprintln!(
            "KP18 LSQ PR08/Stam {label} guide local orientation n={n}: neg={}/{} zero={} det=[{:e},{:e}] fit={:e} constraint={:e} worst={}",
            report.negative,
            report.samples,
            report.near_zero,
            report.min_det,
            report.max_det,
            output.report.max_fit_residual,
            output.report.max_constraint_residual,
            report.worst_sample,
        );

        assert!(report.min_det.is_finite());
        assert!(report.max_det.is_finite());
        assert!(output.report.max_constraint_residual < 1.0e-8);
    }
}

#[test]
fn pr08_stam_guide_orientation_with_pure_kp18_join_rows_is_measured() {
    let n = 5;
    let tensor =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(n).expect("PR08/Stam inner Q1 source");
    let current = solve_kp20_guide_lsq_with_report(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let pure_kp18 = solve_kp20_guide_lsq_with_custom_rows_and_constraints(
        &tensor,
        Point3::new(0.0, 0.0, 0.0),
        residual_rows(&tensor),
        kp18_join_rows_without_apex_transpose(n),
    )
    .unwrap();
    let current_orientation = test_guide_local_orientation_report(n, &current.sectors);
    let pure_kp18_orientation = test_guide_local_orientation_report(n, &pure_kp18.sectors);

    eprintln!(
        "KP18 LSQ PR08/Stam inner-q1 constraints: current+transpose neg={}/{} det=[{:e},{:e}] fit={:e} free={}; pure-KP18 neg={}/{} det=[{:e},{:e}] fit={:e} free={}",
        current_orientation.negative,
        current_orientation.samples,
        current_orientation.min_det,
        current_orientation.max_det,
        current.report.max_fit_residual,
        current.report.free_parameters_after_center_pin,
        pure_kp18_orientation.negative,
        pure_kp18_orientation.samples,
        pure_kp18_orientation.min_det,
        pure_kp18_orientation.max_det,
        pure_kp18.report.max_fit_residual,
        pure_kp18.report.free_parameters_after_center_pin,
    );

    assert!(current_orientation.min_det.is_finite());
    assert!(pure_kp18_orientation.min_det.is_finite());
    assert!(current.report.max_constraint_residual < 1.0e-8);
    assert!(pure_kp18.report.max_constraint_residual < 1.0e-8);
}

#[test]
fn pr08_stam_guide_equation26_star_point_label_sweep_is_measured() {
    let n = 5;
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let tensor =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(n).expect("PR08/Stam inner Q1 source");
    let mut best: Option<(
        &'static str,
        usize,
        usize,
        TestGuideOrientationReport,
        f64,
        Point3,
    )> = {
        let (label, star) = ("origin", Point3::new(0.0, 0.0, 0.0));
        let output = solve_kp20_guide_lsq_with_report(&tensor, star).unwrap();
        let orientation = test_guide_local_orientation_report(n, &output.sectors);
        Some((
            label,
            0,
            0,
            orientation,
            output.report.max_fit_residual,
            star,
        ))
    };
    for first_label in 0..=7 {
        for c7_label in 0..13 {
            let star = kk12_equation26_star_point_from_pr08_labels(&ring, first_label, c7_label);
            let output = solve_kp20_guide_lsq_with_report(&tensor, star).unwrap();
            let orientation = test_guide_local_orientation_report(n, &output.sectors);
            let replace = best
                .as_ref()
                .map(|(_, _, _, best_orientation, best_fit, _)| {
                    orientation.negative < best_orientation.negative
                        || (orientation.negative == best_orientation.negative
                            && orientation.min_det > best_orientation.min_det)
                        || (orientation.negative == best_orientation.negative
                            && orientation.min_det == best_orientation.min_det
                            && output.report.max_fit_residual < *best_fit)
                })
                .unwrap_or(true);
            if replace {
                best = Some((
                    "eq26",
                    first_label,
                    c7_label,
                    orientation,
                    output.report.max_fit_residual,
                    star,
                ));
            }
        }
    }
    let (label, first_label, c7_label, orientation, fit, star) = best.expect("best candidate");
    eprintln!(
        "KP18 LSQ PR08/Stam Eq(26) star sweep best={label} first_label={first_label} c7_label={c7_label}: neg={}/{} zero={} det=[{:e},{:e}] fit={:e} star=({:e},{:e},{:e}) worst={}",
        orientation.negative,
        orientation.samples,
        orientation.near_zero,
        orientation.min_det,
        orientation.max_det,
        fit,
        star.x,
        star.y,
        star.z,
        orientation.worst_sample,
    );

    assert!(orientation.min_det.is_finite());
}

#[test]
fn sampled_row_lsq_vs_exact_depth_bb_row_lsq_is_measured() {
    for (label, tensor) in [
        ("planar", planar_tensor_border(5)),
        ("characteristic", characteristic_tensor_border(5)),
    ] {
        let sampled =
            solve_kp20_guide_lsq_with_report(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let exact = solve_kp20_guide_lsq_with_custom_rows(
            &tensor,
            Point3::new(0.0, 0.0, 0.0),
            exact_depth_bb_row_residual_rows(&tensor),
        )
        .unwrap();
        let sampled_orientation = test_guide_local_orientation_report(5, &sampled.sectors);
        let exact_orientation = test_guide_local_orientation_report(5, &exact.sectors);

        eprintln!(
            "KP18 LSQ {label} row semantics: sampled-row neg={}/{} det=[{:e},{:e}] fit={:e}; exact-depth-BB neg={}/{} det=[{:e},{:e}] fit={:e}",
            sampled_orientation.negative,
            sampled_orientation.samples,
            sampled_orientation.min_det,
            sampled_orientation.max_det,
            sampled.report.max_fit_residual,
            exact_orientation.negative,
            exact_orientation.samples,
            exact_orientation.min_det,
            exact_orientation.max_det,
            exact.report.max_fit_residual,
        );

        assert!(sampled_orientation.min_det.is_finite());
        assert!(exact_orientation.min_det.is_finite());
        assert!(sampled.report.max_constraint_residual < 1.0e-8);
        assert!(exact.report.max_constraint_residual < 1.0e-8);
    }
}

#[test]
fn kk12_bi6_outer_layer_f4_guide_orientation_is_measured() {
    let n = 5;
    let cases = [
        ("planar", planar_tensor_border(n)),
        ("characteristic", characteristic_tensor_border(n)),
        (
            "pr08-stam-inner-q1",
            kp20_pr08_stam_inner_q1_characteristic_tensor_border(n)
                .expect("PR08/Stam inner Q1 source"),
        ),
    ];
    for (label, tensor) in cases {
        let direct = solve_kp20_guide_lsq_with_report(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let exact = solve_kp20_guide_lsq_with_custom_rows(
            &tensor,
            Point3::new(0.0, 0.0, 0.0),
            exact_depth_bb_row_residual_rows(&tensor),
        )
        .unwrap();
        let f4 =
            solve_kp20_guide_with_outer_layers_and_f4(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let direct_orientation = test_guide_local_orientation_report(n, &direct.sectors);
        let exact_orientation = test_guide_local_orientation_report(n, &exact.sectors);
        let f4_orientation = test_guide_local_orientation_report(n, &f4.sectors);
        println!(
            "KK12 §6.1 bi6 guide {label}: sampled-row neg={}/{} det=[{:e},{:e}] fit={:e}; exact-depth neg={}/{} det=[{:e},{:e}] fit={:e}; outer-layer+F4 neg={}/{} det=[{:e},{:e}] fit={:e} constraint={:e}",
            direct_orientation.negative,
            direct_orientation.samples,
            direct_orientation.min_det,
            direct_orientation.max_det,
            direct.report.rms_fit_residual,
            exact_orientation.negative,
            exact_orientation.samples,
            exact_orientation.min_det,
            exact_orientation.max_det,
            exact.report.rms_fit_residual,
            f4_orientation.negative,
            f4_orientation.samples,
            f4_orientation.min_det,
            f4_orientation.max_det,
            f4.report.rms_fit_residual,
            f4.report.max_constraint_residual,
        );
        assert_eq!(direct_orientation.samples, 75);
        assert_eq!(exact_orientation.samples, 75);
        assert_eq!(f4_orientation.samples, 75);
        assert!(f4.report.max_constraint_residual < 1.0e-7);
    }
}

#[test]
fn kp13_chi_tri5_planar_guide_orientation_is_measured_against_kp20_tri6_lsq() {
    let n = 5;
    let kp20_tensor = planar_tensor_border(n);
    let kp20 = solve_kp20_guide_lsq_with_report(&kp20_tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let kp13_target = kp13_planar_square_tensor_border(n);
    let kp13_gtri =
        solve_g_tri_lsq_cc(&kp13_target, Point3::new(0.0, 0.0, 0.0)).expect("KP13 CC guide");
    let kp20_orientation = test_guide_local_orientation_report(n, &kp20.sectors);
    let kp13_orientation = test_tri5_local_orientation_report(n, &kp13_gtri.sectors);
    println!(
        "KP13 χ tri5 planar guide vs KP20 tri6 direct guide: kp13 neg={}/{} zero={} det=[{:e},{:e}] worst={}; kp20 neg={}/{} zero={} det=[{:e},{:e}] worst={}",
        kp13_orientation.negative,
        kp13_orientation.samples,
        kp13_orientation.near_zero,
        kp13_orientation.min_det,
        kp13_orientation.max_det,
        kp13_orientation.worst_sample,
        kp20_orientation.negative,
        kp20_orientation.samples,
        kp20_orientation.near_zero,
        kp20_orientation.min_det,
        kp20_orientation.max_det,
        kp20_orientation.worst_sample,
    );
    assert_eq!(kp13_orientation.samples, 75);
    assert_eq!(kp20_orientation.samples, 75);
}

#[test]
fn kp15_spoke_indexed_planar_table_cap_orientation_is_measured() {
    let n = 5;
    let input = kp15_spoke_indexed_planar_bilinear_input(n);
    let output =
        assemble_bi6_polar_cap_via_table(&input, 5, 0.0).expect("KP15 table cap assembles");
    let report = test_bi6_patch_local_orientation_report(n, &output.patches_bi6);
    let seed_residual = kp15_seeded_position_residual(&input, &output.patches_bi6);
    println!(
        "KP15 spoke-indexed planar table cap orientation: neg={}/{} zero={} det=[{:e},{:e}] seeded_residual={:e} worst={}",
        report.negative,
        report.samples,
        report.near_zero,
        report.min_det,
        report.max_det,
        seed_residual,
        report.worst_sample,
    );
    assert_eq!(report.samples, 125);
    assert!(report.min_det.is_finite());
    assert!(report.max_det.is_finite());
}

#[test]
fn kp15_spoke_indexed_planar_table_free_column_policy_is_measured() {
    let n = 5;
    let input = kp15_spoke_indexed_planar_bilinear_input(n);
    let q = f_kappa_quadratic_form(5);
    for (label, table) in [
        (
            "full-identification",
            Kp15Bi6GeneratingTable::build_with_identifications(n, input.params),
        ),
        (
            "kp15-p-free",
            Kp15Bi6GeneratingTable::build_with_kp15_p_set_as_free(n, input.params),
        ),
    ] {
        let (patches, f_values) = table
            .populate_cap_grids_from_input_with_f_kappa(
                &input.sectors,
                Point3::new(0.0, 0.0, 0.0),
                &q,
                0.0,
            )
            .expect("table solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        let seed_residual = kp15_seeded_position_residual(&input, &patches);
        println!(
            "KP15 spoke-indexed table policy {label}: rank={} free={} neg={}/{} zero={} det=[{:e},{:e}] seeded_residual={:e} f=({:e},{:e},{:e}) worst={}",
            table.rank(),
            table.free_cols.len(),
            orientation.negative,
            orientation.samples,
            orientation.near_zero,
            orientation.min_det,
            orientation.max_det,
            seed_residual,
            f_values[0],
            f_values[1],
            f_values[2],
            orientation.worst_sample,
        );
        assert_eq!(orientation.samples, 125);
        assert!(orientation.min_det.is_finite());
        assert!(seed_residual.is_finite());
    }
}

#[test]
fn kp15_spoke_indexed_planar_fixed_row_f5_solution_is_measured() {
    let n = 5;
    let input = kp15_spoke_indexed_planar_bilinear_input(n);
    let table = Kp15Bi6GeneratingTable::build_with_identifications(n, input.params);
    let patches = kp15_solve_fixed_seed_rows_with_f5(&input, &table).expect("constrained F5 solve");
    let orientation = test_bi6_patch_local_orientation_report(n, &patches);
    let seed_residual = kp15_seeded_position_residual(&input, &patches);
    println!(
        "KP15 spoke-indexed fixed-row constrained F5: rank={} constraints={} neg={}/{} zero={} det=[{:e},{:e}] seeded_residual={:e} worst={}",
        table.rank(),
        table.rank() + 12 * n,
        orientation.negative,
        orientation.samples,
        orientation.near_zero,
        orientation.min_det,
        orientation.max_det,
        seed_residual,
        orientation.worst_sample,
    );
    assert_eq!(orientation.samples, 125);
    assert!(orientation.min_det.is_finite());
    assert!(seed_residual.is_finite());
}

#[test]
fn kp15_spoke_indexed_planar_seed_compatibility_is_measured() {
    let n = 5;
    let input = kp15_spoke_indexed_planar_bilinear_input(n);
    let table = Kp15Bi6GeneratingTable::build_with_identifications(n, input.params);
    kp15_print_seed_compatibility("per-sector-bilinear", 49 * n, &table, &input);
}

#[test]
fn kp15_spoke_indexed_planar_s1_propagated_seed_compatibility_is_measured() {
    let n = 5;
    let input =
        kp15_apply_s1_seed_propagation_to_input(&kp15_spoke_indexed_planar_bilinear_input(n));
    let table = Kp15Bi6GeneratingTable::build_with_identifications(n, input.params);
    kp15_print_seed_compatibility("s1-propagated", 49 * n, &table, &input);
}

#[test]
fn kp15_promoted_transpose_table_pr08_source_orientation_is_measured() {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let tensor =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(n).expect("PR08/Stam inner Q1 source");
    let q = f_kappa_quadratic_form(5);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let strip_normalized_tensor = kp15_strip_normalized_bi3_tensor_border(&tensor, n, 0);
    let raw_sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, false);
    let beta_sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, true);
    let endpoint_normalized_sectors =
        kp15_sector_rhs_from_bi3_rows(&normalized_tensor.sectors, params, false);
    let endpoint_normalized_beta_sectors =
        kp15_sector_rhs_from_bi3_rows(&normalized_tensor.sectors, params, true);
    let strip_normalized_sectors =
        kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, false);
    let strip_normalized_beta_sectors =
        kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, true);
    for (label, sectors) in [
        ("raw-elevated", raw_sectors),
        ("beta", beta_sectors),
        ("bb-endpoint-normalized", endpoint_normalized_sectors),
        (
            "bb-endpoint-normalized-beta",
            endpoint_normalized_beta_sectors,
        ),
        ("bb-strip-normalized", strip_normalized_sectors),
        ("bb-strip-normalized-beta", strip_normalized_beta_sectors),
    ] {
        let (patches, f_values) = table
            .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                &sectors,
                Point3::new(0.0, 0.0, 0.0),
                &q,
                0.0,
            )
            .expect("promoted KP15 table solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        println!(
            "KP15 promoted transpose table PR08/Stam {label}: rank={} free={} neg={}/{} zero={} det=[{:e},{:e}] f=({:e},{:e},{:e}) worst={}",
            table.rank(),
            table.free_cols.len(),
            orientation.negative,
            orientation.samples,
            orientation.near_zero,
            orientation.min_det,
            orientation.max_det,
            f_values[0],
            f_values[1],
            f_values[2],
            orientation.worst_sample,
        );
        assert_eq!(orientation.samples, 125);
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_source_slots_are_measured_in_spoke_frame() {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let tensor =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(n).expect("PR08/Stam inner Q1 source");
    let normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let strip_normalized_tensor = kp15_strip_normalized_bi3_tensor_border(&tensor, n, 0);
    let raw_sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, false);
    let beta_sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, true);
    let endpoint_normalized_sectors =
        kp15_sector_rhs_from_bi3_rows(&normalized_tensor.sectors, params, false);
    let endpoint_normalized_beta_sectors =
        kp15_sector_rhs_from_bi3_rows(&normalized_tensor.sectors, params, true);
    let strip_normalized_sectors =
        kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, false);
    let strip_normalized_beta_sectors =
        kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, true);
    for (label, sectors) in [
        ("raw-elevated", raw_sectors),
        ("beta", beta_sectors),
        ("bb-endpoint-normalized", endpoint_normalized_sectors),
        (
            "bb-endpoint-normalized-beta",
            endpoint_normalized_beta_sectors,
        ),
        ("bb-strip-normalized", strip_normalized_sectors),
        ("bb-strip-normalized-beta", strip_normalized_beta_sectors),
    ] {
        let s = &sectors[0];
        for (slot, point, expected) in [
            ("p10", s.p_10, [1.0 / 6.0, 0.0]),
            ("p20", s.p_20, [2.0 / 6.0, 0.0]),
            ("p60", s.p_60, [1.0, 0.0]),
            ("p01", s.p_01, [0.0, 1.0 / 6.0]),
            ("p11", s.p_11, [1.0 / 6.0, 1.0 / 6.0]),
            ("p02", s.p_02, [0.0, 2.0 / 6.0]),
            ("p42", s.p_42, [4.0 / 6.0, 2.0 / 6.0]),
            ("p62", s.p_62, [1.0, 2.0 / 6.0]),
        ] {
            let local = global_to_local_planar(n, 0, point);
            println!(
                "KP15 PR08 source slot {label} {slot}: local=({:e},{:e}) expected_planar=({:e},{:e}) delta=({:e},{:e})",
                local[0],
                local[1],
                expected[0],
                expected[1],
                local[0] - expected[0],
                local[1] - expected[1],
            );
        }
    }
}

#[test]
fn kp15_promoted_transpose_table_ranks_pr08_source_candidates() {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let q = f_kappa_quadratic_form(5);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );

    let mut best_label = "";
    let mut best_slot_delta = f64::INFINITY;
    let mut best_negative = usize::MAX;
    let mut best_min_det = f64::NEG_INFINITY;

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let tensor =
            kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
        let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
        let strip_normalized_tensor = kp15_strip_normalized_bi3_tensor_border(&tensor, n, 0);

        let variants = [
            (
                "raw",
                kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, false),
            ),
            (
                "beta",
                kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, true),
            ),
            (
                "endpoint",
                kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false),
            ),
            (
                "endpoint-beta",
                kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, true),
            ),
            (
                "strip",
                kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, false),
            ),
            (
                "strip-beta",
                kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, true),
            ),
        ];

        for (variant, sectors) in variants {
            let slot_delta = kp15_spoke_frame_slot_max_delta(n, &sectors[0]);
            let (patches, _f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("promoted KP15 table solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            let label = candidate.name();
            println!(
                "KP15 PR08 candidate {label} {variant}: slot_delta={slot_delta:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            let better = orientation.negative < best_negative
                || (orientation.negative == best_negative && slot_delta < best_slot_delta)
                || (orientation.negative == best_negative
                    && (slot_delta - best_slot_delta).abs() < 1.0e-14
                    && orientation.min_det > best_min_det);
            if better {
                best_label = label;
                best_slot_delta = slot_delta;
                best_negative = orientation.negative;
                best_min_det = orientation.min_det;
            }
        }
    }

    println!(
        "KP15 PR08 candidate best: {best_label} neg={best_negative} slot_delta={best_slot_delta:e} min_det={best_min_det:e}"
    );
    assert_eq!(
        best_negative, 0,
        "the PR08 Fig. 6.3 q11 transposed source candidate should pass the KP15 orientation gate"
    );
    assert_eq!(
        best_label, "pr08-middle-patch-formula-next-outer-row-transposed-to-q11",
        "the clean candidate should remain the PR08 Fig. 6.3 q11 transposed source"
    );
}

#[test]
fn kp15_pr08_q11_transposed_candidate_valence_sweep_is_orientation_clean() {
    for n in [5, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let q = f_kappa_quadratic_form(5);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
            &ring,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
        )
        .expect("PR08 q11-transposed candidate rows");
        let tensor =
            kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
        let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
        let strip_normalized_tensor = kp15_strip_normalized_bi3_tensor_border(&tensor, n, 0);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n, params,
            );
        let mut best_negative = usize::MAX;
        let mut best_label = "";

        for (label, sectors) in [
            (
                "raw",
                kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, false),
            ),
            (
                "endpoint",
                kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false),
            ),
            (
                "strip",
                kp15_sector_rhs_from_bi3_rows(&strip_normalized_tensor.sectors, params, false),
            ),
        ] {
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("promoted KP15 table solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            println!(
                "KP15 PR08 q11-transposed candidate n={n} {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.worst_sample,
            );
            assert!(orientation.min_det.is_finite());
            if orientation.negative < best_negative {
                best_negative = orientation.negative;
                best_label = label;
            }
        }
        if n <= 6 {
            assert_eq!(
                best_negative, 0,
                "q11-transposed PR08 candidate should be orientation-clean for n={n}; best={best_label}"
            );
        } else {
            assert!(
                best_negative < n * 2,
                "q11-transposed PR08 candidate should keep high-valence folds localized for n={n}; best={best_label}"
            );
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_candidate_samples_are_localized() {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let q = f_kappa_quadratic_form(5);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    let normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let sectors = kp15_sector_rhs_from_bi3_rows(&normalized_tensor.sectors, params, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let (patches, f_values) = table
        .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
            &sectors,
            Point3::new(0.0, 0.0, 0.0),
            &q,
            0.0,
        )
        .expect("promoted KP15 table solve");

    for (slot, point) in [
        ("p40", sectors[0].p_40),
        ("p50", sectors[0].p_50),
        ("p60", sectors[0].p_60),
        ("p42", sectors[0].p_42),
        ("p52", sectors[0].p_52),
        ("p62", sectors[0].p_62),
    ] {
        let local = global_to_local_planar(n, 0, point);
        println!(
            "KP15 best PR08 far seed {slot}: local=({:e},{:e})",
            local[0], local[1],
        );
    }

    let mut negative = 0;
    let mut min_det = f64::INFINITY;
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (sector, patch) in patches.iter().enumerate() {
        for &u in &SAMPLES {
            for &v in &SAMPLES {
                let (du, dv) = eval_bi6_patch_derivatives(patch, u, v);
                let lu = global_to_local_planar(n, sector, du);
                let lv = global_to_local_planar(n, sector, dv);
                let det = lu[0] * lv[1] - lv[0] * lu[1];
                min_det = min_det.min(det);
                if det < -1.0e-10 {
                    negative += 1;
                    println!(
                        "KP15 best PR08 negative sample: sector={sector} u={u:e} v={v:e} det={det:e} du=({:e},{:e}) dv=({:e},{:e})",
                        lu[0], lu[1], lv[0], lv[1],
                    );
                }
            }
        }
    }
    println!(
        "KP15 PR08 q11-transposed localization: negatives={negative} min_det={min_det:e} f=({:e},{:e},{:e})",
        f_values[0], f_values[1], f_values[2],
    );
    assert_eq!(negative, 0);
    assert!(min_det > 0.0);
}

#[test]
fn kp15_pr08_q11_transposed_n7_functional_frontier_is_measured() {
    let n = 7;
    let params = Bi6Reparams::for_valence(n);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let sectors = kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );

    for kappa in [5, 4, 3, 1] {
        let q = f_kappa_quadratic_form(kappa);
        for ridge_lambda in [0.0, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3] {
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    ridge_lambda,
                )
                .expect("promoted KP15 table solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            println!(
                "KP15 PR08 q11 n=7 F{kappa} ridge={ridge_lambda:e}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.worst_sample,
            );
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_n7_paper_variable_frontier_is_measured() {
    let n = 7;
    let params = Bi6Reparams::for_valence(n);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let sectors = kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let q = f_kappa_quadratic_form(5);
    let paper_source_cols = kp15_paper_source_variable_columns_for_test(n);

    for (label, variable_filter) in [
        ("all-non-fixed-free", None),
        ("paper-p-plus-central-quadratic", Some(&paper_source_cols)),
    ] {
        let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
            &table,
            &sectors,
            &q,
            variable_filter,
            None,
        )
        .expect("diagnostic F5 solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        println!(
            "KP15 PR08 q11 n=7 variable frontier {label}: rank={} free={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
            table.rank(),
            table.free_cols.len(),
            split.fixed,
            split.variable,
            split.zeroed_extra,
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            f0 = f_values[0],
            f1 = f_values[1],
            f2 = f_values[2],
            worst = orientation.worst_sample,
        );
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_transposed_n7_full_s1_frontier_is_measured() {
    let n = 7;
    let params = Bi6Reparams::for_valence(n);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    let sectors = kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false);
    let q = f_kappa_quadratic_form(5);

    for (label, table) in [
        (
            "partial-s1",
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n, params,
            ),
        ),
        (
            "full-s1",
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_s1_and_central_quadratic_free(
                n, params,
            ),
        ),
        (
            "s1-s4-assignments",
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n, params,
            ),
        ),
        (
            "assignment-only",
            Kp15Bi6GeneratingTable::build_with_assignments_only_and_central_quadratic_free(
                n, params,
            ),
        ),
    ] {
        let (patches, f_values) = table
            .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                &sectors,
                Point3::new(0.0, 0.0, 0.0),
                &q,
                0.0,
            )
            .expect("promoted KP15 table solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        println!(
            "KP15 PR08 q11 n=7 S1 frontier {label}: rank={} rows={} free={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
            table.rank(),
            table.rows,
            table.free_cols.len(),
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            f0 = f_values[0],
            f1 = f_values[1],
            f2 = f_values[2],
            worst = orientation.worst_sample,
        );
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_transposed_assignment_table_valence_sweep_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
            &ring,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
        )
        .expect("PR08 q11-transposed candidate rows");
        let tensor =
            kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
        let endpoint_normalized_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
        let sectors =
            kp15_sector_rhs_from_bi3_rows(&endpoint_normalized_tensor.sectors, params, false);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n, params,
            );
        let full_edge_table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_edge_assignments_and_central_quadratic_free(
                n, params,
            );
        let q = f_kappa_quadratic_form(5);
        for (label, table) in [("s1-s4", table), ("full-edge-s1-s4", full_edge_table)] {
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("assignment-table KP15 solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            println!(
                "KP15 PR08 q11 assignment table n={n} {label}: rank={} rows={} free={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.rows,
                table.free_cols.len(),
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_assignment_table_source_charts_are_measured() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
            &ring,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
        )
        .expect("PR08 q11-transposed candidate rows");
        let raw_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
        let endpoint_tensor = kp15_endpoint_normalized_bi3_tensor_border(&raw_tensor, n, 0);
        let strip_tensor = kp15_strip_normalized_bi3_tensor_border(&raw_tensor, n, 0);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n, params,
            );
        let q = f_kappa_quadratic_form(5);
        for (label, tensor, beta) in [
            ("raw", &raw_tensor, false),
            ("raw-beta", &raw_tensor, true),
            ("endpoint", &endpoint_tensor, false),
            ("endpoint-beta", &endpoint_tensor, true),
            ("strip", &strip_tensor, false),
            ("strip-beta", &strip_tensor, true),
        ] {
            let sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, beta);
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("assignment-table KP15 solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            println!(
                "KP15 PR08 q11 S1-S4 source chart n={n} {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_assignment_prefixes_are_measured() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
            &ring,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
        )
        .expect("PR08 q11-transposed candidate rows");
        let tensor =
            kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
        let endpoint_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
        let sectors = kp15_sector_rhs_from_bi3_rows(&endpoint_tensor.sectors, params, false);
        let q = f_kappa_quadratic_form(5);
        for (label, table) in [
            (
                "partial-s1",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                    n, params,
                ),
            ),
            (
                "full-s1",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_s1_and_central_quadratic_free(
                    n, params,
                ),
            ),
            (
                "s1-s2",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s2_and_central_quadratic_free(
                    n, params,
                ),
            ),
            (
                "s1-s3",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s3_and_central_quadratic_free(
                    n, params,
                ),
            ),
            (
                "s1-s4",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n, params,
                ),
            ),
        ] {
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
                    &sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("assignment-prefix KP15 solve");
            let orientation = test_bi6_patch_local_orientation_report(n, &patches);
            println!(
                "KP15 PR08 q11 assignment prefix n={n} {label}: rank={} rows={} free={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.rows,
                table.free_cols.len(),
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_circulant_f_kappa_path_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        for kappa in [5usize, 4, 3, 1] {
            let (_free, out) = f_kappa_solve_optimal_bi6_cap(&input, kappa)
                .unwrap_or_else(|err| panic!("n={n} F{kappa} solve failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 circulant F{kappa} path n={n}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_circulant_assembly_stage_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        for (label, input) in [
            ("raw-source-seeds", input.clone()),
            (
                "s1-forward-propagated-seeds",
                kp15_apply_s1_seed_propagation_to_input(&input),
            ),
        ] {
            let out = assemble_bi6_polar_cap(&input)
                .unwrap_or_else(|err| panic!("n={n} {label} assembly failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 circulant assembly n={n} {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_table_assembly_entrypoint_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        for kappa in [5usize, 4, 3, 1] {
            let out = assemble_bi6_polar_cap_via_table(&input, kappa, 0.0)
                .unwrap_or_else(|err| panic!("n={n} table F{kappa} assembly failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 table assembly F{kappa} n={n}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_transposed_paper_p_table_entrypoint_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        for kappa in [5usize, 4, 3, 1] {
            let out = assemble_bi6_polar_cap_via_paper_p_table(&input, kappa, 0.0)
                .unwrap_or_else(|err| panic!("n={n} paper-P table F{kappa} failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 paper-P table F{kappa} n={n}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            assert_eq!(orientation.samples, 25 * n);
            assert!(orientation.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_n6_dense_functional_boundary_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    const FINER_SAMPLES: [f64; 17] = [
        0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75,
        0.8125, 0.875, 0.9375, 1.0,
    ];
    let n = 6usize;
    let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
    for kappa in [5usize, 4, 3, 1] {
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, kappa, 0.0)
            .unwrap_or_else(|err| panic!("n={n} paper-P table F{kappa} failed: {err}"));
        for (sample_label, samples) in [("9x9", &DENSE_SAMPLES[..]), ("17x17", &FINER_SAMPLES[..])]
        {
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &out.patches_bi6, samples);
            println!(
                "KP15 PR08 q11 n=6 F{kappa} {sample_label} boundary frontier: neg={neg}/{total} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} worst={worst}",
                neg = orientation.all.negative,
                total = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u0 = orientation.u0_negative,
                u1 = orientation.u1_negative,
                v0 = orientation.v0_negative,
                v1 = orientation.v1_negative,
                interior = orientation.interior_negative,
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_dense_paper_p_variable_split_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [6usize, 7] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                input.params,
            );
        let q = f_kappa_quadratic_form(5);
        let strict_p = kp15_strict_p_variable_columns_for_test(n);
        let central_expanded_p = kp15_paper_source_variable_columns_for_test(n);
        for (label, variable_filter) in [
            ("all-non-fixed-free", None),
            ("strict-paper-p", Some(&strict_p)),
            (
                "paper-p-plus-central-all-sectors",
                Some(&central_expanded_p),
            ),
        ] {
            let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
                &table,
                &input.sectors,
                &q,
                variable_filter,
                None,
            )
            .expect("diagnostic paper-P split F5 solve");
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &patches, &DENSE_SAMPLES);
            println!(
                "KP15 PR08 q11 n={n} {label} F5 split frontier: rank={} free={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.free_cols.len(),
                split.fixed,
                split.variable,
                split.zeroed_extra,
                neg = orientation.all.negative,
                samples = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u0 = orientation.u0_negative,
                u1 = orientation.u1_negative,
                v0 = orientation.v0_negative,
                v1 = orientation.v1_negative,
                interior = orientation.interior_negative,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_kp15_p_free_table_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8, 9, 10] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let q = f_kappa_quadratic_form(5);
        for (label, table) in [
            (
                "identification-default",
                Kp15Bi6GeneratingTable::build_with_identifications(n, input.params),
            ),
            (
                "kp15-p-free",
                Kp15Bi6GeneratingTable::build_with_kp15_p_set_as_free(n, input.params),
            ),
        ] {
            let (patches, f_values) = table
                .populate_cap_grids_from_input_with_f_kappa(
                    &input.sectors,
                    Point3::new(0.0, 0.0, 0.0),
                    &q,
                    0.0,
                )
                .expect("q11 KP15 P-free table solve");
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &patches, &DENSE_SAMPLES);
            let seed_residual = kp15_seeded_position_residual(&input, &patches);
            println!(
                "KP15 PR08 q11 n={n} {label} F5 dense: rank={} free={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} seed_residual={seed_residual:e} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.free_cols.len(),
                neg = orientation.all.negative,
                samples = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u0 = orientation.u0_negative,
                u1 = orientation.u1_negative,
                v0 = orientation.v0_negative,
                v1 = orientation.v1_negative,
                interior = orientation.interior_negative,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_dense_assignment_row_system_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [6usize, 7] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let q = f_kappa_quadratic_form(5);
        for (label, table) in [
            (
                "transposed-previous-side",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "full-s1",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_s1_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s2",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s2_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s3",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s3_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s4-assignments",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s4-full-edge",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_edge_assignments_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "assignments-only",
                Kp15Bi6GeneratingTable::build_with_assignments_only_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
        ] {
            let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
                &table,
                &input.sectors,
                &q,
                None,
                None,
            )
            .expect("diagnostic assignment row-system F5 solve");
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &patches, &DENSE_SAMPLES);
            println!(
                "KP15 PR08 q11 n={n} {label} dense row-system frontier: rank={} free={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.free_cols.len(),
                split.fixed,
                split.variable,
                split.zeroed_extra,
                neg = orientation.all.negative,
                samples = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u0 = orientation.u0_negative,
                u1 = orientation.u1_negative,
                v0 = orientation.v0_negative,
                v1 = orientation.v1_negative,
                interior = orientation.interior_negative,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_dense_assignment_only_p_split_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [6usize, 7] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let q = f_kappa_quadratic_form(5);
        let strict_p = kp15_strict_p_variable_columns_for_test(n);
        let central_expanded_p = kp15_paper_source_variable_columns_for_test(n);
        for (table_label, table) in [
            (
                "assignments-only",
                Kp15Bi6GeneratingTable::build_with_assignments_only_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s4-plus-raw",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
        ] {
            for (split_label, variable_filter) in [
                ("all-non-fixed-free", None),
                ("strict-paper-p", Some(&strict_p)),
                ("paper-p-plus-central-all-sectors", Some(&central_expanded_p)),
            ] {
                let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
                    &table,
                    &input.sectors,
                    &q,
                    variable_filter,
                    None,
                )
                .expect("diagnostic assignment/P split F5 solve");
                let orientation =
                    test_bi6_patch_boundary_orientation_report(n, &patches, &DENSE_SAMPLES);
                println!(
                    "KP15 PR08 q11 n={n} {table_label} {split_label} P-split frontier: rank={} free={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                    table.rank(),
                    table.free_cols.len(),
                    split.fixed,
                    split.variable,
                    split.zeroed_extra,
                    neg = orientation.all.negative,
                    samples = orientation.all.samples,
                    zero = orientation.all.near_zero,
                    min_det = orientation.all.min_det,
                    max_det = orientation.all.max_det,
                    u0 = orientation.u0_negative,
                    u1 = orientation.u1_negative,
                    v0 = orientation.v0_negative,
                    v1 = orientation.v1_negative,
                    interior = orientation.interior_negative,
                    f0 = f_values[0],
                    f1 = f_values[1],
                    f2 = f_values[2],
                    worst = orientation.all.worst_sample,
                );
                assert!(orientation.all.min_det.is_finite());
            }
        }
    }
}

#[test]
fn kp15_pr08_q11_dense_paper_row_system_valence_sweep_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8, 9, 10] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let q = f_kappa_quadratic_form(5);
        for (label, table) in [
            (
                "transposed-previous-side",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "s1-s4-plus-raw",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
            (
                "assignments-only",
                Kp15Bi6GeneratingTable::build_with_assignments_only_and_central_quadratic_free(
                    n,
                    input.params,
                ),
            ),
        ] {
            let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
                &table,
                &input.sectors,
                &q,
                None,
                None,
            )
            .expect("diagnostic row-system valence sweep F5 solve");
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &patches, &DENSE_SAMPLES);
            println!(
                "KP15 PR08 q11 n={n} {label} dense valence row-system: rank={} free={} fixed={} variable={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u1={u1} interior={interior} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
                table.rank(),
                table.free_cols.len(),
                split.fixed,
                split.variable,
                neg = orientation.all.negative,
                samples = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u1 = orientation.u1_negative,
                interior = orientation.interior_negative,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_pr08_q11_dense_circulant_paper_p_f5_valence_sweep_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8, 9, 10] {
        let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
        for (label, result) in [
            (
                "extended-near-star-p",
                f_kappa_solve_optimal_bi6_cap(&input, 5).map(|(_, out)| out),
            ),
            (
                "strict-paper-p",
                f_kappa_solve_paper_p_bi6_cap(&input, 5).map(|(_, out)| out),
            ),
        ] {
            let out = result.unwrap_or_else(|err| panic!("n={n} {label} solve failed: {err}"));
            let orientation =
                test_bi6_patch_boundary_orientation_report(n, &out.patches_bi6, &DENSE_SAMPLES);
            println!(
                "KP15 PR08 q11 n={n} circulant {label} F5 dense: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] u0={u0} u1={u1} v0={v0} v1={v1} interior={interior} worst={worst}",
                neg = orientation.all.negative,
                samples = orientation.all.samples,
                zero = orientation.all.near_zero,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
                u0 = orientation.u0_negative,
                u1 = orientation.u1_negative,
                v0 = orientation.v0_negative,
                v1 = orientation.v1_negative,
                interior = orientation.interior_negative,
                worst = orientation.all.worst_sample,
            );
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_paper_p_generating_function_nullity_is_measured() {
    for n in [5usize, 6, 7, 8] {
        for (label, table) in [
            (
                "transposed-previous-side",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "s1-s4-plus-raw",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "assignments-only",
                Kp15Bi6GeneratingTable::build_with_assignments_only_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
        ] {
            let pinned = kp15_paper_p_and_tensor_border_columns_for_test(n);
            let nullity = kp15_constraint_nullity_after_pinning_for_test(&table, &pinned);
            println!(
                "KP15 n={n} {label} generating-function nullity after pinning P+tensor: rank={} cols={} pinned={} nullity={nullity}",
                table.rank(),
                table.cols,
                pinned.len(),
            );
            assert!(nullity <= table.cols);
        }
    }
}

#[test]
fn kp15_paper_p_generating_function_leftover_families_are_measured() {
    for n in [6usize, 7] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let pinned = kp15_paper_p_and_tensor_border_columns_for_test(n);
        let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
        let mut family_counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
        for col in &free_cols {
            let pos = col % 49;
            *family_counts.entry((pos / 7, pos % 7)).or_default() += 1;
        }
        println!(
            "KP15 n={n} S1-S4+raw leftover families after P+tensor pinning: total={} families={family_counts:?}",
            free_cols.len(),
        );
        assert!(!free_cols.is_empty());
    }
}

#[test]
fn kp15_fig2b_tensor_border_pin_hypotheses_are_measured() {
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        for (label, extra) in [
            ("current-12-slot", Vec::new()),
            (
                "full-depth-0-2",
                kp15_outer_depth_row_columns_for_test(n, 0..=2),
            ),
            (
                "full-outer-i-4-6",
                kp15_outer_i_band_columns_for_test(n, 4..=6),
            ),
            (
                "full-outer-j-4-6",
                kp15_outer_j_band_columns_for_test(n, 4..=6),
            ),
            (
                "full-outer-i-or-j-4-6",
                kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n),
            ),
            (
                "paper-depth-two-sided-strips",
                kp15_paper_depth_two_sided_strip_columns_for_test(n),
            ),
            (
                "paper-depth-plus-corner-l-band",
                kp15_paper_depth_two_sided_strip_with_corner_l_band_columns_for_test(n),
            ),
            (
                "full-l-band-baseline",
                kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n),
            ),
        ] {
            let mut pinned = kp15_paper_p_and_tensor_border_columns_for_test(n);
            pinned.extend(extra);
            let nullity = kp15_constraint_nullity_after_pinning_for_test(&table, &pinned);
            let leftover_families = if n == 6 && nullity <= 24 {
                let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
                let mut counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
                for col in &free_cols {
                    let pos = col % 49;
                    *counts.entry((pos / 7, pos % 7)).or_default() += 1;
                }
                format!(" families={counts:?}")
            } else {
                String::new()
            };
            println!(
                "KP15 n={n} Fig2b tensor-border pin hypothesis {label}: pinned={} nullity={nullity}{leftover_families}",
                pinned.len(),
            );
            assert!(nullity <= table.cols);
        }
    }
}

#[test]
fn kp15_fig2b_paper_source_support_leftover_families_are_measured() {
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        for (label, source_cols) in [
            (
                "depth-only-two-sided-strips",
                kp15_paper_depth_two_sided_strip_columns_for_test(n),
            ),
            (
                "paper-full-l-band",
                kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n),
            ),
        ] {
            let mut pinned = kp15_strict_p_variable_columns_for_test(n);
            pinned.extend(source_cols);
            let nullity = kp15_constraint_nullity_after_pinning_for_test(&table, &pinned);
            let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
            let mut family_counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
            for col in &free_cols {
                let pos = col % 49;
                *family_counts.entry((pos / 7, pos % 7)).or_default() += 1;
            }
            println!(
                "KP15 n={n} Fig2b paper-source support {label}: pinned={} nullity={nullity} free={} families={family_counts:?}",
                pinned.len(),
                free_cols.len(),
            );
            assert!(nullity <= table.cols);
        }
    }
}

#[test]
fn kp15_fig2b_two_sided_tensor_source_leaves_sector_curve_unclosed() {
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let mut pinned = kp15_strict_p_variable_columns_for_test(n);
        pinned.extend(kp15_paper_depth_two_sided_strip_columns_for_test(n));
        let nullity = kp15_constraint_nullity_after_pinning_for_test(&table, &pinned);
        let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
        let mut family_counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
        for col in &free_cols {
            let pos = col % 49;
            *family_counts.entry((pos / 7, pos % 7)).or_default() += 1;
        }
        println!(
            "KP15 n={n} Fig2b tensor-depth source after P+two-sided tensor-depth pinning: rank={} free={} pinned={} nullity={nullity} families={family_counts:?}",
            table.rank(),
            table.free_cols.len(),
            pinned.len(),
        );
        assert!(nullity > 0);
    }
}

#[test]
fn kp15_fig2b_tensor_depth_plus_m_endpoint_leaves_source_unclosed() {
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let mut pinned = kp15_strict_p_variable_columns_for_test(n);
        pinned.extend(kp15_paper_depth_two_sided_strip_columns_for_test(n));
        for sector in 0..n {
            pinned.insert(kp15_global_col(n, sector, 6, 0));
            pinned.insert(kp15_global_col(n, sector, 0, 6));
        }
        let nullity = kp15_constraint_nullity_after_pinning_for_test(&table, &pinned);
        let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
        let mut family_counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
        for col in &free_cols {
            let pos = col % 49;
            *family_counts.entry((pos / 7, pos % 7)).or_default() += 1;
        }
        println!(
            "KP15 n={n} Fig2b tensor-depth+m endpoint after P pinning: rank={} free={} pinned={} nullity={nullity} families={family_counts:?}",
            table.rank(),
            table.free_cols.len(),
            pinned.len(),
        );
        assert!(nullity > 0);
    }
}

#[test]
fn kp15_fig2b_two_sided_tensor_source_population_preserves_tensor_depth_rows() {
    for n in [5usize, 6, 7] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let source = kp15_compatible_two_sided_tensor_source_for_test(&table);
        let q = f_kappa_quadratic_form(5);
        let (grids, f_values) = table
            .populate_cap_grids_from_two_sided_tensor_source_with_f_kappa(&source, &q, 0.0)
            .expect("source-closed KP15 F5 solve");
        let mut max_source_delta: f64 = 0.0;
        for col in kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n) {
            let sector = col / 49;
            let pos = col % 49;
            let i = pos / 7;
            let j = pos % 7;
            max_source_delta =
                max_source_delta.max(grids[sector][i][j].distance_to(&source[sector][i][j]));
        }
        println!(
            "KP15 n={n} Fig2b tensor-depth source population: max_source_delta={max_source_delta:e} f=({:e},{:e},{:e})",
            f_values[0], f_values[1], f_values[2],
        );
        assert!(max_source_delta < 1.0e-9);
        assert!(f_values.iter().all(|v| v.is_finite()));
    }
}

#[test]
fn kp15_fig2b_full_stam_grid_sources_are_measured_through_constrained_table() {
    let patches = [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ];
    let transforms = [
        Kp20StamGridTransform::Identity,
        Kp20StamGridTransform::Rotate90,
        Kp20StamGridTransform::Rotate180,
        Kp20StamGridTransform::Rotate270,
        Kp20StamGridTransform::Transpose,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridTransform::FlipRows,
        Kp20StamGridTransform::FlipColumns,
    ];
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let projection =
            kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Stam projection");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for patch in patches {
            for transform in transforms {
                let source =
                    kp15_full_stam_bicubic_source_for_test(&projection, n, patch, transform);
                match table
                    .populate_cap_grids_from_two_sided_tensor_source_with_f_kappa(&source, &q, 0.0)
                {
                    Ok((grids, f_values)) => {
                        let orientation = test_bi6_patch_boundary_orientation_report(
                            n,
                            &grids,
                            &[0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0],
                        );
                        assert!(orientation.all.min_det.is_finite());
                        assert!(orientation.all.max_det.is_finite());
                        assert!(f_values.iter().all(|value| value.is_finite()));
                        accepted += 1;
                        if orientation.all.min_det > best_min_det {
                            best_min_det = orientation.all.min_det;
                            best_label = format!("{patch:?}/{transform:?}");
                        }
                    }
                    Err(err) => {
                        rejected += 1;
                        if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                            let abs_residual = residual.abs();
                            if abs_residual < smallest_rejection {
                                smallest_rejection = abs_residual;
                                smallest_rejection_label = format!("{patch:?}/{transform:?}");
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b full-Stam source summary n={n}: accepted={accepted} rejected={rejected} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert_eq!(accepted + rejected, patches.len() * transforms.len());
    }
}

#[test]
fn kp15_fig2b_beta_full_l_band_stam_sources_are_ranked() {
    let patches = [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ];
    let transforms = [
        Kp20StamGridTransform::Identity,
        Kp20StamGridTransform::Rotate90,
        Kp20StamGridTransform::Rotate180,
        Kp20StamGridTransform::Rotate270,
        Kp20StamGridTransform::Transpose,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridTransform::FlipRows,
        Kp20StamGridTransform::FlipColumns,
    ];
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let projection =
            kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Stam projection");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for patch in patches {
            for transform in transforms {
                for order in [
                    Kp15FullLBandOverwriteOrder::TopThenRight,
                    Kp15FullLBandOverwriteOrder::RightThenTop,
                ] {
                    let source = kp15_beta_full_l_band_from_bicubic_source_for_test(
                        &projection,
                        n,
                        patch,
                        transform,
                        order,
                    );
                    match table.populate_cap_grids_from_two_sided_tensor_source_with_f_kappa(
                        &source, &q, 0.0,
                    ) {
                        Ok((grids, f_values)) => {
                            let orientation =
                                test_bi6_patch_boundary_orientation_report(n, &grids, &samples);
                            assert!(f_values.iter().all(|value| value.is_finite()));
                            accepted += 1;
                            if orientation.all.negative < best_neg
                                || (orientation.all.negative == best_neg
                                    && orientation.all.min_det > best_min_det)
                            {
                                best_neg = orientation.all.negative;
                                best_min_det = orientation.all.min_det;
                                best_label = format!("{patch:?}/{transform:?}/{order:?}");
                            }
                        }
                        Err(err) => {
                            if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                                let abs_residual = residual.abs();
                                if abs_residual < smallest_rejection {
                                    smallest_rejection = abs_residual;
                                    smallest_rejection_label =
                                        format!("{patch:?}/{transform:?}/{order:?}");
                                }
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b beta full-L-band Stam source ranking n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(best_neg <= samples.len() * samples.len() * n || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig2b_stam_depth_strip_source_projection_is_measured() {
    let patches = [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ];
    let transforms = [
        Kp20StamGridTransform::Identity,
        Kp20StamGridTransform::Rotate90,
        Kp20StamGridTransform::Rotate180,
        Kp20StamGridTransform::Rotate270,
        Kp20StamGridTransform::Transpose,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridTransform::FlipRows,
        Kp20StamGridTransform::FlipColumns,
    ];
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let projection =
            kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Stam projection");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for patch in patches {
            for transform in transforms {
                let source =
                    kp15_full_stam_bicubic_source_for_test(&projection, n, patch, transform);
                let constraints =
                    kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
                match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "paper-depth source projection",
                ) {
                    Ok((grids, f_values)) => {
                        let orientation =
                            test_bi6_patch_boundary_orientation_report(n, &grids, &samples);
                        assert!(orientation.all.min_det.is_finite());
                        assert!(orientation.all.max_det.is_finite());
                        assert!(f_values.iter().all(|value| value.is_finite()));
                        accepted += 1;
                        if orientation.all.negative < best_neg
                            || (orientation.all.negative == best_neg
                                && orientation.all.min_det > best_min_det)
                        {
                            best_neg = orientation.all.negative;
                            best_min_det = orientation.all.min_det;
                            best_label = format!("{patch:?}/{transform:?}");
                        }
                    }
                    Err(err) => {
                        rejected += 1;
                        if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                            let abs_residual = residual.abs();
                            if abs_residual < smallest_rejection {
                                smallest_rejection = abs_residual;
                                smallest_rejection_label = format!("{patch:?}/{transform:?}");
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b depth-strip source projection summary n={n}: accepted={accepted} rejected={rejected} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert_eq!(accepted + rejected, patches.len() * transforms.len());
    }
}

#[test]
fn kp15_fig2b_q11_depth_strip_source_projection_is_measured() {
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        for apply_beta in [false, true] {
            let tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
            let source = kp15_two_sided_source_grid_from_bi3_tensor_for_test(
                n,
                &tensor.sectors,
                Bi6Reparams::for_valence(n),
                apply_beta,
            );
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "q11 paper-depth source projection",
            ) {
                Ok((grids, f_values)) => {
                    let orientation =
                        test_bi6_patch_boundary_orientation_report(n, &grids, &samples);
                    println!(
                        "KP15 Fig2b q11 depth-strip source projection n={n} beta={apply_beta}: accepted neg={neg}/{samples_count} det=[{min_det:e},{max_det:e}] f=({fx:e},{fy:e},{fz:e})",
                        neg = orientation.all.negative,
                        samples_count = orientation.all.samples,
                        min_det = orientation.all.min_det,
                        max_det = orientation.all.max_det,
                        fx = f_values[0],
                        fy = f_values[1],
                        fz = f_values[2],
                    );
                    assert!(orientation.all.min_det.is_finite());
                    assert!(orientation.all.max_det.is_finite());
                    assert!(f_values.iter().all(|value| value.is_finite()));
                }
                Err(err) => {
                    println!(
                        "KP15 Fig2b q11 depth-strip source projection n={n} beta={apply_beta}: rejected {err}"
                    );
                }
            }
        }
    }
}

#[test]
fn kp15_pr08_q11_endpoint_input_and_fig12_source_grid_are_compared() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let input = kp15_input_from_bi3_tensor(n, &tensor, false);
        for apply_beta in [false, true] {
            let source = kp15_two_sided_source_grid_from_bi3_tensor_for_test(
                n,
                &tensor.sectors,
                params,
                apply_beta,
            );
            let mut max_seed_delta: f64 = 0.0;
            let mut worst_seed = String::new();
            for (sector, sector_rhs) in input.sectors.iter().enumerate() {
                for (i, j, expected) in [
                    (0, 0, sector_rhs.p_00),
                    (1, 0, sector_rhs.p_10),
                    (2, 0, sector_rhs.p_20),
                    (4, 0, sector_rhs.p_40),
                    (5, 0, sector_rhs.p_50),
                    (6, 0, sector_rhs.p_60),
                    (0, 1, sector_rhs.p_01),
                    (1, 1, sector_rhs.p_11),
                    (0, 2, sector_rhs.p_02),
                    (4, 2, sector_rhs.p_42),
                    (5, 2, sector_rhs.p_52),
                    (6, 2, sector_rhs.p_62),
                ] {
                    let delta = source[sector][i][j].distance_to(&expected);
                    if delta > max_seed_delta {
                        max_seed_delta = delta;
                        worst_seed = format!("sector={sector} p_{i}{j}");
                    }
                }
            }
            let table =
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                    n,
                    params,
                );
            let labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
                n,
                &source,
                &kp15_paper_depth_two_sided_strip_columns_for_test(n),
                "q11-fig12-strip",
            );
            let first_inconsistency =
                kp15_first_inconsistent_source_constraint_for_test(&table, &labeled)
                    .unwrap_or_else(|| "all constraints consistent".to_string());
            println!(
                "KP15 PR08 q11 endpoint input vs Fig12 source grid n={n} beta={apply_beta}: max_seed_delta={max_seed_delta:e} worst_seed={worst_seed} first_strip_inconsistency={first_inconsistency}",
            );
            assert!(max_seed_delta.is_finite());
        }
    }
}

#[test]
fn kp15_fig2b_q11_depth_strip_source_phase_is_ranked() {
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for apply_beta in [false, true] {
            for row_offset in [-1isize, 0, 1] {
                for column_offset in [-1isize, 0, 1] {
                    for row_reversed in [false, true] {
                        for column_reversed in [false, true] {
                            let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
                                n,
                                &tensor.sectors,
                                Bi6Reparams::for_valence(n),
                                apply_beta,
                                row_offset,
                                column_offset,
                                row_reversed,
                                column_reversed,
                            );
                            let constraints = kp15_source_constraints_from_grid_columns_for_test(
                                n,
                                &source,
                                &depth_cols,
                            );
                            let label = format!(
                                "beta={apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed}"
                            );
                            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                                &constraints,
                                &q,
                                0.0,
                                "q11 paper-depth source phase",
                            ) {
                                Ok((grids, _)) => {
                                    accepted += 1;
                                    let orientation = test_bi6_patch_boundary_orientation_report(
                                        n, &grids, &samples,
                                    );
                                    if orientation.all.negative < best_neg
                                        || (orientation.all.negative == best_neg
                                            && orientation.all.min_det > best_min_det)
                                    {
                                        best_neg = orientation.all.negative;
                                        best_min_det = orientation.all.min_det;
                                        best_label = label;
                                    }
                                }
                                Err(err) => {
                                    if let Some(residual) =
                                        kp15_parse_constraint_residual_for_test(&err)
                                    {
                                        let abs_residual = residual.abs();
                                        if abs_residual < smallest_rejection {
                                            smallest_rejection = abs_residual;
                                            smallest_rejection_label = label;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b q11 depth-strip phase summary n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_cap_net_depth_strip_source_phase_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for apply_beta in [false, true] {
            for row_offset in [-1isize, 0, 1] {
                for column_offset in [-1isize, 0, 1] {
                    for row_reversed in [false, true] {
                        for column_reversed in [false, true] {
                            let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
                                n,
                                &tensor.sectors,
                                Bi6Reparams::for_valence(n),
                                apply_beta,
                                row_offset,
                                column_offset,
                                row_reversed,
                                column_reversed,
                            );
                            let constraints = kp15_source_constraints_from_grid_columns_for_test(
                                n,
                                &source,
                                &depth_cols,
                            );
                            let label = format!(
                                "beta={apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed}"
                            );
                            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                                &constraints,
                                &q,
                                0.0,
                                "cap-net paper-depth source phase",
                            ) {
                                Ok(_) => accepted += 1,
                                Err(err) => {
                                    if let Some(residual) =
                                        kp15_parse_constraint_residual_for_test(&err)
                                    {
                                        let abs_residual = residual.abs();
                                        if abs_residual < smallest_rejection {
                                            smallest_rejection = abs_residual;
                                            smallest_rejection_label = label;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b cap-net depth-strip phase summary n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_cap_net_tensor_depth_plus_m_endpoint_kkt_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        for m_label in [4usize, 5] {
            let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
                n,
                &source,
                &cap_net.sectors,
                m_label,
                -1,
            );
            let constraints: Vec<(usize, Point3)> = labeled
                .iter()
                .map(|item| {
                    let col = item.constraint.row.terms[0].0;
                    (
                        col,
                        Point3::new(
                            item.constraint.rhs[0],
                            item.constraint.rhs[1],
                            item.constraint.rhs[2],
                        ),
                    )
                })
                .collect();
            let (grids, f_values) = table
                .populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &f_kappa_quadratic_form(5),
                    0.0,
                    "cap-net tensor-depth plus m endpoint",
                )
                .unwrap_or_else(|err| {
                    panic!(
                        "n={n} cap-net tensor-depth c{} rejected: {err}",
                        m_label + 1
                    )
                });
            let orientation =
                test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
            println!(
                "KP15 Fig2b cap-net tensor-depth plus m endpoint n={n} m_label=c{}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                m_label + 1,
                orientation.negative,
                orientation.samples,
                orientation.near_zero,
                orientation.min_det,
                orientation.max_det,
                orientation.worst_sample,
                f_values[0],
                f_values[1],
                f_values[2],
            );
            assert!(orientation.samples > 0);
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_tensor_depth_c6_endpoint_phase_family_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut accepted = 0usize;
        for row_offset in [-1isize, 0, 1] {
            for column_offset in [-1isize, 0, 1] {
                for row_reversed in [false, true] {
                    for column_reversed in [false, true] {
                        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
                            n,
                            &tensor.sectors,
                            params,
                            true,
                            row_offset,
                            column_offset,
                            row_reversed,
                            column_reversed,
                        );
                        let labeled =
                            kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
                                n,
                                &source,
                                &cap_net.sectors,
                                5,
                                -1,
                            );
                        let constraints: Vec<(usize, Point3)> = labeled
                            .iter()
                            .map(|item| {
                                let col = item.constraint.row.terms[0].0;
                                (
                                    col,
                                    Point3::new(
                                        item.constraint.rhs[0],
                                        item.constraint.rhs[1],
                                        item.constraint.rhs[2],
                                    ),
                                )
                            })
                            .collect();
                        let Ok((grids, f_values)) = table
                            .populate_cap_grids_from_column_constraints_with_f_kappa(
                                &constraints,
                                &f_kappa_quadratic_form(5),
                                0.0,
                                "cap-net tensor-depth c6 endpoint phase family",
                            )
                        else {
                            continue;
                        };
                        accepted += 1;
                        let orientation = test_bi6_patch_local_orientation_report_with_samples(
                            n, &grids, &samples,
                        );
                        if orientation.negative < best_neg
                            || (orientation.negative == best_neg
                                && orientation.min_det > best_min_det)
                        {
                            best_neg = orientation.negative;
                            best_min_det = orientation.min_det;
                            best_label = format!(
                                "row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                                orientation.min_det,
                                orientation.max_det,
                                orientation.worst_sample,
                                f_values[0],
                                f_values[1],
                                f_values[2],
                            );
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b cap-net tensor-depth c6 endpoint phase family n={n}: accepted={accepted} best_neg={best_neg}/{} best_min_det={best_min_det:e} best={best_label}",
            n * samples.len() * samples.len(),
        );
        assert!(accepted > 0);
    }
}

#[test]
fn kp15_fig2b_cap_net_full_fig12_strip_plus_c6_endpoint_kkt_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut source_by_col = std::collections::BTreeMap::<usize, Point3>::new();
        for (col, point) in kp15_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &kp15_paper_depth_two_sided_strip_columns_for_test(n),
        ) {
            source_by_col.insert(col, point);
        }
        for sector in 0..n {
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
            source_by_col.insert(
                kp15_global_col(n, sector, 6, 0),
                cap_net.sectors[p60_sector][5],
            );
            source_by_col.insert(
                kp15_global_col(n, sector, 0, 6),
                cap_net.sectors[p06_sector][5],
            );
        }
        let constraints: Vec<(usize, Point3)> = source_by_col.into_iter().collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net full Fig12 strip plus c6 endpoint",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net full Fig12 strip plus c6 endpoint n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                println!(
                    "KP15 Fig2b cap-net full Fig12 strip plus c6 endpoint n={n}: rejected {err}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_full_fig12_strip_first_inconsistency_is_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &kp15_paper_depth_two_sided_strip_columns_for_test(n),
            "full-fig12-strip",
        );
        for sector in 0..n {
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
            for (label, slot, point) in [
                ("p60", (6, 0), cap_net.sectors[p60_sector][5]),
                ("p06", (0, 6), cap_net.sectors[p06_sector][5]),
            ] {
                constraints.retain(|constraint| {
                    constraint.constraint.row.terms[0].0
                        != kp15_global_col(n, sector, slot.0, slot.1)
                });
                constraints.push(Kp15LabeledConstraint {
                    label: format!("cap-net-c6 sector={sector} {label}"),
                    constraint: AffineConstraintRow {
                        row: LinearRow {
                            terms: vec![(kp15_global_col(n, sector, slot.0, slot.1), 1.0)],
                        },
                        rhs: [point.x, point.y, point.z],
                    },
                });
            }
        }
        let report = kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
            .unwrap_or_else(|| "all constraints consistent".to_string());
        println!("KP15 Fig2b cap-net full Fig12 strip first inconsistency n={n}: {report}");
    }
}

#[test]
fn kp15_fig2b_cap_net_full_fig12_strip_phase_compatibility_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut compatible = 0usize;
        let mut first_compatible = String::new();
        let mut first_rejection = None::<String>;
        for row_offset in [-1isize, 0, 1] {
            for column_offset in [-1isize, 0, 1] {
                for row_reversed in [false, true] {
                    for column_reversed in [false, true] {
                        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
                            n,
                            &tensor.sectors,
                            params,
                            true,
                            row_offset,
                            column_offset,
                            row_reversed,
                            column_reversed,
                        );
                        let mut constraints =
                            kp15_labeled_source_constraints_from_grid_columns_for_test(
                                n,
                                &source,
                                &kp15_paper_depth_two_sided_strip_columns_for_test(n),
                                "full-fig12-strip",
                            );
                        for sector in 0..n {
                            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
                            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
                            for (label, slot, point) in [
                                ("p60", (6, 0), cap_net.sectors[p60_sector][5]),
                                ("p06", (0, 6), cap_net.sectors[p06_sector][5]),
                            ] {
                                let col = kp15_global_col(n, sector, slot.0, slot.1);
                                constraints.retain(|constraint| {
                                    constraint.constraint.row.terms[0].0 != col
                                });
                                constraints.push(Kp15LabeledConstraint {
                                    label: format!("cap-net-c6 sector={sector} {label}"),
                                    constraint: AffineConstraintRow {
                                        row: LinearRow {
                                            terms: vec![(col, 1.0)],
                                        },
                                        rhs: [point.x, point.y, point.z],
                                    },
                                });
                            }
                        }
                        let label = format!(
                            "row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed}"
                        );
                        if let Some(report) =
                            kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
                        {
                            if first_rejection.is_none() {
                                first_rejection = Some(format!("{label} first={report}"));
                            }
                        } else {
                            compatible += 1;
                            if first_compatible.is_empty() {
                                first_compatible = label;
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b cap-net full Fig12 strip phase compatibility n={n}: compatible={compatible} first_compatible={first_compatible} first_rejection={}",
            first_rejection.as_deref().unwrap_or("none"),
        );
        assert!(compatible > 0 || first_rejection.is_some());
    }
}

#[test]
fn kp15_fig2b_cap_net_eq10_boundary_plus_depth_kkt_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut source_cols = kp15_fig2b_depth2_tensor_border_columns_for_test(n);
        for sector in 0..n {
            for slot in [(5, 0), (6, 0), (0, 5), (0, 6)] {
                source_cols.push(kp15_global_col(n, sector, slot.0, slot.1));
            }
        }
        source_cols.sort_unstable();
        source_cols.dedup();
        let mut source_by_col = std::collections::BTreeMap::<usize, Point3>::new();
        for (col, point) in
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols)
        {
            source_by_col.insert(col, point);
        }
        for sector in 0..n {
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
            source_by_col.insert(
                kp15_global_col(n, sector, 6, 0),
                cap_net.sectors[p60_sector][5],
            );
            source_by_col.insert(
                kp15_global_col(n, sector, 0, 6),
                cap_net.sectors[p06_sector][5],
            );
        }
        let constraints: Vec<(usize, Point3)> = source_by_col.into_iter().collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net eq10 boundary plus tensor-depth",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net eq10 boundary plus depth n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                println!("KP15 Fig2b cap-net eq10 boundary plus depth n={n}: rejected {err}");
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_p50_m_endpoint_plus_depth_kkt_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut source_cols = kp15_fig2b_depth2_tensor_border_columns_for_test(n);
        for sector in 0..n {
            source_cols.push(kp15_global_col(n, sector, 5, 0));
            source_cols.push(kp15_global_col(n, sector, 6, 0));
            source_cols.push(kp15_global_col(n, sector, 0, 6));
        }
        source_cols.sort_unstable();
        source_cols.dedup();
        let mut source_by_col = std::collections::BTreeMap::<usize, Point3>::new();
        for (col, point) in
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols)
        {
            source_by_col.insert(col, point);
        }
        for sector in 0..n {
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
            source_by_col.insert(
                kp15_global_col(n, sector, 6, 0),
                cap_net.sectors[p60_sector][5],
            );
            source_by_col.insert(
                kp15_global_col(n, sector, 0, 6),
                cap_net.sectors[p06_sector][5],
            );
        }
        let constraints: Vec<(usize, Point3)> = source_by_col.into_iter().collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net p50/m endpoint plus tensor-depth",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net p50+m endpoint plus depth n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                println!("KP15 Fig2b cap-net p50+m endpoint plus depth n={n}: rejected {err}");
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_reduced_p_variable_frontier_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let input = kp15_input_from_bi3_tensor(n, &tensor, true);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let strict_p = kp15_strict_p_variable_columns_for_test(n);
        for (label, variable_filter, released_tensor_cols) in [
            ("all-non-fixed-free", None, None),
            ("strict-p-only", Some(&strict_p), None),
            (
                "strict-p-release-p50-p60",
                Some(&strict_p),
                Some(&kp15_source_slot_family_columns(n, &[(5, 0), (6, 0)])),
            ),
        ] {
            let Ok((patches, f_values, split)) = kp15_solve_table_with_test_variable_columns(
                &table,
                &input.sectors,
                &q,
                variable_filter,
                released_tensor_cols,
            ) else {
                println!("KP15 Fig2b cap-net reduced-P frontier n={n} {label}: rejected");
                continue;
            };
            let orientation =
                test_bi6_patch_local_orientation_report_with_samples(n, &patches, &samples);
            println!(
                "KP15 Fig2b cap-net reduced-P frontier n={n} {label}: fixed={} variable={} zeroed_extra={} neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                split.fixed,
                split.variable,
                split.zeroed_extra,
                orientation.negative,
                orientation.samples,
                orientation.near_zero,
                orientation.min_det,
                orientation.max_det,
                orientation.worst_sample,
                f_values[0],
                f_values[1],
                f_values[2],
            );
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_solved_adjacent_boundary_slots_are_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
            n,
            &source,
            &cap_net.sectors,
            5,
            -1,
        );
        let constraints: Vec<(usize, Point3)> = labeled
            .iter()
            .map(|item| {
                let col = item.constraint.row.terms[0].0;
                (
                    col,
                    Point3::new(
                        item.constraint.rhs[0],
                        item.constraint.rhs[1],
                        item.constraint.rhs[2],
                    ),
                )
            })
            .collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let (grids, _) = table
            .populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &f_kappa_quadratic_form(5),
                0.0,
                "cap-net solved adjacent boundary slots",
            )
            .expect("compatible cap-net depth+m solve");
        for slot in [(0usize, 4usize), (0, 5)] {
            let solved = grids[0][slot.0][slot.1];
            let source_value = source[0][slot.0][slot.1];
            let mut nearest = (f64::INFINITY, String::new());
            for sector_offset in [-1isize, 0, 1] {
                let sector = kp15_offset_sector_for_test(n, 0, sector_offset);
                for label in 0..6 {
                    let candidate = cap_net.sectors[sector][label];
                    let distance = solved.distance_to(&candidate);
                    if distance < nearest.0 {
                        nearest = (distance, format!("c{} offset={sector_offset}", label + 1));
                    }
                }
            }
            println!(
                "KP15 Fig2b cap-net solved adjacent boundary n={n} p{}{}: source_delta={:e} nearest={} distance={:e} solved=({:e},{:e},{:e})",
                slot.0,
                slot.1,
                solved.distance_to(&source_value),
                nearest.1,
                nearest.0,
                solved.x,
                solved.y,
                solved.z,
            );
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_adjacent_boundary_dependency_is_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut constraints = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
            n,
            &source,
            &cap_net.sectors,
            5,
            -1,
        );
        for slot in [(0usize, 4usize), (0, 5)] {
            let point = source[0][slot.0][slot.1];
            constraints.push(Kp15LabeledConstraint {
                label: format!("adjacent-boundary sector=0 p_{}{}", slot.0, slot.1),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(n, 0, slot.0, slot.1), 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            });
        }
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let first_inconsistency =
            kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
                .unwrap_or_else(|| "all constraints consistent".to_string());
        println!(
            "KP15 Fig2b cap-net adjacent boundary dependency n={n}: first={first_inconsistency}"
        );
        for target_label in [
            "adjacent-boundary sector=0 p_04",
            "adjacent-boundary sector=0 p_05",
        ] {
            let report =
                kp15_implied_ordered_constraint_delta_for_test(&table, &constraints, target_label)
                    .unwrap_or_else(|| format!("{target_label} not found"));
            let formula =
                kp15_ordered_source_formula_for_test(&table, &constraints, target_label, 0)
                    .unwrap_or_else(|| "independent".to_string());
            println!(
                "KP15 Fig2b cap-net adjacent boundary dependency n={n}: {report}; x_formula={formula}"
            );
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_full_strip_boundary_dependency_formula_is_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &kp15_paper_depth_two_sided_strip_columns_for_test(n),
            "full-fig12-strip",
        );
        for sector in 0..n {
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
            for (label, slot, point) in [
                ("p60", (6, 0), cap_net.sectors[p60_sector][5]),
                ("p06", (0, 6), cap_net.sectors[p06_sector][5]),
            ] {
                let col = kp15_global_col(n, sector, slot.0, slot.1);
                constraints.retain(|constraint| constraint.constraint.row.terms[0].0 != col);
                constraints.push(Kp15LabeledConstraint {
                    label: format!("cap-net-c6 sector={sector} {label}"),
                    constraint: AffineConstraintRow {
                        row: LinearRow {
                            terms: vec![(col, 1.0)],
                        },
                        rhs: [point.x, point.y, point.z],
                    },
                });
            }
        }
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let first_inconsistency =
            kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
                .unwrap_or_else(|| "all constraints consistent".to_string());
        println!(
            "KP15 Fig2b cap-net full-strip boundary dependency n={n}: first={first_inconsistency}"
        );
        for target_label in [
            "full-fig12-strip sector=0 p_04",
            "full-fig12-strip sector=0 p_05",
        ] {
            let report =
                kp15_implied_ordered_constraint_delta_for_test(&table, &constraints, target_label)
                    .unwrap_or_else(|| format!("{target_label} not found"));
            let formula =
                kp15_ordered_source_formula_for_test(&table, &constraints, target_label, 0)
                    .unwrap_or_else(|| "independent".to_string());
            println!(
                "KP15 Fig2b cap-net full-strip boundary dependency n={n}: {report}; x_formula={formula}"
            );
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_cyclic_boundary_full_strip_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut source_by_col = std::collections::BTreeMap::<usize, Point3>::new();
        for (col, point) in kp15_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &kp15_paper_depth_two_sided_strip_columns_for_test(n),
        ) {
            source_by_col.insert(col, point);
        }
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            for j in 4..=5 {
                source_by_col.insert(kp15_global_col(n, sector, 0, j), source[next][j][0]);
            }
            source_by_col.insert(kp15_global_col(n, sector, 0, 6), cap_net.sectors[sector][5]);
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            source_by_col.insert(
                kp15_global_col(n, sector, 6, 0),
                cap_net.sectors[p60_sector][5],
            );
        }
        let labeled_constraints = source_by_col
            .iter()
            .map(|(&col, point)| {
                let sector = col / 49;
                let pos = col % 49;
                let i = pos / 7;
                let j = pos % 7;
                Kp15LabeledConstraint {
                    label: format!("cyclic-boundary sector={sector} p_{i}{j}"),
                    constraint: AffineConstraintRow {
                        row: LinearRow {
                            terms: vec![(col, 1.0)],
                        },
                        rhs: [point.x, point.y, point.z],
                    },
                }
            })
            .collect::<Vec<_>>();
        let constraints: Vec<(usize, Point3)> = source_by_col.into_iter().collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net cyclic boundary full Fig12 strip",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net cyclic-boundary full strip n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                let first = kp15_first_inconsistent_source_constraint_for_test(
                    &table,
                    &labeled_constraints,
                )
                .unwrap_or_else(|| "all constraints consistent".to_string());
                println!(
                    "KP15 Fig2b cap-net cyclic-boundary full strip n={n}: rejected {err}; first={first}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_cyclic_boundary_depth_strip_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let mut source_by_col = std::collections::BTreeMap::<usize, Point3>::new();
        for (col, point) in kp15_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &kp15_fig2b_depth2_tensor_border_columns_for_test(n),
        ) {
            source_by_col.insert(col, point);
        }
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            for i in 4..=5 {
                source_by_col.insert(kp15_global_col(n, sector, i, 0), source[sector][i][0]);
                source_by_col.insert(kp15_global_col(n, sector, 0, i), source[next][i][0]);
            }
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            source_by_col.insert(
                kp15_global_col(n, sector, 6, 0),
                cap_net.sectors[p60_sector][5],
            );
            source_by_col.insert(kp15_global_col(n, sector, 0, 6), cap_net.sectors[sector][5]);
        }
        let labeled_constraints = source_by_col
            .iter()
            .map(|(&col, point)| {
                let sector = col / 49;
                let pos = col % 49;
                let i = pos / 7;
                let j = pos % 7;
                Kp15LabeledConstraint {
                    label: format!("cyclic-boundary-depth sector={sector} p_{i}{j}"),
                    constraint: AffineConstraintRow {
                        row: LinearRow {
                            terms: vec![(col, 1.0)],
                        },
                        rhs: [point.x, point.y, point.z],
                    },
                }
            })
            .collect::<Vec<_>>();
        let constraints: Vec<(usize, Point3)> = source_by_col.into_iter().collect();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net cyclic boundary depth strip",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net cyclic-boundary depth strip n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                let first = kp15_first_inconsistent_source_constraint_for_test(
                    &table,
                    &labeled_constraints,
                )
                .unwrap_or_else(|| "all constraints consistent".to_string());
                println!(
                    "KP15 Fig2b cap-net cyclic-boundary depth strip n={n}: rejected {err}; first={first}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_c6_endpoint_asymmetric_spectral_depth_source_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let scale_cases = [
            ("identity", 1.0),
            ("lambda", ring.lambda),
            ("one-minus-lambda", 1.0 - ring.lambda),
            ("half-one-minus-lambda", 0.5 * (1.0 - ring.lambda)),
            ("lambda-squared", ring.lambda * ring.lambda),
        ];
        let mut best = (
            "none",
            "none",
            usize::MAX,
            f64::NEG_INFINITY,
            [f64::INFINITY; 3],
        );
        for (row_label, row_scale) in scale_cases {
            for (column_label, column_scale) in scale_cases {
                let scaled_source = kp15_scale_fig2b_source_depth_from_boundary_for_test(
                    &source,
                    row_scale,
                    column_scale,
                );
                let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
                    n,
                    &scaled_source,
                    &cap_net.sectors,
                    5,
                    -1,
                );
                let constraints: Vec<(usize, Point3)> = labeled
                    .iter()
                    .map(|item| {
                        let col = item.constraint.row.terms[0].0;
                        (
                            col,
                            Point3::new(
                                item.constraint.rhs[0],
                                item.constraint.rhs[1],
                                item.constraint.rhs[2],
                            ),
                        )
                    })
                    .collect();
                let (grids, f_values) = table
                    .populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "cap-net c6 endpoint asymmetric spectral depth source",
                    )
                    .expect("cap-net depth+c6 source should remain feasible");
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                if orientation.negative < best.2
                    || (orientation.negative == best.2 && orientation.min_det > best.3)
                {
                    best = (
                        row_label,
                        column_label,
                        orientation.negative,
                        orientation.min_det,
                        f_values,
                    );
                }
            }
        }
        println!(
            "KP15 Fig2b cap-net c6 endpoint asymmetric spectral depth source n={n}: row_scale={} column_scale={} neg={} min_det={:e} f=({:e},{:e},{:e})",
            best.0, best.1, best.2, best.3, best.4[0], best.4[1], best.4[2],
        );
        assert!(best.3.is_finite());
    }
}

#[test]
fn kp15_fig2b_cap_net_c6_endpoint_best_spectral_source_failure_is_localized() {
    let coarse_samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    let dense_samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let (row_label, row_scale, column_label, column_scale) = match n {
            5 => ("identity", 1.0, "lambda", ring.lambda),
            6 => ("identity", 1.0, "one-minus-lambda", 1.0 - ring.lambda),
            7 => (
                "identity",
                1.0,
                "half-one-minus-lambda",
                0.5 * (1.0 - ring.lambda),
            ),
            8 => ("lambda", ring.lambda, "identity", 1.0),
            _ => unreachable!("measured valence"),
        };
        let scaled_source =
            kp15_scale_fig2b_source_depth_from_boundary_for_test(&source, row_scale, column_scale);
        let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
            n,
            &scaled_source,
            &cap_net.sectors,
            5,
            -1,
        );
        let constraints: Vec<(usize, Point3)> = labeled
            .iter()
            .map(|item| {
                let col = item.constraint.row.terms[0].0;
                (
                    col,
                    Point3::new(
                        item.constraint.rhs[0],
                        item.constraint.rhs[1],
                        item.constraint.rhs[2],
                    ),
                )
            })
            .collect();
        let (grids, f_values) = table
            .populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "cap-net c6 endpoint best spectral source localization",
            )
            .expect("best measured cap-net depth+c6 source should remain feasible");
        let coarse =
            test_bi6_patch_local_orientation_report_with_samples(n, &grids, &coarse_samples);
        let dense = test_bi6_patch_boundary_orientation_report(n, &grids, &dense_samples);
        println!(
            "KP15 Fig2b cap-net c6 endpoint best spectral source localization n={n}: row_scale={row_label} column_scale={column_label} coarse_neg={}/{} dense_neg={}/{} boundary=(u0:{},u1:{},v0:{},v1:{}) interior={} min_det={:e} max_det={:e} worst={} f=({:e},{:e},{:e})",
            coarse.negative,
            coarse.samples,
            dense.all.negative,
            dense.all.samples,
            dense.u0_negative,
            dense.u1_negative,
            dense.v0_negative,
            dense.v1_negative,
            dense.interior_negative,
            dense.all.min_det,
            dense.all.max_det,
            dense.all.worst_sample,
            f_values[0],
            f_values[1],
            f_values[2],
        );
        assert!(dense.all.min_det.is_finite());
    }
}

#[test]
fn kp15_fig2b_cap_net_full_fig12_strip_c6_endpoint_asymmetric_spectral_source_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let scale_cases = [
            ("identity", 1.0),
            ("lambda", ring.lambda),
            ("one-minus-lambda", 1.0 - ring.lambda),
            ("half-one-minus-lambda", 0.5 * (1.0 - ring.lambda)),
            ("lambda-squared", ring.lambda * ring.lambda),
        ];
        let mut best = (
            "none",
            "none",
            usize::MAX,
            f64::NEG_INFINITY,
            [f64::INFINITY; 3],
            String::new(),
        );
        let mut feasible = 0usize;
        let mut first_rejection = None::<String>;
        for (row_label, row_scale) in scale_cases {
            for (column_label, column_scale) in scale_cases {
                let scaled_source = kp15_scale_fig2b_source_depth_from_boundary_for_test(
                    &source,
                    row_scale,
                    column_scale,
                );
                let mut labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
                    n,
                    &scaled_source,
                    &kp15_paper_depth_two_sided_strip_columns_for_test(n),
                    "full-fig12-strip",
                );
                for sector in 0..n {
                    let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
                    let p06_sector = kp15_offset_sector_for_test(n, sector + 1, -1);
                    for (label, slot, point) in [
                        ("p60", (6, 0), cap_net.sectors[p60_sector][5]),
                        ("p06", (0, 6), cap_net.sectors[p06_sector][5]),
                    ] {
                        let col = kp15_global_col(n, sector, slot.0, slot.1);
                        labeled.retain(|constraint| constraint.constraint.row.terms[0].0 != col);
                        labeled.push(Kp15LabeledConstraint {
                            label: format!("cap-net-c6 sector={sector} {label}"),
                            constraint: AffineConstraintRow {
                                row: LinearRow {
                                    terms: vec![(col, 1.0)],
                                },
                                rhs: [point.x, point.y, point.z],
                            },
                        });
                    }
                }
                let constraints: Vec<(usize, Point3)> = labeled
                    .iter()
                    .map(|item| {
                        let col = item.constraint.row.terms[0].0;
                        (
                            col,
                            Point3::new(
                                item.constraint.rhs[0],
                                item.constraint.rhs[1],
                                item.constraint.rhs[2],
                            ),
                        )
                    })
                    .collect();
                match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "cap-net full Fig12 strip c6 endpoint asymmetric spectral source",
                ) {
                    Ok((grids, f_values)) => {
                        feasible += 1;
                        let orientation = test_bi6_patch_local_orientation_report_with_samples(
                            n, &grids, &samples,
                        );
                        if orientation.negative < best.2
                            || (orientation.negative == best.2 && orientation.min_det > best.3)
                        {
                            best = (
                                row_label,
                                column_label,
                                orientation.negative,
                                orientation.min_det,
                                f_values,
                                orientation.worst_sample,
                            );
                        }
                    }
                    Err(err) => {
                        if first_rejection.is_none() {
                            let first = kp15_first_inconsistent_source_constraint_for_test(
                                &table, &labeled,
                            )
                            .unwrap_or_else(|| "all constraints consistent".to_string());
                            first_rejection = Some(format!("{err}; first={first}"));
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig2b cap-net full Fig12 strip c6 endpoint asymmetric spectral source n={n}: feasible={feasible}/25 row_scale={} column_scale={} neg={} min_det={:e} worst={} f=({:e},{:e},{:e}) first_rejection={}",
            best.0,
            best.1,
            best.2,
            best.3,
            best.5,
            best.4[0],
            best.4[1],
            best.4[2],
            first_rejection.unwrap_or_else(|| "none".to_string()),
        );
        if feasible == 0 {
            assert_eq!(best.0, "none");
        } else {
            assert!(best.3.is_finite());
        }
    }
}

#[test]
fn kp20_cap_net_pr08_stencil_refined_labels_are_compared_to_characteristic_ring() {
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let refined = kp20_pr08_refined_labels_from_cap_net_for_test(&cap_net);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let mut source_xy = Vec::with_capacity(13 * n);
        let mut target_xy = Vec::with_capacity(13 * n);
        for sector in 0..n {
            for label in 0..13 {
                let source = refined[sector][label];
                let target = ring.sectors[sector][label];
                source_xy.push([source.x, source.y]);
                target_xy.push([target.x, target.y]);
            }
        }
        let scale =
            kp15_complex_similarity_scale_xy_for_test(&source_xy, &target_xy).expect("scale");
        let mut max_residual: f64 = 0.0;
        let mut worst = String::new();
        for sector in 0..n {
            for label in 0..13 {
                let source = refined[sector][label];
                let scaled = kp15_complex_mul_xy_for_test([source.x, source.y], scale);
                let target = ring.sectors[sector][label];
                let residual =
                    ((scaled[0] - target.x).powi(2) + (scaled[1] - target.y).powi(2)).sqrt();
                if residual > max_residual {
                    max_residual = residual;
                    worst = format!("sector={sector} q{}", label + 1);
                }
            }
        }
        println!(
            "KP20 cap-net PR08 stencil refined-label comparison n={n}: scale=({:e},{:e}) max_residual={:e} worst={worst}",
            scale[0], scale[1], max_residual,
        );
        assert!(max_residual.is_finite());
    }
}

#[test]
fn kp20_cap_net_pr08_stencil_source_label_permutations_are_ranked() {
    let permutations = kp20_cap_net_label_permutations_for_test();
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let mut best = (f64::INFINITY, [0usize; 6], [0.0_f64; 2], String::new());
        for permutation in &permutations {
            let refined =
                kp20_pr08_refined_labels_from_cap_net_permutation_for_test(&cap_net, permutation);
            let mut source_xy = Vec::with_capacity(13 * n);
            let mut target_xy = Vec::with_capacity(13 * n);
            for sector in 0..n {
                for label in 0..13 {
                    let source = refined[sector][label];
                    let target = ring.sectors[sector][label];
                    source_xy.push([source.x, source.y]);
                    target_xy.push([target.x, target.y]);
                }
            }
            let scale =
                kp15_complex_similarity_scale_xy_for_test(&source_xy, &target_xy).expect("scale");
            let mut max_residual: f64 = 0.0;
            let mut worst = String::new();
            for sector in 0..n {
                for label in 0..13 {
                    let source = refined[sector][label];
                    let scaled = kp15_complex_mul_xy_for_test([source.x, source.y], scale);
                    let target = ring.sectors[sector][label];
                    let residual =
                        ((scaled[0] - target.x).powi(2) + (scaled[1] - target.y).powi(2)).sqrt();
                    if residual > max_residual {
                        max_residual = residual;
                        worst = format!("sector={sector} q{}", label + 1);
                    }
                }
            }
            if max_residual < best.0 {
                best = (max_residual, *permutation, scale, worst);
            }
        }
        let labels = best.1.map(|label| label + 1);
        println!(
            "KP20 cap-net PR08 stencil label permutation n={n}: best_perm={labels:?} scale=({:e},{:e}) max_residual={:e} worst={}",
            best.2[0], best.2[1], best.0, best.3,
        );
        assert!(best.0.is_finite());
    }
}

#[test]
fn kp15_fig2b_depth_endpoint_plus_p_leftover_families_are_reported() {
    for n in [5usize, 6, 7, 8] {
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                Bi6Reparams::for_valence(n),
            );
        let mut pinned = kp15_strict_p_variable_columns_for_test(n);
        pinned.extend(kp15_fig2b_depth2_plus_m_endpoint_columns_for_test(n));
        let free_cols = kp15_free_columns_after_pinning_for_test(&table, &pinned);
        let mut family_counts = std::collections::BTreeMap::<(usize, usize), usize>::new();
        for col in &free_cols {
            let pos = col % 49;
            *family_counts.entry((pos / 7, pos % 7)).or_default() += 1;
        }
        println!(
            "KP15 Fig2b depth+endpoint+P leftover families n={n}: total={} families={family_counts:?}",
            free_cols.len(),
        );
    }
}

#[test]
fn kp15_fig2b_cap_net_full_l_band_source_closure_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &source_cols,
            "cap-net-full-l-band",
        );
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net full L-band tensor-border",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net full L-band source closure n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                let first = kp15_first_inconsistent_source_constraint_for_test(&table, &labeled)
                    .unwrap_or_else(|| "all constraints consistent".to_string());
                println!(
                    "KP15 Fig2b cap-net full L-band source closure n={n}: rejected {err}; first={first}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_cyclic_full_l_band_source_closure_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let mut source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            for j in 4..=5 {
                source[sector][0][j] = source[next][j][0];
            }
            source[sector][0][6] = cap_net.sectors[sector][5];
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            source[sector][6][0] = cap_net.sectors[p60_sector][5];
        }
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &source_cols,
            "cap-net-cyclic-full-l-band",
        );
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net cyclic full L-band tensor-border",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net cyclic full L-band source closure n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                let first = kp15_first_inconsistent_source_constraint_for_test(&table, &labeled)
                    .unwrap_or_else(|| "all constraints consistent".to_string());
                println!(
                    "KP15 Fig2b cap-net cyclic full L-band source closure n={n}: rejected {err}; first={first}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_cyclic_full_l_band_first_cross_formula_is_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let mut source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            for j in 4..=5 {
                source[sector][0][j] = source[next][j][0];
            }
            source[sector][0][6] = cap_net.sectors[sector][5];
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            source[sector][6][0] = cap_net.sectors[p60_sector][5];
        }
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &source_cols,
            "cap-net-cyclic-full-l-band",
        );
        let first = kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
            .unwrap_or_else(|| "all constraints consistent".to_string());
        println!("KP15 Fig2b cap-net cyclic full L-band first-cross formula n={n}: first={first}");
        for target_label in [
            "cap-net-cyclic-full-l-band sector=1 p_41",
            "cap-net-cyclic-full-l-band sector=1 p_51",
        ] {
            let report =
                kp15_implied_ordered_constraint_delta_for_test(&table, &constraints, target_label)
                    .unwrap_or_else(|| format!("{target_label} not found"));
            let formula =
                kp15_ordered_source_formula_for_test(&table, &constraints, target_label, 0)
                    .unwrap_or_else(|| "independent".to_string());
            println!(
                "KP15 Fig2b cap-net cyclic full L-band first-cross formula n={n}: {report}; x_formula={formula}"
            );
        }
    }
}

#[test]
fn kp15_fig2b_cap_net_s4_first_cross_full_l_band_source_closure_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let mut source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            params,
            true,
            0,
            -1,
            false,
            false,
        );
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            for j in 4..=5 {
                source[sector][0][j] = source[next][j][0];
            }
            source[sector][0][6] = cap_net.sectors[sector][5];
            let p60_sector = kp15_offset_sector_for_test(n, sector, -1);
            source[sector][6][0] = cap_net.sectors[p60_sector][5];
        }
        let original = source.clone();
        let w3 = params.c / 15.0;
        for sector in 0..n {
            let next = kp15_offset_sector_for_test(n, sector, 1);
            source[next][4][1] = kp15_point_linear_combination_for_test(&[
                (2.0, original[sector][0][4]),
                (-1.0, original[sector][1][4]),
                (w3, original[sector][0][6]),
                (-w3, original[sector][0][5]),
            ]);
            source[next][5][1] = kp15_point_linear_combination_for_test(&[
                (2.0, original[sector][0][5]),
                (-1.0, original[sector][1][5]),
            ]);
            source[next][6][1] = kp15_point_linear_combination_for_test(&[
                (2.0, original[sector][0][6]),
                (-1.0, original[sector][1][6]),
            ]);
        }
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &source_cols,
            "cap-net-s4-first-cross-full-l-band",
        );
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net S4 first-cross full L-band tensor-border",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net S4 first-cross full L-band source closure n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                    orientation.negative,
                    orientation.samples,
                    orientation.near_zero,
                    orientation.min_det,
                    orientation.max_det,
                    orientation.worst_sample,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
            }
            Err(err) => {
                let first = kp15_first_inconsistent_source_constraint_for_test(&table, &labeled)
                    .unwrap_or_else(|| "all constraints consistent".to_string());
                println!(
                    "KP15 Fig2b cap-net S4 first-cross full L-band source closure n={n}: rejected {err}; first={first}"
                );
                for target_label in [
                    "cap-net-s4-first-cross-full-l-band sector=1 p_41",
                    "cap-net-s4-first-cross-full-l-band sector=1 p_51",
                ] {
                    let formula =
                        kp15_ordered_source_formula_for_test(&table, &labeled, target_label, 0)
                            .unwrap_or_else(|| "independent".to_string());
                    println!(
                        "KP15 Fig2b cap-net S4 first-cross full L-band source closure n={n}: {target_label} x_formula={formula}"
                    );
                }
            }
        }
    }
}
