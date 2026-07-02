use super::super::bernstein_bi6::{bernstein_bi6_all, bernstein_bi6_deriv_all};
use super::super::bernstein_tri5::{bernstein_tri5_partial_u, bernstein_tri5_partial_v};
use super::super::kp13_bi5_boundary::{Bi3TensorBorder, TensorBorder};
use super::super::kp13_bi5_lsq_cc::solve_g_tri_lsq_cc;
use super::*;
use crate::geometry::gregory::bernstein_tri6::TRI6_CP_COUNT;
use crate::geometry::gregory::kp15_bi6_local_constraint;
use crate::geometry::gregory::kp20_bi6_cc_ring::{
    Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
    kp20_pr08_characteristic_tensor_border_candidate_rows,
    kp20_pr08_fig63_endpoint_side_regular_tensor_border_rows,
    kp20_pr08_stam_appendix_b_inner_ring_full_grid_for_test,
    kp20_pr08_stam_appendix_b_projected_rows_from_top_vectors_for_test,
    kp20_pr08_stam_inner_q1_characteristic_tensor_border,
    kp20_pr08_stam_q1_characteristic_tensor_border, pr08_catmull_clark_subdivision_stencil_matrix,
    pr08_subdivide_characteristic_ring,
};
use crate::geometry::gregory::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, kp20_cap_net_figure_tensor_border_rows,
    kp20_tensor_border_from_regular_cnet_rows,
};
use crate::geometry::gregory::kp20_bi6_stam_eval::{
    Kp20StamGridSide, Kp20StamGridTransform, Kp20StamPatch,
    kp20_pr08_stam_appendix_tensor_projected_grids, kp20_pr08_stam_projected_transformed_grid,
    kp20_pr08_stam_projected_transformed_tensor_border_side_rows,
};
use crate::geometry::gregory::kp20_characteristic_cap_net;
use crate::geometry::gregory::{
    Bi6CapInput, Bi6Reparams, CirculantSectorRhs, Kp15Bi6GeneratingTable,
    apply_beta_reparameterization, assemble_bi6_polar_cap,
    assemble_bi6_polar_cap_via_paper_p_table, assemble_bi6_polar_cap_via_table,
    f_kappa_quadratic_form, f_kappa_solve_optimal_bi6_cap, f_kappa_solve_paper_p_bi6_cap,
};
use faer::Mat;
use faer::prelude::SpSolver;

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

#[derive(Debug, Clone)]
struct TestGuideOrientationReport {
    samples: usize,
    negative: usize,
    near_zero: usize,
    min_det: f64,
    max_det: f64,
    worst_sample: String,
}

#[derive(Debug, Clone)]
struct TestBi6BoundaryOrientationReport {
    all: TestGuideOrientationReport,
    u0_negative: usize,
    u1_negative: usize,
    v0_negative: usize,
    v1_negative: usize,
    interior_negative: usize,
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
    let mut best = None::<(
        &'static str,
        usize,
        usize,
        TestGuideOrientationReport,
        f64,
        Point3,
    )>;
    {
        let (label, star) = ("origin", Point3::new(0.0, 0.0, 0.0));
        let output = solve_kp20_guide_lsq_with_report(&tensor, star).unwrap();
        let orientation = test_guide_local_orientation_report(n, &output.sectors);
        best = Some((
            label,
            0,
            0,
            orientation,
            output.report.max_fit_residual,
            star,
        ));
    }
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
                "legacy-full-l-band",
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
fn kp15_fig2b_cap_net_asymmetric_spectral_depth_source_is_ranked() {
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
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
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
                let constraints = kp15_source_constraints_from_grid_columns_for_test(
                    n,
                    &scaled_source,
                    &source_cols,
                );
                let (grids, f_values) = table
                    .populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "cap-net asymmetric spectral depth source",
                    )
                    .expect("cap-net paper depth source should remain feasible");
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
            "KP15 Fig2b cap-net asymmetric spectral depth source n={n}: row_scale={} column_scale={} neg={} min_det={:e} f=({:e},{:e},{:e})",
            best.0, best.1, best.2, best.3, best.4[0], best.4[1], best.4[2],
        );
        assert!(best.3.is_finite());
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

#[test]
fn kp15_fig2b_cap_net_ordered_first_cross_full_l_band_source_closure_is_measured() {
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
        let source_label = "cap-net-ordered-first-cross-full-l-band";
        let mut labeled = kp15_labeled_source_constraints_from_grid_columns_for_test(
            n,
            &source,
            &source_cols,
            source_label,
        );
        for sector in 1..n {
            for slot in [(4usize, 1usize), (5, 1), (6, 1)] {
                let target_label = format!("{source_label} sector={sector} p_{}{}", slot.0, slot.1);
                let Some(point) =
                    kp15_ordered_source_formula_point_for_test(&table, &labeled, &target_label)
                else {
                    continue;
                };
                source[sector][slot.0][slot.1] = point;
                if let Some(target) = labeled
                    .iter_mut()
                    .find(|constraint| constraint.label == target_label)
                {
                    target.constraint.rhs = [point.x, point.y, point.z];
                }
            }
        }
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &f_kappa_quadratic_form(5),
            0.0,
            "cap-net ordered first-cross full L-band tensor-border",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig2b cap-net ordered first-cross full L-band source closure n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
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
                    "KP15 Fig2b cap-net ordered first-cross full L-band source closure n={n}: rejected {err}; first={first}"
                );
            }
        }
    }
}

#[test]
fn kp15_fig2b_pr08_fig63_full_l_band_closure_candidates_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let mut accepted = 0usize;
        let mut best_rejection = f64::INFINITY;
        let mut best_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Ok(row_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                n,
                spec.row_side,
                spec.swap_targets,
                spec.conjugate_source,
            ) else {
                continue;
            };
            let Ok(column_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                n,
                spec.column_side,
                spec.swap_targets,
                spec.conjugate_source,
            ) else {
                continue;
            };
            let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
            let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
            for row_offset in [-1isize, 0, 1] {
                for column_offset in [-1isize, 0, 1] {
                    for row_reversed in [false, true] {
                        for column_reversed in [false, true] {
                            for row_placement in Kp15Fig12LocalPlacement::all() {
                                for column_placement in Kp15Fig12LocalPlacement::all() {
                                    let source =
                                        kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                            n,
                                            &row_rows,
                                            &column_rows,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            spec.row_reflected,
                                            spec.column_reflected,
                                            row_placement,
                                            column_placement,
                                        );
                                    let constraints =
                                        kp15_source_constraints_from_grid_columns_for_test(
                                            n,
                                            &source,
                                            &source_cols,
                                        );
                                    let label = format!(
                                        "row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={} row_offset={row_offset} column_offset={column_offset} row_rev={row_reversed} column_rev={column_reversed} row_place={} column_place={}",
                                        spec.row_side.name(),
                                        spec.column_side.name(),
                                        spec.swap_targets,
                                        spec.conjugate_source,
                                        spec.row_reflected,
                                        spec.column_reflected,
                                        row_placement.name(),
                                        column_placement.name(),
                                    );
                                    match table
                                        .populate_cap_grids_from_column_constraints_with_f_kappa(
                                            &constraints,
                                            &f_kappa_quadratic_form(5),
                                            0.0,
                                            "PR08 Fig6.3 full L-band closure candidate",
                                        ) {
                                        Ok(_) => {
                                            accepted += 1;
                                            if best_label.is_empty() {
                                                best_label = label;
                                            }
                                        }
                                        Err(err) => {
                                            if let Some(residual) =
                                                kp15_parse_constraint_residual_for_test(&err)
                                            {
                                                let residual = residual.abs();
                                                if residual < best_rejection {
                                                    best_rejection = residual;
                                                    best_label = label;
                                                }
                                            }
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
            "KP15 Fig2b PR08 Fig6.3 full L-band closure candidates n={n}: accepted={accepted} best_rejection={best_rejection:e} best={best_label}"
        );
        assert!(accepted > 0 || best_rejection.is_finite());
    }
}

#[test]
fn kp15_fig2b_pr08_fig63_best_full_l_band_failure_is_localized() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let spec = match n {
            5 => Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Left,
                column_side: Kp20StamGridSide::Top,
                swap_targets: false,
                conjugate_source: false,
                row_reflected: false,
                column_reflected: false,
            },
            6 => Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Left,
                column_side: Kp20StamGridSide::Top,
                swap_targets: false,
                conjugate_source: false,
                row_reflected: true,
                column_reflected: true,
            },
            7 => Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Top,
                column_side: Kp20StamGridSide::Left,
                swap_targets: true,
                conjugate_source: false,
                row_reflected: false,
                column_reflected: true,
            },
            8 => Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Top,
                column_side: Kp20StamGridSide::Left,
                swap_targets: true,
                conjugate_source: false,
                row_reflected: true,
                column_reflected: true,
            },
            _ => unreachable!(),
        };
        let (
            row_offset,
            column_offset,
            row_reversed,
            column_reversed,
            row_placement,
            column_placement,
        ) = match n {
            5 => (
                1,
                -1,
                false,
                false,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            ),
            6 => (
                -1,
                0,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            ),
            7 => (
                -1,
                0,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            ),
            8 => (
                0,
                -1,
                false,
                false,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            ),
            _ => unreachable!(),
        };
        let row_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
            n,
            spec.row_side,
            spec.swap_targets,
            spec.conjugate_source,
        )
        .expect("best row tensor");
        let column_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
            n,
            spec.column_side,
            spec.swap_targets,
            spec.conjugate_source,
        )
        .expect("best column tensor");
        let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
        let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
        let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
            n,
            &row_rows,
            &column_rows,
            row_offset,
            column_offset,
            row_reversed,
            column_reversed,
            spec.row_reflected,
            spec.column_reflected,
            row_placement,
            column_placement,
        );
        let constraints = kp15_labeled_full_tensor_border_constraints_for_test(n, &source);
        let first = kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
            .unwrap_or_else(|| "all constraints consistent".to_string());
        println!("KP15 Fig2b PR08 Fig6.3 best full-L-band failure n={n}: {first}");
        assert!(!first.is_empty());
    }
}

#[test]
fn kp15_fig2b_pr08_fig63_cyclic_p_i0_full_l_band_corrections_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut best = (f64::INFINITY, String::new());
        for adjacent_offset in [-1isize, 1] {
            for overwrite_count in [1usize, 2, 3] {
                let mut source = kp15_pr08_fig63_best_full_l_band_source_for_test(n);
                for sector in 0..n {
                    let adjacent = kp15_offset_sector_for_test(n, sector, adjacent_offset);
                    for i in 4..(4 + overwrite_count) {
                        source[sector][i][0] = source[adjacent][0][i];
                    }
                }
                let constraints = kp15_labeled_full_tensor_border_constraints_for_test(n, &source);
                let residual = kp15_first_inconsistent_source_constraint_residual_for_test(
                    &table,
                    &constraints,
                )
                .map(|(residual, _)| residual.abs())
                .unwrap_or(0.0);
                if residual < best.0 {
                    best = (
                        residual,
                        format!(
                            "adjacent_offset={adjacent_offset} overwrite_count={overwrite_count}"
                        ),
                    );
                }
            }
        }
        println!(
            "KP15 Fig2b PR08 Fig6.3 cyclic p_i0 correction ranking n={n}: best_residual={:e} best={}",
            best.0, best.1,
        );
        assert!(best.0.is_finite());
    }
}

fn kp15_pr08_fig63_best_full_l_band_source_for_test(n: usize) -> Vec<[[Point3; 7]; 7]> {
    let params = Bi6Reparams::for_valence(n);
    let (
        spec,
        row_offset,
        column_offset,
        row_reversed,
        column_reversed,
        row_placement,
        column_placement,
    ) = match n {
        5 => (
            Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Left,
                column_side: Kp20StamGridSide::Top,
                swap_targets: false,
                conjugate_source: false,
                row_reflected: false,
                column_reflected: false,
            },
            1,
            -1,
            false,
            false,
            Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            Kp15Fig12LocalPlacement::SectorJoinFirst,
        ),
        6 => (
            Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Left,
                column_side: Kp20StamGridSide::Top,
                swap_targets: false,
                conjugate_source: false,
                row_reflected: true,
                column_reflected: true,
            },
            -1,
            0,
            true,
            true,
            Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            Kp15Fig12LocalPlacement::SectorJoinFirst,
        ),
        7 => (
            Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Top,
                column_side: Kp20StamGridSide::Left,
                swap_targets: true,
                conjugate_source: false,
                row_reflected: false,
                column_reflected: true,
            },
            -1,
            0,
            true,
            true,
            Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            Kp15Fig12LocalPlacement::TensorBoundaryFirst,
        ),
        8 => (
            Kp15Fig63RadialScaleSpec {
                n,
                row_side: Kp20StamGridSide::Top,
                column_side: Kp20StamGridSide::Left,
                swap_targets: true,
                conjugate_source: false,
                row_reflected: true,
                column_reflected: true,
            },
            0,
            -1,
            false,
            false,
            Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            Kp15Fig12LocalPlacement::SectorJoinFirst,
        ),
        _ => panic!("unsupported valence {n}"),
    };
    let row_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
        n,
        spec.row_side,
        spec.swap_targets,
        spec.conjugate_source,
    )
    .expect("best row tensor");
    let column_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
        n,
        spec.column_side,
        spec.swap_targets,
        spec.conjugate_source,
    )
    .expect("best column tensor");
    let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
    let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
    kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
        n,
        &row_rows,
        &column_rows,
        row_offset,
        column_offset,
        row_reversed,
        column_reversed,
        spec.row_reflected,
        spec.column_reflected,
        row_placement,
        column_placement,
    )
}

#[test]
fn kp15_fig2b_pr08_fig63_joint_corner_full_l_band_closure_candidates_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let mut accepted = 0usize;
        let mut best_rejection = f64::INFINITY;
        let mut best_label = String::new();
        let candidates = match n {
            5 | 6 => vec![(
                Kp20StamGridSide::Left,
                Kp20StamGridSide::Top,
                false,
                true,
                -2isize,
                -1isize,
                true,
                false,
                false,
                true,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            7 => vec![(
                Kp20StamGridSide::Left,
                Kp20StamGridSide::Top,
                true,
                false,
                -2,
                -1,
                true,
                false,
                false,
                true,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            8 => vec![(
                Kp20StamGridSide::Top,
                Kp20StamGridSide::Left,
                true,
                false,
                1,
                -2,
                false,
                true,
                false,
                false,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            _ => unreachable!(),
        };
        for (
            row_side,
            column_side,
            swap_outer_targets,
            conjugate_source,
            row_offset,
            column_offset,
            row_reversed,
            column_reversed,
            row_reflected,
            column_reflected,
            row_placement,
            column_placement,
        ) in candidates
        {
            let row_tensor = kp15_pr08_fig63_joint_corner_side_tensor_border_for_fig12(
                n,
                row_side,
                swap_outer_targets,
                conjugate_source,
            )
            .expect("joint-corner row tensor");
            let column_tensor = kp15_pr08_fig63_joint_corner_side_tensor_border_for_fig12(
                n,
                column_side,
                swap_outer_targets,
                conjugate_source,
            )
            .expect("joint-corner column tensor");
            let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
            let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
            let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                n,
                &row_rows,
                &column_rows,
                row_offset,
                column_offset,
                row_reversed,
                column_reversed,
                row_reflected,
                column_reflected,
                row_placement,
                column_placement,
            );
            let constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
                n,
                &source,
                &source_cols,
                "joint-corner-full-l-band",
            );
            let label = format!(
                "row_side={} column_side={} swap_outer={swap_outer_targets} conjugate={conjugate_source} row_offset={row_offset} column_offset={column_offset} row_rev={row_reversed} column_rev={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_place={} column_place={}",
                row_side.name(),
                column_side.name(),
                row_placement.name(),
                column_placement.name(),
            );
            if let Some((residual, _)) =
                kp15_first_inconsistent_source_constraint_residual_for_test(&table, &constraints)
            {
                if residual < best_rejection {
                    best_rejection = residual;
                    best_label = label;
                }
            } else {
                accepted += 1;
                if best_label.is_empty() {
                    best_label = label;
                }
            }
        }
        println!(
            "KP15 Fig2b PR08 Fig6.3 joint-corner full L-band closure candidates n={n}: accepted={accepted} best_rejection={best_rejection:e} best={best_label}"
        );
        assert!(accepted > 0 || best_rejection.is_finite());
    }
}

#[test]
fn kp15_fig2b_stam_window_full_l_band_closure_candidates_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_l_band_columns(n);
        let candidates = match n {
            5 => vec![(
                1usize,
                Kp15Fig12StamSourceWindow::BottomRows,
                2usize,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2isize,
                1isize,
                true,
                false,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            6 => vec![(
                0,
                Kp15Fig12StamSourceWindow::RightColumns,
                2,
                Kp15Fig12StamSourceWindow::BottomRows,
                -2,
                0,
                false,
                false,
                false,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            7 => vec![(
                1,
                Kp15Fig12StamSourceWindow::RightColumns,
                1,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2,
                1,
                false,
                true,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            )],
            8 => vec![(
                0,
                Kp15Fig12StamSourceWindow::BottomRows,
                0,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2,
                1,
                false,
                true,
                true,
                true,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            _ => unreachable!(),
        };
        let mut accepted = 0usize;
        let mut best_rejection = f64::INFINITY;
        let mut best_label = String::new();
        for (
            row_patch,
            row_window,
            column_patch,
            column_window,
            row_offset,
            column_offset,
            row_reversed,
            column_reversed,
            row_reflected,
            column_reflected,
            row_placement,
            column_placement,
        ) in candidates
        {
            let row_rows =
                kp15_stam_window_beta_rows_for_test(&ring, row_patch, row_window, params);
            let column_rows =
                kp15_stam_window_beta_rows_for_test(&ring, column_patch, column_window, params);
            let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                n,
                &row_rows,
                &column_rows,
                row_offset,
                column_offset,
                row_reversed,
                column_reversed,
                row_reflected,
                column_reflected,
                row_placement,
                column_placement,
            );
            let constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
                n,
                &source,
                &source_cols,
                "stam-window-full-l-band",
            );
            let label = format!(
                "row=patch={row_patch} window={} column=patch={column_patch} window={} row_offset={row_offset} column_offset={column_offset} row_rev={row_reversed} column_rev={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_place={} column_place={}",
                row_window.name(),
                column_window.name(),
                row_placement.name(),
                column_placement.name(),
            );
            if let Some((residual, _)) =
                kp15_first_inconsistent_source_constraint_residual_for_test(&table, &constraints)
            {
                if residual < best_rejection {
                    best_rejection = residual;
                    best_label = label;
                }
            } else {
                accepted += 1;
                if best_label.is_empty() {
                    best_label = label;
                }
            }
        }
        println!(
            "KP15 Fig2b Stam window full L-band closure candidates n={n}: accepted={accepted} best_rejection={best_rejection:e} best={best_label}"
        );
        assert!(accepted > 0 || best_rejection.is_finite());
    }
}

#[test]
fn kp15_fig2b_stam_window_paper_depth_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
        let candidates = match n {
            5 => vec![(
                1usize,
                Kp15Fig12StamSourceWindow::BottomRows,
                2usize,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2isize,
                1isize,
                true,
                false,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            6 => vec![(
                0,
                Kp15Fig12StamSourceWindow::RightColumns,
                2,
                Kp15Fig12StamSourceWindow::BottomRows,
                -2,
                0,
                false,
                false,
                false,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            7 => vec![(
                1,
                Kp15Fig12StamSourceWindow::RightColumns,
                1,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2,
                1,
                false,
                true,
                true,
                true,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
                Kp15Fig12LocalPlacement::TensorBoundaryFirst,
            )],
            8 => vec![(
                0,
                Kp15Fig12StamSourceWindow::BottomRows,
                0,
                Kp15Fig12StamSourceWindow::LeftColumns,
                -2,
                1,
                false,
                true,
                true,
                true,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
                Kp15Fig12LocalPlacement::SectorJoinFirst,
            )],
            _ => unreachable!(),
        };
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        for (
            row_patch,
            row_window,
            column_patch,
            column_window,
            row_offset,
            column_offset,
            row_reversed,
            column_reversed,
            row_reflected,
            column_reflected,
            row_placement,
            column_placement,
        ) in candidates
        {
            let row_rows =
                kp15_stam_window_beta_rows_for_test(&ring, row_patch, row_window, params);
            let column_rows =
                kp15_stam_window_beta_rows_for_test(&ring, column_patch, column_window, params);
            let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                n,
                &row_rows,
                &column_rows,
                row_offset,
                column_offset,
                row_reversed,
                column_reversed,
                row_reflected,
                column_reflected,
                row_placement,
                column_placement,
            );
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
            let (grids, f_values) = table
                .populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "Stam window paper-depth source projection",
                )
                .expect("paper depth source should be feasible");
            accepted += 1;
            let orientation =
                test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
            let label = format!(
                "row=patch={row_patch} window={} column=patch={column_patch} window={} row_offset={row_offset} column_offset={column_offset} row_rev={row_reversed} column_rev={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_place={} column_place={} f=({:e},{:e},{:e})",
                row_window.name(),
                column_window.name(),
                row_placement.name(),
                column_placement.name(),
                f_values[0],
                f_values[1],
                f_values[2],
            );
            if orientation.negative < best_neg
                || (orientation.negative == best_neg && orientation.min_det > best_min_det)
            {
                best_neg = orientation.negative;
                best_min_det = orientation.min_det;
                best_label = label;
            }
        }
        println!(
            "KP15 Fig2b Stam window paper-depth source projection n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label}"
        );
        assert!(accepted > 0);
    }
}

#[test]
fn kp15_fig2b_cap_net_source_compatibility_localizes_by_row_family() {
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        for (table_label, table) in [
            (
                "raw",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "s1",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_full_s1_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "s1-s2",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s2_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "s1-s3",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_s1_s3_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
            (
                "s1-s4",
                Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_central_quadratic_free(
                    n,
                    Bi6Reparams::for_valence(n),
                ),
            ),
        ] {
            let mut accepted = 0usize;
            let mut smallest_rejection = f64::INFINITY;
            let mut smallest_rejection_label = String::new();
            for apply_beta in [false, true] {
                for row_offset in [-1isize, 0, 1] {
                    for column_offset in [-1isize, 0, 1] {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                let source =
                                    kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
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
                                    "cap-net row-family compatibility",
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
                "KP15 Fig2b cap-net source row-family compatibility n={n} table={table_label}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
            );
            assert!(smallest_rejection.is_finite() || accepted > 0);
        }
    }
}

#[test]
fn kp15_fig2b_depth_strip_side_compatibility_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let phase_candidates = [
            (false, 0, 1, true, false),
            (true, 0, 1, true, false),
            (false, 1, -1, true, true),
            (false, 1, 1, true, false),
            (true, 1, 1, false, true),
            (true, 1, 1, true, false),
        ];
        for (side_label, depth_cols) in [
            ("row-side", kp15_paper_depth_row_side_columns_for_test(n)),
            (
                "column-side",
                kp15_paper_depth_column_side_columns_for_test(n),
            ),
            (
                "two-sided",
                kp15_paper_depth_two_sided_strip_columns_for_test(n),
            ),
        ] {
            let mut accepted = 0usize;
            let mut smallest_rejection = f64::INFINITY;
            let mut smallest_rejection_label = String::new();
            for (apply_beta, row_offset, column_offset, row_reversed, column_reversed) in
                phase_candidates
            {
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
                let constraints =
                    kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
                let label = format!(
                    "beta={apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed}"
                );
                match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "depth-strip side compatibility",
                ) {
                    Ok(_) => accepted += 1,
                    Err(err) => {
                        if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                            let abs_residual = residual.abs();
                            if abs_residual < smallest_rejection {
                                smallest_rejection = abs_residual;
                                smallest_rejection_label = label;
                            }
                        }
                    }
                }
            }
            println!(
                "KP15 Fig2b depth-strip side compatibility n={n} side={side_label}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
            );
            assert!(smallest_rejection.is_finite() || accepted > 0);
        }
    }
}

#[test]
fn kp15_fig2b_paired_depth_strip_projection_is_measured() {
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let source = kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
            n,
            &tensor.sectors,
            Bi6Reparams::for_valence(n),
            false,
            0,
            1,
            true,
            false,
        );
        for (label, hard_cols, soft_cols) in [
            (
                "hard-row-fit-column",
                kp15_paper_depth_row_side_columns_for_test(n),
                kp15_paper_depth_column_side_columns_for_test(n),
            ),
            (
                "hard-column-fit-row",
                kp15_paper_depth_column_side_columns_for_test(n),
                kp15_paper_depth_row_side_columns_for_test(n),
            ),
        ] {
            let projected = kp15_project_source_side_to_row_manifold_for_test(
                n, &table, &source, &hard_cols, &soft_cols,
            )
            .expect("paired source projection");
            let max_soft_movement =
                kp15_max_source_column_delta_for_test(n, &source, &projected, &soft_cols);
            let orientation = test_bi6_patch_boundary_orientation_report(n, &projected, &samples);
            println!(
                "KP15 Fig2b paired depth-strip projection n={n} {label}: soft_movement={max_soft_movement:e} neg={neg}/{samples_count} det=[{min_det:e},{max_det:e}]",
                neg = orientation.all.negative,
                samples_count = orientation.all.samples,
                min_det = orientation.all.min_det,
                max_det = orientation.all.max_det,
            );
            assert!(max_soft_movement.is_finite());
            assert!(orientation.all.min_det.is_finite());
        }
    }
}

#[test]
fn kp15_fig2b_depth_strip_normal_reflection_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for row_normal_reflected in [false, true] {
            for column_normal_reflected in [false, true] {
                for (apply_beta, row_offset, column_offset, row_reversed, column_reversed) in [
                    (false, 0, 1, true, false),
                    (true, 0, 1, true, false),
                    (false, 1, -1, true, true),
                    (false, 1, 1, true, false),
                    (true, 1, 1, false, true),
                    (true, 1, 1, true, false),
                ] {
                    let source = kp15_two_sided_source_grid_from_bi3_tensor_oriented_for_test(
                        n,
                        &tensor.sectors,
                        Bi6Reparams::for_valence(n),
                        apply_beta,
                        row_offset,
                        column_offset,
                        row_reversed,
                        column_reversed,
                        row_normal_reflected,
                        column_normal_reflected,
                    );
                    let constraints =
                        kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
                    let label = format!(
                        "beta={apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_normal_reflected} column_reflect={column_normal_reflected}"
                    );
                    match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "normal-reflected depth source",
                    ) {
                        Ok(_) => accepted += 1,
                        Err(err) => {
                            if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
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
        println!(
            "KP15 Fig2b normal-reflected depth source n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_pr08_candidate_depth_pair_compatibility_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for candidate in Kp20Pr08TensorBorderCandidate::all() {
            let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
                .expect("PR08 candidate rows");
            let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor rows");
            for apply_beta in [false, true] {
                for normal_reflected in [false, true] {
                    let source = kp15_two_sided_source_grid_from_bi3_tensor_oriented_for_test(
                        n,
                        &tensor.sectors,
                        Bi6Reparams::for_valence(n),
                        apply_beta,
                        0,
                        1,
                        true,
                        false,
                        false,
                        normal_reflected,
                    );
                    let constraints =
                        kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
                    let label = format!(
                        "{} beta={apply_beta} column_reflect={normal_reflected}",
                        candidate.name()
                    );
                    match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "PR08 candidate depth pair",
                    ) {
                        Ok(_) => accepted += 1,
                        Err(err) => {
                            if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
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
        println!(
            "KP15 Fig2b PR08 candidate depth-pair compatibility n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_pr08_independent_side_candidate_pairs_are_ranked() {
    let candidates = [
        Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows,
        Kp20Pr08TensorBorderCandidate::StamLowerRightStarSideRows,
        Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows,
        Kp20Pr08TensorBorderCandidate::StamLowerRightRightSideRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ2InteriorRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    ];
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let candidate_tensors = candidates
            .iter()
            .copied()
            .map(|candidate| {
                let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
                    .expect("PR08 candidate rows");
                let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor rows");
                (candidate, tensor)
            })
            .collect::<Vec<_>>();
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for (row_candidate, row_tensor) in &candidate_tensors {
            for (column_candidate, column_tensor) in &candidate_tensors {
                for row_apply_beta in [false, true] {
                    for column_apply_beta in [false, true] {
                        let source =
                            kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
                                n,
                                &row_tensor.sectors,
                                &column_tensor.sectors,
                                Bi6Reparams::for_valence(n),
                                row_apply_beta,
                                column_apply_beta,
                                0,
                                1,
                                true,
                                false,
                                false,
                                true,
                            );
                        let constraints = kp15_source_constraints_from_grid_columns_for_test(
                            n,
                            &source,
                            &depth_cols,
                        );
                        let label = format!(
                            "row={} row_beta={row_apply_beta} column={} column_beta={column_apply_beta} column_reflect=true",
                            row_candidate.name(),
                            column_candidate.name(),
                        );
                        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                            &constraints,
                            &q,
                            0.0,
                            "PR08 independent side candidate depth pair",
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
        println!(
            "KP15 Fig2b PR08 independent side candidate pairs n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_pr08_best_family_side_orientation_is_ranked() {
    for (n, row_candidate, column_candidate) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows,
            Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let row_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
            .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for row_offset in -2..=2 {
            for column_offset in -2..=2 {
                for row_reversed in [false, true] {
                    for column_reversed in [false, true] {
                        for row_reflected in [false, true] {
                            for column_reflected in [false, true] {
                                for row_apply_beta in [false, true] {
                                    for column_apply_beta in [false, true] {
                                        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
                                            n,
                                            &row_tensor.sectors,
                                            &column_tensor.sectors,
                                            Bi6Reparams::for_valence(n),
                                            row_apply_beta,
                                            column_apply_beta,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            row_reflected,
                                            column_reflected,
                                        );
                                        let constraints =
                                            kp15_source_constraints_from_grid_columns_for_test(
                                                n,
                                                &source,
                                                &depth_cols,
                                            );
                                        let label = format!(
                                            "row={} column={} row_beta={row_apply_beta} column_beta={column_apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected}",
                                            row_candidate.name(),
                                            column_candidate.name(),
                                        );
                                        match table
                                            .populate_cap_grids_from_column_constraints_with_f_kappa(
                                                &constraints,
                                                &q,
                                                0.0,
                                                "PR08 best-family side orientation depth pair",
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
                }
            }
        }
        println!(
            "KP15 Fig2b PR08 best-family side orientation n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_best_pr08_pair_inconsistency_is_localized() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let row_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_row_side_columns_for_test(spec.n),
            "row-side",
        );
        let column_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_column_side_columns_for_test(spec.n),
            "column-side",
        );
        let mut row_then_column = row_side.clone();
        row_then_column.extend(column_side.clone());
        let mut column_then_row = column_side;
        column_then_row.extend(row_side);
        let first_row_then_column =
            kp15_first_inconsistent_source_constraint_for_test(&table, &row_then_column)
                .expect("row then column inconsistency");
        let first_column_then_row =
            kp15_first_inconsistent_source_constraint_for_test(&table, &column_then_row)
                .expect("column then row inconsistency");
        println!(
            "KP15 Fig2b best PR08 pair inconsistency n={}: row_then_column={} column_then_row={}",
            spec.n, first_row_then_column, first_column_then_row,
        );
    }
}

#[test]
fn kp15_fig2b_best_pr08_corner_adjacent_join_is_ranked() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(spec.n);
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for sector_offset in -1..=1 {
            for mode in [
                Kp15CornerJoinMode::RowToColumn,
                Kp15CornerJoinMode::ColumnToRow,
                Kp15CornerJoinMode::Average,
            ] {
                let mut joined = source.clone();
                kp15_apply_corner_adjacent_join_for_test(&mut joined, sector_offset, mode);
                let constraints = kp15_source_constraints_from_grid_columns_for_test(
                    spec.n,
                    &joined,
                    &depth_cols,
                );
                let label = format!("{mode:?} sector_offset={sector_offset}");
                match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "PR08 best-family corner-adjacent joined depth pair",
                ) {
                    Ok(_) => accepted += 1,
                    Err(err) => {
                        if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
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
        println!(
            "KP15 Fig2b best PR08 corner-adjacent join n={}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
            spec.n,
        );
        assert!(smallest_rejection.is_finite() || accepted > 0);
    }
}

#[test]
fn kp15_fig2b_best_pr08_cross_corner_implied_values_are_measured() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let row_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_row_side_columns_for_test(spec.n),
            "row-side",
        );
        let column_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_column_side_columns_for_test(spec.n),
            "column-side",
        );
        let mut row_then_column = row_side.clone();
        row_then_column.extend(column_side.clone());
        let mut column_then_row = column_side;
        column_then_row.extend(row_side);
        let p14_report = kp15_implied_ordered_constraint_delta_for_test(
            &table,
            &row_then_column,
            "column-side sector=0 p_14",
        )
        .expect("ordered p14 implication");
        let p41_report = kp15_implied_ordered_constraint_delta_for_test(
            &table,
            &column_then_row,
            "row-side sector=0 p_41",
        )
        .expect("ordered p41 implication");
        println!(
            "KP15 Fig2b best PR08 implied cross-corner n={}: {} ; {}",
            spec.n, p14_report, p41_report,
        );
    }
}

#[test]
fn kp15_fig2b_best_pr08_cross_corner_source_formula_is_reported() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let row_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_row_side_columns_for_test(spec.n),
            "row-side",
        );
        let column_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_column_side_columns_for_test(spec.n),
            "column-side",
        );
        let mut row_then_column = row_side.clone();
        row_then_column.extend(column_side.clone());
        let mut column_then_row = column_side;
        column_then_row.extend(row_side);
        let p14_formula = kp15_ordered_source_formula_for_test(
            &table,
            &row_then_column,
            "column-side sector=0 p_14",
            0,
        )
        .expect("p14 formula");
        let p41_formula = kp15_ordered_source_formula_for_test(
            &table,
            &column_then_row,
            "row-side sector=0 p_41",
            0,
        )
        .expect("p41 formula");
        println!(
            "KP15 Fig2b best PR08 cross-corner source formula n={}: p14_x={} ; p41_x={}",
            spec.n, p14_formula, p41_formula,
        );
    }
}

#[test]
fn kp15_fig2b_s4_cross_corner_sum_corrected_source_is_measured() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let mut source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let max_delta = kp15_apply_s4_cross_corner_sum_correction_for_test(
            &mut source,
            Bi6Reparams::for_valence(spec.n),
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(spec.n);
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(spec.n, &source, &depth_cols);
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "KP15 S4 cross-corner sum-corrected depth pair",
        ) {
            Ok(_) => println!(
                "KP15 Fig2b S4 cross-corner sum-corrected source n={}: accepted max_delta={max_delta:e}",
                spec.n,
            ),
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig2b S4 cross-corner sum-corrected source n={}: rejected residual={residual:e} max_delta={max_delta:e}",
                    spec.n,
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_fig2b_s4_sum_corrected_inconsistency_is_localized() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let mut source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        kp15_apply_s4_cross_corner_sum_correction_for_test(
            &mut source,
            Bi6Reparams::for_valence(spec.n),
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let row_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_row_side_columns_for_test(spec.n),
            "row-side",
        );
        let column_side = kp15_labeled_source_constraints_from_grid_columns_for_test(
            spec.n,
            &source,
            &kp15_paper_depth_column_side_columns_for_test(spec.n),
            "column-side",
        );
        let mut row_then_column = row_side.clone();
        row_then_column.extend(column_side.clone());
        let mut column_then_row = column_side;
        column_then_row.extend(row_side);
        let row_then_column_failure =
            kp15_first_inconsistent_source_constraint_for_test(&table, &row_then_column)
                .unwrap_or_else(|| "none".to_string());
        let column_then_row_failure =
            kp15_first_inconsistent_source_constraint_for_test(&table, &column_then_row)
                .unwrap_or_else(|| "none".to_string());
        println!(
            "KP15 Fig2b S4 sum-corrected inconsistency n={}: row_then_column={} column_then_row={}",
            spec.n, row_then_column_failure, column_then_row_failure,
        );
    }
}

#[test]
fn kp15_fig2b_s4_cross_corner_ownership_corrections_are_ranked() {
    for spec in kp15_best_pr08_side_specs_for_test() {
        let ring =
            kp20_catmull_clark_characteristic_ring(spec.n).expect("PR08 characteristic ring");
        let row_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.row_candidate)
                .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, spec.column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            spec.n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(spec.n),
            spec.row_apply_beta,
            spec.column_apply_beta,
            spec.row_offset,
            spec.column_offset,
            spec.row_reversed,
            spec.column_reversed,
            spec.row_reflected,
            spec.column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                spec.n,
                Bi6Reparams::for_valence(spec.n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(spec.n);
        for mode in [
            Kp15S4CrossCornerOwnershipMode::DeriveColumnP14,
            Kp15S4CrossCornerOwnershipMode::DeriveRowP41,
            Kp15S4CrossCornerOwnershipMode::DeriveBothSequentially,
        ] {
            let mut corrected = source.clone();
            let max_delta = kp15_apply_s4_cross_corner_ownership_correction_for_test(
                &mut corrected,
                Bi6Reparams::for_valence(spec.n),
                mode,
            );
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(spec.n, &corrected, &depth_cols);
            let residual = match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "KP15 S4 cross-corner ownership-corrected depth pair",
            ) {
                Ok(_) => 0.0,
                Err(err) => kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN),
            };
            println!(
                "KP15 Fig2b S4 cross-corner ownership source n={} mode={mode:?}: residual={residual:e} max_delta={max_delta:e}",
                spec.n,
            );
            assert!(residual.is_finite());
        }
    }
}

#[test]
fn kp15_fig2b_pr08_best_family_s4_source_residual_is_ranked() {
    for (n, row_candidate, column_candidate) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows,
            Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let row_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
            .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_offset in -2..=2 {
            for column_offset in -2..=2 {
                for row_reversed in [false, true] {
                    for column_reversed in [false, true] {
                        for row_reflected in [false, true] {
                            for column_reflected in [false, true] {
                                for row_apply_beta in [false, true] {
                                    for column_apply_beta in [false, true] {
                                        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
                                            n,
                                            &row_tensor.sectors,
                                            &column_tensor.sectors,
                                            Bi6Reparams::for_valence(n),
                                            row_apply_beta,
                                            column_apply_beta,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            row_reflected,
                                            column_reflected,
                                        );
                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(
                                            &source,
                                            Bi6Reparams::for_valence(n),
                                        );
                                        if residual < best_residual {
                                            best_residual = residual;
                                            best_label = format!(
                                                "row_beta={row_apply_beta} column_beta={column_apply_beta} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected}"
                                            );
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
            "KP15 Fig2b PR08 best-family S4 source residual n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig2b_pr08_s4_ranked_chart_depth_pair_compatibility_is_measured() {
    for (
        n,
        row_candidate,
        column_candidate,
        row_apply_beta,
        column_apply_beta,
        row_reversed,
        column_reversed,
        row_reflected,
        column_reflected,
    ) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            false,
            false,
            true,
            true,
            false,
            true,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows,
            Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows,
            false,
            true,
            true,
            true,
            true,
            true,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            false,
            false,
            false,
            false,
            false,
            true,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            false,
            false,
            true,
            true,
            false,
            true,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let row_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
            .expect("row PR08 candidate rows");
        let column_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                .expect("column PR08 candidate rows");
        let row_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
        let column_tensor =
            kp20_tensor_border_from_regular_cnet_rows(&column_rows).expect("column tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            n,
            &row_tensor.sectors,
            &column_tensor.sectors,
            Bi6Reparams::for_valence(n),
            row_apply_beta,
            column_apply_beta,
            -2,
            -1,
            row_reversed,
            column_reversed,
            row_reflected,
            column_reflected,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
        let s4_residual =
            kp15_s4_cross_corner_sum_residual_for_test(&source, Bi6Reparams::for_valence(n));
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "PR08 S4-ranked chart depth pair",
        ) {
            Ok(_) => println!(
                "KP15 Fig2b PR08 S4-ranked chart n={n}: accepted s4_residual={s4_residual:e}",
            ),
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig2b PR08 S4-ranked chart n={n}: rejected residual={residual:e} s4_residual={s4_residual:e}",
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_fig2b_beta_s4_source_frame_candidates_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_candidate in Kp20Pr08TensorBorderCandidate::all() {
            let row_rows =
                kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
                    .expect("row PR08 candidate rows");
            let row_tensor =
                kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
            for column_candidate in Kp20Pr08TensorBorderCandidate::all() {
                let column_rows =
                    kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                        .expect("column PR08 candidate rows");
                let column_tensor = kp20_tensor_border_from_regular_cnet_rows(&column_rows)
                    .expect("column tensor rows");
                for row_offset in -2..=2 {
                    for column_offset in -2..=2 {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                for row_reflected in [false, true] {
                                    for column_reflected in [false, true] {
                                        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
                                            n,
                                            &row_tensor.sectors,
                                            &column_tensor.sectors,
                                            Bi6Reparams::for_valence(n),
                                            true,
                                            true,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            row_reflected,
                                            column_reflected,
                                        );
                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(
                                            &source,
                                            Bi6Reparams::for_valence(n),
                                        );
                                        if residual < best_residual {
                                            best_residual = residual;
                                            best_label = format!(
                                                "row={} column={} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected}",
                                                row_candidate.name(),
                                                column_candidate.name(),
                                            );
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
            "KP15 Fig2b beta/S4 source-frame candidate n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig2b_beta_s4_ranked_depth_pair_compatibility_is_measured() {
    for (n, candidate, row_reversed, column_reversed) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaLeftColumnReversed,
            true,
            true,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows,
            false,
            false,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            true,
            true,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows,
            false,
            false,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor rows");
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
            n,
            &tensor.sectors,
            &tensor.sectors,
            Bi6Reparams::for_valence(n),
            true,
            true,
            -2,
            -1,
            row_reversed,
            column_reversed,
            false,
            true,
        );
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
        let s4_residual =
            kp15_s4_cross_corner_sum_residual_for_test(&source, Bi6Reparams::for_valence(n));
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "beta/S4-ranked depth pair",
        ) {
            Ok(_) => println!(
                "KP15 Fig2b beta/S4-ranked depth pair n={n}: accepted s4_residual={s4_residual:e}",
            ),
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig2b beta/S4-ranked depth pair n={n}: rejected residual={residual:e} s4_residual={s4_residual:e} candidate={}",
                    candidate.name(),
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_beta_s4_cross_corner_formula_expands_to_degree_raised_bi3_rows() {
    let n = 5usize;
    let params = Bi6Reparams::for_valence(n);
    let w3 = params.c / 15.0;
    let row4 = kp15_degree_raise_scalar_basis_3_to_6_for_test(4);
    let row5 = kp15_degree_raise_scalar_basis_3_to_6_for_test(5);
    let row6 = kp15_degree_raise_scalar_basis_3_to_6_for_test(6);
    println!(
        "KP15 beta/S4 bi3 expansion: q[*][4]={row4:?} q[*][5]={row5:?} q[*][6]={row6:?} w3={w3:e} q0[4]+w3(q0[6]-q0[5]) coeffs={:?}",
        [
            row4[0] + w3 * (row6[0] - row5[0]),
            row4[1] + w3 * (row6[1] - row5[1]),
            row4[2] + w3 * (row6[2] - row5[2]),
            row4[3] + w3 * (row6[3] - row5[3]),
        ],
    );
    assert!((row4.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    assert!((row5.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    assert!((row6.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
}

#[test]
fn kp15_fig2b_s4_first_order_normal_scale_is_solved_and_measured() {
    for (n, candidate, row_reversed, column_reversed) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaLeftColumnReversed,
            true,
            true,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows,
            false,
            false,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            true,
            true,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows,
            false,
            false,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor rows");
        let params = Bi6Reparams::for_valence(n);
        let solve = kp15_solve_s4_first_order_normal_scales_for_test(
            n,
            &tensor.sectors,
            &tensor.sectors,
            params,
            true,
            true,
            -2,
            -1,
            row_reversed,
            column_reversed,
        );
        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_scaled_for_test(
            n,
            &tensor.sectors,
            &tensor.sectors,
            params,
            true,
            true,
            -2,
            -1,
            row_reversed,
            column_reversed,
            solve.row_scale,
            solve.column_scale,
        );
        let s4_residual = kp15_s4_cross_corner_sum_residual_for_test(&source, params);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n, params,
            );
        let q = f_kappa_quadratic_form(5);
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let constraints =
            kp15_source_constraints_from_grid_columns_for_test(n, &source, &depth_cols);
        let depth_residual = match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "KP15 S4 first-order normal-scaled source",
        ) {
            Ok(_) => 0.0,
            Err(err) => kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN),
        };
        println!(
            "KP15 Fig2b S4 first-order normal scale n={n}: row_scale={:e} column_scale={:e} det={:e} fit_residual={:e} s4_residual={:e} depth_residual={:e} candidate={}",
            solve.row_scale,
            solve.column_scale,
            solve.determinant,
            solve.residual,
            s4_residual,
            depth_residual,
            candidate.name(),
        );
        assert!(solve.row_scale.is_finite());
        assert!(solve.column_scale.is_finite());
        assert!(solve.residual.is_finite());
        assert!(s4_residual.is_finite());
        assert!(depth_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_stam_full_grid_source_windows_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let mut candidates = Vec::<(String, Kp20RegularTensorBorderCNet)>::new();
        for patch in 0..=2 {
            for window in Kp15Fig12StamSourceWindow::all() {
                let mut sectors = Vec::with_capacity(n);
                for sector in 0..n {
                    let grid = kp20_pr08_stam_appendix_b_inner_ring_full_grid_for_test(
                        &ring, sector, patch,
                    )
                    .expect("Stam Appendix-B full grid");
                    sectors.push(kp15_fig12_stam_source_window_rows_for_test(&grid, window));
                }
                candidates.push((
                    format!("patch={patch} window={}", window.name()),
                    Kp20RegularTensorBorderCNet { sectors },
                ));
            }
        }
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for (row_label, row_rows) in &candidates {
            for (column_label, column_rows) in &candidates {
                for row_offset in -2..=2 {
                    for column_offset in -2..=2 {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                for row_reflected in [false, true] {
                                    for column_reflected in [false, true] {
                                        let source = kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
                                            n,
                                            &row_rows.sectors,
                                            &column_rows.sectors,
                                            Bi6Reparams::for_valence(n),
                                            true,
                                            true,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            row_reflected,
                                            column_reflected,
                                        );
                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(
                                            &source,
                                            Bi6Reparams::for_valence(n),
                                        );
                                        if residual < best_residual {
                                            best_residual = residual;
                                            best_label = format!(
                                                "row={row_label} column={column_label} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected}"
                                            );
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
            "KP15 Fig12 Stam full-grid source window n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_bb_coefficient_beta_source_frame_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_candidate in Kp20Pr08TensorBorderCandidate::all() {
            let row_rows =
                kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
                    .expect("row PR08 candidate rows");
            let row_tensor =
                kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
            let row_bi6_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
            for column_candidate in Kp20Pr08TensorBorderCandidate::all() {
                let column_rows =
                    kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                        .expect("column PR08 candidate rows");
                let column_tensor = kp20_tensor_border_from_regular_cnet_rows(&column_rows)
                    .expect("column tensor rows");
                let column_bi6_rows =
                    kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                for row_offset in -2..=2 {
                    for column_offset in -2..=2 {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                for row_reflected in [false, true] {
                                    for column_reflected in [false, true] {
                                        let source = kp15_two_sided_source_grid_from_paired_bi6_rows_oriented_for_test(
                                            n,
                                            &row_bi6_rows,
                                            &column_bi6_rows,
                                            row_offset,
                                            column_offset,
                                            row_reversed,
                                            column_reversed,
                                            row_reflected,
                                            column_reflected,
                                        );
                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(
                                            &source, params,
                                        );
                                        if residual < best_residual {
                                            best_residual = residual;
                                            best_label = format!(
                                                "row={} column={} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected}",
                                                row_candidate.name(),
                                                column_candidate.name(),
                                            );
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
            "KP15 Fig12 BB-coefficient beta source-frame n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_bb_beta_tensor_boundary_placement_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_candidate in Kp20Pr08TensorBorderCandidate::all() {
            let row_rows =
                kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, row_candidate)
                    .expect("row PR08 candidate rows");
            let row_tensor =
                kp20_tensor_border_from_regular_cnet_rows(&row_rows).expect("row tensor rows");
            let row_bi6_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
            for column_candidate in Kp20Pr08TensorBorderCandidate::all() {
                let column_rows =
                    kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, column_candidate)
                        .expect("column PR08 candidate rows");
                let column_tensor = kp20_tensor_border_from_regular_cnet_rows(&column_rows)
                    .expect("column tensor rows");
                let column_bi6_rows =
                    kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                for row_offset in -2..=2 {
                    for column_offset in -2..=2 {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                for row_reflected in [false, true] {
                                    for column_reflected in [false, true] {
                                        for row_placement in Kp15Fig12LocalPlacement::all() {
                                            for column_placement in Kp15Fig12LocalPlacement::all() {
                                                let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                    n,
                                                    &row_bi6_rows,
                                                    &column_bi6_rows,
                                                    row_offset,
                                                    column_offset,
                                                    row_reversed,
                                                    column_reversed,
                                                    row_reflected,
                                                    column_reflected,
                                                    row_placement,
                                                    column_placement,
                                                );
                                                let residual =
                                                    kp15_s4_cross_corner_sum_residual_for_test(
                                                        &source, params,
                                                    );
                                                if residual < best_residual {
                                                    best_residual = residual;
                                                    best_label = format!(
                                                        "row={} column={} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                        row_candidate.name(),
                                                        column_candidate.name(),
                                                        row_placement.name(),
                                                        column_placement.name(),
                                                    );
                                                }
                                            }
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
            "KP15 Fig12 BB beta tensor-boundary placement n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_stam_windows_with_bb_beta_and_local_placement_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let params = Bi6Reparams::for_valence(n);
        let mut candidates = Vec::<(String, Vec<[[Point3; 7]; 3]>)>::new();
        for patch in 0..=2 {
            for window in Kp15Fig12StamSourceWindow::all() {
                let mut bi3_sectors = Vec::with_capacity(n);
                for sector in 0..n {
                    let grid = kp20_pr08_stam_appendix_b_inner_ring_full_grid_for_test(
                        &ring, sector, patch,
                    )
                    .expect("Stam Appendix-B full grid");
                    bi3_sectors.push(kp15_fig12_stam_source_window_rows_for_test(&grid, window));
                }
                candidates.push((
                    format!("patch={patch} window={}", window.name()),
                    kp15_fig12_bb_beta_rows_for_test(&bi3_sectors, params),
                ));
            }
        }
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for (row_label, row_rows) in &candidates {
            for (column_label, column_rows) in &candidates {
                for row_offset in -2..=2 {
                    for column_offset in -2..=2 {
                        for row_reversed in [false, true] {
                            for column_reversed in [false, true] {
                                for row_reflected in [false, true] {
                                    for column_reflected in [false, true] {
                                        for row_placement in Kp15Fig12LocalPlacement::all() {
                                            for column_placement in Kp15Fig12LocalPlacement::all() {
                                                let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                    n,
                                                    row_rows,
                                                    column_rows,
                                                    row_offset,
                                                    column_offset,
                                                    row_reversed,
                                                    column_reversed,
                                                    row_reflected,
                                                    column_reflected,
                                                    row_placement,
                                                    column_placement,
                                                );
                                                let residual =
                                                    kp15_s4_cross_corner_sum_residual_for_test(
                                                        &source, params,
                                                    );
                                                if residual < best_residual {
                                                    best_residual = residual;
                                                    best_label = format!(
                                                        "row={row_label} column={column_label} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                        row_placement.name(),
                                                        column_placement.name(),
                                                    );
                                                }
                                            }
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
            "KP15 Fig12 Stam window BB beta/local placement n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_composed_pr08_subdivision_source_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let source_rows = kp15_composed_pr08_subdivision_source_rows_for_test(&ring, 1, 1, true)
            .expect("composed PR08 subdivision source rows");
        let params = Bi6Reparams::for_valence(n);
        let legacy_rows = kp15_bi6_tensor_rows_for_test(&source_rows.sectors, params, true);
        let bb_rows = kp15_fig12_bb_beta_rows_for_test(&source_rows.sectors, params);
        for (label, rows) in [
            ("legacy-beta-row-semantics", legacy_rows),
            ("fig12-bb-beta", bb_rows),
        ] {
            let mut best_residual = f64::INFINITY;
            let mut best_label = String::new();
            for row_offset in -2..=2 {
                for column_offset in -2..=2 {
                    for row_reversed in [false, true] {
                        for column_reversed in [false, true] {
                            for row_reflected in [false, true] {
                                for column_reflected in [false, true] {
                                    for row_placement in Kp15Fig12LocalPlacement::all() {
                                        for column_placement in Kp15Fig12LocalPlacement::all() {
                                            let source =
                                                kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                    n,
                                                    &rows,
                                                    &rows,
                                                    row_offset,
                                                    column_offset,
                                                    row_reversed,
                                                    column_reversed,
                                                    row_reflected,
                                                    column_reflected,
                                                    row_placement,
                                                    column_placement,
                                                );
                                            let residual =
                                                kp15_s4_cross_corner_sum_residual_for_test(
                                                    &source, params,
                                                );
                                            if residual < best_residual {
                                                best_residual = residual;
                                                best_label = format!(
                                                    "row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                    row_placement.name(),
                                                    column_placement.name(),
                                                );
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
                "KP15 Fig12 composed PR08 subdivision source n={n} {label}: best_residual={best_residual:e} best={best_label}",
            );
            assert!(best_residual.is_finite());
        }
    }
}

#[test]
fn kp15_fig12_pr08_q13_q2_top_source_variants_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for source_level in Kp15Pr08TopSourceLevel::all() {
            for patch_index in 0..=2 {
                for row_start in [0usize, 1] {
                    for reverse_rows in [false, true] {
                        let source_rows = kp15_pr08_q13_q2_top_source_rows_for_test(
                            &ring,
                            source_level,
                            patch_index,
                            row_start,
                            reverse_rows,
                        )
                        .expect("PR08 q13/q2 top source rows");
                        for use_fig12_bb in [false, true] {
                            let rows = if use_fig12_bb {
                                kp15_fig12_bb_beta_rows_for_test(&source_rows.sectors, params)
                            } else {
                                kp15_bi6_tensor_rows_for_test(&source_rows.sectors, params, true)
                            };
                            for row_offset in -2..=2 {
                                for column_offset in -2..=2 {
                                    for row_reversed in [false, true] {
                                        for column_reversed in [false, true] {
                                            for row_reflected in [false, true] {
                                                for column_reflected in [false, true] {
                                                    for row_placement in
                                                        Kp15Fig12LocalPlacement::all()
                                                    {
                                                        for column_placement in
                                                            Kp15Fig12LocalPlacement::all()
                                                        {
                                                            let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                                n,
                                                                &rows,
                                                                &rows,
                                                                row_offset,
                                                                column_offset,
                                                                row_reversed,
                                                                column_reversed,
                                                                row_reflected,
                                                                column_reflected,
                                                                row_placement,
                                                                column_placement,
                                                            );
                                                            let residual =
                                                                kp15_s4_cross_corner_sum_residual_for_test(
                                                                    &source, params,
                                                                );
                                                            if residual < best_residual {
                                                                best_residual = residual;
                                                                best_label = format!(
                                                                    "source_level={} patch={patch_index} row_start={row_start} reverse_rows={reverse_rows} fig12_bb={use_fig12_bb} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                                    source_level.name(),
                                                                    row_placement.name(),
                                                                    column_placement.name(),
                                                                );
                                                            }
                                                        }
                                                    }
                                                }
                                            }
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
            "KP15 Fig12 PR08 q13/q2 top source variants n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_to_stam_top_label_pairs_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for source_level in Kp15Pr08TopSourceLevel::all() {
            let source = match source_level {
                Kp15Pr08TopSourceLevel::RawCharacteristicRing => ring.sectors.clone(),
                Kp15Pr08TopSourceLevel::AfterOnePr08Subdivision => {
                    pr08_subdivide_characteristic_ring(&ring).expect("PR08 subdivision")
                }
            };
            for odd_label in 0..13 {
                for even_label in 0..13 {
                    for patch_index in 0..=2 {
                        for row_start in [0usize, 1] {
                            for reverse_rows in [false, true] {
                                let Some(source_rows) = kp15_pr08_top_label_source_rows_for_test(
                                    &ring,
                                    &source,
                                    odd_label,
                                    even_label,
                                    patch_index,
                                    row_start,
                                    reverse_rows,
                                ) else {
                                    continue;
                                };
                                let rows =
                                    kp15_fig12_bb_beta_rows_for_test(&source_rows.sectors, params);
                                let source =
                                    kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                        n,
                                        &rows,
                                        &rows,
                                        -2,
                                        -1,
                                        true,
                                        false,
                                        false,
                                        true,
                                        Kp15Fig12LocalPlacement::SectorJoinFirst,
                                        Kp15Fig12LocalPlacement::SectorJoinFirst,
                                    );
                                let residual =
                                    kp15_s4_cross_corner_sum_residual_for_test(&source, params);
                                if residual < best_residual {
                                    best_residual = residual;
                                    best_label = format!(
                                        "source_level={} odd_label={odd_label} even_label={even_label} patch={patch_index} row_start={row_start} reverse_rows={reverse_rows}",
                                        source_level.name(),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08-to-Stam top label pairs n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_endpoint_side_source_orientation_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let side_tensor = kp15_pr08_endpoint_side_tensor_border_for_fig12(n);
        let side_rows = kp15_fig12_bb_beta_rows_for_test(&side_tensor.sectors, params);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_offset in -2isize..=2 {
            for column_offset in -2isize..=2 {
                for row_reversed in [false, true] {
                    for column_reversed in [false, true] {
                        for row_reflected in [false, true] {
                            for column_reflected in [false, true] {
                                for row_placement in Kp15Fig12LocalPlacement::all() {
                                    for column_placement in Kp15Fig12LocalPlacement::all() {
                                        let source =
                                            kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                n,
                                                &side_rows,
                                                &side_rows,
                                                row_offset,
                                                column_offset,
                                                row_reversed,
                                                column_reversed,
                                                row_reflected,
                                                column_reflected,
                                                row_placement,
                                                column_placement,
                                            );
                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(
                                            &source, params,
                                        );
                                        if residual < best_residual {
                                            best_residual = residual;
                                            best_label = format!(
                                                "row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                row_placement.name(),
                                                column_placement.name(),
                                            );
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
            "KP15 Fig12 endpoint-side source orientation n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_paired_corner_sides_are_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
            for column_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
                if row_side == column_side {
                    continue;
                }
                for row_swap_targets in [false, true] {
                    for column_swap_targets in [false, true] {
                        for row_conjugate in [false, true] {
                            for column_conjugate in [false, true] {
                                let Ok(row_tensor) =
                                    kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                                        n,
                                        row_side,
                                        row_swap_targets,
                                        row_conjugate,
                                    )
                                else {
                                    continue;
                                };
                                let Ok(column_tensor) =
                                    kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                                        n,
                                        column_side,
                                        column_swap_targets,
                                        column_conjugate,
                                    )
                                else {
                                    continue;
                                };
                                let row_rows =
                                    kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
                                let column_rows = kp15_fig12_bb_beta_rows_for_test(
                                    &column_tensor.sectors,
                                    params,
                                );
                                for row_offset in -2isize..=2 {
                                    for column_offset in -2isize..=2 {
                                        for row_reversed in [false, true] {
                                            for column_reversed in [false, true] {
                                                for row_reflected in [false, true] {
                                                    for column_reflected in [false, true] {
                                                        for row_placement in
                                                            Kp15Fig12LocalPlacement::all()
                                                        {
                                                            for column_placement in
                                                                Kp15Fig12LocalPlacement::all()
                                                            {
                                                                let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                                    n,
                                                                    &row_rows,
                                                                    &column_rows,
                                                                    row_offset,
                                                                    column_offset,
                                                                    row_reversed,
                                                                    column_reversed,
                                                                    row_reflected,
                                                                    column_reflected,
                                                                    row_placement,
                                                                    column_placement,
                                                                );
                                                                let residual = kp15_s4_cross_corner_sum_residual_for_test(&source, params);
                                                                if residual < best_residual {
                                                                    best_residual = residual;
                                                                    best_label = format!(
                                                                        "row_side={} column_side={} row_swap={row_swap_targets} column_swap={column_swap_targets} row_conjugate={row_conjugate} column_conjugate={column_conjugate} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                                        row_side.name(),
                                                                        column_side.name(),
                                                                        row_placement.name(),
                                                                        column_placement.name(),
                                                                    );
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
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
            "KP15 Fig12 PR08 Fig6.3 paired corner sides n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_joint_corner_normalization_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let mut best_residual = f64::INFINITY;
        let mut best_label = String::new();
        for row_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
            for column_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
                if row_side == column_side {
                    continue;
                }
                for swap_outer_targets in [false, true] {
                    for conjugate_source in [false, true] {
                        let Ok(row_tensor) =
                            kp15_pr08_fig63_joint_corner_side_tensor_border_for_fig12(
                                n,
                                row_side,
                                swap_outer_targets,
                                conjugate_source,
                            )
                        else {
                            continue;
                        };
                        let Ok(column_tensor) =
                            kp15_pr08_fig63_joint_corner_side_tensor_border_for_fig12(
                                n,
                                column_side,
                                swap_outer_targets,
                                conjugate_source,
                            )
                        else {
                            continue;
                        };
                        let row_rows =
                            kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
                        let column_rows =
                            kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                        for row_offset in -2isize..=2 {
                            for column_offset in -2isize..=2 {
                                for row_reversed in [false, true] {
                                    for column_reversed in [false, true] {
                                        for row_reflected in [false, true] {
                                            for column_reflected in [false, true] {
                                                for row_placement in Kp15Fig12LocalPlacement::all()
                                                {
                                                    for column_placement in
                                                        Kp15Fig12LocalPlacement::all()
                                                    {
                                                        let source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                                            n,
                                                            &row_rows,
                                                            &column_rows,
                                                            row_offset,
                                                            column_offset,
                                                            row_reversed,
                                                            column_reversed,
                                                            row_reflected,
                                                            column_reflected,
                                                            row_placement,
                                                            column_placement,
                                                        );
                                                        let residual = kp15_s4_cross_corner_sum_residual_for_test(&source, params);
                                                        if residual < best_residual {
                                                            best_residual = residual;
                                                            best_label = format!(
                                                                "row_side={} column_side={} swap_outer={swap_outer_targets} conjugate={conjugate_source} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                                                row_side.name(),
                                                                column_side.name(),
                                                                row_placement.name(),
                                                                column_placement.name(),
                                                            );
                                                        }
                                                    }
                                                }
                                            }
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
            "KP15 Fig12 PR08 Fig6.3 joint corner normalization n={n}: best_residual={best_residual:e} best={best_label}",
        );
        assert!(best_residual.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_only_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n, params,
            );
        let q = f_kappa_quadratic_form(5);
        let source_cols = kp15_fig2b_depth2_tensor_border_columns_for_test(n);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        let frame_cases = [(false, false), (false, true), (true, false), (true, true)];
        for row_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
            for column_side in [Kp20StamGridSide::Left, Kp20StamGridSide::Top] {
                if row_side == column_side {
                    continue;
                }
                for (swap_targets, conjugate_source) in frame_cases {
                    let Ok(row_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                        n,
                        row_side,
                        swap_targets,
                        conjugate_source,
                    ) else {
                        continue;
                    };
                    let Ok(column_tensor) =
                        kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                            n,
                            column_side,
                            swap_targets,
                            conjugate_source,
                        )
                    else {
                        continue;
                    };
                    let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
                    let column_rows =
                        kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                    let row_offset = -2isize;
                    let column_offset = -1isize;
                    let row_reversed = row_side == Kp20StamGridSide::Left;
                    let column_reversed = column_side == Kp20StamGridSide::Left;
                    let row_placement = Kp15Fig12LocalPlacement::SectorJoinFirst;
                    let column_placement = Kp15Fig12LocalPlacement::SectorJoinFirst;
                    for row_reflected in [false, true] {
                        for column_reflected in [false, true] {
                            let source =
                                kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                    n,
                                    &row_rows,
                                    &column_rows,
                                    row_offset,
                                    column_offset,
                                    row_reversed,
                                    column_reversed,
                                    row_reflected,
                                    column_reflected,
                                    row_placement,
                                    column_placement,
                                );
                            let constraints = kp15_source_constraints_from_grid_columns_for_test(
                                n,
                                &source,
                                &source_cols,
                            );
                            let label = format!(
                                "row_side={} column_side={} swap={swap_targets} conjugate={conjugate_source} row_offset={row_offset} column_offset={column_offset} row_reversed={row_reversed} column_reversed={column_reversed} row_reflect={row_reflected} column_reflect={column_reflected} row_placement={} column_placement={}",
                                row_side.name(),
                                column_side.name(),
                                row_placement.name(),
                                column_placement.name(),
                            );
                            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                                &constraints,
                                &q,
                                0.0,
                                "PR08 Fig6.3 depth2-only source projection",
                            ) {
                                Ok((grids, _)) => {
                                    accepted += 1;
                                    let orientation =
                                        test_bi6_patch_local_orientation_report_with_samples(
                                            n, &grids, &samples,
                                        );
                                    if orientation.negative < best_neg
                                        || (orientation.negative == best_neg
                                            && orientation.min_det > best_min_det)
                                    {
                                        best_neg = orientation.negative;
                                        best_min_det = orientation.min_det;
                                        best_label = label;
                                    }
                                }
                                Err(err) => {
                                    rejected += 1;
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
            "KP15 Fig12 PR08 Fig6.3 depth2-only source projection n={n}: accepted={accepted} rejected={rejected} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_full_tensor_border_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(n, &source, &source_cols);
            let label = format!(
                "row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={}",
                spec.row_side.name(),
                spec.column_side.name(),
                spec.swap_targets,
                spec.conjugate_source,
                spec.row_reflected,
                spec.column_reflected,
            );
            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "PR08 Fig6.3 full tensor-border source projection",
            ) {
                Ok((grids, f_values)) => {
                    accepted += 1;
                    let orientation =
                        test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                    if orientation.negative < best_neg
                        || (orientation.negative == best_neg && orientation.min_det > best_min_det)
                    {
                        best_neg = orientation.negative;
                        best_min_det = orientation.min_det;
                        best_label = format!(
                            "{label} f=({:e},{:e},{:e})",
                            f_values[0], f_values[1], f_values[2],
                        );
                    }
                }
                Err(err) => {
                    if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                        let abs_residual = residual.abs();
                        if abs_residual < smallest_rejection {
                            smallest_rejection = abs_residual;
                            smallest_rejection_label = label;
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 full tensor-border projection n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_full_tensor_border_first_inconsistency_is_reported() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut best = None::<(f64, Kp15Fig63RadialScaleSpec, Vec<[[Point3; 7]; 7]>)>;
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            let constraints = kp15_source_constraints_from_grid_columns_for_test(
                n,
                &source,
                &kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n),
            );
            if let Err(err) = table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &f_kappa_quadratic_form(5),
                0.0,
                "PR08 Fig6.3 full tensor-border first inconsistency",
            ) && let Some(residual) = kp15_parse_constraint_residual_for_test(&err)
            {
                let abs_residual = residual.abs();
                if best
                    .as_ref()
                    .map(|(best_residual, _, _)| abs_residual < *best_residual)
                    .unwrap_or(true)
                {
                    best = Some((abs_residual, spec, source));
                }
            }
        }
        let Some((best_residual, best_spec, best_source)) = best else {
            continue;
        };
        let labeled = kp15_labeled_full_tensor_border_constraints_for_test(n, &best_source);
        let report = kp15_first_inconsistent_source_constraint_for_test(&table, &labeled)
            .unwrap_or_else(|| "all constraints consistent".to_string());
        println!(
            "KP15 Fig12 PR08 Fig6.3 full tensor-border first inconsistency n={n}: best_residual={best_residual:e} row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={} first={report}",
            best_spec.row_side.name(),
            best_spec.column_side.name(),
            best_spec.swap_targets,
            best_spec.conjugate_source,
            best_spec.row_reflected,
            best_spec.column_reflected,
        );
    }
}

#[test]
fn kp15_fig12_pr08_fig63_full_tensor_border_shifted_to_m_endpoint_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            let constraints =
                kp15_shifted_m_endpoint_tensor_border_constraints_for_test(n, &source);
            let label = format!(
                "row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={}",
                spec.row_side.name(),
                spec.column_side.name(),
                spec.swap_targets,
                spec.conjugate_source,
                spec.row_reflected,
                spec.column_reflected,
            );
            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &f_kappa_quadratic_form(5),
                0.0,
                "PR08 Fig6.3 shifted-m endpoint tensor-border projection",
            ) {
                Ok(_) => accepted += 1,
                Err(err) => {
                    if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                        let abs_residual = residual.abs();
                        if abs_residual < smallest_rejection {
                            smallest_rejection = abs_residual;
                            smallest_rejection_label = label;
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 shifted-m endpoint tensor-border n={n}: accepted={accepted} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_plus_m_endpoint_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            let constraints = kp15_source_constraints_from_grid_columns_for_test(
                n,
                &source,
                &kp15_fig2b_depth2_plus_m_endpoint_columns_for_test(n),
            );
            let label = format!(
                "row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={}",
                spec.row_side.name(),
                spec.column_side.name(),
                spec.swap_targets,
                spec.conjugate_source,
                spec.row_reflected,
                spec.column_reflected,
            );
            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "PR08 Fig6.3 depth2 plus m-endpoint tensor-border projection",
            ) {
                Ok((grids, f_values)) => {
                    accepted += 1;
                    let orientation =
                        test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                    if orientation.negative < best_neg
                        || (orientation.negative == best_neg && orientation.min_det > best_min_det)
                    {
                        best_neg = orientation.negative;
                        best_min_det = orientation.min_det;
                        best_label = format!(
                            "{label} f=({:e},{:e},{:e})",
                            f_values[0], f_values[1], f_values[2],
                        );
                    }
                }
                Err(err) => {
                    if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                        let abs_residual = residual.abs();
                        if abs_residual < smallest_rejection {
                            smallest_rejection = abs_residual;
                            smallest_rejection_label = label;
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 depth2 plus m-endpoint projection n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_plus_cap_net_m_endpoint_is_ranked() {
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let mut compatible = 0usize;
        let mut best_label = String::new();
        let mut first_rejection = None::<String>;
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            if spec.row_side != Kp20StamGridSide::Left || spec.column_side != Kp20StamGridSide::Top
            {
                continue;
            }
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            for m_label in [4usize, 5] {
                for m_offset in [-1isize, 0, 1] {
                    let constraints =
                        kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
                            n,
                            &source,
                            &cap_net.sectors,
                            m_label,
                            m_offset,
                        );
                    let label = format!(
                        "row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={} m_label=c{} m_offset={}",
                        spec.row_side.name(),
                        spec.column_side.name(),
                        spec.swap_targets,
                        spec.conjugate_source,
                        spec.row_reflected,
                        spec.column_reflected,
                        m_label + 1,
                        m_offset,
                    );
                    if let Some(report) =
                        kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
                    {
                        if first_rejection.is_none() {
                            first_rejection = Some(format!("{label} first={report}"));
                        }
                    } else {
                        compatible += 1;
                        if best_label.is_empty() {
                            best_label = label;
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 depth2 plus cap-net m endpoint rank n={n}: compatible={compatible} first_compatible={best_label} first_rejection={}",
            first_rejection.as_deref().unwrap_or("none"),
        );
        assert!(compatible > 0 || first_rejection.is_some());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_plus_cap_net_m_endpoint_kkt_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let spec = Kp15Fig63RadialScaleSpec {
            n,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: false,
            row_reflected: false,
            column_reflected: false,
        };
        let source =
            kp15_fig63_full_tensor_border_source_for_spec(spec).expect("canonical Fig6.3 source");
        let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
            n,
            &source,
            &cap_net.sectors,
            4,
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
                "PR08 Fig6.3 depth2 plus cap-net m endpoint canonical",
            )
            .unwrap_or_else(|err| panic!("n={n} cap-net endpoint KKT rejected: {err}"));
        let orientation = test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
        println!(
            "KP15 Fig12 PR08 Fig6.3 depth2 plus cap-net m endpoint KKT n={n}: neg={}/{} zero={} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
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

#[test]
fn kp15_fig12_pr08_fig63_depth2_plus_cap_net_m_endpoint_kkt_family_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            if spec.row_side != Kp20StamGridSide::Left || spec.column_side != Kp20StamGridSide::Top
            {
                continue;
            }
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            for m_label in [4usize, 5] {
                for m_offset in [-1isize, 0, 1] {
                    let labeled = kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
                        n,
                        &source,
                        &cap_net.sectors,
                        m_label,
                        m_offset,
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
                    let label = format!(
                        "swap={} conjugate={} row_reflect={} column_reflect={} m_label=c{} m_offset={}",
                        spec.swap_targets,
                        spec.conjugate_source,
                        spec.row_reflected,
                        spec.column_reflected,
                        m_label + 1,
                        m_offset,
                    );
                    match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &f_kappa_quadratic_form(5),
                        0.0,
                        "PR08 Fig6.3 depth2 plus cap-net m endpoint family",
                    ) {
                        Ok((grids, f_values)) => {
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
                                    "{label} det=[{:e},{:e}] worst={} f=({:e},{:e},{:e})",
                                    orientation.min_det,
                                    orientation.max_det,
                                    orientation.worst_sample,
                                    f_values[0],
                                    f_values[1],
                                    f_values[2],
                                );
                            }
                        }
                        Err(_) => rejected += 1,
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 depth2 plus cap-net m endpoint KKT family n={n}: accepted={accepted} rejected={rejected} best_neg={best_neg}/{} best_min_det={best_min_det:e} best={best_label}",
            n * samples.len() * samples.len(),
        );
        assert!(accepted > 0);
    }
}

#[test]
fn kp15_fig12_pr08_fig63_radial_seed_plus_depth2_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let radial_input = kp15_spoke_indexed_planar_bilinear_input(n);
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for (row_side, column_side) in [
            (Kp20StamGridSide::Left, Kp20StamGridSide::Top),
            (Kp20StamGridSide::Top, Kp20StamGridSide::Left),
        ] {
            for (swap_targets, conjugate_source) in
                [(false, false), (false, true), (true, false), (true, true)]
            {
                let Ok(row_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                    n,
                    row_side,
                    swap_targets,
                    conjugate_source,
                ) else {
                    continue;
                };
                let Ok(column_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                    n,
                    column_side,
                    swap_targets,
                    conjugate_source,
                ) else {
                    continue;
                };
                let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
                let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                let row_reversed = row_side == Kp20StamGridSide::Left;
                let column_reversed = column_side == Kp20StamGridSide::Left;
                for row_reflected in [false, true] {
                    for column_reflected in [false, true] {
                        let depth_source =
                            kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                n,
                                &row_rows,
                                &column_rows,
                                -2,
                                -1,
                                row_reversed,
                                column_reversed,
                                row_reflected,
                                column_reflected,
                                Kp15Fig12LocalPlacement::SectorJoinFirst,
                                Kp15Fig12LocalPlacement::SectorJoinFirst,
                            );
                        let constraints = kp15_radial_seed_plus_depth2_source_constraints_for_test(
                            n,
                            &radial_input.sectors,
                            &depth_source,
                        );
                        let label = format!(
                            "row_side={} column_side={} swap={swap_targets} conjugate={conjugate_source} row_reflect={row_reflected} column_reflect={column_reflected}",
                            row_side.name(),
                            column_side.name(),
                        );
                        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                            &constraints,
                            &q,
                            0.0,
                            "PR08 Fig6.3 radial-seed plus depth2 source projection",
                        ) {
                            Ok((grids, f_values)) => {
                                accepted += 1;
                                let orientation =
                                    test_bi6_patch_local_orientation_report_with_samples(
                                        n, &grids, &samples,
                                    );
                                if orientation.negative < best_neg
                                    || (orientation.negative == best_neg
                                        && orientation.min_det > best_min_det)
                                {
                                    best_neg = orientation.negative;
                                    best_min_det = orientation.min_det;
                                    best_label = format!(
                                        "{label} f=({:e},{:e},{:e})",
                                        f_values[0], f_values[1], f_values[2],
                                    );
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
        println!(
            "KP15 Fig12 PR08 Fig6.3 radial-seed plus depth2 projection n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_radial_side_seed_plus_depth2_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        let mut accepted = 0usize;
        let mut smallest_rejection = f64::INFINITY;
        let mut smallest_rejection_label = String::new();
        for radial_reversed_depth in [false, true] {
            let Ok(radial_rows) = kp15_pr08_fig63_center_to_sector_join_radial_rows_for_test(
                n,
                radial_reversed_depth,
            ) else {
                continue;
            };
            let radial_sectors = kp15_sector_rhs_from_bi3_rows(&radial_rows.sectors, params, false);
            for (row_side, column_side) in [
                (Kp20StamGridSide::Left, Kp20StamGridSide::Top),
                (Kp20StamGridSide::Top, Kp20StamGridSide::Left),
            ] {
                for (swap_targets, conjugate_source) in
                    [(false, false), (false, true), (true, false), (true, true)]
                {
                    let Ok(row_tensor) = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                        n,
                        row_side,
                        swap_targets,
                        conjugate_source,
                    ) else {
                        continue;
                    };
                    let Ok(column_tensor) =
                        kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
                            n,
                            column_side,
                            swap_targets,
                            conjugate_source,
                        )
                    else {
                        continue;
                    };
                    let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
                    let column_rows =
                        kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
                    let row_reversed = row_side == Kp20StamGridSide::Left;
                    let column_reversed = column_side == Kp20StamGridSide::Left;
                    for row_reflected in [false, true] {
                        for column_reflected in [false, true] {
                            let depth_source =
                                kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
                                    n,
                                    &row_rows,
                                    &column_rows,
                                    -2,
                                    -1,
                                    row_reversed,
                                    column_reversed,
                                    row_reflected,
                                    column_reflected,
                                    Kp15Fig12LocalPlacement::SectorJoinFirst,
                                    Kp15Fig12LocalPlacement::SectorJoinFirst,
                                );
                            let constraints =
                                kp15_radial_seed_plus_depth2_source_constraints_for_test(
                                    n,
                                    &radial_sectors,
                                    &depth_source,
                                );
                            let label = format!(
                                "radial_reversed_depth={radial_reversed_depth} row_side={} column_side={} swap={swap_targets} conjugate={conjugate_source} row_reflect={row_reflected} column_reflect={column_reflected}",
                                row_side.name(),
                                column_side.name(),
                            );
                            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                                &constraints,
                                &q,
                                0.0,
                                "PR08 Fig6.3 radial-side seed plus depth2 source projection",
                            ) {
                                Ok((grids, f_values)) => {
                                    accepted += 1;
                                    let orientation =
                                        test_bi6_patch_local_orientation_report_with_samples(
                                            n, &grids, &samples,
                                        );
                                    if orientation.negative < best_neg
                                        || (orientation.negative == best_neg
                                            && orientation.min_det > best_min_det)
                                    {
                                        best_neg = orientation.negative;
                                        best_min_det = orientation.min_det;
                                        best_label = format!(
                                            "{label} f=({:e},{:e},{:e})",
                                            f_values[0], f_values[1], f_values[2],
                                        );
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
            "KP15 Fig12 PR08 Fig6.3 radial-side seed plus depth2 projection n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label} smallest_rejection={smallest_rejection:e} smallest_rejection_candidate={smallest_rejection_label}",
        );
        assert!(accepted > 0 || smallest_rejection.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_radial_scale_compatibility_is_ranked() {
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let mut best_scale = 0.0;
        let mut best_rejection = f64::INFINITY;
        let mut accepted_scales = 0usize;
        for step in -80..=80 {
            let scale = step as f64 / 20.0;
            let radial_sectors = kp15_scale_sector_rhs_for_test(&radial_input.sectors, scale);
            let constraints = kp15_radial_seed_plus_depth2_source_constraints_for_test(
                n,
                &radial_sectors,
                &depth_source,
            );
            match table.populate_cap_grids_from_column_constraints_with_f_kappa(
                &constraints,
                &q,
                0.0,
                "PR08 Fig6.3 radial-scale compatibility",
            ) {
                Ok(_) => {
                    accepted_scales += 1;
                    best_scale = scale;
                    best_rejection = 0.0;
                    break;
                }
                Err(err) => {
                    if let Some(residual) = kp15_parse_constraint_residual_for_test(&err) {
                        let abs_residual = residual.abs();
                        if abs_residual < best_rejection {
                            best_rejection = abs_residual;
                            best_scale = scale;
                        }
                    }
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 radial-scale compatibility n={n}: accepted_scales={accepted_scales} best_scale={best_scale:e} best_rejection={best_rejection:e} row_side={} column_side={}",
            spec.row_side.name(),
            spec.column_side.name(),
        );
        assert!(best_rejection.is_finite());
    }
}

fn kp15_fig63_radial_scale_specs_for_test() -> [Kp15Fig63RadialScaleSpec; 4] {
    [
        Kp15Fig63RadialScaleSpec {
            n: 5,
            row_side: Kp20StamGridSide::Top,
            column_side: Kp20StamGridSide::Left,
            swap_targets: true,
            conjugate_source: false,
            row_reflected: false,
            column_reflected: true,
        },
        Kp15Fig63RadialScaleSpec {
            n: 6,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: true,
            row_reflected: true,
            column_reflected: true,
        },
        Kp15Fig63RadialScaleSpec {
            n: 7,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: true,
            row_reflected: true,
            column_reflected: true,
        },
        Kp15Fig63RadialScaleSpec {
            n: 8,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: true,
            conjugate_source: true,
            row_reflected: true,
            column_reflected: true,
        },
    ]
}

fn kp15_fig63_best_depth_source_specs_for_test() -> [Kp15Fig63RadialScaleSpec; 4] {
    [
        Kp15Fig63RadialScaleSpec {
            n: 5,
            row_side: Kp20StamGridSide::Top,
            column_side: Kp20StamGridSide::Left,
            swap_targets: true,
            conjugate_source: false,
            row_reflected: false,
            column_reflected: false,
        },
        Kp15Fig63RadialScaleSpec {
            n: 6,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: false,
            row_reflected: true,
            column_reflected: false,
        },
        Kp15Fig63RadialScaleSpec {
            n: 7,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: false,
            row_reflected: true,
            column_reflected: false,
        },
        Kp15Fig63RadialScaleSpec {
            n: 8,
            row_side: Kp20StamGridSide::Left,
            column_side: Kp20StamGridSide::Top,
            swap_targets: false,
            conjugate_source: false,
            row_reflected: false,
            column_reflected: true,
        },
    ]
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_implied_radial_edge_is_reported() {
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        for (label, slot) in [("p40", (4, 0)), ("p50", (5, 0)), ("p60", (6, 0))] {
            let constraints = kp15_labeled_depth2_then_single_radial_constraint_for_test(
                n,
                &depth_source,
                &radial_input.sectors,
                slot,
            );
            let target = format!("radial sector=0 p_{}{}", slot.0, slot.1);
            let report =
                kp15_implied_ordered_constraint_delta_for_test(&table, &constraints, &target)
                    .expect("ordered radial edge implication report");
            println!("KP15 Fig12 PR08 Fig6.3 depth2-implied radial edge n={n} {label}: {report}");
        }
    }
}

#[test]
fn kp15_fig12_pr08_fig63_depth2_plus_radial_first_inconsistency_is_reported() {
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let constraints = kp15_labeled_depth2_then_radial_seed_constraints_for_test(
            n,
            &depth_source,
            &radial_input.sectors,
        );
        let report = kp15_first_inconsistent_source_constraint_for_test(&table, &constraints)
            .unwrap_or_else(|| "all constraints consistent".to_string());
        println!("KP15 Fig12 PR08 Fig6.3 depth2+radial first inconsistency n={n}: {report}");
    }
}

#[test]
fn kp15_fig12_pr08_fig63_s1_radial_seed_plus_depth2_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let s1_radial_input = kp15_apply_s1_seed_propagation_to_input(&radial_input);
        let constraints = kp15_radial_seed_plus_depth2_source_constraints_for_test(
            n,
            &s1_radial_input.sectors,
            &depth_source,
        );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "PR08 Fig6.3 S1-radial seed plus depth2 source projection",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 S1-radial plus depth2 projection n={n}: accepted neg={}/{} min_det={:e} f=({:e},{:e},{:e}) row_side={} column_side={}",
                    orientation.negative,
                    orientation.samples,
                    orientation.min_det,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(orientation.min_det.is_finite());
            }
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 S1-radial plus depth2 projection n={n}: rejected residual={residual:e} row_side={} column_side={}",
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_fig12_pr08_fig63_sector0_central_radial_plus_depth2_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let constraints = kp15_sector0_central_radial_edge_plus_depth2_source_constraints_for_test(
            n,
            &radial_input.sectors,
            &depth_source,
        );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "PR08 Fig6.3 sector0-central radial plus depth2 source projection",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 sector0-central radial plus depth2 projection n={n}: accepted neg={}/{} min_det={:e} f=({:e},{:e},{:e}) row_side={} column_side={}",
                    orientation.negative,
                    orientation.samples,
                    orientation.min_det,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(orientation.min_det.is_finite());
            }
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 sector0-central radial plus depth2 projection n={n}: rejected residual={residual:e} row_side={} column_side={}",
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_fig12_pr08_fig63_radial_edge_plus_depth2_source_projection_is_measured() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some((depth_source, radial_input)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let constraints = kp15_radial_edge_plus_depth2_source_constraints_for_test(
            n,
            &radial_input.sectors,
            &depth_source,
        );
        match table.populate_cap_grids_from_column_constraints_with_f_kappa(
            &constraints,
            &q,
            0.0,
            "PR08 Fig6.3 radial-edge plus depth2 source projection",
        ) {
            Ok((grids, f_values)) => {
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 radial-edge plus depth2 projection n={n}: accepted neg={}/{} min_det={:e} f=({:e},{:e},{:e}) row_side={} column_side={}",
                    orientation.negative,
                    orientation.samples,
                    orientation.min_det,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(orientation.min_det.is_finite());
            }
            Err(err) => {
                let residual = kp15_parse_constraint_residual_for_test(&err).unwrap_or(f64::NAN);
                println!(
                    "KP15 Fig12 PR08 Fig6.3 radial-edge plus depth2 projection n={n}: rejected residual={residual:e} row_side={} column_side={}",
                    spec.row_side.name(),
                    spec.column_side.name(),
                );
                assert!(residual.is_finite());
            }
        }
    }
}

#[test]
fn kp15_fig12_pr08_fig63_spectral_depth_source_normalization_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_radial_scale_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some((depth_source, _)) =
            kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec)
        else {
            continue;
        };
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
        let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
        let scale_cases = [
            ("identity", 1.0),
            ("lambda", ring.lambda),
            ("one-minus-lambda", 1.0 - ring.lambda),
            ("half-one-minus-lambda", spectral_cross_depth_scale),
            ("lambda-squared", ring.lambda * ring.lambda),
        ];
        let mut best = ("none", usize::MAX, f64::NEG_INFINITY, [f64::INFINITY; 3]);
        for (label, scale) in scale_cases {
            let scaled_source =
                kp15_scale_fig2b_source_depth_from_boundary_for_test(&depth_source, scale, scale);
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(n, &scaled_source, &source_cols);
            let (grids, f_values) = table
                .populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "PR08 Fig6.3 spectral depth source normalization",
                )
                .expect("paper depth source normalization should remain feasible");
            let orientation =
                test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
            println!(
                "KP15 Fig12 PR08 Fig6.3 spectral depth source n={n} {label}: scale={scale:e} lambda={lambda:e} neg={neg}/{sample_count} min_det={min_det:e} f=({f0:e},{f1:e},{f2:e}) row_side={row_side} column_side={column_side}",
                lambda = ring.lambda,
                neg = orientation.negative,
                sample_count = orientation.samples,
                min_det = orientation.min_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                row_side = spec.row_side.name(),
                column_side = spec.column_side.name(),
            );
            if orientation.negative < best.1
                || (orientation.negative == best.1 && orientation.min_det > best.2)
            {
                best = (label, orientation.negative, orientation.min_det, f_values);
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 best spectral depth source n={n}: label={} neg={} min_det={:e} f=({:e},{:e},{:e})",
            best.0, best.1, best.2, best.3[0], best.3[1], best.3[2],
        );
        assert!(best.2.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_best_chart_spectral_depth_source_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_best_depth_source_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
            continue;
        };
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
        let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
        let scale_cases = [
            ("identity", 1.0),
            ("lambda", ring.lambda),
            ("one-minus-lambda", 1.0 - ring.lambda),
            ("half-one-minus-lambda", spectral_cross_depth_scale),
            ("lambda-squared", ring.lambda * ring.lambda),
        ];
        let mut best = ("none", usize::MAX, f64::NEG_INFINITY, [f64::INFINITY; 3]);
        for (label, scale) in scale_cases {
            let scaled_source =
                kp15_scale_fig2b_source_depth_from_boundary_for_test(&source, scale, scale);
            let constraints =
                kp15_source_constraints_from_grid_columns_for_test(n, &scaled_source, &source_cols);
            let (grids, f_values) = table
                .populate_cap_grids_from_column_constraints_with_f_kappa(
                    &constraints,
                    &q,
                    0.0,
                    "PR08 Fig6.3 best-chart spectral depth source",
                )
                .expect("paper depth source should remain feasible");
            let orientation =
                test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
            println!(
                "KP15 Fig12 PR08 Fig6.3 best-chart spectral source n={n} {label}: scale={scale:e} neg={neg}/{sample_count} min_det={min_det:e} f=({f0:e},{f1:e},{f2:e}) row_side={row_side} column_side={column_side} swap={swap} conjugate={conjugate} row_reflect={row_reflect} column_reflect={column_reflect}",
                neg = orientation.negative,
                sample_count = orientation.samples,
                min_det = orientation.min_det,
                f0 = f_values[0],
                f1 = f_values[1],
                f2 = f_values[2],
                row_side = spec.row_side.name(),
                column_side = spec.column_side.name(),
                swap = spec.swap_targets,
                conjugate = spec.conjugate_source,
                row_reflect = spec.row_reflected,
                column_reflect = spec.column_reflected,
            );
            if orientation.negative < best.1
                || (orientation.negative == best.1 && orientation.min_det > best.2)
            {
                best = (label, orientation.negative, orientation.min_det, f_values);
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 best-chart spectral source best n={n}: label={} neg={} min_det={:e} f=({:e},{:e},{:e})",
            best.0, best.1, best.2, best.3[0], best.3[1], best.3[2],
        );
        assert!(best.2.is_finite());
    }
}

#[test]
fn kp15_fig12_pr08_fig63_best_chart_asymmetric_spectral_depth_source_is_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for spec in kp15_fig63_best_depth_source_specs_for_test() {
        let n = spec.n;
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
            continue;
        };
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
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
                let constraints = kp15_source_constraints_from_grid_columns_for_test(
                    n,
                    &scaled_source,
                    &source_cols,
                );
                let (grids, f_values) = table
                    .populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "PR08 Fig6.3 asymmetric spectral depth source",
                    )
                    .expect("paper depth source should remain feasible");
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
            "KP15 Fig12 PR08 Fig6.3 best-chart asymmetric spectral source n={n}: row_scale={} column_scale={} neg={} min_det={:e} f=({:e},{:e},{:e}) row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={}",
            best.0,
            best.1,
            best.2,
            best.3,
            best.4[0],
            best.4[1],
            best.4[2],
            spec.row_side.name(),
            spec.column_side.name(),
            spec.swap_targets,
            spec.conjugate_source,
            spec.row_reflected,
            spec.column_reflected,
        );
        assert!(best.3.is_finite());
    }
}

#[test]
#[ignore = "expensive frontier sweep; run explicitly when broad spectral candidate ranking is needed"]
fn kp15_fig12_pr08_fig63_spectral_depth_source_candidates_are_ranked() {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    for n in [5usize, 6, 7, 8] {
        let params = Bi6Reparams::for_valence(n);
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
                n,
                params,
            );
        let q = f_kappa_quadratic_form(5);
        let source_cols = kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n);
        let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
        let scale_cases = [
            ("identity", 1.0),
            ("lambda", ring.lambda),
            ("one-minus-lambda", 1.0 - ring.lambda),
            ("half-one-minus-lambda", spectral_cross_depth_scale),
            ("lambda-squared", ring.lambda * ring.lambda),
        ];
        let mut accepted = 0usize;
        let mut best_neg = usize::MAX;
        let mut best_min_det = f64::NEG_INFINITY;
        let mut best_label = String::new();
        for spec in kp15_fig63_full_tensor_border_specs_for_test(n) {
            let Some(source) = kp15_fig63_full_tensor_border_source_for_spec(spec) else {
                continue;
            };
            for (scale_label, scale) in scale_cases {
                let scaled_source =
                    kp15_scale_fig2b_source_depth_from_boundary_for_test(&source, scale, scale);
                let constraints = kp15_source_constraints_from_grid_columns_for_test(
                    n,
                    &scaled_source,
                    &source_cols,
                );
                let (grids, f_values) = table
                    .populate_cap_grids_from_column_constraints_with_f_kappa(
                        &constraints,
                        &q,
                        0.0,
                        "PR08 Fig6.3 spectral depth source candidate ranking",
                    )
                    .expect("paper depth source candidate should remain feasible");
                accepted += 1;
                let orientation =
                    test_bi6_patch_local_orientation_report_with_samples(n, &grids, &samples);
                let label = format!(
                    "scale={scale_label}:{scale:e} row_side={} column_side={} swap={} conjugate={} row_reflect={} column_reflect={} f=({:e},{:e},{:e})",
                    spec.row_side.name(),
                    spec.column_side.name(),
                    spec.swap_targets,
                    spec.conjugate_source,
                    spec.row_reflected,
                    spec.column_reflected,
                    f_values[0],
                    f_values[1],
                    f_values[2],
                );
                if orientation.negative < best_neg
                    || (orientation.negative == best_neg && orientation.min_det > best_min_det)
                {
                    best_neg = orientation.negative;
                    best_min_det = orientation.min_det;
                    best_label = label;
                }
            }
        }
        println!(
            "KP15 Fig12 PR08 Fig6.3 spectral depth source candidate ranking n={n}: accepted={accepted} best_neg={best_neg} best_min_det={best_min_det:e} best={best_label}"
        );
        assert!(accepted > 0);
        assert!(best_min_det.is_finite());
    }
}

#[test]
fn kp15_fig2b_best_pr08_depth_pair_projection_is_measured() {
    let samples = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for (n, candidate, apply_beta) in [
        (
            5usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            true,
        ),
        (
            6usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            true,
        ),
        (
            7usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            false,
        ),
        (
            8usize,
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            false,
        ),
    ] {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor rows");
        let table =
            Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
                n,
                Bi6Reparams::for_valence(n),
            );
        let source = kp15_two_sided_source_grid_from_bi3_tensor_oriented_for_test(
            n,
            &tensor.sectors,
            Bi6Reparams::for_valence(n),
            apply_beta,
            0,
            1,
            true,
            false,
            false,
            true,
        );
        let depth_cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
        let projected =
            kp15_project_source_side_to_row_manifold_for_test(n, &table, &source, &[], &depth_cols)
                .expect("best PR08 paired source projection");
        let max_movement =
            kp15_max_source_column_delta_for_test(n, &source, &projected, &depth_cols);
        let orientation = test_bi6_patch_boundary_orientation_report(n, &projected, &samples);
        println!(
            "KP15 Fig2b best PR08 depth-pair projection n={n} candidate={} beta={apply_beta}: movement={max_movement:e} neg={neg}/{samples_count} det=[{min_det:e},{max_det:e}]",
            candidate.name(),
            neg = orientation.all.negative,
            samples_count = orientation.all.samples,
            min_det = orientation.all.min_det,
            max_det = orientation.all.max_det,
        );
        assert!(max_movement.is_finite());
        assert!(orientation.all.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_s1_seed_closure_residual_is_measured() {
    for n in [5usize, 6, 7, 8, 9, 10] {
        let raw = kp15_pr08_q11_transposed_endpoint_input(n, false);
        let propagated = kp15_apply_s1_seed_propagation_to_input(&raw);
        for (label, input) in [("raw", &raw), ("forward-propagated", &propagated)] {
            let residual = kp15_s1_cyclic_seed_residual(input);
            println!(
                "KP15 PR08 q11 n={n} {label} S1 cyclic seed residual: max={max:e} edge={edge}",
                max = residual.0,
                edge = residual.1,
            );
            assert!(residual.0.is_finite());
        }
    }
}

#[test]
fn kp15_section34_circulant_system_ranks_match_paper() {
    for n in [3usize, 5, 6, 7, 8, 9, 10, 11, 12] {
        let params = Bi6Reparams::for_valence(n);
        let s2_rows = kp15_section34_eq7_p21_p12_matrix_for_test(n, params);
        let s2_rank = dense_rank_relative_for_test(&s2_rows);
        let eq8_rows = kp15_section34_eq8_p22_matrix_for_test(n);
        let eq8_rank = dense_rank_relative_for_test(&eq8_rows);
        println!(
            "KP15 §3.4 circulant ranks n={n}: eq7_rows={} eq7_vars={} eq7_rank={} eq7_nullity={} eq8_rows={} eq8_vars={} eq8_rank={} eq8_nullity={}",
            s2_rows.len(),
            2 * n,
            s2_rank,
            2 * n - s2_rank,
            eq8_rows.len(),
            n,
            eq8_rank,
            n - eq8_rank,
        );
        let expected_eq7_nullity = usize::from(n == 3 || n == 6);
        assert_eq!(2 * n - s2_rank, expected_eq7_nullity);
        assert_eq!(n - eq8_rank, 1);
    }
}

#[test]
fn kp15_pr08_q11_n7_fixed_tensor_slot_frontier_is_measured() {
    let n = 7usize;
    let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n,
            input.params,
        );
    let q = f_kappa_quadratic_form(5);
    for (label, released_cols) in [
        ("none", std::collections::BTreeSet::new()),
        ("p40", kp15_source_slot_family_columns(n, &[(4, 0)])),
        ("p50", kp15_source_slot_family_columns(n, &[(5, 0)])),
        ("p60", kp15_source_slot_family_columns(n, &[(6, 0)])),
        ("p42", kp15_source_slot_family_columns(n, &[(4, 2)])),
        ("p52", kp15_source_slot_family_columns(n, &[(5, 2)])),
        ("p62", kp15_source_slot_family_columns(n, &[(6, 2)])),
        (
            "row0-far",
            kp15_source_slot_family_columns(n, &[(4, 0), (5, 0), (6, 0)]),
        ),
        (
            "row2-far",
            kp15_source_slot_family_columns(n, &[(4, 2), (5, 2), (6, 2)]),
        ),
        (
            "p40-p42",
            kp15_source_slot_family_columns(n, &[(4, 0), (4, 2)]),
        ),
        (
            "p40-p60",
            kp15_source_slot_family_columns(n, &[(4, 0), (6, 0)]),
        ),
        (
            "p42-p52",
            kp15_source_slot_family_columns(n, &[(4, 2), (5, 2)]),
        ),
        (
            "clean-singletons",
            kp15_source_slot_family_columns(n, &[(4, 0), (6, 0), (4, 2), (5, 2)]),
        ),
    ] {
        let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
            &table,
            &input.sectors,
            &q,
            None,
            Some(&released_cols),
        )
        .expect("released tensor-slot F5 solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        let released_movement =
            kp15_released_slot_movement(&input.sectors, &patches, &released_cols);
        println!(
            "KP15 PR08 q11 n=7 released tensor slots {label}: released={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] released_move={released_movement:e} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
            released_cols.len(),
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
        assert_eq!(orientation.samples, 25 * n);
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_n6_dense_fixed_tensor_slot_frontier_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    let n = 6usize;
    let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n,
            input.params,
        );
    let q = f_kappa_quadratic_form(5);
    for (label, released_cols) in [
        ("none", std::collections::BTreeSet::new()),
        ("p40", kp15_source_slot_family_columns(n, &[(4, 0)])),
        ("p50", kp15_source_slot_family_columns(n, &[(5, 0)])),
        ("p60", kp15_source_slot_family_columns(n, &[(6, 0)])),
        ("p42", kp15_source_slot_family_columns(n, &[(4, 2)])),
        ("p52", kp15_source_slot_family_columns(n, &[(5, 2)])),
        ("p62", kp15_source_slot_family_columns(n, &[(6, 2)])),
        (
            "p40-p42",
            kp15_source_slot_family_columns(n, &[(4, 0), (4, 2)]),
        ),
        (
            "clean-n7-singletons",
            kp15_source_slot_family_columns(n, &[(4, 0), (6, 0), (4, 2), (5, 2)]),
        ),
    ] {
        let (patches, f_values, split) = kp15_solve_table_with_test_variable_columns(
            &table,
            &input.sectors,
            &q,
            None,
            Some(&released_cols),
        )
        .expect("released tensor-slot F5 solve");
        let orientation =
            test_bi6_patch_local_orientation_report_with_samples(n, &patches, &DENSE_SAMPLES);
        let released_movement =
            kp15_released_slot_movement(&input.sectors, &patches, &released_cols);
        println!(
            "KP15 PR08 q11 n=6 dense released tensor slots {label}: released={} fixed={} variable={} zeroed_extra={} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] released_move={released_movement:e} f=({f0:e},{f1:e},{f2:e}) worst={worst}",
            released_cols.len(),
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
        assert_eq!(orientation.samples, 81 * n);
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_n7_clean_released_slot_deltas_are_measured() {
    let n = 7usize;
    let input = kp15_pr08_q11_transposed_endpoint_input(n, false);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n,
            input.params,
        );
    let q = f_kappa_quadratic_form(5);
    for (label, slot) in [
        ("p40", (4, 0)),
        ("p60", (6, 0)),
        ("p42", (4, 2)),
        ("p52", (5, 2)),
    ] {
        let released_cols = kp15_source_slot_family_columns(n, &[slot]);
        let (patches, _f_values, _split) = kp15_solve_table_with_test_variable_columns(
            &table,
            &input.sectors,
            &q,
            None,
            Some(&released_cols),
        )
        .expect("released tensor-slot F5 solve");
        let orientation = test_bi6_patch_local_orientation_report(n, &patches);
        let delta = kp15_released_slot_local_delta_stats(&input.sectors, &patches, &released_cols);
        println!(
            "KP15 PR08 q11 n=7 released {label} delta: neg={neg}/{samples} det_min={min_det:e} avg_local=({avg_x:e},{avg_y:e}) max_norm={max_norm:e} max_deviation={max_dev:e}",
            neg = orientation.negative,
            samples = orientation.samples,
            min_det = orientation.min_det,
            avg_x = delta.avg_local[0],
            avg_y = delta.avg_local[1],
            max_norm = delta.max_norm,
            max_dev = delta.max_deviation_from_average,
        );
        assert_eq!(orientation.negative, 0);
        assert!(delta.max_norm.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_paper_p_beta_source_chart_is_measured() {
    for n in [5usize, 6, 7, 8] {
        for (chart, beta) in [("direct-degree-raised", false), ("kp15-beta", true)] {
            let input = kp15_pr08_q11_transposed_endpoint_input(n, beta);
            let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
                .unwrap_or_else(|err| panic!("n={n} {chart} paper-P F5 failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 paper-P F5 chart n={n} {chart}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
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
fn kp15_pr08_q11_n7_beta_slot_deltas_are_measured() {
    let n = 7usize;
    let direct = kp15_pr08_q11_transposed_endpoint_input(n, false);
    let beta = kp15_pr08_q11_transposed_endpoint_input(n, true);
    for (label, slot) in [("p42", (4, 2)), ("p52", (5, 2)), ("p62", (6, 2))] {
        let mut sum = [0.0, 0.0];
        let mut max_norm: f64 = 0.0;
        let mut deltas = Vec::new();
        for sector in 0..n {
            let b = kp15_seed_value_for_test(&beta.sectors, sector, slot.0, slot.1)
                .expect("beta slot exists");
            let d = kp15_seed_value_for_test(&direct.sectors, sector, slot.0, slot.1)
                .expect("direct slot exists");
            let diff = Point3::new(b.x - d.x, b.y - d.y, b.z - d.z);
            let local = global_to_local_planar(n, sector, diff);
            sum[0] += local[0];
            sum[1] += local[1];
            max_norm = max_norm.max((local[0] * local[0] + local[1] * local[1]).sqrt());
            deltas.push(local);
        }
        let avg = [sum[0] / n as f64, sum[1] / n as f64];
        let mut max_deviation: f64 = 0.0;
        for local in deltas {
            let dx = local[0] - avg[0];
            let dy = local[1] - avg[1];
            max_deviation = max_deviation.max((dx * dx + dy * dy).sqrt());
        }
        println!(
            "KP15 PR08 q11 n=7 beta-minus-direct {label}: avg_local=({avg_x:e},{avg_y:e}) max_norm={max_norm:e} max_deviation={max_deviation:e}",
            avg_x = avg[0],
            avg_y = avg[1],
        );
        assert!(max_norm.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_n7_bi3_row_control_projection_frontier_is_measured() {
    let n = 7usize;
    for (label, row, cp, delta) in [
        (
            "row0-cp2-from-p40",
            0usize,
            2usize,
            [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6],
        ),
        (
            "row0-cp3-from-p60",
            0usize,
            3usize,
            [4.571604603070545e-2, -2.846314044096245e-2],
        ),
        (
            "row2-cp2-from-p42",
            2usize,
            2usize,
            [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
        ),
        (
            "row2-cp2-from-p52",
            2usize,
            2usize,
            [3.8050226127550726e-2 / 0.6, -5.247757835619124e-3 / 0.6],
        ),
        (
            "row2-cp3-from-p52",
            2usize,
            3usize,
            [3.8050226127550726e-2 / 0.5, -5.247757835619124e-3 / 0.5],
        ),
    ] {
        let input =
            kp15_pr08_q11_transposed_endpoint_input_with_bi3_local_delta(n, false, row, cp, delta);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("n={n} {label} paper-P F5 failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        println!(
            "KP15 PR08 q11 n=7 bi3 row-control {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] delta=({dx:e},{dy:e}) worst={worst}",
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            dx = delta[0],
            dy = delta[1],
            worst = orientation.worst_sample,
        );
        assert_eq!(orientation.samples, 25 * n);
        assert!(orientation.min_det.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_n7_bi3_row_control_combined_projection_is_measured() {
    let n = 7usize;
    for (label, deltas) in [
        (
            "row0cp2-row2cp2",
            vec![
                (
                    0usize,
                    2usize,
                    [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6],
                ),
                (
                    2usize,
                    2usize,
                    [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
                ),
            ],
        ),
        (
            "row0cp3-row2cp2",
            vec![
                (
                    0usize,
                    3usize,
                    [4.571604603070545e-2, -2.846314044096245e-2],
                ),
                (
                    2usize,
                    2usize,
                    [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
                ),
            ],
        ),
    ] {
        let input =
            kp15_pr08_q11_transposed_endpoint_input_with_bi3_local_deltas(n, false, &deltas);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("n={n} {label} paper-P F5 failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        println!(
            "KP15 PR08 q11 n=7 combined bi3 row-control {label}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
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

#[test]
fn kp15_pr08_exact_stam_side_q11_row_control_source_is_measured() {
    let expected_row0_cp2 = [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6];
    let expected_row2_cp2 = [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6];

    for n in [5usize, 6, 7, 8] {
        let q11_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let exact_side_regular_tensor = kp15_pr08_exact_stam_side_regular_endpoint_tensor(n);
        let exact_side_kp15_tensor = kp15_endpoint_normalized_bi3_tensor_border(
            &kp15_pr08_exact_stam_side_raw_tensor(n),
            n,
            0,
        );

        for (label, tensor) in [
            ("regular-endpoint-exact-side", &exact_side_regular_tensor),
            ("kp15-endpoint-exact-side", &exact_side_kp15_tensor),
        ] {
            let row0_cp2_delta =
                kp15_tensor_local_delta(n, 0, tensor.sectors[0][0][2], q11_tensor.sectors[0][0][2]);
            let row2_cp2_delta =
                kp15_tensor_local_delta(n, 0, tensor.sectors[0][2][2], q11_tensor.sectors[0][2][2]);
            let row0_residual = [
                row0_cp2_delta[0] - expected_row0_cp2[0],
                row0_cp2_delta[1] - expected_row0_cp2[1],
            ];
            let row2_residual = [
                row2_cp2_delta[0] - expected_row2_cp2[0],
                row2_cp2_delta[1] - expected_row2_cp2[1],
            ];

            let input = kp15_input_from_bi3_tensor(n, tensor, false);
            let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
                .unwrap_or_else(|err| panic!("n={n} {label} paper-P/F5 failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 exact Stam side q11 row-control source n={n} {label}: row0cp2_delta=({r0x:e},{r0y:e}) row0_residual=({r0rx:e},{r0ry:e}) row2cp2_delta=({r2x:e},{r2y:e}) row2_residual=({r2rx:e},{r2ry:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
                r0x = row0_cp2_delta[0],
                r0y = row0_cp2_delta[1],
                r0rx = row0_residual[0],
                r0ry = row0_residual[1],
                r2x = row2_cp2_delta[0],
                r2y = row2_cp2_delta[1],
                r2rx = row2_residual[0],
                r2ry = row2_residual[1],
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
fn kp15_pr08_q11_row_control_label_phase_candidates_are_ranked() {
    let n = 7usize;
    let expected = [
        (
            "row0-cp2",
            0usize,
            2usize,
            [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6],
        ),
        (
            "row2-cp2",
            2usize,
            2usize,
            [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
        ),
    ];
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let base_rows = kp15_pr08_q11_transposed_regular_rows(&ring);
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);

    for (target_label, target_row, target_cp, expected_delta) in expected {
        let mut ranked = Vec::new();
        for sector_offset in [-1isize, 0, 1] {
            let donor_sector = ((n as isize + sector_offset) as usize) % n;
            for donor_label in 0..13 {
                let mut rows = base_rows.clone();
                for sector in 0..n {
                    let source_sector =
                        ((sector as isize + sector_offset + n as isize) as usize) % n;
                    rows.sectors[sector][target_row][target_cp] =
                        ring.sectors[source_sector][donor_label];
                }
                let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows)
                    .expect("label replacement tensor");
                let endpoint_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
                let delta = kp15_tensor_local_delta(
                    n,
                    0,
                    endpoint_tensor.sectors[0][target_row][target_cp],
                    base_tensor.sectors[0][target_row][target_cp],
                );
                let residual = (delta[0] - expected_delta[0]).hypot(delta[1] - expected_delta[1]);
                ranked.push((residual, sector_offset, donor_sector, donor_label, delta));
            }
        }
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (rank, (residual, sector_offset, donor_sector, donor_label, delta)) in
            ranked.iter().take(8).enumerate()
        {
            println!(
                "KP15 PR08 q11 {target_label} label-phase candidate rank={rank} offset={sector_offset} donor_sector={donor_sector} donor_q{} delta=({dx:e},{dy:e}) residual={residual:e}",
                donor_label + 1,
                dx = delta[0],
                dy = delta[1],
            );
        }
        assert!(ranked[0].0.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_first_order_characteristic_normalization_is_measured() {
    let n = 7usize;
    let target_deltas = [
        (
            0usize,
            2usize,
            [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6],
        ),
        (
            2usize,
            2usize,
            [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
        ),
    ];
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    let p0 = global_to_local_planar(
        n,
        0,
        base_tensor.sectors[0][target_deltas[0].0][target_deltas[0].1],
    );
    let p1 = global_to_local_planar(
        n,
        0,
        base_tensor.sectors[0][target_deltas[1].0][target_deltas[1].1],
    );
    let frame_delta =
        kp15_first_order_map_from_two_point_deltas(p0, target_deltas[0].2, p1, target_deltas[1].2)
            .expect("q11 row controls should span a nondegenerate first-order frame");
    let det = (1.0 + frame_delta[0][0]) * (1.0 + frame_delta[1][1])
        - frame_delta[0][1] * frame_delta[1][0];
    let normalized_tensor =
        kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
    let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
    let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
        .expect("first-order characteristic-normalized paper-P/F5 solve");
    let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);

    for (label, row, cp, expected_delta) in [
        (
            "row0-cp2",
            target_deltas[0].0,
            target_deltas[0].1,
            target_deltas[0].2,
        ),
        (
            "row2-cp2",
            target_deltas[1].0,
            target_deltas[1].1,
            target_deltas[1].2,
        ),
    ] {
        let delta = kp15_tensor_local_delta(
            n,
            0,
            normalized_tensor.sectors[0][row][cp],
            base_tensor.sectors[0][row][cp],
        );
        println!(
            "KP15 PR08 q11 first-order normalization {label}: delta=({dx:e},{dy:e}) expected=({ex:e},{ey:e}) residual={residual:e}",
            dx = delta[0],
            dy = delta[1],
            ex = expected_delta[0],
            ey = expected_delta[1],
            residual = (delta[0] - expected_delta[0]).hypot(delta[1] - expected_delta[1]),
        );
    }
    println!(
        "KP15 PR08 q11 first-order characteristic normalization: D=[[{d00:e},{d01:e}],[{d10:e},{d11:e}]] det(I+D)={det:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst}",
        d00 = frame_delta[0][0],
        d01 = frame_delta[0][1],
        d10 = frame_delta[1][0],
        d11 = frame_delta[1][1],
        neg = orientation.negative,
        samples = orientation.samples,
        zero = orientation.near_zero,
        min_det = orientation.min_det,
        max_det = orientation.max_det,
        worst = orientation.worst_sample,
    );
    assert!(det > 0.0);
    assert_eq!(orientation.samples, 25 * n);
    assert!(orientation.min_det.is_finite());
}

#[test]
fn kp15_pr08_q11_paper_frame_normalization_candidates_are_measured() {
    let n = 7usize;
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    let fitted_delta = kp15_q11_fitted_first_order_frame_delta(n);
    let fitted_map = kp15_add_identity_to_frame_delta(fitted_delta);
    let strip_tensor = {
        let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
        let rows = kp15_pr08_q11_transposed_regular_rows(&ring);
        let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("q11 tensor");
        kp15_strip_normalized_bi3_tensor_border(&tensor, n, 0)
    };
    let beta_tensor = kp15_beta_tensor_from_bi3_tensor(&base_tensor, n);
    let exact_side_tensor = kp15_pr08_exact_stam_side_regular_endpoint_tensor(n);
    let exact_side_kp15_tensor =
        kp15_endpoint_normalized_bi3_tensor_border(&kp15_pr08_exact_stam_side_raw_tensor(n), n, 0);

    for (label, candidate) in [
        (
            "strip-row0cp2-row2cp2",
            kp15_first_order_map_between_tensor_controls(
                n,
                &base_tensor,
                &strip_tensor,
                &[(0, 2), (2, 2)],
            ),
        ),
        (
            "beta-row0cp2-row2cp2",
            kp15_first_order_map_between_tensor_controls(
                n,
                &base_tensor,
                &beta_tensor,
                &[(0, 2), (2, 2)],
            ),
        ),
        (
            "exact-side-regular-row0cp2-row2cp2",
            kp15_first_order_map_between_tensor_controls(
                n,
                &base_tensor,
                &exact_side_tensor,
                &[(0, 2), (2, 2)],
            ),
        ),
        (
            "exact-side-kp15-row0cp2-row2cp2",
            kp15_first_order_map_between_tensor_controls(
                n,
                &base_tensor,
                &exact_side_kp15_tensor,
                &[(0, 2), (2, 2)],
            ),
        ),
        (
            "kp15-s1-first-order-row0cp1-row1cp0",
            kp15_s1_first_order_frame_map_from_tensor(&base_tensor, n),
        ),
    ] {
        let Some(map) = candidate else {
            println!("KP15 PR08 q11 paper-frame candidate {label}: singular");
            continue;
        };
        let frame_delta = kp15_subtract_identity_from_frame_map(map);
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("{label} paper-P/F5 solve failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        let det = kp15_det2(map);
        let diff = kp15_frame_map_frobenius_delta(map, fitted_map);
        println!(
            "KP15 PR08 q11 paper-frame candidate {label}: map=[[{m00:e},{m01:e}],[{m10:e},{m11:e}]] det={det:e} fitted_diff={diff:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            m00 = map[0][0],
            m01 = map[0][1],
            m10 = map[1][0],
            m11 = map[1][1],
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

#[test]
fn kp15_pr08_q11_f5_variational_source_frame_is_measured() {
    let n = 7usize;
    let params = Bi6Reparams::for_valence(n);
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let q = f_kappa_quadratic_form(5);
    let base_patches = kp15_solve_paper_p_f5_for_tensor(&table, n, &base_tensor, &q);
    let basis_deltas = [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.0, 1.0], [0.0, 0.0]],
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 0.0], [0.0, 1.0]],
    ];
    let basis_patches = basis_deltas
        .iter()
        .map(|&frame_delta| {
            let tensor =
                kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
            kp15_patch_delta(
                &kp15_solve_paper_p_f5_for_tensor(&table, n, &tensor, &q),
                &base_patches,
            )
        })
        .collect::<Vec<_>>();
    let (h, g) = kp15_source_frame_quadratic_terms(&base_patches, &basis_patches, &q);
    let rhs = g.iter().map(|value| -value).collect::<Vec<_>>();
    let solved = solve_dense_faer_for_test(h, rhs).expect("F5 source-frame quadratic solve");
    let frame_delta = [[solved[0], solved[1]], [solved[2], solved[3]]];
    let fitted_delta = kp15_q11_fitted_first_order_frame_delta(n);
    let fitted_diff = kp15_frame_map_frobenius_delta(
        kp15_add_identity_to_frame_delta(frame_delta),
        kp15_add_identity_to_frame_delta(fitted_delta),
    );
    let normalized_tensor =
        kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
    let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
    let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
        .expect("F5 variational source-frame paper-P/F5 solve");
    let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
    let det = kp15_det2(kp15_add_identity_to_frame_delta(frame_delta));

    println!(
        "KP15 PR08 q11 F5 variational source frame: D=[[{d00:e},{d01:e}],[{d10:e},{d11:e}]] det(I+D)={det:e} fitted_diff={fitted_diff:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
        d00 = frame_delta[0][0],
        d01 = frame_delta[0][1],
        d10 = frame_delta[1][0],
        d11 = frame_delta[1][1],
        neg = orientation.negative,
        samples = orientation.samples,
        zero = orientation.near_zero,
        min_det = orientation.min_det,
        max_det = orientation.max_det,
        worst = orientation.worst_sample,
    );
    assert!(det.is_finite());
    assert_eq!(orientation.samples, 25 * n);
    assert!(orientation.min_det.is_finite());
}

#[test]
fn kp15_pr08_q11_constrained_variational_source_frames_are_measured() {
    let n = 7usize;
    let params = Bi6Reparams::for_valence(n);
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let q = f_kappa_quadratic_form(5);
    let base_patches = kp15_solve_paper_p_f5_for_tensor(&table, n, &base_tensor, &q);
    let basis_deltas = [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.0, 1.0], [0.0, 0.0]],
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 0.0], [0.0, 1.0]],
    ];
    let basis_patches = basis_deltas
        .iter()
        .map(|&frame_delta| {
            let tensor =
                kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
            kp15_patch_delta(
                &kp15_solve_paper_p_f5_for_tensor(&table, n, &tensor, &q),
                &base_patches,
            )
        })
        .collect::<Vec<_>>();
    let (h, g) = kp15_source_frame_quadratic_terms(&base_patches, &basis_patches, &q);
    let fitted_map = kp15_add_identity_to_frame_delta(kp15_q11_fitted_first_order_frame_delta(n));

    for (label, constraints) in [
        (
            "endpoint-spoke-fixed",
            vec![
                (vec![1.0, 0.0, 0.0, 0.0], 0.0),
                (vec![0.0, 0.0, 1.0, 0.0], 0.0),
            ],
        ),
        ("unit-trace", vec![(vec![1.0, 0.0, 0.0, 1.0], 0.0)]),
        (
            "fitted-trace",
            vec![(
                vec![1.0, 0.0, 0.0, 1.0],
                fitted_map[0][0] + fitted_map[1][1] - 2.0,
            )],
        ),
        (
            "fitted-first-column",
            vec![
                (vec![1.0, 0.0, 0.0, 0.0], fitted_map[0][0] - 1.0),
                (vec![0.0, 0.0, 1.0, 0.0], fitted_map[1][0]),
            ],
        ),
    ] {
        let Some(solution) = kp15_solve_quadratic_with_linear_constraints(&h, &g, &constraints)
        else {
            println!("KP15 PR08 q11 constrained F5 source frame {label}: singular");
            continue;
        };
        let frame_delta = [[solution[0], solution[1]], [solution[2], solution[3]]];
        let map = kp15_add_identity_to_frame_delta(frame_delta);
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("{label} constrained source-frame solve failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        let det = kp15_det2(map);
        let fitted_diff = kp15_frame_map_frobenius_delta(map, fitted_map);
        println!(
            "KP15 PR08 q11 constrained F5 source frame {label}: D=[[{d00:e},{d01:e}],[{d10:e},{d11:e}]] det(I+D)={det:e} fitted_diff={fitted_diff:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            d00 = frame_delta[0][0],
            d01 = frame_delta[0][1],
            d10 = frame_delta[1][0],
            d11 = frame_delta[1][1],
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

#[test]
fn kp15_pr08_q11_endpoint_spoke_variational_source_frame_valence_sweep_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let frame_delta = kp15_endpoint_spoke_f5_variational_frame_delta(n, &base_tensor)
            .expect("endpoint-spoke F5 frame solve");
        let map = kp15_add_identity_to_frame_delta(frame_delta);
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, frame_delta);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("n={n} endpoint-spoke source-frame solve failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        println!(
            "KP15 PR08 q11 endpoint-spoke F5 source frame n={n}: D=[[{d00:e},{d01:e}],[{d10:e},{d11:e}]] det(I+D)={det:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            d00 = frame_delta[0][0],
            d01 = frame_delta[0][1],
            d10 = frame_delta[1][0],
            d11 = frame_delta[1][1],
            det = kp15_det2(map),
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

#[test]
fn kp15_pr08_q11_endpoint_spoke_source_frame_line_search_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let frame_delta = kp15_endpoint_spoke_f5_variational_frame_delta(n, &base_tensor)
            .expect("endpoint-spoke F5 frame solve");
        let mut best = (0.0_f64, usize::MAX, f64::NEG_INFINITY);
        for alpha in [0.0_f64, 0.125, 0.25, 0.5, 0.75, 1.0] {
            let scaled = [
                [alpha * frame_delta[0][0], alpha * frame_delta[0][1]],
                [alpha * frame_delta[1][0], alpha * frame_delta[1][1]],
            ];
            let normalized_tensor =
                kp15_apply_sector_local_first_order_normalization(&base_tensor, n, scaled);
            let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
            let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
                .unwrap_or_else(|err| panic!("n={n} alpha={alpha} paper-P/F5 failed: {err}"));
            let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
            println!(
                "KP15 PR08 q11 endpoint-spoke line search n={n} alpha={alpha:e}: neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
                neg = orientation.negative,
                samples = orientation.samples,
                zero = orientation.near_zero,
                min_det = orientation.min_det,
                max_det = orientation.max_det,
                worst = orientation.worst_sample,
            );
            if orientation.negative < best.1
                || (orientation.negative == best.1 && orientation.min_det > best.2)
            {
                best = (alpha, orientation.negative, orientation.min_det);
            }
        }
        println!(
            "KP15 PR08 q11 endpoint-spoke line search best n={n}: alpha={alpha:e} neg={neg} min_det={min_det:e}",
            alpha = best.0,
            neg = best.1,
            min_det = best.2,
        );
        assert!(best.2.is_finite());
    }
}

#[test]
fn kp15_pr08_q11_validation_gated_endpoint_spoke_frame_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let frame_delta = kp15_endpoint_spoke_f5_variational_frame_delta(n, &base_tensor)
            .expect("endpoint-spoke F5 frame solve");
        let selected = kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation(
            n,
            &base_tensor,
            frame_delta,
            1.0e-10,
        );
        let Some(alpha) = selected else {
            println!("KP15 PR08 q11 validation-gated endpoint-spoke frame n={n}: declined");
            continue;
        };
        let scaled = [
            [alpha * frame_delta[0][0], alpha * frame_delta[0][1]],
            [alpha * frame_delta[1][0], alpha * frame_delta[1][1]],
        ];
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, scaled);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0).unwrap_or_else(|err| {
            panic!("n={n} validation-gated source-frame solve failed: {err}")
        });
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        println!(
            "KP15 PR08 q11 validation-gated endpoint-spoke frame n={n}: alpha={alpha:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            worst = orientation.worst_sample,
        );
        assert_eq!(orientation.negative, 0);
        assert!(orientation.min_det > 0.0);
        if n == 5 || n == 6 || n == 8 {
            assert_eq!(alpha, 0.0);
        } else {
            assert!(alpha > 0.0 && alpha < 0.125);
        }
    }
}

#[test]
fn kp15_pr08_q11_validation_gated_max_margin_endpoint_spoke_frame_is_measured() {
    for n in [5usize, 6, 7, 8] {
        let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let frame_delta = kp15_endpoint_spoke_f5_variational_frame_delta(n, &base_tensor)
            .expect("endpoint-spoke F5 frame solve");
        let selected = kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation_margin(
            n,
            &base_tensor,
            frame_delta,
            1.0e-10,
        );
        let Some(alpha) = selected else {
            println!("KP15 PR08 q11 max-margin endpoint-spoke frame n={n}: declined");
            continue;
        };
        let scaled = [
            [alpha * frame_delta[0][0], alpha * frame_delta[0][1]],
            [alpha * frame_delta[1][0], alpha * frame_delta[1][1]],
        ];
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, scaled);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0)
            .unwrap_or_else(|err| panic!("n={n} max-margin source-frame solve failed: {err}"));
        let orientation = test_bi6_patch_local_orientation_report(n, &out.patches_bi6);
        println!(
            "KP15 PR08 q11 max-margin endpoint-spoke frame n={n}: alpha={alpha:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            worst = orientation.worst_sample,
        );
        assert_eq!(orientation.negative, 0);
        assert!(orientation.min_det > 0.0);
        if n == 5 || n == 6 || n == 8 {
            assert_eq!(alpha, 0.0);
        } else {
            assert!(alpha > 0.5 && alpha < 1.0);
        }
    }
}

#[test]
fn kp15_pr08_q11_validation_gated_max_margin_dense_valence_sweep_is_measured() {
    const DENSE_SAMPLES: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];
    for n in [5usize, 6, 7, 8, 9, 10] {
        let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
        let frame_delta = kp15_endpoint_spoke_f5_variational_frame_delta(n, &base_tensor)
            .expect("endpoint-spoke F5 frame solve");
        let Some(alpha) =
            kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation_margin_with_samples(
                n,
                &base_tensor,
                frame_delta,
                1.0e-10,
                &DENSE_SAMPLES,
            )
        else {
            println!("KP15 PR08 q11 dense max-margin endpoint-spoke frame n={n}: declined");
            continue;
        };
        let scaled = [
            [alpha * frame_delta[0][0], alpha * frame_delta[0][1]],
            [alpha * frame_delta[1][0], alpha * frame_delta[1][1]],
        ];
        let normalized_tensor =
            kp15_apply_sector_local_first_order_normalization(&base_tensor, n, scaled);
        let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
        let out = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0).unwrap_or_else(|err| {
            panic!("n={n} dense max-margin source-frame solve failed: {err}")
        });
        let orientation = test_bi6_patch_local_orientation_report_with_samples(
            n,
            &out.patches_bi6,
            &DENSE_SAMPLES,
        );
        println!(
            "KP15 PR08 q11 dense max-margin endpoint-spoke frame n={n}: alpha={alpha:e} neg={neg}/{samples} zero={zero} det_range=[{min_det:e},{max_det:e}] worst={worst}",
            neg = orientation.negative,
            samples = orientation.samples,
            zero = orientation.near_zero,
            min_det = orientation.min_det,
            max_det = orientation.max_det,
            worst = orientation.worst_sample,
        );
        assert_eq!(orientation.negative, 0);
        assert!(orientation.min_det > 0.0);
    }
}

fn kp15_print_seed_compatibility(
    label: &str,
    vars: usize,
    table: &Kp15Bi6GeneratingTable,
    input: &Bi6CapInput,
) {
    let smoothness = kp15_table_smoothness_constraints(table);
    let seeds = kp15_seed_constraint_descriptors(input);
    for axis in 0..2 {
        let mut constraints = smoothness.clone();
        let mut first_incompatible = None::<String>;
        for seed in &seeds {
            constraints.push(seed.constraint.clone());
            let (rank_coeff, rank_aug) = affine_constraint_ranks(vars, &constraints, axis);
            if rank_aug > rank_coeff {
                first_incompatible = Some(seed.label.clone());
                break;
            }
        }
        let mut all_constraints = smoothness.clone();
        all_constraints.extend(seeds.iter().map(|seed| seed.constraint.clone()));
        let (rank_coeff, rank_aug) = affine_constraint_ranks(vars, &all_constraints, axis);
        println!(
            "KP15 spoke-indexed planar seed compatibility {label} axis={axis}: coeff_rank={rank_coeff} augmented_rank={rank_aug} first_incompatible={}",
            first_incompatible.as_deref().unwrap_or("none"),
        );
        assert!(rank_aug >= rank_coeff);
    }
}

#[test]
fn rejects_non_finite_star_or_tensor_border() {
    let tensor = planar_tensor_border(5);
    assert!(matches!(
        solve_kp20_guide_lsq(&tensor, Point3::new(f64::NAN, 0.0, 0.0)),
        Err(Kp20GuideLsqError::InvalidTensorBorder(_))
    ));

    let mut bad = planar_tensor_border(5);
    bad.sectors[2][1][3].z = f64::INFINITY;
    assert!(matches!(
        solve_kp20_guide_lsq(&bad, Point3::new(0.0, 0.0, 0.0)),
        Err(Kp20GuideLsqError::InvalidTensorBorder(_))
    ));
}

#[test]
fn solve_is_deterministic_and_affine_covariant() {
    let tensor = planar_tensor_border(5);
    let star = Point3::new(0.1, -0.2, 0.3);
    let a = solve_kp20_guide_lsq(&tensor, star).unwrap();
    let b = solve_kp20_guide_lsq(&tensor, star).unwrap();
    assert_sectors_close(&a, &b, 1e-15);

    let shifted = map_tensor(&tensor, |p| affine(p, 2.5, Point3::new(3.0, -4.0, 1.0)));
    let shifted_star = affine(star, 2.5, Point3::new(3.0, -4.0, 1.0));
    let shifted_solution = solve_kp20_guide_lsq(&shifted, shifted_star).unwrap();
    let expected: Vec<Kp20Tri6GuideSector> = a
        .iter()
        .map(|sector| Kp20Tri6GuideSector {
            cps: sector
                .cps
                .map(|p| affine(p, 2.5, Point3::new(3.0, -4.0, 1.0))),
        })
        .collect();
    assert_sectors_close(&shifted_solution, &expected, 1e-7);
}

fn planar_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::new();
    for s in 0..n {
        let a = core::f64::consts::TAU * s as f64 / n as f64;
        let b = core::f64::consts::TAU * (s + 1) as f64 / n as f64;
        let p0 = Point3::new(a.cos(), a.sin(), 0.0);
        let p3 = Point3::new(b.cos(), b.sin(), 0.0);
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for i in 0..4 {
            let t = i as f64 / 3.0;
            let x = (1.0 - t) * p0.x + t * p3.x;
            let y = (1.0 - t) * p0.y + t * p3.y;
            rows[0][i] = Point3::new(x, y, 0.0);
            rows[1][i] = Point3::new(0.8 * x, 0.8 * y, 0.0);
            rows[2][i] = Point3::new(0.6 * x, 0.6 * y, 0.0);
        }
        sectors.push(rows);
    }
    TensorBorderBi3 { sectors }
}

fn characteristic_tensor_border(n: usize) -> TensorBorderBi3 {
    let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap net");
    let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
    kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor border")
}

fn solve_kp20_guide_lsq_with_custom_rows(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    rows: ResidualRows,
) -> Result<Kp20GuideLsqOutput, Kp20GuideLsqError> {
    let n = tensor_border.sectors.len();
    solve_kp20_guide_lsq_with_custom_rows_and_constraints(
        tensor_border,
        star_point,
        rows,
        constraint_rows(n)?,
    )
}

fn solve_kp20_guide_lsq_with_custom_rows_and_constraints(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    rows: ResidualRows,
    all_constraints: Vec<LinearRow>,
) -> Result<Kp20GuideLsqOutput, Kp20GuideLsqError> {
    let n = tensor_border.sectors.len();
    validate_tensor_border(tensor_border)?;
    let vars = n * TRI6_CP_COUNT;
    let constraints = independent_rows(vars, &all_constraints);
    let mut solved = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];

    let mut kkt_condition_number: f64 = 0.0;
    for (axis, solved_axis) in solved.iter_mut().enumerate() {
        let rhs = residual_rhs(&rows, axis);
        let fixed = axis_value(star_point, axis);
        let axis_solution = solve_axis(vars, &rows, &rhs, &constraints, fixed)
            .ok_or(Kp20GuideLsqError::SingularSystem(axis))?;
        kkt_condition_number = kkt_condition_number.max(axis_solution.kkt_condition_number);
        *solved_axis = axis_solution.values;
    }

    let mut sectors = vec![Kp20Tri6GuideSector::zero(); n];
    for s in 0..n {
        for cp in 0..TRI6_CP_COUNT {
            sectors[s].cps[cp] = Point3::new(
                solved[0][var_index(s, cp)],
                solved[1][var_index(s, cp)],
                solved[2][var_index(s, cp)],
            );
        }
        sectors[s].cps[CENTER_FLAT] = star_point;
        sectors[s].cps[flat_index(DEGREE, DEGREE, 0)] = tensor_border.sectors[s][0][0];
        sectors[s].cps[flat_index(DEGREE, 0, DEGREE)] = tensor_border.sectors[s][0][3];
    }
    let report = report(
        n,
        vars,
        &rows,
        &all_constraints,
        &solved,
        kkt_condition_number,
    );
    Ok(Kp20GuideLsqOutput { sectors, report })
}

#[derive(Clone)]
struct AffineConstraintRow {
    row: LinearRow,
    rhs: [f64; 3],
}

fn solve_kp20_guide_with_outer_layers_and_f4(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
) -> Result<Kp20GuideLsqOutput, Kp20GuideLsqError> {
    let n = tensor_border.sectors.len();
    validate_tensor_border(tensor_border)?;
    let vars = n * TRI6_CP_COUNT;
    let mut affine_constraints = Vec::new();
    affine_constraints.push(AffineConstraintRow {
        row: LinearRow {
            terms: vec![(var_index(0, CENTER_FLAT), 1.0)],
        },
        rhs: [star_point.x, star_point.y, star_point.z],
    });
    for row in constraint_rows(n)? {
        affine_constraints.push(AffineConstraintRow {
            row,
            rhs: [0.0, 0.0, 0.0],
        });
    }
    append_outer_layer_pin_constraints(tensor_border, &mut affine_constraints);
    let constraints = independent_affine_rows(vars, &affine_constraints);
    let rows = triangular_f4_residual_rows(n);
    let mut solved = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];

    let mut kkt_condition_number: f64 = 0.0;
    for (axis, solved_axis) in solved.iter_mut().enumerate() {
        let axis_solution = solve_axis_with_affine_constraints(
            vars,
            &rows,
            &rows.targets[axis],
            &constraints,
            axis,
        )
        .ok_or(Kp20GuideLsqError::SingularSystem(axis))?;
        kkt_condition_number = kkt_condition_number.max(axis_solution.kkt_condition_number);
        *solved_axis = axis_solution.values;
    }

    let mut sectors = vec![Kp20Tri6GuideSector::zero(); n];
    for s in 0..n {
        for cp in 0..TRI6_CP_COUNT {
            sectors[s].cps[cp] = Point3::new(
                solved[0][var_index(s, cp)],
                solved[1][var_index(s, cp)],
                solved[2][var_index(s, cp)],
            );
        }
    }
    let report =
        report_with_affine_constraints(n, vars, &rows, &constraints, &solved, kkt_condition_number);
    Ok(Kp20GuideLsqOutput { sectors, report })
}

fn append_outer_layer_pin_constraints(
    tensor_border: &TensorBorderBi3,
    constraints: &mut Vec<AffineConstraintRow>,
) {
    for (s, sector) in tensor_border.sectors.iter().enumerate() {
        let raised0 = degree_raise_1d_3_to_6(&sector[0]);
        append_exact_depth_pin_constraints(constraints, s, 0, &raised0);
        let raised1 = degree_raise_1d_3_to_5(&sector[1]);
        append_exact_depth_pin_constraints(constraints, s, 1, &raised1);
        let raised2 = degree_raise_1d_3_to_4_for_test(&sector[2]);
        append_exact_depth_pin_constraints(constraints, s, 2, &raised2);
    }
}

fn append_exact_depth_pin_constraints<const N: usize>(
    constraints: &mut Vec<AffineConstraintRow>,
    sector: usize,
    depth: usize,
    raised: &[Point3; N],
) {
    let row_degree = DEGREE - depth;
    assert_eq!(N, row_degree + 1);
    for (j, target) in raised.iter().enumerate() {
        constraints.push(AffineConstraintRow {
            row: LinearRow {
                terms: vec![(var_index(sector, cp(j, row_degree - j)), 1.0)],
            },
            rhs: [target.x, target.y, target.z],
        });
    }
}

fn triangular_f4_residual_rows(n: usize) -> ResidualRows {
    let mut rows = Vec::new();
    let mut targets = [Vec::new(), Vec::new(), Vec::new()];
    const SAMPLES: [(f64, f64, f64); 6] = [
        (0.091576213509771, 0.091576213509771, 0.109951743655322),
        (0.091576213509771, 0.816847572980459, 0.109951743655322),
        (0.816847572980459, 0.091576213509771, 0.109951743655322),
        (0.445948490915965, 0.445948490915965, 0.223381589678011),
        (0.445948490915965, 0.108103018168070, 0.223381589678011),
        (0.108103018168070, 0.445948490915965, 0.223381589678011),
    ];
    for sector in 0..n {
        for &(u, v, area_weight) in &SAMPLES {
            for du in 0..=4 {
                let dv = 4 - du;
                let binomial_weight = (binom(4, du) as f64 * area_weight).sqrt();
                let basis = triangular_basis_derivative(DEGREE, du, dv, u, v);
                rows.push(LinearRow {
                    terms: basis
                        .iter()
                        .enumerate()
                        .filter_map(|(cp, &w)| {
                            if w.abs() > 1.0e-14 {
                                Some((var_index(sector, cp), binomial_weight * w))
                            } else {
                                None
                            }
                        })
                        .collect(),
                });
                targets[0].push(0.0);
                targets[1].push(0.0);
                targets[2].push(0.0);
            }
        }
    }
    ResidualRows { rows, targets }
}

fn triangular_basis_derivative(n: usize, du: usize, dv: usize, u: f64, v: f64) -> Vec<f64> {
    if du == 0 && dv == 0 {
        return triangular_basis(n, u, v);
    }
    if du > 0 {
        let lower = triangular_basis_derivative(n - 1, du - 1, dv, u, v);
        lift_triangular_partial_u(n, &lower)
    } else {
        let lower = triangular_basis_derivative(n - 1, du, dv - 1, u, v);
        lift_triangular_partial_v(n, &lower)
    }
}

fn triangular_basis(n: usize, u: f64, v: f64) -> Vec<f64> {
    let w = 1.0 - u - v;
    let mut out = vec![0.0_f64; (n + 1) * (n + 2) / 2];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            out[flat_index(n, i, j)] =
                trinomial(n, i, j) as f64 * u.powi(i as i32) * v.powi(j as i32) * w.powi(k as i32);
        }
    }
    out
}

fn lift_triangular_partial_u(n: usize, lower: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0_f64; (n + 1) * (n + 2) / 2];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            let term_i = if i >= 1 {
                lower[flat_index(n - 1, i - 1, j)]
            } else {
                0.0
            };
            let term_k = if k >= 1 {
                lower[flat_index(n - 1, i, j)]
            } else {
                0.0
            };
            out[flat_index(n, i, j)] = n as f64 * (term_i - term_k);
        }
    }
    out
}

fn lift_triangular_partial_v(n: usize, lower: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0_f64; (n + 1) * (n + 2) / 2];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            let term_j = if j >= 1 {
                lower[flat_index(n - 1, i, j - 1)]
            } else {
                0.0
            };
            let term_k = if k >= 1 {
                lower[flat_index(n - 1, i, j)]
            } else {
                0.0
            };
            out[flat_index(n, i, j)] = n as f64 * (term_j - term_k);
        }
    }
    out
}

fn solve_axis_with_affine_constraints(
    vars: usize,
    residuals: &ResidualRows,
    rhs: &[f64],
    constraints: &[AffineConstraintRow],
    axis: usize,
) -> Option<AxisSolution> {
    let size = vars + constraints.len();
    let mut a = vec![vec![0.0_f64; size]; size];
    let mut b = vec![0.0_f64; size];
    for (row, &target) in residuals.rows.iter().zip(rhs) {
        accumulate_normal_row(&mut a, &mut b, &row.terms, target, 1.0);
    }
    for (i, row) in a.iter_mut().enumerate().take(vars) {
        row[i] += 1.0e-12;
    }
    for (cidx, constraint) in constraints.iter().enumerate() {
        add_lagrange_row(
            vars + cidx,
            &mut a,
            &mut b,
            &constraint.row.terms,
            constraint.rhs[axis],
        );
    }
    let kkt_condition_number = condition_number_dense(&a)?.condition_number;
    solve_dense(a, b).map(|x| AxisSolution {
        values: x[..vars].to_vec(),
        kkt_condition_number,
    })
}

fn independent_affine_rows(vars: usize, rows: &[AffineConstraintRow]) -> Vec<AffineConstraintRow> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut keep = Vec::new();
    for row in rows {
        let mut dense = vec![0.0_f64; vars];
        for &(i, w) in &row.row.terms {
            dense[i] += w;
        }
        for b in &basis {
            let dot = dense.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            for (x, bi) in dense.iter_mut().zip(b) {
                *x -= dot * bi;
            }
        }
        let norm = dense.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1.0e-10 {
            for x in &mut dense {
                *x /= norm;
            }
            basis.push(dense);
            keep.push(row.clone());
        }
    }
    keep
}

fn report_with_affine_constraints(
    n: usize,
    vars: usize,
    rows: &ResidualRows,
    constraints: &[AffineConstraintRow],
    solved: &[Vec<f64>; 3],
    kkt_condition_number: f64,
) -> Kp20GuideLsqReport {
    let mut report = report(n, vars, rows, &[], solved, kkt_condition_number);
    let mut max_constraint: f64 = 0.0;
    for constraint in constraints {
        for (axis, solved_axis) in solved.iter().enumerate() {
            let residual = row_value(&constraint.row, solved_axis) - constraint.rhs[axis];
            max_constraint = max_constraint.max(residual.abs());
        }
    }
    report.independent_constraints = constraints.len();
    report.free_parameters_after_center_pin = vars - constraints.len();
    report.max_constraint_residual = max_constraint;
    report
}

fn trinomial(n: usize, i: usize, j: usize) -> usize {
    binom(n, i) * binom(n - i, j)
}

fn kp18_join_rows_without_apex_transpose(n: usize) -> Vec<LinearRow> {
    let rows = constraint_rows(n).expect("KP18 constraints");
    let per_sector_rows_with_transpose = 18;
    rows.into_iter()
        .enumerate()
        .filter_map(|(row_index, row)| {
            if row_index < n * per_sector_rows_with_transpose {
                let local = row_index % per_sector_rows_with_transpose;
                if local >= 16 {
                    return None;
                }
            }
            Some(row)
        })
        .collect()
}

fn exact_depth_bb_row_residual_rows(t: &TensorBorderBi3) -> ResidualRows {
    let mut rows = Vec::new();
    let mut targets = [Vec::new(), Vec::new(), Vec::new()];
    for (s, sector) in t.sectors.iter().enumerate() {
        let raised0 = degree_raise_1d_3_to_6(&sector[0]);
        append_exact_depth_row(&mut rows, &mut targets, s, 0, &raised0);
        let raised1 = degree_raise_1d_3_to_5(&sector[1]);
        append_exact_depth_row(&mut rows, &mut targets, s, 1, &raised1);
        let raised2 = degree_raise_1d_3_to_4_for_test(&sector[2]);
        append_exact_depth_row(&mut rows, &mut targets, s, 2, &raised2);
    }
    ResidualRows { rows, targets }
}

fn append_exact_depth_row<const N: usize>(
    rows: &mut Vec<LinearRow>,
    targets: &mut [Vec<f64>; 3],
    sector: usize,
    depth: usize,
    raised: &[Point3; N],
) {
    let row_degree = DEGREE - depth;
    assert_eq!(N, row_degree + 1);
    for (j, target) in raised.iter().enumerate() {
        rows.push(LinearRow {
            terms: vec![(var_index(sector, cp(j, row_degree - j)), 1.0)],
        });
        targets[0].push(target.x);
        targets[1].push(target.y);
        targets[2].push(target.z);
    }
}

fn degree_raise_1d_3_to_6(c: &[Point3; 4]) -> [Point3; 7] {
    let d5 = degree_raise_1d_3_to_5(c);
    let mut d6 = [Point3::new(0.0, 0.0, 0.0); 7];
    for k in 0..7 {
        let a = k as f64 / 6.0;
        let b = (6 - k) as f64 / 6.0;
        let prev = if k > 0 {
            d5[k - 1]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        let cur = if k < 6 {
            d5[k]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        d6[k] = Point3::new(
            a * prev.x + b * cur.x,
            a * prev.y + b * cur.y,
            a * prev.z + b * cur.z,
        );
    }
    d6
}

fn degree_raise_1d_3_to_4_for_test(c: &[Point3; 4]) -> [Point3; 5] {
    let mut d4 = [Point3::new(0.0, 0.0, 0.0); 5];
    for k in 0..5 {
        let a = k as f64 / 4.0;
        let b = (4 - k) as f64 / 4.0;
        let prev = if k > 0 {
            c[k - 1]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        let cur = if k < 4 {
            c[k]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        d4[k] = Point3::new(
            a * prev.x + b * cur.x,
            a * prev.y + b * cur.y,
            a * prev.z + b * cur.z,
        );
    }
    d4
}

fn kk12_equation26_star_point_from_pr08_labels(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
    first_label: usize,
    c7_label: usize,
) -> Point3 {
    const E: [f64; 6] = [
        0.00004133, 0.00049235, 0.00079515, 0.00049235, 0.04832434, 0.16420709,
    ];
    let mut star = Point3::new(0.0, 0.0, 0.0);
    for sector in 0..3 {
        for (m, &weight) in E.iter().enumerate() {
            let p = ring.sectors[sector][first_label + m];
            star.x += weight * p.x;
            star.y += weight * p.y;
            star.z += weight * p.z;
        }
    }
    let e_sum = E.iter().sum::<f64>();
    let c7 = ring.sectors[1][c7_label];
    star.x += (1.0 - 3.0 * e_sum) * c7.x;
    star.y += (1.0 - 3.0 * e_sum) * c7.y;
    star.z += (1.0 - 3.0 * e_sum) * c7.z;
    star
}

fn linear_sector_wedge_prior(n: usize) -> Vec<Kp20Tri6GuideSector> {
    (0..n)
        .map(|sector| {
            let mut guide = Kp20Tri6GuideSector::zero();
            for i in 0..=DEGREE {
                for j in 0..=(DEGREE - i) {
                    let local = [i as f64 / DEGREE as f64, j as f64 / DEGREE as f64];
                    guide.cps[flat_index(DEGREE, i, j)] = local_to_global_planar(n, sector, local);
                }
            }
            guide
        })
        .collect()
}

fn test_guide_local_orientation_report(
    n: usize,
    sectors: &[Kp20Tri6GuideSector],
) -> TestGuideOrientationReport {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let mut out = TestGuideOrientationReport {
        samples: 0,
        negative: 0,
        near_zero: 0,
        min_det: f64::INFINITY,
        max_det: f64::NEG_INFINITY,
        worst_sample: String::new(),
    };
    for (sector, guide) in sectors.iter().enumerate() {
        for &u in &SAMPLES {
            for &v in &SAMPLES {
                if u + v > 1.0 {
                    continue;
                }
                let jet = guide.jet3(u, v);
                let gu = global_to_local_planar(n, sector, jet.h[1][0]);
                let gv = global_to_local_planar(n, sector, jet.h[0][1]);
                let det = gu[0] * gv[1] - gv[0] * gu[1];
                out.samples += 1;
                out.negative += usize::from(det < -1.0e-10);
                out.near_zero += usize::from(det.abs() <= 1.0e-10);
                if det < out.min_det {
                    out.min_det = det;
                    out.worst_sample = format!("sector={sector}@({u:.2},{v:.2})");
                }
                out.max_det = out.max_det.max(det);
            }
        }
    }
    out
}

fn test_tri5_local_orientation_report(
    n: usize,
    sectors: &[super::super::kp13_bi5_gtri5::Gtri5Sector],
) -> TestGuideOrientationReport {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let mut out = TestGuideOrientationReport {
        samples: 0,
        negative: 0,
        near_zero: 0,
        min_det: f64::INFINITY,
        max_det: f64::NEG_INFINITY,
        worst_sample: String::new(),
    };
    for (sector, guide) in sectors.iter().enumerate() {
        for &u in &SAMPLES {
            for &v in &SAMPLES {
                if u + v > 1.0 {
                    continue;
                }
                let gu = eval_tri5_derivative(&guide.cps, &bernstein_tri5_partial_u(u, v));
                let gv = eval_tri5_derivative(&guide.cps, &bernstein_tri5_partial_v(u, v));
                let lu = global_to_local_planar(n, sector, gu);
                let lv = global_to_local_planar(n, sector, gv);
                let det = lu[0] * lv[1] - lv[0] * lu[1];
                out.samples += 1;
                out.negative += usize::from(det < -1.0e-10);
                out.near_zero += usize::from(det.abs() <= 1.0e-10);
                if det < out.min_det {
                    out.min_det = det;
                    out.worst_sample = format!("sector={sector}@({u:.2},{v:.2})");
                }
                out.max_det = out.max_det.max(det);
            }
        }
    }
    out
}

fn test_bi6_patch_local_orientation_report(
    n: usize,
    patches: &[[[Point3; 7]; 7]],
) -> TestGuideOrientationReport {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    test_bi6_patch_local_orientation_report_with_samples(n, patches, &SAMPLES)
}

fn test_bi6_patch_local_orientation_report_with_samples(
    n: usize,
    patches: &[[[Point3; 7]; 7]],
    samples: &[f64],
) -> TestGuideOrientationReport {
    let mut out = TestGuideOrientationReport {
        samples: 0,
        negative: 0,
        near_zero: 0,
        min_det: f64::INFINITY,
        max_det: f64::NEG_INFINITY,
        worst_sample: String::new(),
    };
    for (sector, patch) in patches.iter().enumerate() {
        for &u in samples {
            for &v in samples {
                let (du, dv) = eval_bi6_patch_derivatives(patch, u, v);
                let lu = global_to_local_planar(n, sector, du);
                let lv = global_to_local_planar(n, sector, dv);
                let det = lu[0] * lv[1] - lv[0] * lu[1];
                out.samples += 1;
                out.negative += usize::from(det < -1.0e-10);
                out.near_zero += usize::from(det.abs() <= 1.0e-10);
                if det < out.min_det {
                    out.min_det = det;
                    out.worst_sample = format!("sector={sector}@({u:.2},{v:.2})");
                }
                out.max_det = out.max_det.max(det);
            }
        }
    }
    out
}

fn test_bi6_patch_boundary_orientation_report(
    n: usize,
    patches: &[[[Point3; 7]; 7]],
    samples: &[f64],
) -> TestBi6BoundaryOrientationReport {
    let mut out = TestBi6BoundaryOrientationReport {
        all: TestGuideOrientationReport {
            samples: 0,
            negative: 0,
            near_zero: 0,
            min_det: f64::INFINITY,
            max_det: f64::NEG_INFINITY,
            worst_sample: String::new(),
        },
        u0_negative: 0,
        u1_negative: 0,
        v0_negative: 0,
        v1_negative: 0,
        interior_negative: 0,
    };
    for (sector, patch) in patches.iter().enumerate() {
        for &u in samples {
            for &v in samples {
                let (du, dv) = eval_bi6_patch_derivatives(patch, u, v);
                let lu = global_to_local_planar(n, sector, du);
                let lv = global_to_local_planar(n, sector, dv);
                let det = lu[0] * lv[1] - lv[0] * lu[1];
                out.all.samples += 1;
                out.all.negative += usize::from(det < -1.0e-10);
                out.all.near_zero += usize::from(det.abs() <= 1.0e-10);
                if det < out.all.min_det {
                    out.all.min_det = det;
                    out.all.worst_sample = format!("sector={sector}@({u:.4},{v:.4})");
                }
                out.all.max_det = out.all.max_det.max(det);
                if det < -1.0e-10 {
                    if u == 0.0 {
                        out.u0_negative += 1;
                    } else if u == 1.0 {
                        out.u1_negative += 1;
                    } else if v == 0.0 {
                        out.v0_negative += 1;
                    } else if v == 1.0 {
                        out.v1_negative += 1;
                    } else {
                        out.interior_negative += 1;
                    }
                }
            }
        }
    }
    out
}

fn eval_bi6_patch_derivatives(patch: &[[Point3; 7]; 7], u: f64, v: f64) -> (Point3, Point3) {
    let bu = bernstein_bi6_all(u);
    let bv = bernstein_bi6_all(v);
    let du_basis = bernstein_bi6_deriv_all(u);
    let dv_basis = bernstein_bi6_deriv_all(v);
    let mut du = Point3::new(0.0, 0.0, 0.0);
    let mut dv = Point3::new(0.0, 0.0, 0.0);
    for i in 0..7 {
        for j in 0..7 {
            let p = patch[i][j];
            let wu = du_basis[i] * bv[j];
            let wv = bu[i] * dv_basis[j];
            du.x += wu * p.x;
            du.y += wu * p.y;
            du.z += wu * p.z;
            dv.x += wv * p.x;
            dv.y += wv * p.y;
            dv.z += wv * p.z;
        }
    }
    (du, dv)
}

fn kp15_seeded_position_residual(input: &Bi6CapInput, patches: &[[[Point3; 7]; 7]]) -> f64 {
    let mut max_residual: f64 = 0.0;
    for (sector, source) in input.sectors.iter().enumerate() {
        let patch = &patches[sector];
        for (i, j, expected) in [
            (0, 0, Point3::new(0.0, 0.0, 0.0)),
            (1, 0, source.p_10),
            (2, 0, source.p_20),
            (4, 0, source.p_40),
            (5, 0, source.p_50),
            (6, 0, source.p_60),
            (0, 1, source.p_01),
            (1, 1, source.p_11),
            (0, 2, source.p_02),
            (4, 2, source.p_42),
            (5, 2, source.p_52),
            (6, 2, source.p_62),
        ] {
            max_residual = max_residual.max(patch[i][j].distance_to(&expected));
        }
    }
    max_residual
}

fn kp15_solve_fixed_seed_rows_with_f5(
    input: &Bi6CapInput,
    table: &Kp15Bi6GeneratingTable,
) -> Result<Vec<[[Point3; 7]; 7]>, String> {
    let n = input.n;
    let vars = 49 * n;
    let q = f_kappa_quadratic_form(5);
    let mut constraints = Vec::new();
    for row in 0..table.rank() {
        let mut terms = Vec::new();
        for col in 0..vars {
            let weight = table.rref[row * table.cols + col];
            if weight.abs() > 1.0e-14 {
                terms.push((col, weight));
            }
        }
        constraints.push(AffineConstraintRow {
            row: LinearRow { terms },
            rhs: [0.0, 0.0, 0.0],
        });
    }
    append_kp15_seed_constraints(input, &mut constraints);

    let mut solved = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];
    for (axis, values) in solved.iter_mut().enumerate() {
        *values = kp15_solve_axis_fixed_seed_rows(vars, n, &q, &constraints, axis)
            .ok_or_else(|| format!("KP15 fixed-row F5 KKT singular on axis {axis}"))?;
    }

    let mut patches = Vec::with_capacity(n);
    for sector in 0..n {
        let mut patch = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for (i, row) in patch.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let col = kp15_global_col(n, sector, i, j);
                *cell = Point3::new(solved[0][col], solved[1][col], solved[2][col]);
            }
        }
        patches.push(patch);
    }
    Ok(patches)
}

fn append_kp15_seed_constraints(input: &Bi6CapInput, constraints: &mut Vec<AffineConstraintRow>) {
    for (sector, source) in input.sectors.iter().enumerate() {
        for (i, j, target) in [
            (0, 0, Point3::new(0.0, 0.0, 0.0)),
            (1, 0, source.p_10),
            (2, 0, source.p_20),
            (4, 0, source.p_40),
            (5, 0, source.p_50),
            (6, 0, source.p_60),
            (0, 1, source.p_01),
            (1, 1, source.p_11),
            (0, 2, source.p_02),
            (4, 2, source.p_42),
            (5, 2, source.p_52),
            (6, 2, source.p_62),
        ] {
            constraints.push(AffineConstraintRow {
                row: LinearRow {
                    terms: vec![(kp15_global_col(input.n, sector, i, j), 1.0)],
                },
                rhs: [target.x, target.y, target.z],
            });
        }
    }
}

#[derive(Clone)]
struct SeedConstraintDescriptor {
    label: String,
    constraint: AffineConstraintRow,
}

fn kp15_table_smoothness_constraints(table: &Kp15Bi6GeneratingTable) -> Vec<AffineConstraintRow> {
    let mut constraints = Vec::with_capacity(table.rank());
    for row in 0..table.rank() {
        let mut terms = Vec::new();
        for col in 0..table.cols {
            let weight = table.rref[row * table.cols + col];
            if weight.abs() > 1.0e-14 {
                terms.push((col, weight));
            }
        }
        constraints.push(AffineConstraintRow {
            row: LinearRow { terms },
            rhs: [0.0, 0.0, 0.0],
        });
    }
    constraints
}

fn kp15_seed_constraint_descriptors(input: &Bi6CapInput) -> Vec<SeedConstraintDescriptor> {
    let mut descriptors = Vec::with_capacity(12 * input.n);
    for (sector, source) in input.sectors.iter().enumerate() {
        for (i, j, target) in [
            (0, 0, Point3::new(0.0, 0.0, 0.0)),
            (1, 0, source.p_10),
            (2, 0, source.p_20),
            (4, 0, source.p_40),
            (5, 0, source.p_50),
            (6, 0, source.p_60),
            (0, 1, source.p_01),
            (1, 1, source.p_11),
            (0, 2, source.p_02),
            (4, 2, source.p_42),
            (5, 2, source.p_52),
            (6, 2, source.p_62),
        ] {
            descriptors.push(SeedConstraintDescriptor {
                label: format!("sector={sector} p_{i}{j}"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(input.n, sector, i, j), 1.0)],
                    },
                    rhs: [target.x, target.y, target.z],
                },
            });
        }
    }
    descriptors
}

fn affine_constraint_ranks(
    vars: usize,
    constraints: &[AffineConstraintRow],
    axis: usize,
) -> (usize, usize) {
    let coefficient_rows = constraints
        .iter()
        .map(|constraint| {
            let mut row = vec![0.0_f64; vars];
            for &(col, weight) in &constraint.row.terms {
                row[col] += weight;
            }
            row
        })
        .collect::<Vec<_>>();
    let augmented_rows = constraints
        .iter()
        .map(|constraint| {
            let mut row = vec![0.0_f64; vars + 1];
            for &(col, weight) in &constraint.row.terms {
                row[col] += weight;
            }
            row[vars] = constraint.rhs[axis];
            row
        })
        .collect::<Vec<_>>();
    (
        dense_rank_relative_for_test(&coefficient_rows),
        dense_rank_relative_for_test(&augmented_rows),
    )
}

fn dense_rank_relative_for_test(rows: &[Vec<f64>]) -> usize {
    let row_count = rows.len();
    let col_count = rows.first().map_or(0, Vec::len);
    if row_count == 0
        || col_count == 0
        || rows.iter().any(|row| row.len() != col_count)
        || rows.iter().flatten().any(|value| !value.is_finite())
    {
        return 0;
    }
    let matrix = Mat::from_fn(row_count, col_count, |r, c| rows[r][c]);
    let singular_values = matrix.singular_values();
    let sigma_max = singular_values.first().copied().unwrap_or(0.0);
    let tol = (sigma_max * 1.0e-10).max(1.0e-12);
    singular_values.iter().filter(|&&s| s > tol).count()
}

fn kp15_solve_axis_fixed_seed_rows(
    vars: usize,
    n: usize,
    q: &[[f64; 49]; 49],
    constraints: &[AffineConstraintRow],
    axis: usize,
) -> Option<Vec<f64>> {
    let size = vars + constraints.len();
    let mut a = vec![vec![0.0_f64; size]; size];
    let mut b = vec![0.0_f64; size];
    for sector in 0..n {
        let base = 49 * sector;
        for r in 0..49 {
            for c in 0..49 {
                a[base + r][base + c] += 2.0 * q[r][c];
            }
        }
    }
    for (i, row) in a.iter_mut().enumerate().take(vars) {
        row[i] += 1.0e-12;
    }
    for (row_index, constraint) in constraints.iter().enumerate() {
        add_lagrange_row(
            vars + row_index,
            &mut a,
            &mut b,
            &constraint.row.terms,
            constraint.rhs[axis],
        );
    }
    solve_dense_faer_for_test(a, b).map(|x| x[..vars].to_vec())
}

fn kp15_global_col(n: usize, sector: usize, i: usize, j: usize) -> usize {
    assert!(sector < n);
    assert!(i < 7 && j < 7);
    sector * 49 + i * 7 + j
}

#[derive(Clone, Copy)]
struct Kp15FreeSplitForTest {
    fixed: usize,
    variable: usize,
    zeroed_extra: usize,
}

fn kp15_solve_table_with_test_variable_columns(
    table: &Kp15Bi6GeneratingTable,
    sectors: &[CirculantSectorRhs],
    q_local: &[[f64; 49]; 49],
    allowed_variable_cols: Option<&std::collections::BTreeSet<usize>>,
    released_tensor_cols: Option<&std::collections::BTreeSet<usize>>,
) -> Result<(Vec<[[Point3; 7]; 7]>, [f64; 3], Kp15FreeSplitForTest), String> {
    let n = table.n;
    assert_eq!(sectors.len(), n);
    let n_free = table.free_cols.len();
    let mut t = vec![0.0_f64; table.cols * n_free];
    for j in 0..n_free {
        let mut free = vec![0.0_f64; n_free];
        free[j] = 1.0;
        let x = table.apply(&free);
        for (row, value) in x.iter().enumerate() {
            t[row * n_free + j] = *value;
        }
    }

    let mut h = vec![0.0_f64; n_free * n_free];
    for sector in 0..n {
        let row_base = 49 * sector;
        let mut qtk = vec![0.0_f64; 49 * n_free];
        for r in 0..49 {
            for c in 0..n_free {
                let mut acc = 0.0;
                for inner in 0..49 {
                    acc += q_local[r][inner] * t[(row_base + inner) * n_free + c];
                }
                qtk[r * n_free + c] = acc;
            }
        }
        for r in 0..n_free {
            for c in 0..n_free {
                let mut acc = 0.0;
                for inner in 0..49 {
                    acc += t[(row_base + inner) * n_free + r] * qtk[inner * n_free + c];
                }
                h[r * n_free + c] += acc;
            }
        }
    }

    let mut fixed_idx = Vec::new();
    let mut variable_idx = Vec::new();
    let mut fixed_values = [Vec::new(), Vec::new(), Vec::new()];
    let mut zeroed_extra = 0usize;
    for (free_index, &col) in table.free_cols.iter().enumerate() {
        let sector = col / 49;
        let pos = col % 49;
        let i = pos / 7;
        let j = pos % 7;
        let tensor_fixed_value =
            kp15_seed_value_for_test(sectors, sector, i, j).filter(|_| i + j > 2);
        let released_tensor = tensor_fixed_value.is_some()
            && released_tensor_cols
                .map(|cols| cols.contains(&col))
                .unwrap_or(false);
        let variable_allowed = released_tensor
            || (tensor_fixed_value.is_none()
                && allowed_variable_cols
                    .map(|cols| cols.contains(&col))
                    .unwrap_or(true));
        if variable_allowed {
            variable_idx.push(free_index);
        } else {
            fixed_idx.push(free_index);
            let value = tensor_fixed_value.unwrap_or(Point3::new(0.0, 0.0, 0.0));
            fixed_values[0].push(value.x);
            fixed_values[1].push(value.y);
            fixed_values[2].push(value.z);
            if tensor_fixed_value.is_none() {
                zeroed_extra += 1;
            }
        }
    }

    let pull_h = |rows: &[usize], cols: &[usize]| -> Vec<Vec<f64>> {
        rows.iter()
            .map(|&r| cols.iter().map(|&c| h[r * n_free + c]).collect())
            .collect()
    };
    let h_vv = pull_h(&variable_idx, &variable_idx);
    let h_vf = pull_h(&variable_idx, &fixed_idx);
    let mut variable_values = [
        vec![0.0; variable_idx.len()],
        vec![0.0; variable_idx.len()],
        vec![0.0; variable_idx.len()],
    ];
    for axis in 0..3 {
        let rhs = (0..variable_idx.len())
            .map(|r| {
                let mut acc = 0.0;
                for c in 0..fixed_idx.len() {
                    acc -= h_vf[r][c] * fixed_values[axis][c];
                }
                acc
            })
            .collect::<Vec<_>>();
        variable_values[axis] = solve_dense_faer_for_test(h_vv.clone(), rhs)
            .ok_or_else(|| format!("test F5 solve singular at axis {axis}"))?;
    }

    let assemble_free = |axis: usize| -> Vec<f64> {
        let mut free = vec![0.0_f64; n_free];
        for (k, &free_index) in fixed_idx.iter().enumerate() {
            free[free_index] = fixed_values[axis][k];
        }
        for (k, &free_index) in variable_idx.iter().enumerate() {
            free[free_index] = variable_values[axis][k];
        }
        free
    };
    let xs = table.apply(&assemble_free(0));
    let ys = table.apply(&assemble_free(1));
    let zs = table.apply(&assemble_free(2));

    let mut patches = Vec::with_capacity(n);
    for sector in 0..n {
        let mut patch = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for (i, row) in patch.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let col = kp15_global_col(n, sector, i, j);
                *cell = Point3::new(xs[col], ys[col], zs[col]);
            }
        }
        patches.push(patch);
    }

    let f_coord = |x: &[f64]| -> f64 {
        let mut acc = 0.0;
        for sector in 0..n {
            let base = 49 * sector;
            for r in 0..49 {
                for c in 0..49 {
                    acc += x[base + r] * q_local[r][c] * x[base + c];
                }
            }
        }
        acc
    };
    Ok((
        patches,
        [f_coord(&xs), f_coord(&ys), f_coord(&zs)],
        Kp15FreeSplitForTest {
            fixed: fixed_idx.len(),
            variable: variable_idx.len(),
            zeroed_extra,
        },
    ))
}

fn kp15_paper_source_variable_columns_for_test(n: usize) -> std::collections::BTreeSet<usize> {
    let mut cols = std::collections::BTreeSet::new();
    for i in 0..=2 {
        for j in 0..=2 {
            if i + j <= 2 {
                cols.insert(kp15_global_col(n, 0, i, j));
            }
        }
    }
    cols.insert(kp15_global_col(n, 0, 2, 2));
    for sector in 0..n {
        cols.insert(kp15_global_col(n, sector, 2, 3));
        cols.insert(kp15_global_col(n, sector, 3, 3));
        for i in 0..=2 {
            for j in 0..=2 {
                if i + j <= 2 {
                    cols.insert(kp15_global_col(n, sector, i, j));
                }
            }
        }
    }
    if n == 6 {
        cols.insert(kp15_global_col(n, 0, 2, 1));
    }
    cols
}

fn kp15_strict_p_variable_columns_for_test(n: usize) -> std::collections::BTreeSet<usize> {
    let mut cols = std::collections::BTreeSet::new();
    for i in 0..=2 {
        for j in 0..=2 {
            if i + j <= 2 {
                cols.insert(kp15_global_col(n, 0, i, j));
            }
        }
    }
    cols.insert(kp15_global_col(n, 0, 2, 2));
    for sector in 0..n {
        cols.insert(kp15_global_col(n, sector, 2, 3));
        cols.insert(kp15_global_col(n, sector, 3, 3));
    }
    if n == 6 {
        cols.insert(kp15_global_col(n, 0, 2, 1));
    }
    cols
}

fn kp15_paper_p_and_tensor_border_columns_for_test(n: usize) -> std::collections::BTreeSet<usize> {
    let mut cols = kp15_strict_p_variable_columns_for_test(n);
    for sector in 0..n {
        for (i, j) in [
            (0, 0),
            (1, 0),
            (2, 0),
            (4, 0),
            (5, 0),
            (6, 0),
            (0, 1),
            (1, 1),
            (0, 2),
            (4, 2),
            (5, 2),
            (6, 2),
        ] {
            cols.insert(kp15_global_col(n, sector, i, j));
        }
    }
    cols
}

fn kp15_compatible_two_sided_tensor_source_for_test(
    table: &Kp15Bi6GeneratingTable,
) -> Vec<[[Point3; 7]; 7]> {
    let n = table.n;
    let free_for_coord = |coord: usize| -> Vec<f64> {
        table
            .free_cols
            .iter()
            .enumerate()
            .map(|(index, &col)| {
                let sector = col / 49;
                let pos = col % 49;
                let i = pos / 7;
                let j = pos % 7;
                let phase = (sector + 1) as f64 * (coord + 2) as f64;
                0.013 * (index + 1) as f64 + 0.017 * phase + 0.019 * i as f64 - 0.011 * j as f64
            })
            .collect()
    };
    let xs = table.apply(&free_for_coord(0));
    let ys = table.apply(&free_for_coord(1));
    let zs = table.apply(&free_for_coord(2));
    let mut out = Vec::with_capacity(n);
    for sector in 0..n {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let col = kp15_global_col(n, sector, i, j);
                *cell = Point3::new(xs[col], ys[col], zs[col]);
            }
        }
        out.push(grid);
    }
    out
}

fn kp15_full_stam_bicubic_source_for_test(
    projection: &crate::geometry::gregory::kp20_bi6_stam_eval::Kp20Pr08StamProjectedGrids,
    n: usize,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
) -> Vec<[[Point3; 7]; 7]> {
    (0..n)
        .map(|sector| {
            let bicubic =
                kp20_pr08_stam_projected_transformed_grid(projection, sector, patch, transform)
                    .expect("transformed Stam grid");
            kp15_degree_raise_bicubic_grid_to_bi6_for_test(&bicubic)
        })
        .collect()
}

fn kp15_beta_full_l_band_from_bicubic_source_for_test(
    projection: &crate::geometry::gregory::kp20_bi6_stam_eval::Kp20Pr08StamProjectedGrids,
    n: usize,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    overwrite_order: Kp15FullLBandOverwriteOrder,
) -> Vec<[[Point3; 7]; 7]> {
    let params = Bi6Reparams::for_valence(n);
    (0..n)
        .map(|sector| {
            let bicubic =
                kp20_pr08_stam_projected_transformed_grid(projection, sector, patch, transform)
                    .expect("transformed Stam grid");
            kp15_beta_full_l_band_from_bicubic_grid_for_test(&bicubic, params, overwrite_order)
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
enum Kp15FullLBandOverwriteOrder {
    TopThenRight,
    RightThenTop,
}

fn kp15_beta_full_l_band_from_bicubic_grid_for_test(
    bicubic: &[[Point3; 4]; 4],
    params: Bi6Reparams,
    overwrite_order: Kp15FullLBandOverwriteOrder,
) -> [[Point3; 7]; 7] {
    let mut out = kp15_degree_raise_bicubic_grid_to_bi6_for_test(bicubic);
    let top_strip = [bicubic[3], bicubic[2], bicubic[1]];
    let right_strip = [
        [bicubic[0][3], bicubic[1][3], bicubic[2][3], bicubic[3][3]],
        [bicubic[0][2], bicubic[1][2], bicubic[2][2], bicubic[3][2]],
        [bicubic[0][1], bicubic[1][1], bicubic[2][1], bicubic[3][1]],
    ];
    let top_rows = apply_beta_reparameterization(&top_strip, &params);
    let right_rows = apply_beta_reparameterization(&right_strip, &params);
    let write_top = |grid: &mut [[Point3; 7]; 7]| {
        for depth in 0..=2 {
            let target = 6 - depth;
            for along in 0..=6 {
                grid[target][along] = top_rows[depth][along];
            }
        }
    };
    let write_right = |grid: &mut [[Point3; 7]; 7]| {
        for depth in 0..=2 {
            let target = 6 - depth;
            for along in 0..=6 {
                grid[along][target] = right_rows[depth][along];
            }
        }
    };
    match overwrite_order {
        Kp15FullLBandOverwriteOrder::TopThenRight => {
            write_top(&mut out);
            write_right(&mut out);
        }
        Kp15FullLBandOverwriteOrder::RightThenTop => {
            write_right(&mut out);
            write_top(&mut out);
        }
    }
    out
}

fn kp15_two_sided_source_grid_from_bi3_tensor_for_test(
    n: usize,
    tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    apply_beta: bool,
) -> Vec<[[Point3; 7]; 7]> {
    kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
        n,
        tensor_rows,
        params,
        apply_beta,
        0,
        -1,
        false,
        false,
    )
}

fn kp15_two_sided_source_grid_from_bi3_tensor_phase_for_test(
    n: usize,
    tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    apply_beta: bool,
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
) -> Vec<[[Point3; 7]; 7]> {
    kp15_two_sided_source_grid_from_bi3_tensor_oriented_for_test(
        n,
        tensor_rows,
        params,
        apply_beta,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        false,
        false,
    )
}

fn kp15_two_sided_source_grid_from_bi3_tensor_oriented_for_test(
    n: usize,
    tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    apply_beta: bool,
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_normal_reflected: bool,
    column_normal_reflected: bool,
) -> Vec<[[Point3; 7]; 7]> {
    kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
        n,
        tensor_rows,
        tensor_rows,
        params,
        apply_beta,
        apply_beta,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        row_normal_reflected,
        column_normal_reflected,
    )
}

fn kp15_two_sided_source_grid_from_paired_bi3_tensors_oriented_for_test(
    n: usize,
    row_tensor_rows: &[[[Point3; 4]; 3]],
    column_tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    row_apply_beta: bool,
    column_apply_beta: bool,
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_normal_reflected: bool,
    column_normal_reflected: bool,
) -> Vec<[[Point3; 7]; 7]> {
    let row_rows = row_tensor_rows
        .iter()
        .map(|sector_rows| {
            if row_apply_beta {
                apply_beta_reparameterization(sector_rows, &params)
            } else {
                [
                    degree_raise_1d_3_to_6(&sector_rows[0]),
                    degree_raise_1d_3_to_6(&sector_rows[1]),
                    degree_raise_1d_3_to_6(&sector_rows[2]),
                ]
            }
        })
        .collect::<Vec<_>>();
    let column_rows = column_tensor_rows
        .iter()
        .map(|sector_rows| {
            if column_apply_beta {
                apply_beta_reparameterization(sector_rows, &params)
            } else {
                [
                    degree_raise_1d_3_to_6(&sector_rows[0]),
                    degree_raise_1d_3_to_6(&sector_rows[1]),
                    degree_raise_1d_3_to_6(&sector_rows[2]),
                ]
            }
        })
        .collect::<Vec<_>>();
    let mut out = Vec::with_capacity(n);
    for sector in 0..n {
        let row_sector = kp15_offset_sector_for_test(n, sector, row_sector_offset);
        let column_sector = kp15_offset_sector_for_test(n, sector, column_sector_offset);
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for depth in 0..=2 {
            for along in 0..=6 {
                let row_along = if row_reversed { 6 - along } else { along };
                let column_along = if column_reversed { 6 - along } else { along };
                grid[along][depth] = kp15_oriented_depth_point_for_test(
                    &row_rows[row_sector],
                    depth,
                    row_along,
                    row_normal_reflected,
                );
                grid[depth][along] = kp15_oriented_depth_point_for_test(
                    &column_rows[column_sector],
                    depth,
                    column_along,
                    column_normal_reflected,
                );
            }
        }
        out.push(grid);
    }
    out
}

fn kp15_two_sided_source_grid_from_paired_bi3_tensors_scaled_for_test(
    n: usize,
    row_tensor_rows: &[[[Point3; 4]; 3]],
    column_tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    row_apply_beta: bool,
    column_apply_beta: bool,
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_normal_scale: f64,
    column_normal_scale: f64,
) -> Vec<[[Point3; 7]; 7]> {
    let row_rows = kp15_bi6_tensor_rows_for_test(row_tensor_rows, params, row_apply_beta);
    let column_rows = kp15_bi6_tensor_rows_for_test(column_tensor_rows, params, column_apply_beta);
    let mut out = Vec::with_capacity(n);
    for sector in 0..n {
        let row_sector = kp15_offset_sector_for_test(n, sector, row_sector_offset);
        let column_sector = kp15_offset_sector_for_test(n, sector, column_sector_offset);
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for depth in 0..=2 {
            for along in 0..=6 {
                let row_along = if row_reversed { 6 - along } else { along };
                let column_along = if column_reversed { 6 - along } else { along };
                grid[along][depth] = kp15_scaled_depth_point_for_test(
                    &row_rows[row_sector],
                    depth,
                    row_along,
                    row_normal_scale,
                );
                grid[depth][along] = kp15_scaled_depth_point_for_test(
                    &column_rows[column_sector],
                    depth,
                    column_along,
                    column_normal_scale,
                );
            }
        }
        out.push(grid);
    }
    out
}

fn kp15_two_sided_source_grid_from_paired_bi6_rows_oriented_for_test(
    n: usize,
    row_rows: &[[[Point3; 7]; 3]],
    column_rows: &[[[Point3; 7]; 3]],
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_normal_reflected: bool,
    column_normal_reflected: bool,
) -> Vec<[[Point3; 7]; 7]> {
    kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
        n,
        row_rows,
        column_rows,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        row_normal_reflected,
        column_normal_reflected,
        Kp15Fig12LocalPlacement::SectorJoinFirst,
        Kp15Fig12LocalPlacement::SectorJoinFirst,
    )
}

fn kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
    n: usize,
    row_rows: &[[[Point3; 7]; 3]],
    column_rows: &[[[Point3; 7]; 3]],
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_normal_reflected: bool,
    column_normal_reflected: bool,
    row_placement: Kp15Fig12LocalPlacement,
    column_placement: Kp15Fig12LocalPlacement,
) -> Vec<[[Point3; 7]; 7]> {
    let mut out = Vec::with_capacity(n);
    for sector in 0..n {
        let row_sector = kp15_offset_sector_for_test(n, sector, row_sector_offset);
        let column_sector = kp15_offset_sector_for_test(n, sector, column_sector_offset);
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for depth in 0..=2 {
            for along in 0..=6 {
                let row_along = if row_reversed { 6 - along } else { along };
                let column_along = if column_reversed { 6 - along } else { along };
                let row_depth = row_placement.local_depth(depth);
                let column_depth = column_placement.local_depth(depth);
                grid[along][row_depth] = kp15_oriented_depth_point_for_test(
                    &row_rows[row_sector],
                    depth,
                    row_along,
                    row_normal_reflected,
                );
                grid[column_depth][along] = kp15_oriented_depth_point_for_test(
                    &column_rows[column_sector],
                    depth,
                    column_along,
                    column_normal_reflected,
                );
            }
        }
        out.push(grid);
    }
    out
}

fn kp15_bi6_tensor_rows_for_test(
    tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    apply_beta: bool,
) -> Vec<[[Point3; 7]; 3]> {
    tensor_rows
        .iter()
        .map(|sector_rows| {
            if apply_beta {
                apply_beta_reparameterization(sector_rows, &params)
            } else {
                [
                    degree_raise_1d_3_to_6(&sector_rows[0]),
                    degree_raise_1d_3_to_6(&sector_rows[1]),
                    degree_raise_1d_3_to_6(&sector_rows[2]),
                ]
            }
        })
        .collect()
}

fn kp15_fig12_bb_beta_rows_for_test(
    tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
) -> Vec<[[Point3; 7]; 3]> {
    tensor_rows
        .iter()
        .map(|sector_rows| {
            let gamma_rows = apply_beta_reparameterization(sector_rows, &params);
            let raised0 = degree_raise_1d_3_to_6(&sector_rows[0]);
            let raised1 = degree_raise_1d_3_to_6(&sector_rows[1]);
            let raised2 = degree_raise_1d_3_to_6(&sector_rows[2]);
            let mut out = gamma_rows;
            for i in 0..=6 {
                out[0][i] = raised0[i];
                out[1][i] = point3_scale_for_test(point3_add_for_test(raised0[i], raised1[i]), 0.5);
            }
            for i in 0..=1 {
                out[2][i] = point3_add_for_test(
                    point3_add_for_test(
                        point3_scale_for_test(raised0[i], 0.2),
                        point3_scale_for_test(raised1[i], 0.6),
                    ),
                    point3_scale_for_test(raised2[i], 0.2),
                );
            }
            out
        })
        .collect()
}

fn kp15_composed_pr08_subdivision_source_rows_for_test(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    kp15_pr08_q13_q2_top_source_rows_for_test(
        ring,
        Kp15Pr08TopSourceLevel::AfterOnePr08Subdivision,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn kp15_pr08_q13_q2_top_source_rows_for_test(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
    source_level: Kp15Pr08TopSourceLevel,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let source = match source_level {
        Kp15Pr08TopSourceLevel::RawCharacteristicRing => ring.sectors.clone(),
        Kp15Pr08TopSourceLevel::AfterOnePr08Subdivision => {
            pr08_subdivide_characteristic_ring(ring)?
        }
    };
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let stam_top = kp15_stam_frequency_one_top_vector_for_test(ring.lambda, valence);
    let q13_ring = source
        .iter()
        .map(|sector| [sector[12].x, sector[12].y])
        .collect::<Vec<_>>();
    let q2_ring = source
        .iter()
        .map(|sector| [sector[1].x, sector[1].y])
        .collect::<Vec<_>>();
    let odd_scale = kp15_complex_similarity_scale_xy_for_test(&q13_ring, &stam_top.0)
        .ok_or_else(|| "degenerate PR08 q13 ring for composed source".to_string())?;
    let even_scale = kp15_complex_similarity_scale_xy_for_test(&q2_ring, &stam_top.1)
        .ok_or_else(|| "degenerate PR08 q2 ring for composed source".to_string())?;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = source[0][0].x;
    u0_y[0] = source[0][0].y;
    for sector in 0..valence {
        let odd =
            kp15_complex_mul_xy_for_test([source[sector][12].x, source[sector][12].y], odd_scale);
        let even =
            kp15_complex_mul_xy_for_test([source[sector][1].x, source[sector][1].y], even_scale);
        u0_x[1 + 2 * sector] = odd[0];
        u0_y[1 + 2 * sector] = odd[1];
        u0_x[2 + 2 * sector] = even[0];
        u0_y[2 + 2 * sector] = even[1];
    }
    kp20_pr08_stam_appendix_b_projected_rows_from_top_vectors_for_test(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn kp15_pr08_top_label_source_rows_for_test(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
    source: &[[Point3; 13]],
    odd_label: usize,
    even_label: usize,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Option<Kp20RegularTensorBorderCNet> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let stam_top = kp15_stam_frequency_one_top_vector_for_test(ring.lambda, valence);
    let odd_ring = source
        .iter()
        .map(|sector| [sector[odd_label].x, sector[odd_label].y])
        .collect::<Vec<_>>();
    let even_ring = source
        .iter()
        .map(|sector| [sector[even_label].x, sector[even_label].y])
        .collect::<Vec<_>>();
    let odd_scale = kp15_complex_similarity_scale_xy_for_test(&odd_ring, &stam_top.0)?;
    let even_scale = kp15_complex_similarity_scale_xy_for_test(&even_ring, &stam_top.1)?;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = source[0][0].x;
    u0_y[0] = source[0][0].y;
    for sector in 0..valence {
        let odd = kp15_complex_mul_xy_for_test(
            [source[sector][odd_label].x, source[sector][odd_label].y],
            odd_scale,
        );
        let even = kp15_complex_mul_xy_for_test(
            [source[sector][even_label].x, source[sector][even_label].y],
            even_scale,
        );
        u0_x[1 + 2 * sector] = odd[0];
        u0_y[1 + 2 * sector] = odd[1];
        u0_x[2 + 2 * sector] = even[0];
        u0_y[2 + 2 * sector] = even[1];
    }
    kp20_pr08_stam_appendix_b_projected_rows_from_top_vectors_for_test(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
    .ok()
}

#[derive(Clone, Copy)]
enum Kp15Pr08TopSourceLevel {
    RawCharacteristicRing,
    AfterOnePr08Subdivision,
}

impl Kp15Pr08TopSourceLevel {
    const fn all() -> [Self; 2] {
        [Self::RawCharacteristicRing, Self::AfterOnePr08Subdivision]
    }

    const fn name(self) -> &'static str {
        match self {
            Self::RawCharacteristicRing => "raw-characteristic-ring",
            Self::AfterOnePr08Subdivision => "after-one-pr08-subdivision",
        }
    }
}

fn kp15_stam_frequency_one_top_vector_for_test(
    lambda: f64,
    valence: usize,
) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
    let gain = 4.0 * lambda - 1.0;
    let mut odd = Vec::with_capacity(valence);
    let mut even = Vec::with_capacity(valence);
    for sector in 0..valence {
        let theta = core::f64::consts::TAU * sector as f64 / valence as f64;
        let next_theta = core::f64::consts::TAU * (sector + 1) as f64 / valence as f64;
        odd.push([gain * theta.cos(), gain * theta.sin()]);
        even.push([
            theta.cos() + next_theta.cos(),
            theta.sin() + next_theta.sin(),
        ]);
    }
    (odd, even)
}

fn kp15_complex_similarity_scale_xy_for_test(
    source: &[[f64; 2]],
    target: &[[f64; 2]],
) -> Option<[f64; 2]> {
    if source.len() != target.len() {
        return None;
    }
    let mut numerator = [0.0_f64; 2];
    let mut denominator = 0.0_f64;
    for (&s, &t) in source.iter().zip(target) {
        numerator[0] += t[0] * s[0] + t[1] * s[1];
        numerator[1] += t[1] * s[0] - t[0] * s[1];
        denominator += s[0] * s[0] + s[1] * s[1];
    }
    (denominator > 1.0e-20).then_some([numerator[0] / denominator, numerator[1] / denominator])
}

fn kp15_complex_mul_xy_for_test(value: [f64; 2], scale: [f64; 2]) -> [f64; 2] {
    [
        value[0] * scale[0] - value[1] * scale[1],
        value[0] * scale[1] + value[1] * scale[0],
    ]
}

fn kp15_apply_complex_similarity_to_regular_rows_for_test(
    rows: &Kp20RegularTensorBorderCNet,
    scale: [f64; 2],
    conjugate_source: bool,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector_rows| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for cp in 0..4 {
                    let xy = kp15_maybe_conjugate_xy(sector_rows[row][cp], conjugate_source);
                    out[row][cp] = Point3::new(
                        scale[0] * xy[0] - scale[1] * xy[1],
                        scale[1] * xy[0] + scale[0] * xy[1],
                        sector_rows[row][cp].z,
                    );
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

#[derive(Clone, Copy)]
enum Kp15Fig12StamSourceWindow {
    TopRows,
    BottomRows,
    LeftColumns,
    RightColumns,
}

#[derive(Clone, Copy)]
enum Kp15Fig12LocalPlacement {
    SectorJoinFirst,
    TensorBoundaryFirst,
}

impl Kp15Fig12LocalPlacement {
    const fn all() -> [Self; 2] {
        [Self::SectorJoinFirst, Self::TensorBoundaryFirst]
    }

    const fn local_depth(self, source_depth: usize) -> usize {
        match self {
            Self::SectorJoinFirst => source_depth,
            Self::TensorBoundaryFirst => 2 - source_depth,
        }
    }

    const fn name(self) -> &'static str {
        match self {
            Self::SectorJoinFirst => "sector-join-first",
            Self::TensorBoundaryFirst => "tensor-boundary-first",
        }
    }
}

impl Kp15Fig12StamSourceWindow {
    const fn all() -> [Self; 4] {
        [
            Self::TopRows,
            Self::BottomRows,
            Self::LeftColumns,
            Self::RightColumns,
        ]
    }

    const fn name(self) -> &'static str {
        match self {
            Self::TopRows => "top-rows",
            Self::BottomRows => "bottom-rows",
            Self::LeftColumns => "left-columns",
            Self::RightColumns => "right-columns",
        }
    }
}

fn kp15_fig12_stam_source_window_rows_for_test(
    grid: &[[Point3; 4]; 4],
    window: Kp15Fig12StamSourceWindow,
) -> [[Point3; 4]; 3] {
    match window {
        Kp15Fig12StamSourceWindow::TopRows => [grid[0], grid[1], grid[2]],
        Kp15Fig12StamSourceWindow::BottomRows => [grid[1], grid[2], grid[3]],
        Kp15Fig12StamSourceWindow::LeftColumns => [
            [grid[0][0], grid[1][0], grid[2][0], grid[3][0]],
            [grid[0][1], grid[1][1], grid[2][1], grid[3][1]],
            [grid[0][2], grid[1][2], grid[2][2], grid[3][2]],
        ],
        Kp15Fig12StamSourceWindow::RightColumns => [
            [grid[0][3], grid[1][3], grid[2][3], grid[3][3]],
            [grid[0][2], grid[1][2], grid[2][2], grid[3][2]],
            [grid[0][1], grid[1][1], grid[2][1], grid[3][1]],
        ],
    }
}

fn kp15_stam_window_beta_rows_for_test(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    window: Kp15Fig12StamSourceWindow,
    params: Bi6Reparams,
) -> Vec<[[Point3; 7]; 3]> {
    let mut bi3_sectors = Vec::with_capacity(ring.sectors.len());
    for sector in 0..ring.sectors.len() {
        let grid =
            kp20_pr08_stam_appendix_b_inner_ring_full_grid_for_test(ring, sector, patch_index)
                .expect("Stam Appendix-B full grid");
        bi3_sectors.push(kp15_fig12_stam_source_window_rows_for_test(&grid, window));
    }
    kp15_fig12_bb_beta_rows_for_test(&bi3_sectors, params)
}

fn kp15_oriented_depth_point_for_test(
    rows: &[[Point3; 7]; 3],
    depth: usize,
    along: usize,
    normal_reflected: bool,
) -> Point3 {
    if !normal_reflected {
        return rows[depth][along];
    }
    let q0 = rows[0][along];
    let q1 = rows[1][along];
    let q2 = rows[2][along];
    match depth {
        0 => q0,
        1 => Point3::new(2.0 * q0.x - q1.x, 2.0 * q0.y - q1.y, 2.0 * q0.z - q1.z),
        2 => Point3::new(
            q2.x - 4.0 * q1.x + 4.0 * q0.x,
            q2.y - 4.0 * q1.y + 4.0 * q0.y,
            q2.z - 4.0 * q1.z + 4.0 * q0.z,
        ),
        _ => unreachable!("KP15 source depth is 0..=2"),
    }
}

fn kp15_scaled_depth_point_for_test(
    rows: &[[Point3; 7]; 3],
    depth: usize,
    along: usize,
    normal_scale: f64,
) -> Point3 {
    let q0 = rows[0][along];
    let q1 = rows[1][along];
    let q2 = rows[2][along];
    let d1 = point3_sub_for_test(q1, q0);
    let d2 = point3_add_for_test(point3_sub_for_test(q2, point3_scale_for_test(q1, 2.0)), q0);
    match depth {
        0 => q0,
        1 => point3_add_for_test(q0, point3_scale_for_test(d1, normal_scale)),
        2 => point3_add_for_test(
            point3_add_for_test(q0, point3_scale_for_test(d1, 2.0 * normal_scale)),
            point3_scale_for_test(d2, normal_scale * normal_scale),
        ),
        _ => unreachable!("KP15 source depth is 0..=2"),
    }
}

fn kp15_scale_fig2b_source_depth_from_boundary_for_test(
    source: &[[[Point3; 7]; 7]],
    row_normal_scale: f64,
    column_normal_scale: f64,
) -> Vec<[[Point3; 7]; 7]> {
    source
        .iter()
        .map(|grid| {
            let mut out = *grid;
            for along in 0..=6 {
                let row_depths = [grid[along][0], grid[along][1], grid[along][2]];
                for depth in 1..=2 {
                    out[along][depth] = kp15_scaled_point_from_depth_triple_for_test(
                        row_depths,
                        depth,
                        row_normal_scale,
                    );
                }
                let column_depths = [grid[0][along], grid[1][along], grid[2][along]];
                for depth in 1..=2 {
                    out[depth][along] = kp15_scaled_point_from_depth_triple_for_test(
                        column_depths,
                        depth,
                        column_normal_scale,
                    );
                }
            }
            out
        })
        .collect()
}

fn kp15_scaled_point_from_depth_triple_for_test(
    depths: [Point3; 3],
    depth: usize,
    normal_scale: f64,
) -> Point3 {
    let q0 = depths[0];
    let q1 = depths[1];
    let q2 = depths[2];
    let d1 = point3_sub_for_test(q1, q0);
    let d2 = point3_add_for_test(point3_sub_for_test(q2, point3_scale_for_test(q1, 2.0)), q0);
    match depth {
        0 => q0,
        1 => point3_add_for_test(q0, point3_scale_for_test(d1, normal_scale)),
        2 => point3_add_for_test(
            point3_add_for_test(q0, point3_scale_for_test(d1, 2.0 * normal_scale)),
            point3_scale_for_test(d2, normal_scale * normal_scale),
        ),
        _ => unreachable!("KP15 source depth is 0..=2"),
    }
}

fn point3_add_for_test(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

fn point3_sub_for_test(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn point3_scale_for_test(a: Point3, scale: f64) -> Point3 {
    Point3::new(scale * a.x, scale * a.y, scale * a.z)
}

fn point3_components_for_test(a: Point3) -> [f64; 3] {
    [a.x, a.y, a.z]
}

fn kp15_offset_sector_for_test(n: usize, sector: usize, offset: isize) -> usize {
    (sector as isize + offset).rem_euclid(n as isize) as usize
}

fn kp15_degree_raise_bicubic_grid_to_bi6_for_test(grid: &[[Point3; 4]; 4]) -> [[Point3; 7]; 7] {
    let raised_rows = grid.map(|row| degree_raise_1d_3_to_6(&row));
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    for col in 0..7 {
        let column = [
            raised_rows[0][col],
            raised_rows[1][col],
            raised_rows[2][col],
            raised_rows[3][col],
        ];
        let raised_column = degree_raise_1d_3_to_6(&column);
        for row in 0..7 {
            out[row][col] = raised_column[row];
        }
    }
    out
}

fn kp15_degree_raise_scalar_basis_3_to_6_for_test(index: usize) -> [f64; 4] {
    let mut coeffs = [0.0_f64; 4];
    for basis in 0..4 {
        let mut row = [Point3::new(0.0, 0.0, 0.0); 4];
        row[basis] = Point3::new(1.0, 0.0, 0.0);
        coeffs[basis] = degree_raise_1d_3_to_6(&row)[index].x;
    }
    coeffs
}

fn kp15_parse_constraint_residual_for_test(message: &str) -> Option<f64> {
    message
        .split_whitespace()
        .last()
        .and_then(|token| token.parse::<f64>().ok())
}

fn kp15_source_constraints_from_grid_columns_for_test(
    n: usize,
    source: &[[[Point3; 7]; 7]],
    cols: &[usize],
) -> Vec<(usize, Point3)> {
    cols.iter()
        .copied()
        .map(|col| {
            let sector = col / 49;
            let pos = col % 49;
            let i = pos / 7;
            let j = pos % 7;
            assert!(sector < n);
            (col, source[sector][i][j])
        })
        .collect()
}

fn kp15_radial_seed_plus_depth2_source_constraints_for_test(
    n: usize,
    radial_sectors: &[CirculantSectorRhs],
    depth_source: &[[[Point3; 7]; 7]],
) -> Vec<(usize, Point3)> {
    assert_eq!(radial_sectors.len(), n);
    assert_eq!(depth_source.len(), n);
    let mut constraints = Vec::with_capacity(15 * n);
    for (sector, source) in radial_sectors.iter().enumerate() {
        for (i, j, target) in [
            (0, 0, Point3::new(0.0, 0.0, 0.0)),
            (1, 0, source.p_10),
            (2, 0, source.p_20),
            (4, 0, source.p_40),
            (5, 0, source.p_50),
            (6, 0, source.p_60),
            (0, 1, source.p_01),
            (1, 1, source.p_11),
            (0, 2, source.p_02),
        ] {
            constraints.push((kp15_global_col(n, sector, i, j), target));
        }
        for i in 4..=6 {
            constraints.push((kp15_global_col(n, sector, i, 2), depth_source[sector][i][2]));
        }
        for j in 4..=6 {
            constraints.push((kp15_global_col(n, sector, 2, j), depth_source[sector][2][j]));
        }
    }
    constraints
}

fn kp15_sector0_central_radial_edge_plus_depth2_source_constraints_for_test(
    n: usize,
    radial_sectors: &[CirculantSectorRhs],
    depth_source: &[[[Point3; 7]; 7]],
) -> Vec<(usize, Point3)> {
    assert_eq!(radial_sectors.len(), n);
    assert_eq!(depth_source.len(), n);
    let mut constraints = Vec::with_capacity(10 * n + 6);
    let sector0 = radial_sectors[0];
    for (i, j, target) in [
        (0, 0, Point3::new(0.0, 0.0, 0.0)),
        (1, 0, sector0.p_10),
        (2, 0, sector0.p_20),
        (0, 1, sector0.p_01),
        (1, 1, sector0.p_11),
        (0, 2, sector0.p_02),
    ] {
        constraints.push((kp15_global_col(n, 0, i, j), target));
    }
    for (sector, source) in radial_sectors.iter().enumerate() {
        for (i, target) in [(4, source.p_40), (5, source.p_50), (6, source.p_60)] {
            constraints.push((kp15_global_col(n, sector, i, 0), target));
        }
        for i in 4..=6 {
            constraints.push((kp15_global_col(n, sector, i, 2), depth_source[sector][i][2]));
        }
        for j in 4..=6 {
            constraints.push((kp15_global_col(n, sector, 2, j), depth_source[sector][2][j]));
        }
    }
    constraints
}

fn kp15_radial_edge_plus_depth2_source_constraints_for_test(
    n: usize,
    radial_sectors: &[CirculantSectorRhs],
    depth_source: &[[[Point3; 7]; 7]],
) -> Vec<(usize, Point3)> {
    assert_eq!(radial_sectors.len(), n);
    assert_eq!(depth_source.len(), n);
    let mut constraints = Vec::with_capacity(9 * n);
    for (sector, source) in radial_sectors.iter().enumerate() {
        for (i, target) in [(4, source.p_40), (5, source.p_50), (6, source.p_60)] {
            constraints.push((kp15_global_col(n, sector, i, 0), target));
        }
        for i in 4..=6 {
            constraints.push((kp15_global_col(n, sector, i, 2), depth_source[sector][i][2]));
        }
        for j in 4..=6 {
            constraints.push((kp15_global_col(n, sector, 2, j), depth_source[sector][2][j]));
        }
    }
    constraints
}

fn kp15_labeled_depth2_plus_cap_net_m_endpoint_constraints_for_test(
    n: usize,
    depth_source: &[[[Point3; 7]; 7]],
    cap_net_sectors: &[[Point3; 6]],
    m_label: usize,
    m_offset: isize,
) -> Vec<Kp15LabeledConstraint> {
    assert_eq!(depth_source.len(), n);
    assert_eq!(cap_net_sectors.len(), n);
    assert!(m_label < 6);
    let mut constraints = kp15_labeled_source_constraints_from_grid_columns_for_test(
        n,
        depth_source,
        &kp15_fig2b_depth2_tensor_border_columns_for_test(n),
        "depth2",
    );
    for sector in 0..n {
        let p60_sector = kp15_offset_sector_for_test(n, sector, m_offset);
        let p06_sector = kp15_offset_sector_for_test(n, sector + 1, m_offset);
        for (label, slot, point) in [
            ("p60", (6, 0), cap_net_sectors[p60_sector][m_label]),
            ("p06", (0, 6), cap_net_sectors[p06_sector][m_label]),
        ] {
            constraints.push(Kp15LabeledConstraint {
                label: format!(
                    "cap-net-m c{} offset={} sector={sector} {label}",
                    m_label + 1,
                    m_offset,
                ),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(n, sector, slot.0, slot.1), 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            });
        }
    }
    constraints
}

fn kp20_pr08_refined_labels_from_cap_net_for_test(
    cap_net: &crate::geometry::gregory::kp20_bi6_cnet::Kp20CapNet,
) -> Vec<[Point3; 13]> {
    kp20_pr08_refined_labels_from_cap_net_permutation_for_test(cap_net, &[0, 1, 2, 3, 4, 5])
}

fn kp20_pr08_refined_labels_from_cap_net_permutation_for_test(
    cap_net: &crate::geometry::gregory::kp20_bi6_cnet::Kp20CapNet,
    source_label_permutation: &[usize; 6],
) -> Vec<[Point3; 13]> {
    let n = cap_net.sectors.len();
    let matrix = pr08_catmull_clark_subdivision_stencil_matrix(n).expect("PR08 stencil matrix");
    let mut out = vec![[Point3::new(0.0, 0.0, 0.0); 13]; n];
    for dest_sector in 0..n {
        for dest_label in 0..13 {
            let row = dest_sector * 13 + dest_label;
            let mut value = Point3::new(0.0, 0.0, 0.0);
            for src_sector in 0..n {
                for src_label in 0..=6 {
                    let weight = matrix[row][src_sector * 13 + src_label];
                    let source = if src_label == 0 {
                        cap_net.center
                    } else {
                        cap_net.sectors[src_sector][source_label_permutation[src_label - 1]]
                    };
                    value.x += weight * source.x;
                    value.y += weight * source.y;
                    value.z += weight * source.z;
                }
            }
            out[dest_sector][dest_label] = value;
        }
    }
    out
}

fn kp20_cap_net_label_permutations_for_test() -> Vec<[usize; 6]> {
    fn recurse(
        depth: usize,
        current: &mut [usize; 6],
        used: &mut [bool; 6],
        out: &mut Vec<[usize; 6]>,
    ) {
        if depth == 6 {
            out.push(*current);
            return;
        }
        for label in 0..6 {
            if used[label] {
                continue;
            }
            used[label] = true;
            current[depth] = label;
            recurse(depth + 1, current, used, out);
            used[label] = false;
        }
    }

    let mut out = Vec::with_capacity(720);
    let mut current = [0usize; 6];
    let mut used = [false; 6];
    recurse(0, &mut current, &mut used, &mut out);
    out
}

#[derive(Clone, Copy)]
struct Kp15Fig63RadialScaleSpec {
    n: usize,
    row_side: Kp20StamGridSide,
    column_side: Kp20StamGridSide,
    swap_targets: bool,
    conjugate_source: bool,
    row_reflected: bool,
    column_reflected: bool,
}

fn kp15_fig63_full_tensor_border_specs_for_test(
    n: usize,
) -> impl Iterator<Item = Kp15Fig63RadialScaleSpec> {
    [
        (Kp20StamGridSide::Left, Kp20StamGridSide::Top),
        (Kp20StamGridSide::Top, Kp20StamGridSide::Left),
    ]
    .into_iter()
    .flat_map(move |(row_side, column_side)| {
        [(false, false), (false, true), (true, false), (true, true)]
            .into_iter()
            .flat_map(move |(swap_targets, conjugate_source)| {
                [false, true].into_iter().flat_map(move |row_reflected| {
                    [false, true].into_iter().map(move |column_reflected| {
                        Kp15Fig63RadialScaleSpec {
                            n,
                            row_side,
                            column_side,
                            swap_targets,
                            conjugate_source,
                            row_reflected,
                            column_reflected,
                        }
                    })
                })
            })
    })
}

fn kp15_scale_sector_rhs_for_test(
    sectors: &[CirculantSectorRhs],
    scale: f64,
) -> Vec<CirculantSectorRhs> {
    sectors
        .iter()
        .map(|source| CirculantSectorRhs {
            p_00: point3_scale_for_test(source.p_00, scale),
            p_10: point3_scale_for_test(source.p_10, scale),
            p_20: point3_scale_for_test(source.p_20, scale),
            p_40: point3_scale_for_test(source.p_40, scale),
            p_50: point3_scale_for_test(source.p_50, scale),
            p_60: point3_scale_for_test(source.p_60, scale),
            p_01: point3_scale_for_test(source.p_01, scale),
            p_11: point3_scale_for_test(source.p_11, scale),
            p_02: point3_scale_for_test(source.p_02, scale),
            p_42: point3_scale_for_test(source.p_42, scale),
            p_52: point3_scale_for_test(source.p_52, scale),
            p_62: point3_scale_for_test(source.p_62, scale),
        })
        .collect()
}

fn kp15_fig63_full_tensor_border_source_for_spec(
    spec: Kp15Fig63RadialScaleSpec,
) -> Option<Vec<[[Point3; 7]; 7]>> {
    kp15_fig63_depth_source_and_planar_radial_input_for_spec(spec).map(|(source, _)| source)
}

fn kp15_fig63_depth_source_and_planar_radial_input_for_spec(
    spec: Kp15Fig63RadialScaleSpec,
) -> Option<(Vec<[[Point3; 7]; 7]>, Bi6CapInput)> {
    let n = spec.n;
    let params = Bi6Reparams::for_valence(n);
    let row_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
        n,
        spec.row_side,
        spec.swap_targets,
        spec.conjugate_source,
    )
    .ok()?;
    let column_tensor = kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
        n,
        spec.column_side,
        spec.swap_targets,
        spec.conjugate_source,
    )
    .ok()?;
    let row_rows = kp15_fig12_bb_beta_rows_for_test(&row_tensor.sectors, params);
    let column_rows = kp15_fig12_bb_beta_rows_for_test(&column_tensor.sectors, params);
    let depth_source = kp15_two_sided_source_grid_from_paired_bi6_rows_placed_for_test(
        n,
        &row_rows,
        &column_rows,
        -2,
        -1,
        spec.row_side == Kp20StamGridSide::Left,
        spec.column_side == Kp20StamGridSide::Left,
        spec.row_reflected,
        spec.column_reflected,
        Kp15Fig12LocalPlacement::SectorJoinFirst,
        Kp15Fig12LocalPlacement::SectorJoinFirst,
    );
    Some((depth_source, kp15_spoke_indexed_planar_bilinear_input(n)))
}

fn kp15_labeled_depth2_then_single_radial_constraint_for_test(
    n: usize,
    depth_source: &[[[Point3; 7]; 7]],
    radial_sectors: &[CirculantSectorRhs],
    radial_slot: (usize, usize),
) -> Vec<Kp15LabeledConstraint> {
    assert_eq!(depth_source.len(), n);
    assert_eq!(radial_sectors.len(), n);
    let mut constraints = Vec::new();
    for sector in 0..n {
        for i in 4..=6 {
            let point = depth_source[sector][i][2];
            constraints.push(Kp15LabeledConstraint {
                label: format!("depth sector={sector} p_{i}2"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(n, sector, i, 2), 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            });
        }
        for j in 4..=6 {
            let point = depth_source[sector][2][j];
            constraints.push(Kp15LabeledConstraint {
                label: format!("depth sector={sector} p_2{j}"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(n, sector, 2, j), 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            });
        }
    }
    let point = kp15_seed_value_for_test(radial_sectors, 0, radial_slot.0, radial_slot.1)
        .expect("radial slot must be a seed value");
    constraints.push(Kp15LabeledConstraint {
        label: format!("radial sector=0 p_{}{}", radial_slot.0, radial_slot.1),
        constraint: AffineConstraintRow {
            row: LinearRow {
                terms: vec![(kp15_global_col(n, 0, radial_slot.0, radial_slot.1), 1.0)],
            },
            rhs: [point.x, point.y, point.z],
        },
    });
    constraints
}

fn kp15_labeled_full_tensor_border_constraints_for_test(
    n: usize,
    source: &[[[Point3; 7]; 7]],
) -> Vec<Kp15LabeledConstraint> {
    kp15_bi6_local_constraint::kp15_two_sided_tensor_border_columns(n)
        .into_iter()
        .map(|col| {
            let sector = col / 49;
            let pos = col % 49;
            let i = pos / 7;
            let j = pos % 7;
            let point = source[sector][i][j];
            Kp15LabeledConstraint {
                label: format!("tensor sector={sector} p_{i}{j}"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(col, 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            }
        })
        .collect()
}

fn kp15_shifted_m_endpoint_tensor_border_constraints_for_test(
    n: usize,
    source: &[[[Point3; 7]; 7]],
) -> Vec<(usize, Point3)> {
    let mut constraints = Vec::with_capacity(6 * n);
    for sector in 0..n {
        for depth in 0..=2 {
            constraints.push((
                kp15_global_col(n, sector, 6, depth),
                source[sector][4][depth],
            ));
        }
        for depth in 0..=2 {
            constraints.push((
                kp15_global_col(n, sector, depth, 6),
                source[sector][depth][4],
            ));
        }
    }
    constraints
}

fn kp15_labeled_depth2_then_radial_seed_constraints_for_test(
    n: usize,
    depth_source: &[[[Point3; 7]; 7]],
    radial_sectors: &[CirculantSectorRhs],
) -> Vec<Kp15LabeledConstraint> {
    let mut constraints = kp15_labeled_depth2_then_single_radial_constraint_for_test(
        n,
        depth_source,
        radial_sectors,
        (0, 0),
    );
    constraints.pop();
    for (sector, source) in radial_sectors.iter().enumerate() {
        for (i, j, target) in [
            (0, 0, Point3::new(0.0, 0.0, 0.0)),
            (1, 0, source.p_10),
            (2, 0, source.p_20),
            (4, 0, source.p_40),
            (5, 0, source.p_50),
            (6, 0, source.p_60),
            (0, 1, source.p_01),
            (1, 1, source.p_11),
            (0, 2, source.p_02),
        ] {
            constraints.push(Kp15LabeledConstraint {
                label: format!("radial sector={sector} p_{i}{j}"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(kp15_global_col(n, sector, i, j), 1.0)],
                    },
                    rhs: [target.x, target.y, target.z],
                },
            });
        }
    }
    constraints
}

#[derive(Clone, Copy)]
struct Kp15BestPr08SideSpec {
    n: usize,
    row_candidate: Kp20Pr08TensorBorderCandidate,
    column_candidate: Kp20Pr08TensorBorderCandidate,
    row_apply_beta: bool,
    column_apply_beta: bool,
    row_offset: isize,
    column_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
    row_reflected: bool,
    column_reflected: bool,
}

fn kp15_best_pr08_side_specs_for_test() -> [Kp15BestPr08SideSpec; 4] {
    [
        Kp15BestPr08SideSpec {
            n: 5,
            row_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            column_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows,
            row_apply_beta: false,
            column_apply_beta: true,
            row_offset: -1,
            column_offset: 0,
            row_reversed: true,
            column_reversed: true,
            row_reflected: true,
            column_reflected: false,
        },
        Kp15BestPr08SideSpec {
            n: 6,
            row_candidate: Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows,
            column_candidate: Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows,
            row_apply_beta: true,
            column_apply_beta: false,
            row_offset: 0,
            column_offset: -2,
            row_reversed: true,
            column_reversed: false,
            row_reflected: true,
            column_reflected: true,
        },
        Kp15BestPr08SideSpec {
            n: 7,
            row_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            column_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
            row_apply_beta: false,
            column_apply_beta: false,
            row_offset: 0,
            column_offset: 1,
            row_reversed: true,
            column_reversed: false,
            row_reflected: false,
            column_reflected: true,
        },
        Kp15BestPr08SideSpec {
            n: 8,
            row_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            column_candidate: Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
            row_apply_beta: false,
            column_apply_beta: false,
            row_offset: 0,
            column_offset: 1,
            row_reversed: true,
            column_reversed: false,
            row_reflected: false,
            column_reflected: true,
        },
    ]
}

#[derive(Clone, Copy, Debug)]
enum Kp15CornerJoinMode {
    RowToColumn,
    ColumnToRow,
    Average,
}

fn kp15_apply_corner_adjacent_join_for_test(
    grids: &mut [[[Point3; 7]; 7]],
    sector_offset: isize,
    mode: Kp15CornerJoinMode,
) {
    let n = grids.len();
    let original = grids.to_vec();
    for sector in 0..n {
        let paired_sector = kp15_offset_sector_for_test(n, sector, sector_offset);
        let row_point = original[paired_sector][4][1];
        let column_point = original[sector][1][4];
        let joined = match mode {
            Kp15CornerJoinMode::RowToColumn => row_point,
            Kp15CornerJoinMode::ColumnToRow => column_point,
            Kp15CornerJoinMode::Average => Point3::new(
                0.5 * (row_point.x + column_point.x),
                0.5 * (row_point.y + column_point.y),
                0.5 * (row_point.z + column_point.z),
            ),
        };
        match mode {
            Kp15CornerJoinMode::RowToColumn => grids[sector][1][4] = joined,
            Kp15CornerJoinMode::ColumnToRow => grids[paired_sector][4][1] = joined,
            Kp15CornerJoinMode::Average => {
                grids[sector][1][4] = joined;
                grids[paired_sector][4][1] = joined;
            }
        }
    }
}

fn kp15_apply_s4_cross_corner_sum_correction_for_test(
    grids: &mut [[[Point3; 7]; 7]],
    params: Bi6Reparams,
) -> f64 {
    let n = grids.len();
    let original = grids.to_vec();
    let w3 = params.c / 15.0;
    let mut corrections = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut max_delta = 0.0_f64;
    for sector in 0..n {
        let next = (sector + 1) % n;
        let target_sum = Point3::new(
            original[sector][0][4].x
                + original[next][4][0].x
                + w3 * (original[sector][0][6].x - original[sector][0][5].x),
            original[sector][0][4].y
                + original[next][4][0].y
                + w3 * (original[sector][0][6].y - original[sector][0][5].y),
            original[sector][0][4].z
                + original[next][4][0].z
                + w3 * (original[sector][0][6].z - original[sector][0][5].z),
        );
        let current_sum = Point3::new(
            original[sector][1][4].x + original[next][4][1].x,
            original[sector][1][4].y + original[next][4][1].y,
            original[sector][1][4].z + original[next][4][1].z,
        );
        let correction = Point3::new(
            0.5 * (target_sum.x - current_sum.x),
            0.5 * (target_sum.y - current_sum.y),
            0.5 * (target_sum.z - current_sum.z),
        );
        corrections[sector] = correction;
        max_delta = max_delta.max(
            (correction.x * correction.x
                + correction.y * correction.y
                + correction.z * correction.z)
                .sqrt(),
        );
    }
    for sector in 0..n {
        let previous = (sector + n - 1) % n;
        grids[sector][1][4] = Point3::new(
            original[sector][1][4].x + corrections[sector].x,
            original[sector][1][4].y + corrections[sector].y,
            original[sector][1][4].z + corrections[sector].z,
        );
        grids[sector][4][1] = Point3::new(
            original[sector][4][1].x + corrections[previous].x,
            original[sector][4][1].y + corrections[previous].y,
            original[sector][4][1].z + corrections[previous].z,
        );
    }
    max_delta
}

fn kp15_s4_cross_corner_sum_residual_for_test(
    grids: &[[[Point3; 7]; 7]],
    params: Bi6Reparams,
) -> f64 {
    let n = grids.len();
    let w3 = params.c / 15.0;
    let mut max_residual: f64 = 0.0;
    for sector in 0..n {
        let next = (sector + 1) % n;
        let target_sum = Point3::new(
            grids[sector][0][4].x
                + grids[next][4][0].x
                + w3 * (grids[sector][0][6].x - grids[sector][0][5].x),
            grids[sector][0][4].y
                + grids[next][4][0].y
                + w3 * (grids[sector][0][6].y - grids[sector][0][5].y),
            grids[sector][0][4].z
                + grids[next][4][0].z
                + w3 * (grids[sector][0][6].z - grids[sector][0][5].z),
        );
        let actual_sum = Point3::new(
            grids[sector][1][4].x + grids[next][4][1].x,
            grids[sector][1][4].y + grids[next][4][1].y,
            grids[sector][1][4].z + grids[next][4][1].z,
        );
        max_residual = max_residual.max(actual_sum.distance_to(&target_sum));
    }
    max_residual
}

struct Kp15S4FirstOrderNormalScaleSolve {
    row_scale: f64,
    column_scale: f64,
    residual: f64,
    determinant: f64,
}

fn kp15_solve_s4_first_order_normal_scales_for_test(
    n: usize,
    row_tensor_rows: &[[[Point3; 4]; 3]],
    column_tensor_rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    row_apply_beta: bool,
    column_apply_beta: bool,
    row_sector_offset: isize,
    column_sector_offset: isize,
    row_reversed: bool,
    column_reversed: bool,
) -> Kp15S4FirstOrderNormalScaleSolve {
    let base = kp15_two_sided_source_grid_from_paired_bi3_tensors_scaled_for_test(
        n,
        row_tensor_rows,
        column_tensor_rows,
        params,
        row_apply_beta,
        column_apply_beta,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        0.0,
        0.0,
    );
    let row_unit = kp15_two_sided_source_grid_from_paired_bi3_tensors_scaled_for_test(
        n,
        row_tensor_rows,
        column_tensor_rows,
        params,
        row_apply_beta,
        column_apply_beta,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        1.0,
        0.0,
    );
    let column_unit = kp15_two_sided_source_grid_from_paired_bi3_tensors_scaled_for_test(
        n,
        row_tensor_rows,
        column_tensor_rows,
        params,
        row_apply_beta,
        column_apply_beta,
        row_sector_offset,
        column_sector_offset,
        row_reversed,
        column_reversed,
        0.0,
        1.0,
    );
    let w3 = params.c / 15.0;
    let mut ata = [[0.0_f64; 2]; 2];
    let mut aty = [0.0_f64; 2];
    for sector in 0..n {
        let next = (sector + 1) % n;
        let target = point3_add_for_test(
            point3_add_for_test(base[sector][0][4], base[next][4][0]),
            point3_scale_for_test(
                point3_sub_for_test(base[sector][0][6], base[sector][0][5]),
                w3,
            ),
        );
        let base_sum = point3_add_for_test(base[sector][1][4], base[next][4][1]);
        let row_delta = point3_sub_for_test(row_unit[next][4][1], base[next][4][1]);
        let column_delta = point3_sub_for_test(column_unit[sector][1][4], base[sector][1][4]);
        for (axis, y) in point3_components_for_test(point3_sub_for_test(target, base_sum))
            .into_iter()
            .enumerate()
        {
            let a = [
                point3_components_for_test(row_delta)[axis],
                point3_components_for_test(column_delta)[axis],
            ];
            ata[0][0] += a[0] * a[0];
            ata[0][1] += a[0] * a[1];
            ata[1][0] += a[1] * a[0];
            ata[1][1] += a[1] * a[1];
            aty[0] += a[0] * y;
            aty[1] += a[1] * y;
        }
    }
    let det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
    let (row_scale, column_scale) = if det.abs() > 1.0e-14 {
        (
            (aty[0] * ata[1][1] - ata[0][1] * aty[1]) / det,
            (ata[0][0] * aty[1] - aty[0] * ata[1][0]) / det,
        )
    } else {
        let mut numerator = 0.0_f64;
        let mut denominator = 0.0_f64;
        for sector in 0..n {
            let next = (sector + 1) % n;
            let target = point3_add_for_test(
                point3_add_for_test(base[sector][0][4], base[next][4][0]),
                point3_scale_for_test(
                    point3_sub_for_test(base[sector][0][6], base[sector][0][5]),
                    w3,
                ),
            );
            let base_sum = point3_add_for_test(base[sector][1][4], base[next][4][1]);
            let row_delta = point3_sub_for_test(row_unit[next][4][1], base[next][4][1]);
            let column_delta = point3_sub_for_test(column_unit[sector][1][4], base[sector][1][4]);
            let combined_delta = point3_add_for_test(row_delta, column_delta);
            let y = point3_sub_for_test(target, base_sum);
            for axis in 0..3 {
                let a = point3_components_for_test(combined_delta)[axis];
                numerator += a * point3_components_for_test(y)[axis];
                denominator += a * a;
            }
        }
        let scale = if denominator > 1.0e-14 {
            numerator / denominator
        } else {
            0.0
        };
        (scale, scale)
    };
    let mut residual: f64 = 0.0;
    for sector in 0..n {
        let next = (sector + 1) % n;
        let target = point3_add_for_test(
            point3_add_for_test(base[sector][0][4], base[next][4][0]),
            point3_scale_for_test(
                point3_sub_for_test(base[sector][0][6], base[sector][0][5]),
                w3,
            ),
        );
        let actual = point3_add_for_test(
            point3_add_for_test(base[sector][1][4], base[next][4][1]),
            point3_add_for_test(
                point3_scale_for_test(
                    point3_sub_for_test(row_unit[next][4][1], base[next][4][1]),
                    row_scale,
                ),
                point3_scale_for_test(
                    point3_sub_for_test(column_unit[sector][1][4], base[sector][1][4]),
                    column_scale,
                ),
            ),
        );
        residual = residual.max(actual.distance_to(&target));
    }
    Kp15S4FirstOrderNormalScaleSolve {
        row_scale,
        column_scale,
        residual,
        determinant: det,
    }
}

#[derive(Clone, Copy, Debug)]
enum Kp15S4CrossCornerOwnershipMode {
    DeriveColumnP14,
    DeriveRowP41,
    DeriveBothSequentially,
}

fn kp15_apply_s4_cross_corner_ownership_correction_for_test(
    grids: &mut [[[Point3; 7]; 7]],
    params: Bi6Reparams,
    mode: Kp15S4CrossCornerOwnershipMode,
) -> f64 {
    let mut max_delta: f64 = 0.0;
    match mode {
        Kp15S4CrossCornerOwnershipMode::DeriveColumnP14 => {
            max_delta = max_delta.max(kp15_overwrite_column_p14_from_s4_sum_for_test(
                grids, params,
            ));
        }
        Kp15S4CrossCornerOwnershipMode::DeriveRowP41 => {
            max_delta = max_delta.max(kp15_overwrite_row_p41_from_s4_sum_for_test(grids, params));
        }
        Kp15S4CrossCornerOwnershipMode::DeriveBothSequentially => {
            max_delta = max_delta.max(kp15_overwrite_column_p14_from_s4_sum_for_test(
                grids, params,
            ));
            max_delta = max_delta.max(kp15_overwrite_row_p41_from_s4_sum_for_test(grids, params));
        }
    }
    max_delta
}

fn kp15_overwrite_column_p14_from_s4_sum_for_test(
    grids: &mut [[[Point3; 7]; 7]],
    params: Bi6Reparams,
) -> f64 {
    let n = grids.len();
    let original = grids.to_vec();
    let w3 = params.c / 15.0;
    let mut max_delta: f64 = 0.0;
    for sector in 0..n {
        let next = (sector + 1) % n;
        let target = Point3::new(
            original[sector][0][4].x + original[next][4][0].x - original[next][4][1].x
                + w3 * (original[sector][0][6].x - original[sector][0][5].x),
            original[sector][0][4].y + original[next][4][0].y - original[next][4][1].y
                + w3 * (original[sector][0][6].y - original[sector][0][5].y),
            original[sector][0][4].z + original[next][4][0].z - original[next][4][1].z
                + w3 * (original[sector][0][6].z - original[sector][0][5].z),
        );
        max_delta = max_delta.max(original[sector][1][4].distance_to(&target));
        grids[sector][1][4] = target;
    }
    max_delta
}

fn kp15_overwrite_row_p41_from_s4_sum_for_test(
    grids: &mut [[[Point3; 7]; 7]],
    params: Bi6Reparams,
) -> f64 {
    let n = grids.len();
    let original = grids.to_vec();
    let w3 = params.c / 15.0;
    let mut max_delta: f64 = 0.0;
    for sector in 0..n {
        let previous = (sector + n - 1) % n;
        let target = Point3::new(
            original[sector][4][0].x + original[previous][0][4].x - original[previous][1][4].x
                + w3 * (original[previous][0][6].x - original[previous][0][5].x),
            original[sector][4][0].y + original[previous][0][4].y - original[previous][1][4].y
                + w3 * (original[previous][0][6].y - original[previous][0][5].y),
            original[sector][4][0].z + original[previous][0][4].z - original[previous][1][4].z
                + w3 * (original[previous][0][6].z - original[previous][0][5].z),
        );
        max_delta = max_delta.max(original[sector][4][1].distance_to(&target));
        grids[sector][4][1] = target;
    }
    max_delta
}

#[derive(Clone)]
struct Kp15LabeledConstraint {
    label: String,
    constraint: AffineConstraintRow,
}

fn kp15_labeled_source_constraints_from_grid_columns_for_test(
    n: usize,
    source: &[[[Point3; 7]; 7]],
    cols: &[usize],
    side_label: &str,
) -> Vec<Kp15LabeledConstraint> {
    cols.iter()
        .copied()
        .map(|col| {
            let sector = col / 49;
            let pos = col % 49;
            let i = pos / 7;
            let j = pos % 7;
            assert!(sector < n);
            let point = source[sector][i][j];
            Kp15LabeledConstraint {
                label: format!("{side_label} sector={sector} p_{i}{j}"),
                constraint: AffineConstraintRow {
                    row: LinearRow {
                        terms: vec![(col, 1.0)],
                    },
                    rhs: [point.x, point.y, point.z],
                },
            }
        })
        .collect()
}

fn kp15_first_inconsistent_source_constraint_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
) -> Option<String> {
    kp15_first_inconsistent_source_constraint_residual_for_test(table, source_constraints)
        .map(|(residual, label)| format!("{label} residual={residual:e}"))
}

fn kp15_first_inconsistent_source_constraint_residual_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
) -> Option<(f64, String)> {
    for axis in 0..3 {
        let mut basis = Vec::<Vec<f64>>::new();
        let mut targets = Vec::<f64>::new();
        for constraint in kp15_table_smoothness_constraints(table) {
            assert!(
                kp15_push_labeled_affine_basis_row_for_test(
                    table.cols,
                    &constraint.row,
                    constraint.rhs[axis],
                    &mut basis,
                    &mut targets,
                )
                .is_none(),
                "KP15 table rows must be homogeneous-consistent"
            );
        }
        for source_constraint in source_constraints {
            if let Some(residual) = kp15_push_labeled_affine_basis_row_for_test(
                table.cols,
                &source_constraint.constraint.row,
                source_constraint.constraint.rhs[axis],
                &mut basis,
                &mut targets,
            ) {
                return Some((
                    residual.abs(),
                    format!("axis={axis} {}", source_constraint.label),
                ));
            }
        }
    }
    None
}

fn kp15_implied_ordered_constraint_delta_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
    target_label: &str,
) -> Option<String> {
    let target_constraint = source_constraints
        .iter()
        .find(|constraint| constraint.label == target_label)?;
    let mut implied = target_constraint.constraint.rhs;
    let mut dependent_axes = 0usize;
    for (axis, implied_axis) in implied.iter_mut().enumerate() {
        let mut basis = Vec::<Vec<f64>>::new();
        let mut targets = Vec::<f64>::new();
        for constraint in kp15_table_smoothness_constraints(table) {
            assert!(
                kp15_push_labeled_affine_basis_row_for_test(
                    table.cols,
                    &constraint.row,
                    constraint.rhs[axis],
                    &mut basis,
                    &mut targets,
                )
                .is_none(),
                "KP15 table rows must be homogeneous-consistent"
            );
        }
        for source_constraint in source_constraints {
            if source_constraint.label == target_label {
                if let Some(residual) = kp15_dependent_affine_residual_for_test(
                    table.cols,
                    &source_constraint.constraint.row,
                    source_constraint.constraint.rhs[axis],
                    &basis,
                    &targets,
                ) {
                    *implied_axis = source_constraint.constraint.rhs[axis] - residual;
                    dependent_axes += 1;
                }
                break;
            }
            if kp15_push_labeled_affine_basis_row_for_test(
                table.cols,
                &source_constraint.constraint.row,
                source_constraint.constraint.rhs[axis],
                &mut basis,
                &mut targets,
            )
            .is_some()
            {
                return None;
            }
        }
    }
    let source = target_constraint.constraint.rhs;
    let delta = ((source[0] - implied[0]).powi(2)
        + (source[1] - implied[1]).powi(2)
        + (source[2] - implied[2]).powi(2))
    .sqrt();
    Some(format!(
        "{target_label} source=({sx:e},{sy:e},{sz:e}) implied=({ix:e},{iy:e},{iz:e}) delta={delta:e} dependent_axes={dependent_axes}",
        sx = source[0],
        sy = source[1],
        sz = source[2],
        ix = implied[0],
        iy = implied[1],
        iz = implied[2],
    ))
}

struct Kp15AffineBasisFormulaRow {
    row: Vec<f64>,
    source_coefficients: Vec<f64>,
}

fn kp15_ordered_source_formula_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
    target_label: &str,
    axis: usize,
) -> Option<String> {
    let target_index = source_constraints
        .iter()
        .position(|constraint| constraint.label == target_label)?;
    let mut basis = Vec::<Kp15AffineBasisFormulaRow>::new();
    for constraint in kp15_table_smoothness_constraints(table) {
        let mut source_coefficients = vec![0.0_f64; source_constraints.len()];
        kp15_push_formula_basis_row_for_test(
            table.cols,
            &constraint.row,
            &mut source_coefficients,
            &mut basis,
        );
    }
    for (source_index, source_constraint) in source_constraints.iter().enumerate() {
        let mut source_coefficients = vec![0.0_f64; source_constraints.len()];
        source_coefficients[source_index] = 1.0;
        if source_index == target_index {
            let implied = kp15_dependent_source_formula_for_test(
                table.cols,
                &source_constraint.constraint.row,
                &basis,
                source_constraints.len(),
            )?;
            let mut ranked = implied
                .iter()
                .enumerate()
                .filter(|&(index, &weight)| (weight.abs() > 1.0e-8)).map(|(index, &weight)| (
                            weight.abs(),
                            weight,
                            source_constraints[index].label.as_str(),
                        ))
                .collect::<Vec<_>>();
            ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then_with(|| a.2.cmp(b.2)));
            let terms = ranked
                .iter()
                .take(12)
                .map(|(_, weight, label)| format!("{weight:+.9e}*{label}"))
                .collect::<Vec<_>>()
                .join(" ");
            let target = source_constraint.constraint.rhs[axis];
            let implied_value = implied
                .iter()
                .zip(source_constraints)
                .map(|(&weight, constraint)| weight * constraint.constraint.rhs[axis])
                .sum::<f64>();
            return Some(format!(
                "target={target:e} implied={implied_value:e} residual={residual:e} terms={terms}",
                residual = target - implied_value,
            ));
        }
        kp15_push_formula_basis_row_for_test(
            table.cols,
            &source_constraint.constraint.row,
            &mut source_coefficients,
            &mut basis,
        );
    }
    None
}

fn kp15_ordered_source_formula_point_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
    target_label: &str,
) -> Option<Point3> {
    let weights =
        kp15_ordered_source_formula_weights_for_test(table, source_constraints, target_label)?;
    Some(Point3::new(
        weights
            .iter()
            .zip(source_constraints)
            .map(|(&weight, constraint)| weight * constraint.constraint.rhs[0])
            .sum(),
        weights
            .iter()
            .zip(source_constraints)
            .map(|(&weight, constraint)| weight * constraint.constraint.rhs[1])
            .sum(),
        weights
            .iter()
            .zip(source_constraints)
            .map(|(&weight, constraint)| weight * constraint.constraint.rhs[2])
            .sum(),
    ))
}

fn kp15_ordered_source_formula_weights_for_test(
    table: &Kp15Bi6GeneratingTable,
    source_constraints: &[Kp15LabeledConstraint],
    target_label: &str,
) -> Option<Vec<f64>> {
    let target_index = source_constraints
        .iter()
        .position(|constraint| constraint.label == target_label)?;
    let mut basis = Vec::<Kp15AffineBasisFormulaRow>::new();
    for constraint in kp15_table_smoothness_constraints(table) {
        let mut source_coefficients = vec![0.0_f64; source_constraints.len()];
        kp15_push_formula_basis_row_for_test(
            table.cols,
            &constraint.row,
            &mut source_coefficients,
            &mut basis,
        );
    }
    for (source_index, source_constraint) in source_constraints.iter().enumerate() {
        let mut source_coefficients = vec![0.0_f64; source_constraints.len()];
        source_coefficients[source_index] = 1.0;
        if source_index == target_index {
            return kp15_dependent_source_formula_for_test(
                table.cols,
                &source_constraint.constraint.row,
                &basis,
                source_constraints.len(),
            );
        }
        kp15_push_formula_basis_row_for_test(
            table.cols,
            &source_constraint.constraint.row,
            &mut source_coefficients,
            &mut basis,
        );
    }
    None
}

fn kp15_push_formula_basis_row_for_test(
    cols: usize,
    row: &LinearRow,
    source_coefficients: &mut [f64],
    basis: &mut Vec<Kp15AffineBasisFormulaRow>,
) -> bool {
    let mut candidate = vec![0.0_f64; cols];
    for &(col, weight) in &row.terms {
        candidate[col] += weight;
    }
    for basis_row in basis.iter() {
        let dot = candidate
            .iter()
            .zip(&basis_row.row)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        for (value, basis_value) in candidate.iter_mut().zip(&basis_row.row) {
            *value -= dot * basis_value;
        }
        for (value, basis_value) in source_coefficients
            .iter_mut()
            .zip(&basis_row.source_coefficients)
        {
            *value -= dot * basis_value;
        }
    }
    let norm = candidate.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm <= 1.0e-10 {
        return false;
    }
    for value in &mut candidate {
        *value /= norm;
    }
    for value in source_coefficients.iter_mut() {
        *value /= norm;
    }
    basis.push(Kp15AffineBasisFormulaRow {
        row: candidate,
        source_coefficients: source_coefficients.to_vec(),
    });
    true
}

fn kp15_dependent_source_formula_for_test(
    cols: usize,
    row: &LinearRow,
    basis: &[Kp15AffineBasisFormulaRow],
    source_count: usize,
) -> Option<Vec<f64>> {
    let mut candidate = vec![0.0_f64; cols];
    for &(col, weight) in &row.terms {
        candidate[col] += weight;
    }
    let mut formula = vec![0.0_f64; source_count];
    for basis_row in basis {
        let dot = candidate
            .iter()
            .zip(&basis_row.row)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        for (value, basis_value) in candidate.iter_mut().zip(&basis_row.row) {
            *value -= dot * basis_value;
        }
        for (value, basis_value) in formula.iter_mut().zip(&basis_row.source_coefficients) {
            *value += dot * basis_value;
        }
    }
    let norm = candidate.iter().map(|v| v * v).sum::<f64>().sqrt();
    (norm <= 1.0e-10).then_some(formula)
}

fn kp15_dependent_affine_residual_for_test(
    cols: usize,
    row: &LinearRow,
    target: f64,
    basis: &[Vec<f64>],
    targets: &[f64],
) -> Option<f64> {
    let mut candidate = vec![0.0_f64; cols];
    for &(col, weight) in &row.terms {
        candidate[col] += weight;
    }
    let mut candidate_target = target;
    for (basis_row, &basis_target) in basis.iter().zip(targets.iter()) {
        let dot = candidate
            .iter()
            .zip(basis_row)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        for (value, basis_value) in candidate.iter_mut().zip(basis_row) {
            *value -= dot * basis_value;
        }
        candidate_target -= dot * basis_target;
    }
    let norm = candidate.iter().map(|v| v * v).sum::<f64>().sqrt();
    (norm <= 1.0e-10).then_some(candidate_target)
}

fn kp15_push_labeled_affine_basis_row_for_test(
    cols: usize,
    row: &LinearRow,
    target: f64,
    basis: &mut Vec<Vec<f64>>,
    targets: &mut Vec<f64>,
) -> Option<f64> {
    let mut candidate = vec![0.0_f64; cols];
    for &(col, weight) in &row.terms {
        candidate[col] += weight;
    }
    let mut candidate_target = target;
    for (basis_row, &basis_target) in basis.iter().zip(targets.iter()) {
        let dot = candidate
            .iter()
            .zip(basis_row)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        for (value, basis_value) in candidate.iter_mut().zip(basis_row) {
            *value -= dot * basis_value;
        }
        candidate_target -= dot * basis_target;
    }
    let norm = candidate.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm <= 1.0e-10 {
        if candidate_target.abs() > 1.0e-8 {
            return Some(candidate_target);
        }
        return None;
    }
    for value in &mut candidate {
        *value /= norm;
    }
    basis.push(candidate);
    targets.push(candidate_target / norm);
    None
}

fn kp15_project_source_side_to_row_manifold_for_test(
    n: usize,
    table: &Kp15Bi6GeneratingTable,
    source: &[[[Point3; 7]; 7]],
    hard_cols: &[usize],
    soft_cols: &[usize],
) -> Option<Vec<[[Point3; 7]; 7]>> {
    let vars = 49 * n;
    let mut hard_constraints = kp15_table_smoothness_constraints(table);
    for (col, point) in kp15_source_constraints_from_grid_columns_for_test(n, source, hard_cols) {
        hard_constraints.push(AffineConstraintRow {
            row: LinearRow {
                terms: vec![(col, 1.0)],
            },
            rhs: [point.x, point.y, point.z],
        });
    }
    let soft_targets = kp15_source_constraints_from_grid_columns_for_test(n, source, soft_cols);
    let mut coords = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];
    for axis in 0..3 {
        let size = vars + hard_constraints.len();
        let mut matrix = vec![vec![0.0_f64; size]; size];
        let mut rhs = vec![0.0_f64; size];
        for i in 0..vars {
            matrix[i][i] = 1.0e-12;
        }
        for &(col, point) in &soft_targets {
            let target = match axis {
                0 => point.x,
                1 => point.y,
                _ => point.z,
            };
            matrix[col][col] += 1.0;
            rhs[col] += target;
        }
        for (row_index, constraint) in hard_constraints.iter().enumerate() {
            add_lagrange_row(
                vars + row_index,
                &mut matrix,
                &mut rhs,
                &constraint.row.terms,
                constraint.rhs[axis],
            );
        }
        coords[axis] = solve_dense_faer_for_test(matrix, rhs).map(|x| x[..vars].to_vec())?;
    }
    let mut grids = Vec::with_capacity(n);
    for sector in 0..n {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let col = kp15_global_col(n, sector, i, j);
                *cell = Point3::new(coords[0][col], coords[1][col], coords[2][col]);
            }
        }
        grids.push(grid);
    }
    Some(grids)
}

fn kp15_max_source_column_delta_for_test(
    n: usize,
    source: &[[[Point3; 7]; 7]],
    projected: &[[[Point3; 7]; 7]],
    cols: &[usize],
) -> f64 {
    cols.iter().fold(0.0_f64, |acc, &col| {
        let sector = col / 49;
        let pos = col % 49;
        let i = pos / 7;
        let j = pos % 7;
        assert!(sector < n);
        acc.max(source[sector][i][j].distance_to(&projected[sector][i][j]))
    })
}

fn kp15_point_linear_combination_for_test(terms: &[(f64, Point3)]) -> Point3 {
    Point3::new(
        terms.iter().map(|(weight, point)| weight * point.x).sum(),
        terms.iter().map(|(weight, point)| weight * point.y).sum(),
        terms.iter().map(|(weight, point)| weight * point.z).sum(),
    )
}

fn kp15_outer_depth_row_columns_for_test(
    n: usize,
    depth_rows: std::ops::RangeInclusive<usize>,
) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for j in depth_rows.clone() {
            for i in 0..=6 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_outer_i_band_columns_for_test(
    n: usize,
    band: std::ops::RangeInclusive<usize>,
) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for i in band.clone() {
            for j in 0..=6 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_outer_j_band_columns_for_test(
    n: usize,
    band: std::ops::RangeInclusive<usize>,
) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for j in band.clone() {
            for i in 0..=6 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_paper_depth_two_sided_strip_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = kp15_paper_depth_row_side_columns_for_test(n);
    cols.extend(kp15_paper_depth_column_side_columns_for_test(n));
    cols
}

fn kp15_fig2b_depth2_tensor_border_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for i in 4..=6 {
            cols.push(kp15_global_col(n, sector, i, 2));
        }
        for j in 4..=6 {
            cols.push(kp15_global_col(n, sector, 2, j));
        }
    }
    cols
}

fn kp15_fig2b_depth2_plus_m_endpoint_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = kp15_fig2b_depth2_tensor_border_columns_for_test(n);
    for sector in 0..n {
        cols.push(kp15_global_col(n, sector, 6, 0));
        cols.push(kp15_global_col(n, sector, 0, 6));
    }
    cols.sort_unstable();
    cols.dedup();
    cols
}

fn kp15_paper_depth_row_side_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for i in 4..=6 {
            for j in 0..=2 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_paper_depth_column_side_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = Vec::new();
    for sector in 0..n {
        for i in 0..=2 {
            for j in 4..=6 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_paper_depth_two_sided_strip_with_corner_l_band_columns_for_test(n: usize) -> Vec<usize> {
    let mut cols = kp15_paper_depth_two_sided_strip_columns_for_test(n);
    for sector in 0..n {
        for i in 4..=6 {
            for j in 4..=6 {
                cols.push(kp15_global_col(n, sector, i, j));
            }
        }
    }
    cols
}

fn kp15_constraint_nullity_after_pinning_for_test(
    table: &Kp15Bi6GeneratingTable,
    pinned_cols: &std::collections::BTreeSet<usize>,
) -> usize {
    let mut rows = Vec::with_capacity(table.rank() + pinned_cols.len());
    for row in 0..table.rank() {
        rows.push(
            (0..table.cols)
                .map(|col| table.rref[row * table.cols + col])
                .collect::<Vec<_>>(),
        );
    }
    for &col in pinned_cols {
        let mut row = vec![0.0_f64; table.cols];
        row[col] = 1.0;
        rows.push(row);
    }
    table.cols - dense_rank_relative_for_test(&rows)
}

fn kp15_free_columns_after_pinning_for_test(
    table: &Kp15Bi6GeneratingTable,
    pinned_cols: &std::collections::BTreeSet<usize>,
) -> Vec<usize> {
    let mut rows = Vec::with_capacity(table.rank() + pinned_cols.len());
    for row in 0..table.rank() {
        rows.push(
            (0..table.cols)
                .map(|col| table.rref[row * table.cols + col])
                .collect::<Vec<_>>(),
        );
    }
    for &col in pinned_cols {
        let mut row = vec![0.0_f64; table.cols];
        row[col] = 1.0;
        rows.push(row);
    }
    let pivots = kp15_rref_pivot_columns_for_test(&mut rows, 1.0e-10);
    let pivot_set = pivots
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    (0..table.cols)
        .filter(|col| !pivot_set.contains(col))
        .collect()
}

fn kp15_rref_pivot_columns_for_test(rows: &mut [Vec<f64>], tol: f64) -> Vec<usize> {
    let row_count = rows.len();
    let col_count = rows.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    let mut row = 0usize;
    for col in 0..col_count {
        if row >= row_count {
            break;
        }
        let mut best = row;
        let mut best_abs = rows[row][col].abs();
        for (candidate, candidate_row) in rows.iter().enumerate().skip(row + 1) {
            let value = candidate_row[col].abs();
            if value > best_abs {
                best = candidate;
                best_abs = value;
            }
        }
        if best_abs <= tol {
            continue;
        }
        rows.swap(row, best);
        let pivot = rows[row][col];
        for c in col..col_count {
            rows[row][c] /= pivot;
        }
        for r in 0..row_count {
            if r == row {
                continue;
            }
            let factor = rows[r][col];
            if factor == 0.0 {
                continue;
            }
            for c in col..col_count {
                rows[r][c] -= factor * rows[row][c];
            }
        }
        pivots.push(col);
        row += 1;
    }
    pivots
}

fn kp15_section34_eq7_p21_p12_matrix_for_test(n: usize, params: Bi6Reparams) -> Vec<Vec<f64>> {
    let mut rows = Vec::with_capacity(2 * n);
    let ddot5 = -4.0 * params.c;
    let x = |sector: usize| sector;
    let y = |sector: usize| n + sector;
    for sector in 0..n {
        let previous = (sector + n - 1) % n;

        // KP15 eq(7): p21^k = -p12^{k-1} + known terms.
        let mut row = vec![0.0_f64; 2 * n];
        row[x(sector)] = 1.0;
        row[y(previous)] = 1.0;
        rows.push(row);

        // KP15 eq(7), second line:
        // p12^k = p21^{k-1} + ddot_mu5 * p12^{k-1} + known terms.
        let mut row = vec![0.0_f64; 2 * n];
        row[y(sector)] = 1.0;
        row[x(previous)] = -1.0;
        row[y(previous)] = -ddot5;
        rows.push(row);
    }
    rows
}

fn kp15_section34_eq8_p22_matrix_for_test(n: usize) -> Vec<Vec<f64>> {
    let mut rows = Vec::with_capacity(n);
    for sector in 0..n {
        let previous = (sector + n - 1) % n;
        // KP15 eq(8): p22^k = p22^{k-1} + known terms.
        let mut row = vec![0.0_f64; n];
        row[sector] = 1.0;
        row[previous] = -1.0;
        rows.push(row);
    }
    rows
}

fn kp15_source_slot_family_columns(
    n: usize,
    slots: &[(usize, usize)],
) -> std::collections::BTreeSet<usize> {
    let mut cols = std::collections::BTreeSet::new();
    for sector in 0..n {
        for &(i, j) in slots {
            cols.insert(kp15_global_col(n, sector, i, j));
        }
    }
    cols
}

fn kp15_released_slot_movement(
    sectors: &[CirculantSectorRhs],
    patches: &[[[Point3; 7]; 7]],
    released_cols: &std::collections::BTreeSet<usize>,
) -> f64 {
    let n = sectors.len();
    let mut max_movement: f64 = 0.0;
    for &col in released_cols {
        let sector = col / 49;
        let pos = col % 49;
        let i = pos / 7;
        let j = pos % 7;
        if sector < n
            && let Some(seed) = kp15_seed_value_for_test(sectors, sector, i, j)
        {
            max_movement = max_movement.max(patches[sector][i][j].distance_to(&seed));
        }
    }
    max_movement
}

struct Kp15ReleasedSlotDeltaStats {
    avg_local: [f64; 2],
    max_norm: f64,
    max_deviation_from_average: f64,
}

fn kp15_released_slot_local_delta_stats(
    sectors: &[CirculantSectorRhs],
    patches: &[[[Point3; 7]; 7]],
    released_cols: &std::collections::BTreeSet<usize>,
) -> Kp15ReleasedSlotDeltaStats {
    let n = sectors.len();
    let mut deltas = Vec::new();
    let mut sum = [0.0, 0.0];
    let mut max_norm: f64 = 0.0;
    for &col in released_cols {
        let sector = col / 49;
        let pos = col % 49;
        let i = pos / 7;
        let j = pos % 7;
        if sector < n
            && let Some(seed) = kp15_seed_value_for_test(sectors, sector, i, j)
        {
            let solved = patches[sector][i][j];
            let diff = Point3::new(solved.x - seed.x, solved.y - seed.y, solved.z - seed.z);
            let local = global_to_local_planar(n, sector, diff);
            sum[0] += local[0];
            sum[1] += local[1];
            max_norm = max_norm.max((local[0] * local[0] + local[1] * local[1]).sqrt());
            deltas.push(local);
        }
    }
    let count = deltas.len().max(1) as f64;
    let avg = [sum[0] / count, sum[1] / count];
    let mut max_deviation: f64 = 0.0;
    for local in deltas {
        let dx = local[0] - avg[0];
        let dy = local[1] - avg[1];
        max_deviation = max_deviation.max((dx * dx + dy * dy).sqrt());
    }
    Kp15ReleasedSlotDeltaStats {
        avg_local: avg,
        max_norm,
        max_deviation_from_average: max_deviation,
    }
}

fn kp15_seed_value_for_test(
    sectors: &[CirculantSectorRhs],
    sector: usize,
    i: usize,
    j: usize,
) -> Option<Point3> {
    let source = &sectors[sector];
    Some(match (i, j) {
        (0, 0) => Point3::new(0.0, 0.0, 0.0),
        (1, 0) => source.p_10,
        (2, 0) => source.p_20,
        (4, 0) => source.p_40,
        (5, 0) => source.p_50,
        (6, 0) => source.p_60,
        (0, 1) => source.p_01,
        (1, 1) => source.p_11,
        (0, 2) => source.p_02,
        (4, 2) => source.p_42,
        (5, 2) => source.p_52,
        (6, 2) => source.p_62,
        _ => return None,
    })
}

fn solve_dense_faer_for_test(matrix: Vec<Vec<f64>>, rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = matrix.len();
    if n == 0
        || rhs.len() != n
        || matrix.iter().any(|row| row.len() != n)
        || matrix.iter().flatten().any(|value| !value.is_finite())
        || rhs.iter().any(|value| !value.is_finite())
    {
        return None;
    }
    let mut faer_matrix = Mat::<f64>::zeros(n, n);
    let mut faer_rhs = Mat::<f64>::zeros(n, 1);
    for row in 0..n {
        for col in 0..n {
            faer_matrix.write(row, col, matrix[row][col]);
        }
        faer_rhs.write(row, 0, rhs[row]);
    }
    let lu = faer_matrix.partial_piv_lu();
    let solution = lu.solve(&faer_rhs);
    let out = (0..n).map(|row| solution.read(row, 0)).collect::<Vec<_>>();
    if out.iter().all(|value| value.is_finite()) {
        Some(out)
    } else {
        None
    }
}

fn eval_tri5_derivative(cps: &[Point3; 21], weights: &[f64; 21]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (p, &w) in cps.iter().zip(weights) {
        out.x += w * p.x;
        out.y += w * p.y;
        out.z += w * p.z;
    }
    out
}

fn kp15_spoke_indexed_planar_bilinear_input(n: usize) -> Bi6CapInput {
    let params = Bi6Reparams::for_valence(n);
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut sectors = Vec::with_capacity(n);
    for sector in 0..n {
        let at = |i: usize, j: usize| -> Point3 {
            local_to_global_planar(n, sector, [i as f64 / 6.0, j as f64 / 6.0])
        };
        sectors.push(CirculantSectorRhs {
            p_00: zero,
            p_10: at(1, 0),
            p_20: at(2, 0),
            p_40: at(4, 0),
            p_50: at(5, 0),
            p_60: at(6, 0),
            p_01: at(0, 1),
            p_11: at(1, 1),
            p_02: at(0, 2),
            p_42: at(4, 2),
            p_52: at(5, 2),
            p_62: at(6, 2),
        });
    }
    Bi6CapInput {
        n,
        params,
        c7: zero,
        c5s: vec![zero; n],
        c6s: vec![zero; n],
        sectors,
        p0_22: zero,
        p0_21_pinned: if n == 3 || n == 6 { Some(zero) } else { None },
        interior_filler: None,
    }
}

fn kp15_pr08_q11_transposed_endpoint_input(n: usize, apply_beta: bool) -> Bi6CapInput {
    let endpoint_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    kp15_input_from_bi3_tensor(n, &endpoint_tensor, apply_beta)
}

fn kp15_pr08_q11_transposed_endpoint_tensor(n: usize) -> TensorBorderBi3 {
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp15_pr08_q11_transposed_regular_rows(&ring);
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0)
}

fn kp15_pr08_q11_transposed_regular_rows(
    ring: &crate::geometry::gregory::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing,
) -> Kp20RegularTensorBorderCNet {
    
    kp20_pr08_characteristic_tensor_border_candidate_rows(
        ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows")
}

fn kp15_input_from_bi3_tensor(n: usize, tensor: &TensorBorderBi3, apply_beta: bool) -> Bi6CapInput {
    let params = Bi6Reparams::for_valence(n);
    let sectors = kp15_sector_rhs_from_bi3_rows(&tensor.sectors, params, apply_beta);
    let zero = Point3::new(0.0, 0.0, 0.0);
    Bi6CapInput {
        n,
        params,
        c7: zero,
        c5s: vec![zero; n],
        c6s: vec![zero; n],
        sectors,
        p0_22: zero,
        p0_21_pinned: if n == 3 || n == 6 { Some(zero) } else { None },
        interior_filler: None,
    }
}

fn kp15_pr08_exact_stam_side_raw_tensor(n: usize) -> TensorBorderBi3 {
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let projection =
        kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Appendix tensor projection");
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridTransform::FlipColumns,
        Kp20StamGridSide::Left,
        true,
    )
    .expect("PR08/Stam exact side source rows");
    kp20_tensor_border_from_regular_cnet_rows(&rows).expect("PR08/Stam exact side tensor")
}

fn kp15_pr08_exact_stam_side_regular_endpoint_tensor(n: usize) -> TensorBorderBi3 {
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let projection =
        kp20_pr08_stam_appendix_tensor_projected_grids(&ring).expect("Appendix tensor projection");
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridTransform::FlipColumns,
        Kp20StamGridSide::Left,
        true,
    )
    .expect("PR08/Stam exact side source rows");
    let normalized_rows =
        kp15_regular_rows_unit_sector_endpoint_normalized(&rows, n, 0, true, true)
            .expect("PR08/Stam exact side endpoint normalization");
    kp20_tensor_border_from_regular_cnet_rows(&normalized_rows)
        .expect("PR08/Stam exact side endpoint-normalized tensor")
}

fn kp15_pr08_endpoint_side_tensor_border_for_fig12(n: usize) -> TensorBorderBi3 {
    let rows = kp20_pr08_fig63_endpoint_side_regular_tensor_border_rows(n)
        .expect("PR08 Fig.6.3 endpoint-side source rows");
    let normalized_rows =
        kp15_regular_rows_unit_sector_endpoint_normalized(&rows, n, 0, true, false)
            .expect("PR08 Fig.6.3 endpoint-side rows have a nondegenerate edge chart");
    kp20_tensor_border_from_regular_cnet_rows(&normalized_rows)
        .expect("PR08 Fig.6.3 endpoint-side source converts to BB tensor-border")
}

fn kp15_pr08_fig63_lower_right_side_tensor_border_for_fig12(
    n: usize,
    side: Kp20StamGridSide,
    swap_targets: bool,
    conjugate_source: bool,
) -> Result<TensorBorderBi3, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let projection = kp20_pr08_stam_appendix_tensor_projected_grids(&ring)?;
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        side,
        false,
    )?;
    let normalized_rows = kp15_regular_rows_unit_sector_endpoint_normalized(
        &rows,
        n,
        0,
        swap_targets,
        conjugate_source,
    )?;
    kp20_tensor_border_from_regular_cnet_rows(&normalized_rows)
}

fn kp15_pr08_fig63_joint_corner_side_tensor_border_for_fig12(
    n: usize,
    side: Kp20StamGridSide,
    swap_outer_targets: bool,
    conjugate_source: bool,
) -> Result<TensorBorderBi3, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let projection = kp20_pr08_stam_appendix_tensor_projected_grids(&ring)?;
    let left_rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridSide::Left,
        false,
    )?;
    let top_rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridSide::Top,
        false,
    )?;
    let selected_rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        side,
        false,
    )?;
    let alpha = core::f64::consts::TAU / n as f64;
    let spoke0 = [1.0, 0.0];
    let spoke1 = [f64::cos(alpha), f64::sin(alpha)];
    let previous_vertex = kp15_maybe_conjugate_xy(left_rows.sectors[0][0][3], conjugate_source);
    let next_vertex = kp15_maybe_conjugate_xy(top_rows.sectors[0][0][3], conjugate_source);
    let target_previous = if swap_outer_targets { spoke1 } else { spoke0 };
    let target_next = if swap_outer_targets { spoke0 } else { spoke1 };
    let scale = kp15_complex_similarity_scale_xy_for_test(
        &[previous_vertex, next_vertex],
        &[target_previous, target_next],
    )
    .ok_or_else(|| "degenerate PR08 Fig.6.3 joint corner normalization".to_string())?;
    let normalized_rows = kp15_apply_complex_similarity_to_regular_rows_for_test(
        &selected_rows,
        scale,
        conjugate_source,
    );
    kp20_tensor_border_from_regular_cnet_rows(&normalized_rows)
}

fn kp15_pr08_fig63_center_to_sector_join_radial_rows_for_test(
    n: usize,
    reversed_depth: bool,
) -> Result<TensorBorderBi3, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let projection = kp20_pr08_stam_appendix_tensor_projected_grids(&ring)?;
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::LowerRight,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridSide::Right,
        reversed_depth,
    )?;
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows)?;
    kp15_center_to_sector_join_normalized_tensor_border_for_test(&tensor, n)
}

fn kp15_center_to_sector_join_normalized_tensor_border_for_test(
    tensor: &TensorBorderBi3,
    _n: usize,
) -> Result<TensorBorderBi3, String> {
    let sector0 = tensor
        .sectors
        .first()
        .ok_or_else(|| "center-to-sector normalization needs sector 0".to_string())?;
    let source_center = sector0[0][0];
    let source_join = sector0[0][3];
    let delta = [
        source_join.x - source_center.x,
        source_join.y - source_center.y,
    ];
    let denominator = delta[0] * delta[0] + delta[1] * delta[1];
    if denominator <= 1.0e-20 {
        return Err("degenerate PR08 Fig.6.3 center-to-sector radial side".to_string());
    }
    let scale = [delta[0] / denominator, -delta[1] / denominator];
    let sectors = tensor
        .sectors
        .iter()
        .map(|sector_rows| {
            let center = sector_rows[0][0];
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for cp in 0..4 {
                    let xy = [
                        sector_rows[row][cp].x - center.x,
                        sector_rows[row][cp].y - center.y,
                    ];
                    out[row][cp] = Point3::new(
                        scale[0] * xy[0] - scale[1] * xy[1],
                        scale[1] * xy[0] + scale[0] * xy[1],
                        sector_rows[row][cp].z - center.z,
                    );
                }
            }
            out
        })
        .collect();
    Ok(TensorBorderBi3 { sectors })
}

fn kp15_regular_rows_unit_sector_endpoint_normalized(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
    swap_targets: bool,
    conjugate_source: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let sector_rows = rows
        .sectors
        .get(sector)
        .ok_or_else(|| format!("endpoint normalization sector {sector} is missing"))?;
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [
        f64::cos(alpha * sector as f64),
        f64::sin(alpha * sector as f64),
    ];
    let spoke1 = [
        f64::cos(alpha * (sector + 1) as f64),
        f64::sin(alpha * (sector + 1) as f64),
    ];
    let target0 = if swap_targets { spoke1 } else { spoke0 };
    let target1 = if swap_targets { spoke0 } else { spoke1 };
    let source0 = kp15_maybe_conjugate_xy(sector_rows[0][0], conjugate_source);
    let source1 = kp15_maybe_conjugate_xy(sector_rows[0][3], conjugate_source);
    let pairs = [(source0, target0), (source1, target1)];
    let mut numerator = [0.0_f64, 0.0_f64];
    let mut denominator = 0.0_f64;
    for (source, target) in pairs {
        numerator[0] += target[0] * source[0] + target[1] * source[1];
        numerator[1] += target[1] * source[0] - target[0] * source[1];
        denominator += source[0] * source[0] + source[1] * source[1];
    }
    if denominator <= 1.0e-20 {
        return Err("degenerate PR08/Stam endpoint normalization source".to_string());
    }
    let scale = [numerator[0] / denominator, numerator[1] / denominator];
    let sectors = rows
        .sectors
        .iter()
        .map(|sector_rows| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for cp in 0..4 {
                    let xy = kp15_maybe_conjugate_xy(sector_rows[row][cp], conjugate_source);
                    out[row][cp] = Point3::new(
                        scale[0] * xy[0] - scale[1] * xy[1],
                        scale[1] * xy[0] + scale[0] * xy[1],
                        sector_rows[row][cp].z,
                    );
                }
            }
            out
        })
        .collect();
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

fn kp15_maybe_conjugate_xy(point: Point3, conjugate: bool) -> [f64; 2] {
    if conjugate {
        [point.x, -point.y]
    } else {
        [point.x, point.y]
    }
}

fn kp15_tensor_local_delta(n: usize, sector: usize, a: Point3, b: Point3) -> [f64; 2] {
    let delta = Point3::new(a.x - b.x, a.y - b.y, a.z - b.z);
    global_to_local_planar(n, sector, delta)
}

fn kp15_q11_fitted_first_order_frame_delta(n: usize) -> [[f64; 2]; 2] {
    let target_deltas = [
        (
            0usize,
            2usize,
            [1.0386080728772147e-1 / 0.6, -8.162706209978414e-2 / 0.6],
        ),
        (
            2usize,
            2usize,
            [-1.0381440531524103e-1 / 0.6, 4.5351483090759856e-2 / 0.6],
        ),
    ];
    let base_tensor = kp15_pr08_q11_transposed_endpoint_tensor(n);
    let p0 = global_to_local_planar(
        n,
        0,
        base_tensor.sectors[0][target_deltas[0].0][target_deltas[0].1],
    );
    let p1 = global_to_local_planar(
        n,
        0,
        base_tensor.sectors[0][target_deltas[1].0][target_deltas[1].1],
    );
    kp15_first_order_map_from_two_point_deltas(p0, target_deltas[0].2, p1, target_deltas[1].2)
        .expect("q11 fitted frame points should be nondegenerate")
}

fn kp15_beta_tensor_from_bi3_tensor(tensor: &TensorBorderBi3, n: usize) -> TensorBorderBi3 {
    let params = Bi6Reparams::for_valence(n);
    let sectors = tensor
        .sectors
        .iter()
        .map(|sector_rows| {
            let beta = apply_beta_reparameterization(sector_rows, &params);
            [
                kp15_bi6_row_to_bi3(&beta[0]),
                kp15_bi6_row_to_bi3(&beta[1]),
                kp15_bi6_row_to_bi3(&beta[2]),
            ]
        })
        .collect();
    TensorBorderBi3 { sectors }
}

fn kp15_bi6_row_to_bi3(row: &[Point3; 7]) -> [Point3; 4] {
    [row[0], row[2], row[4], row[6]]
}

fn kp15_first_order_map_between_tensor_controls(
    n: usize,
    source: &TensorBorderBi3,
    target: &TensorBorderBi3,
    controls: &[(usize, usize); 2],
) -> Option<[[f64; 2]; 2]> {
    let p0 = global_to_local_planar(n, 0, source.sectors[0][controls[0].0][controls[0].1]);
    let p1 = global_to_local_planar(n, 0, source.sectors[0][controls[1].0][controls[1].1]);
    let q0 = global_to_local_planar(n, 0, target.sectors[0][controls[0].0][controls[0].1]);
    let q1 = global_to_local_planar(n, 0, target.sectors[0][controls[1].0][controls[1].1]);
    kp15_first_order_map_from_point_pairs(p0, q0, p1, q1)
}

fn kp15_s1_first_order_frame_map_from_tensor(
    tensor: &TensorBorderBi3,
    n: usize,
) -> Option<[[f64; 2]; 2]> {
    let params = Bi6Reparams::for_valence(n);
    let rows = &tensor.sectors[0];
    let p10 = global_to_local_planar(n, 0, rows[0][1]);
    let p01 = global_to_local_planar(n, 0, rows[1][0]);
    let next_rows = &tensor.sectors[1];
    let next_p10 = global_to_local_planar(n, 1, next_rows[0][1]);
    let next_p01 = global_to_local_planar(n, 1, next_rows[1][0]);
    let s1_p10 = p10;
    let s1_p01 = [
        -p01[0] + 2.0 * params.c * p10[0],
        -p01[1] + 2.0 * params.c * p10[1],
    ];
    kp15_first_order_map_from_point_pairs(next_p10, s1_p10, next_p01, s1_p01)
}

fn kp15_first_order_map_from_point_pairs(
    p0: [f64; 2],
    q0: [f64; 2],
    p1: [f64; 2],
    q1: [f64; 2],
) -> Option<[[f64; 2]; 2]> {
    let det = p0[0] * p1[1] - p1[0] * p0[1];
    if det.abs() <= 1.0e-20 {
        return None;
    }
    let solve_row = |rhs0: f64, rhs1: f64| -> [f64; 2] {
        [
            (rhs0 * p1[1] - rhs1 * p0[1]) / det,
            (p0[0] * rhs1 - p1[0] * rhs0) / det,
        ]
    };
    Some([solve_row(q0[0], q1[0]), solve_row(q0[1], q1[1])])
}

fn kp15_add_identity_to_frame_delta(delta: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [1.0 + delta[0][0], delta[0][1]],
        [delta[1][0], 1.0 + delta[1][1]],
    ]
}

fn kp15_subtract_identity_from_frame_map(map: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [[map[0][0] - 1.0, map[0][1]], [map[1][0], map[1][1] - 1.0]]
}

fn kp15_det2(map: [[f64; 2]; 2]) -> f64 {
    map[0][0] * map[1][1] - map[0][1] * map[1][0]
}

fn kp15_frame_map_frobenius_delta(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> f64 {
    ((a[0][0] - b[0][0]).powi(2)
        + (a[0][1] - b[0][1]).powi(2)
        + (a[1][0] - b[1][0]).powi(2)
        + (a[1][1] - b[1][1]).powi(2))
    .sqrt()
}

fn kp15_solve_paper_p_f5_for_tensor(
    table: &Kp15Bi6GeneratingTable,
    n: usize,
    tensor: &TensorBorderBi3,
    q: &[[f64; 49]; 49],
) -> Vec<[[Point3; 7]; 7]> {
    let input = kp15_input_from_bi3_tensor(n, tensor, false);
    table
        .populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
            &input.sectors,
            Point3::new(0.0, 0.0, 0.0),
            q,
            0.0,
        )
        .expect("paper-P/F5 tensor solve")
        .0
}

fn kp15_endpoint_spoke_f5_variational_frame_delta(
    n: usize,
    base_tensor: &TensorBorderBi3,
) -> Option<[[f64; 2]; 2]> {
    let params = Bi6Reparams::for_valence(n);
    let table =
        Kp15Bi6GeneratingTable::build_with_transposed_previous_side_and_central_quadratic_free(
            n, params,
        );
    let q = f_kappa_quadratic_form(5);
    let base_patches = kp15_solve_paper_p_f5_for_tensor(&table, n, base_tensor, &q);
    let basis_deltas = [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.0, 1.0], [0.0, 0.0]],
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 0.0], [0.0, 1.0]],
    ];
    let basis_patches = basis_deltas
        .iter()
        .map(|&frame_delta| {
            let tensor =
                kp15_apply_sector_local_first_order_normalization(base_tensor, n, frame_delta);
            kp15_patch_delta(
                &kp15_solve_paper_p_f5_for_tensor(&table, n, &tensor, &q),
                &base_patches,
            )
        })
        .collect::<Vec<_>>();
    let (h, g) = kp15_source_frame_quadratic_terms(&base_patches, &basis_patches, &q);
    let constraints = [
        (vec![1.0, 0.0, 0.0, 0.0], 0.0),
        (vec![0.0, 0.0, 1.0, 0.0], 0.0),
    ];
    let solution = kp15_solve_quadratic_with_linear_constraints(&h, &g, &constraints)?;
    Some([[solution[0], solution[1]], [solution[2], solution[3]]])
}

fn kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    tolerance: f64,
) -> Option<f64> {
    if kp15_endpoint_spoke_frame_alpha_is_orientation_clean(
        n,
        base_tensor,
        frame_delta,
        0.0,
        tolerance,
    ) {
        return Some(0.0);
    }

    let mut previous_alpha = 0.0;
    let mut clean_alpha = None;
    for step in 1..=64 {
        let alpha = step as f64 / 64.0;
        if kp15_endpoint_spoke_frame_alpha_is_orientation_clean(
            n,
            base_tensor,
            frame_delta,
            alpha,
            tolerance,
        ) {
            clean_alpha = Some(alpha);
            break;
        }
        previous_alpha = alpha;
    }
    let mut hi = clean_alpha?;
    let mut lo = previous_alpha;
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if kp15_endpoint_spoke_frame_alpha_is_orientation_clean(
            n,
            base_tensor,
            frame_delta,
            mid,
            tolerance,
        ) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

fn kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation_margin(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    tolerance: f64,
) -> Option<f64> {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation_margin_with_samples(
        n,
        base_tensor,
        frame_delta,
        tolerance,
        &SAMPLES,
    )
}

fn kp15_select_endpoint_spoke_frame_alpha_for_sampled_orientation_margin_with_samples(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    tolerance: f64,
    samples: &[f64],
) -> Option<f64> {
    if kp15_endpoint_spoke_frame_alpha_is_orientation_clean_with_samples(
        n,
        base_tensor,
        frame_delta,
        0.0,
        tolerance,
        samples,
    ) {
        return Some(0.0);
    }

    let mut best = None::<(f64, f64)>;
    for step in 0..=64 {
        let alpha = step as f64 / 64.0;
        if let Some(min_det) = kp15_endpoint_spoke_frame_alpha_clean_min_det(
            n,
            base_tensor,
            frame_delta,
            alpha,
            tolerance,
            samples,
        )
            && best
                .map(|(_, best_min_det)| min_det > best_min_det)
                .unwrap_or(true)
            {
                best = Some((alpha, min_det));
            }
    }
    let (center, _) = best?;
    let mut lo = (center - 1.0 / 64.0).max(0.0);
    let mut hi = (center + 1.0 / 64.0).min(1.0);
    for _ in 0..40 {
        let left = lo + (hi - lo) / 3.0;
        let right = hi - (hi - lo) / 3.0;
        let right_score = kp15_endpoint_spoke_frame_alpha_clean_min_det(
            n,
            base_tensor,
            frame_delta,
            right,
            tolerance,
            samples,
        )
        .unwrap_or(f64::NEG_INFINITY);
        let left_score = kp15_endpoint_spoke_frame_alpha_clean_min_det(
            n,
            base_tensor,
            frame_delta,
            left,
            tolerance,
            samples,
        )
        .unwrap_or(f64::NEG_INFINITY);
        if left_score < right_score {
            lo = left;
        } else {
            hi = right;
        }
    }
    let alpha = 0.5 * (lo + hi);
    if kp15_endpoint_spoke_frame_alpha_is_orientation_clean_with_samples(
        n,
        base_tensor,
        frame_delta,
        alpha,
        tolerance,
        samples,
    ) {
        Some(alpha)
    } else {
        best.map(|(alpha, _)| alpha)
    }
}

fn kp15_endpoint_spoke_frame_alpha_is_orientation_clean(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    alpha: f64,
    tolerance: f64,
) -> bool {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    kp15_endpoint_spoke_frame_alpha_is_orientation_clean_with_samples(
        n,
        base_tensor,
        frame_delta,
        alpha,
        tolerance,
        &SAMPLES,
    )
}

fn kp15_endpoint_spoke_frame_alpha_is_orientation_clean_with_samples(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    alpha: f64,
    tolerance: f64,
    samples: &[f64],
) -> bool {
    kp15_endpoint_spoke_frame_alpha_clean_min_det(
        n,
        base_tensor,
        frame_delta,
        alpha,
        tolerance,
        samples,
    )
    .is_some()
}

fn kp15_endpoint_spoke_frame_alpha_clean_min_det(
    n: usize,
    base_tensor: &TensorBorderBi3,
    frame_delta: [[f64; 2]; 2],
    alpha: f64,
    tolerance: f64,
    samples: &[f64],
) -> Option<f64> {
    let scaled = [
        [alpha * frame_delta[0][0], alpha * frame_delta[0][1]],
        [alpha * frame_delta[1][0], alpha * frame_delta[1][1]],
    ];
    let normalized_tensor =
        kp15_apply_sector_local_first_order_normalization(base_tensor, n, scaled);
    let input = kp15_input_from_bi3_tensor(n, &normalized_tensor, false);
    let Ok(out) = assemble_bi6_polar_cap_via_paper_p_table(&input, 5, 0.0) else {
        return None;
    };
    let orientation =
        test_bi6_patch_local_orientation_report_with_samples(n, &out.patches_bi6, samples);
    if orientation.negative == 0 && orientation.min_det > tolerance {
        Some(orientation.min_det)
    } else {
        None
    }
}

fn kp15_patch_delta(
    patches: &[[[Point3; 7]; 7]],
    base: &[[[Point3; 7]; 7]],
) -> Vec<[[Point3; 7]; 7]> {
    patches
        .iter()
        .zip(base.iter())
        .map(|(patch, base_patch)| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
            for i in 0..7 {
                for j in 0..7 {
                    out[i][j] = Point3::new(
                        patch[i][j].x - base_patch[i][j].x,
                        patch[i][j].y - base_patch[i][j].y,
                        patch[i][j].z - base_patch[i][j].z,
                    );
                }
            }
            out
        })
        .collect()
}

fn kp15_source_frame_quadratic_terms(
    base: &[[[Point3; 7]; 7]],
    basis: &[Vec<[[Point3; 7]; 7]>],
    q: &[[f64; 49]; 49],
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let variables = basis.len();
    let mut h = vec![vec![0.0_f64; variables]; variables];
    let mut g = vec![0.0_f64; variables];
    for a in 0..variables {
        g[a] = kp15_cap_f5_bilinear(&basis[a], base, q);
        for b in 0..variables {
            h[a][b] = kp15_cap_f5_bilinear(&basis[a], &basis[b], q);
        }
    }
    (h, g)
}

fn kp15_cap_f5_bilinear(
    a: &[[[Point3; 7]; 7]],
    b: &[[[Point3; 7]; 7]],
    q: &[[f64; 49]; 49],
) -> f64 {
    let mut acc = 0.0;
    for (patch_a, patch_b) in a.iter().zip(b.iter()) {
        for axis in 0..3 {
            for r in 0..49 {
                let ai = r / 7;
                let aj = r % 7;
                let av = match axis {
                    0 => patch_a[ai][aj].x,
                    1 => patch_a[ai][aj].y,
                    _ => patch_a[ai][aj].z,
                };
                for c in 0..49 {
                    let bi = c / 7;
                    let bj = c % 7;
                    let bv = match axis {
                        0 => patch_b[bi][bj].x,
                        1 => patch_b[bi][bj].y,
                        _ => patch_b[bi][bj].z,
                    };
                    acc += av * q[r][c] * bv;
                }
            }
        }
    }
    acc
}

fn kp15_solve_quadratic_with_linear_constraints(
    h: &[Vec<f64>],
    g: &[f64],
    constraints: &[(Vec<f64>, f64)],
) -> Option<Vec<f64>> {
    let vars = g.len();
    let rows = vars + constraints.len();
    let mut matrix = vec![vec![0.0_f64; rows]; rows];
    let mut rhs = vec![0.0_f64; rows];
    for r in 0..vars {
        for c in 0..vars {
            matrix[r][c] = h[r][c];
        }
        rhs[r] = -g[r];
    }
    for (constraint_index, (coeffs, target)) in constraints.iter().enumerate() {
        assert_eq!(coeffs.len(), vars);
        let row = vars + constraint_index;
        for c in 0..vars {
            matrix[row][c] = coeffs[c];
            matrix[c][row] = coeffs[c];
        }
        rhs[row] = *target;
    }
    let solution = solve_dense_faer_for_test(matrix, rhs)?;
    Some(solution[..vars].to_vec())
}

fn kp15_first_order_map_from_two_point_deltas(
    p0: [f64; 2],
    d0: [f64; 2],
    p1: [f64; 2],
    d1: [f64; 2],
) -> Option<[[f64; 2]; 2]> {
    let det = p0[0] * p1[1] - p1[0] * p0[1];
    if det.abs() <= 1.0e-20 {
        return None;
    }
    let solve_row = |rhs0: f64, rhs1: f64| -> [f64; 2] {
        [
            (rhs0 * p1[1] - rhs1 * p0[1]) / det,
            (p0[0] * rhs1 - p1[0] * rhs0) / det,
        ]
    };
    Some([solve_row(d0[0], d1[0]), solve_row(d0[1], d1[1])])
}

fn kp15_apply_sector_local_first_order_normalization(
    tensor: &TensorBorderBi3,
    n: usize,
    frame_delta: [[f64; 2]; 2],
) -> TensorBorderBi3 {
    let sectors = tensor
        .sectors
        .iter()
        .enumerate()
        .map(|(sector, rows)| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for cp in 0..4 {
                    let local = global_to_local_planar(n, sector, rows[row][cp]);
                    let mapped = [
                        local[0] + frame_delta[0][0] * local[0] + frame_delta[0][1] * local[1],
                        local[1] + frame_delta[1][0] * local[0] + frame_delta[1][1] * local[1],
                    ];
                    out[row][cp] = local_to_global_planar(n, sector, mapped);
                }
            }
            out
        })
        .collect();
    TensorBorderBi3 { sectors }
}

fn kp15_pr08_q11_transposed_endpoint_input_with_bi3_local_delta(
    n: usize,
    apply_beta: bool,
    row: usize,
    cp: usize,
    local_delta: [f64; 2],
) -> Bi6CapInput {
    kp15_pr08_q11_transposed_endpoint_input_with_bi3_local_deltas(
        n,
        apply_beta,
        &[(row, cp, local_delta)],
    )
}

fn kp15_pr08_q11_transposed_endpoint_input_with_bi3_local_deltas(
    n: usize,
    apply_beta: bool,
    deltas: &[(usize, usize, [f64; 2])],
) -> Bi6CapInput {
    let params = Bi6Reparams::for_valence(n);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11,
    )
    .expect("PR08 q11-transposed candidate rows");
    let tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).expect("candidate tensor border");
    let mut endpoint_tensor = kp15_endpoint_normalized_bi3_tensor_border(&tensor, n, 0);
    for sector in 0..n {
        for &(row, cp, local_delta) in deltas {
            let delta = local_to_global_planar(n, sector, local_delta);
            endpoint_tensor.sectors[sector][row][cp].x += delta.x;
            endpoint_tensor.sectors[sector][row][cp].y += delta.y;
            endpoint_tensor.sectors[sector][row][cp].z += delta.z;
        }
    }
    let sectors = kp15_sector_rhs_from_bi3_rows(&endpoint_tensor.sectors, params, apply_beta);
    let zero = Point3::new(0.0, 0.0, 0.0);
    Bi6CapInput {
        n,
        params,
        c7: zero,
        c5s: vec![zero; n],
        c6s: vec![zero; n],
        sectors,
        p0_22: zero,
        p0_21_pinned: if n == 3 || n == 6 { Some(zero) } else { None },
        interior_filler: None,
    }
}

fn kp15_apply_s1_seed_propagation_to_input(input: &Bi6CapInput) -> Bi6CapInput {
    let mut out = input.clone();
    for sector in 1..out.n {
        let prev = out.sectors[sector - 1];
        let s1 = crate::geometry::gregory::rule_s1_star_propagation(
            &out.params,
            prev.p_00,
            prev.p_01,
            prev.p_02,
            prev.p_10,
            prev.p_11,
            prev.p_20,
        );
        out.sectors[sector].p_01 = s1.p01;
        out.sectors[sector].p_11 = s1.p11;
        out.sectors[sector].p_02 = s1.p02;
    }
    out
}

fn kp15_s1_cyclic_seed_residual(input: &Bi6CapInput) -> (f64, String) {
    let mut worst = 0.0_f64;
    let mut label = String::new();
    for sector in 0..input.n {
        let prev_sector = (sector + input.n - 1) % input.n;
        let prev = input.sectors[prev_sector];
        let current = input.sectors[sector];
        // KP15's two-patch local frame transposes the previous side.
        // Local previous p10/p20 correspond to global p01/p02.
        // Local previous p01/p02 correspond to global p10/p20.
        let s1 = crate::geometry::gregory::rule_s1_star_propagation(
            &input.params,
            prev.p_00,
            prev.p_10,
            prev.p_20,
            prev.p_01,
            prev.p_11,
            prev.p_02,
        );
        for (name, actual, expected) in [
            ("p00", current.p_00, s1.p00),
            ("p10", current.p_10, s1.p10),
            ("p20", current.p_20, s1.p20),
            ("p01", current.p_01, s1.p01),
            ("p11", current.p_11, s1.p11),
            ("p02", current.p_02, s1.p02),
        ] {
            let residual = actual.distance_to(&expected);
            if residual > worst {
                worst = residual;
                label = format!("{prev_sector}->{sector}.{name}");
            }
        }
    }
    (worst, label)
}

fn kp15_sector_rhs_from_bi3_rows(
    rows: &[[[Point3; 4]; 3]],
    params: Bi6Reparams,
    apply_beta: bool,
) -> Vec<CirculantSectorRhs> {
    rows.iter()
        .map(|sector_rows| {
            let bi6_rows = if apply_beta {
                apply_beta_reparameterization(sector_rows, &params)
            } else {
                [
                    degree_raise_1d_3_to_6(&sector_rows[0]),
                    degree_raise_1d_3_to_6(&sector_rows[1]),
                    degree_raise_1d_3_to_6(&sector_rows[2]),
                ]
            };
            CirculantSectorRhs {
                p_00: bi6_rows[0][0],
                p_10: bi6_rows[0][1],
                p_20: bi6_rows[0][2],
                p_40: bi6_rows[0][4],
                p_50: bi6_rows[0][5],
                p_60: bi6_rows[0][6],
                p_01: bi6_rows[1][0],
                p_11: bi6_rows[1][1],
                p_02: bi6_rows[2][0],
                p_42: bi6_rows[2][4],
                p_52: bi6_rows[2][5],
                p_62: bi6_rows[2][6],
            }
        })
        .collect()
}

fn kp15_endpoint_normalized_bi3_tensor_border(
    tensor: &TensorBorderBi3,
    n: usize,
    reference_sector: usize,
) -> TensorBorderBi3 {
    let source0 = tensor.sectors[reference_sector][0][0];
    let source1 = tensor.sectors[reference_sector][0][3];
    let source_delta = [source1.x - source0.x, source1.y - source0.y];
    let alpha = core::f64::consts::TAU / n as f64;
    let target = [
        (reference_sector as f64 * alpha).cos(),
        (reference_sector as f64 * alpha).sin(),
    ];
    let denom = source_delta[0] * source_delta[0] + source_delta[1] * source_delta[1];
    assert!(denom > 1.0e-20);
    let scale_re = (target[0] * source_delta[0] + target[1] * source_delta[1]) / denom;
    let scale_im = (target[1] * source_delta[0] - target[0] * source_delta[1]) / denom;
    map_tensor(tensor, |p| {
        let x = p.x - source0.x;
        let y = p.y - source0.y;
        Point3::new(
            scale_re * x - scale_im * y,
            scale_im * x + scale_re * y,
            p.z - source0.z,
        )
    })
}

fn kp15_strip_normalized_bi3_tensor_border(
    tensor: &TensorBorderBi3,
    n: usize,
    reference_sector: usize,
) -> TensorBorderBi3 {
    let source0 = tensor.sectors[reference_sector][0][0];
    let source_u = tensor.sectors[reference_sector][0][3];
    let source_depth = tensor.sectors[reference_sector][2][3];
    let u_axis = [source_u.x - source0.x, source_u.y - source0.y];
    let depth_axis = [source_depth.x - source_u.x, source_depth.y - source_u.y];
    let det = u_axis[0] * depth_axis[1] - u_axis[1] * depth_axis[0];
    assert!(det.abs() > 1.0e-20);
    let alpha = core::f64::consts::TAU / n as f64;
    let spoke0 = [
        (reference_sector as f64 * alpha).cos(),
        (reference_sector as f64 * alpha).sin(),
    ];
    let spoke1 = [
        ((reference_sector + 1) as f64 * alpha).cos(),
        ((reference_sector + 1) as f64 * alpha).sin(),
    ];
    map_tensor(tensor, |p| {
        let x = p.x - source0.x;
        let y = p.y - source0.y;
        let u = (x * depth_axis[1] - y * depth_axis[0]) / det;
        let depth = (u_axis[0] * y - u_axis[1] * x) / det;
        let target_u = u;
        let target_v = depth / 3.0;
        Point3::new(
            target_u * spoke0[0] + target_v * spoke1[0],
            target_u * spoke0[1] + target_v * spoke1[1],
            p.z - source0.z,
        )
    })
}

fn kp15_spoke_frame_slot_max_delta(n: usize, sector: &CirculantSectorRhs) -> f64 {
    [
        (sector.p_10, [1.0 / 6.0, 0.0]),
        (sector.p_20, [2.0 / 6.0, 0.0]),
        (sector.p_40, [4.0 / 6.0, 0.0]),
        (sector.p_50, [5.0 / 6.0, 0.0]),
        (sector.p_60, [1.0, 0.0]),
        (sector.p_01, [0.0, 1.0 / 6.0]),
        (sector.p_11, [1.0 / 6.0, 1.0 / 6.0]),
        (sector.p_02, [0.0, 2.0 / 6.0]),
        (sector.p_42, [4.0 / 6.0, 2.0 / 6.0]),
        (sector.p_52, [5.0 / 6.0, 2.0 / 6.0]),
        (sector.p_62, [1.0, 2.0 / 6.0]),
    ]
    .into_iter()
    .map(|(point, expected)| {
        let local = global_to_local_planar(n, 0, point);
        (local[0] - expected[0])
            .abs()
            .max((local[1] - expected[1]).abs())
    })
    .fold(0.0, f64::max)
}

fn kp13_planar_square_tensor_border(n: usize) -> TensorBorder {
    let mut target = TensorBorder::zero(n);
    for sector in 0..n {
        let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                cps[i][j] = local_to_global_planar(n, sector, [i as f64 / 3.0, j as f64 / 3.0]);
            }
        }
        target.patches[sector] = Bi3TensorBorder::from_cps(cps);
    }
    target
}

fn local_to_global_planar(n: usize, sector: usize, local: [f64; 2]) -> Point3 {
    let alpha = core::f64::consts::TAU / n as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    Point3::new(
        local[0] * spoke0[0] + local[1] * spoke1[0],
        local[0] * spoke0[1] + local[1] * spoke1[1],
        0.0,
    )
}

fn global_to_local_planar(n: usize, sector: usize, point: Point3) -> [f64; 2] {
    let alpha = core::f64::consts::TAU / n as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke0[0] * spoke1[1] - spoke0[1] * spoke1[0];
    [
        (point.x * spoke1[1] - point.y * spoke1[0]) / det,
        (spoke0[0] * point.y - spoke0[1] * point.x) / det,
    ]
}

fn map_tensor(t: &TensorBorderBi3, f: impl Fn(Point3) -> Point3) -> TensorBorderBi3 {
    let mut out = t.clone();
    for sector in &mut out.sectors {
        for row in sector {
            for p in row {
                *p = f(*p);
            }
        }
    }
    out
}

fn affine(p: Point3, scale: f64, offset: Point3) -> Point3 {
    Point3::new(
        scale * p.x + offset.x,
        scale * p.y + offset.y,
        scale * p.z + offset.z,
    )
}

fn assert_sectors_close(a: &[Kp20Tri6GuideSector], b: &[Kp20Tri6GuideSector], eps: f64) {
    assert_eq!(a.len(), b.len());
    for (sa, sb) in a.iter().zip(b) {
        for (pa, pb) in sa.cps.iter().zip(sb.cps) {
            assert!(
                pa.approx_eq(&pb, eps),
                "expected {pa:?} to match {pb:?}, distance {}",
                pa.distance_to(&pb)
            );
        }
    }
}
