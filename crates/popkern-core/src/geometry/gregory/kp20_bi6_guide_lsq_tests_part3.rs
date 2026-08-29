use super::*;
#[allow(unused_imports)]
use crate::geometry::gregory::kp15_bi6_local_constraint;

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
