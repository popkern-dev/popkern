use super::*;
#[allow(unused_imports)]
use crate::geometry::gregory::kp15_bi6_local_constraint;

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
        let baseline_rows = kp15_bi6_tensor_rows_for_test(&source_rows.sectors, params, true);
        let bb_rows = kp15_fig12_bb_beta_rows_for_test(&source_rows.sectors, params);
        for (label, rows) in [
            ("baseline-beta-row-semantics", baseline_rows),
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
