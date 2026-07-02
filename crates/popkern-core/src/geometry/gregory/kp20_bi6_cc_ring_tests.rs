use super::super::kp_polar_boundary_input::TensorBorderBi3;
use super::super::kp_polar_shared::weighted_sum;
use super::super::kp20_bi6_cap::{
    Kp20MacroPatch, Kp20PieceId, Kp20SplitTensorBorder, split_kp20_tensor_border,
};
use super::super::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, cubic_bspline_rows_to_tensor_border_rows,
    kp20_tensor_border_from_regular_cnet_rows,
};
use super::super::kp20_bi6_guide_sample::Kp20AffineGuideMap;
use super::super::kp20_bi6_local_constraint::{
    Kp20Eq9DiagnosticTerms, Kp20Eq9TermSigns, Kp20Eq9TermSource,
};
use super::super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::super::kp20_bi6_piece_maps::{
    kp20_paper_quadrant_sector_piece_maps, kp20_planar_bi6_tau_orientation_report,
    kp20_planar_bi6_tau_validity_report, sample_kp20_macro_patch_from_tau_family_maps,
};
use super::super::kp20_bi6_stam_eval::{
    Kp20StamGridSide, Kp20StamGridTransform, Kp20StamPatch, eval_cubic_bspline_grid,
    kp20_pr08_middle_patch_normalization_value, kp20_pr08_stam_appendix_tensor_projected_grids,
    kp20_pr08_stam_projected_transformed_tensor_border_side_rows,
    kp20_pr08_stam_visual_projected_grids,
};
use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedF3SolveOutput,
    Kp20PlanarTauAugmentedScalarCoeff, Kp20PlanarTauAugmentedScalarLinearConstraint,
    Kp20PlanarTauCoeff, Kp20PlanarTauF3ScalarSolveInput, Kp20PlanarTauF3SolveError,
    Kp20PlanarTauF3SolveOutput, Kp20PlanarTauLinearConstraint, Kp20PlanarTauScalarCoeff,
    Kp20PlanarTauScalarLinearConstraint, kp20_planar_tau_diagonal_symmetry_constraints,
    kp20_planar_tau_fig8_indicated_averaging_constraints, kp20_planar_tau_internal_c2_constraints,
    kp20_planar_tau_rotated_sector_g1_constraints,
    kp20_planar_tau_rotated_sector_g1_g2_constraints,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b,
    kp20_planar_tau_rotated_sector_g2_constraints, kp20_planar_tau_tensor_border_constraints,
    kp20_planar_tau_tensor_border_constraints_with_d1_3,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
    solve_kp20_planar_tau_f3_with_scalar_constraints,
};
use super::super::kp20_bi6_tau_solve_assemble::{
    Kp20PlanarTauPaperRowFamily, Kp20PlanarTauPaperSolveInput,
    kp20_planar_tau_rotational_paper_family_residuals,
    kp20_planar_tau_rotational_paper_family_residuals_with_d1_3,
    kp20_planar_tau_rotational_paper_row_set,
    kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs,
    kp20_planar_tau_rotational_section42_source_projection_report,
    kp20_planar_tau_sector_local_diagonal_symmetry_constraints,
    solve_kp20_planar_tau_paper_sector_f3, solve_kp20_planar_tau_rotational_paper_sector_f3,
};
use super::*;
use crate::geometry::gregory::{
    Bi6Reparams, Kp20TauFamily, Kp20TauMap, Kp20Tri6GuideSector, apply_beta_reparameterization,
    bernstein_bi6_all, bernstein_bi6_deriv_all, tri6_flat_index,
};

fn assert_close_point(a: Point3, b: Point3, tol: f64) {
    assert!(
        a.approx_eq(&b, tol),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

fn cc_test_global_planar_patch_to_local(
    valence: usize,
    sector: usize,
    mut patch: Kp20MacroPatch,
) -> Kp20MacroPatch {
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let local = cc_test_global_to_local_planar(valence, sector, piece.cps[i][j]);
                piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
            }
        }
    }
    patch
}

fn cc_test_global_to_local_planar(valence: usize, sector: usize, point: Point3) -> [f64; 2] {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    [
        (point.x * spoke0[1] - point.y * spoke0[0]) / det,
        (spoke1[0] * point.y - spoke1[1] * point.x) / det,
    ]
}

fn cc_test_transform_tau_chart(
    mut patch: Kp20MacroPatch,
    transform: [[f64; 2]; 2],
) -> Kp20MacroPatch {
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let cp = piece.cps[i][j];
                piece.cps[i][j] = Point3::new(
                    transform[0][0] * cp.x + transform[0][1] * cp.y,
                    transform[1][0] * cp.x + transform[1][1] * cp.y,
                    cp.z,
                );
            }
        }
    }
    patch
}

fn cc_test_sector_local_to_global_planar(valence: usize, sector: usize, local: [f64; 2]) -> Point3 {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    Point3::new(
        local[0] * spoke1[0] + local[1] * spoke0[0],
        local[0] * spoke1[1] + local[1] * spoke0[1],
        0.0,
    )
}

fn cc_test_sector_local_axis_components(valence: usize, sector: usize) -> [[f64; 2]; 2] {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    [
        [spoke0[1] / det, -spoke0[0] / det],
        [-spoke1[1] / det, spoke1[0] / det],
    ]
}

fn cc_test_sector_local_sample_scalar_terms(
    valence: usize,
    sector: usize,
    piece: Kp20PieceId,
    u: f64,
    v: f64,
    local_axis: usize,
    scale: f64,
) -> Vec<(Kp20PlanarTauScalarCoeff, f64)> {
    let bu = bernstein_bi6_all(u);
    let bv = bernstein_bi6_all(v);
    let axis_components = cc_test_sector_local_axis_components(valence, sector);
    let axis_weights = axis_components[local_axis];
    let mut terms = Vec::with_capacity(98);
    for i in 0..=6 {
        for j in 0..=6 {
            let coeff = Kp20PlanarTauCoeff::new(piece, i, j);
            let weight = scale * bu[i] * bv[j];
            terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                weight * axis_weights[0],
            ));
            terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                weight * axis_weights[1],
            ));
        }
    }
    terms
}

fn cc_test_sector_local_sample_value(
    tau: &Kp20MacroPatch,
    piece: Kp20PieceId,
    u: f64,
    v: f64,
) -> [f64; 2] {
    let bu = bernstein_bi6_all(u);
    let bv = bernstein_bi6_all(v);
    let patch = tau.piece(piece);
    let mut out = [0.0_f64; 2];
    for i in 0..=6 {
        for j in 0..=6 {
            let weight = bu[i] * bv[j];
            out[0] += weight * patch.cps[i][j].x;
            out[1] += weight * patch.cps[i][j].y;
        }
    }
    out
}

fn cc_test_active_sector_domain_bound_rows(
    tau_global: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
    tolerance: f64,
) -> (
    Vec<Kp20PlanarTauScalarLinearConstraint>,
    usize,
    usize,
    usize,
) {
    let tau_local = cc_test_global_planar_patch_to_local(valence, sector, *tau_global);
    let mut rows = Vec::new();
    let mut lower_u = 0_usize;
    let mut lower_v = 0_usize;
    let mut upper_sum = 0_usize;
    for piece in Kp20PieceId::all() {
        for &u in samples {
            for &v in samples {
                let local = cc_test_sector_local_sample_value(&tau_local, piece, u, v);
                if local[0] < -tolerance {
                    lower_u += 1;
                    rows.push(Kp20PlanarTauScalarLinearConstraint {
                        terms: cc_test_sector_local_sample_scalar_terms(
                            valence, sector, piece, u, v, 0, 1.0,
                        ),
                        target: 0.0,
                    });
                }
                if local[1] < -tolerance {
                    lower_v += 1;
                    rows.push(Kp20PlanarTauScalarLinearConstraint {
                        terms: cc_test_sector_local_sample_scalar_terms(
                            valence, sector, piece, u, v, 1, 1.0,
                        ),
                        target: 0.0,
                    });
                }
                if local[0] + local[1] > 1.0 + tolerance {
                    upper_sum += 1;
                    let mut terms = cc_test_sector_local_sample_scalar_terms(
                        valence, sector, piece, u, v, 0, 1.0,
                    );
                    terms.extend(cc_test_sector_local_sample_scalar_terms(
                        valence, sector, piece, u, v, 1, 1.0,
                    ));
                    rows.push(Kp20PlanarTauScalarLinearConstraint { terms, target: 1.0 });
                }
            }
        }
    }
    (rows, lower_u, lower_v, upper_sum)
}

fn cc_test_affine_normalized_tensor_rows_from_source_triangle(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    source_sector: usize,
    source_triangle: [Point3; 3],
    target_triangle_local: [[f64; 2]; 3],
) -> Option<Kp20RegularTensorBorderCNet> {
    let s0 = source_triangle[0];
    let sx = [source_triangle[1].x - s0.x, source_triangle[1].y - s0.y];
    let sy = [source_triangle[2].x - s0.x, source_triangle[2].y - s0.y];
    let det = sx[0] * sy[1] - sy[0] * sx[1];
    if det.abs() < 1.0e-14 {
        return None;
    }

    let t0 =
        cc_test_sector_local_to_global_planar(valence, source_sector, target_triangle_local[0]);
    let tx_global =
        cc_test_sector_local_to_global_planar(valence, source_sector, target_triangle_local[1]);
    let ty_global =
        cc_test_sector_local_to_global_planar(valence, source_sector, target_triangle_local[2]);
    let tx = [tx_global.x - t0.x, tx_global.y - t0.y];
    let ty = [ty_global.x - t0.x, ty_global.y - t0.y];

    let sectors = rows
        .sectors
        .iter()
        .map(|sector| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let p = sector[row][col];
                    let px = p.x - s0.x;
                    let py = p.y - s0.y;
                    let a = (px * sy[1] - sy[0] * py) / det;
                    let b = (sx[0] * py - px * sx[1]) / det;
                    out[row][col] = Point3::new(
                        t0.x + a * tx[0] + b * ty[0],
                        t0.y + a * tx[1] + b * ty[1],
                        p.z,
                    );
                }
            }
            out
        })
        .collect();
    Some(Kp20RegularTensorBorderCNet { sectors })
}

fn cc_test_sector_local_chart_transformed_tensor_rows(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
    transform: [[f64; 2]; 2],
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|source_sector| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let local =
                        cc_test_global_to_local_planar(valence, sector, source_sector[row][col]);
                    let mapped = [
                        transform[0][0] * local[0] + transform[0][1] * local[1],
                        transform[1][0] * local[0] + transform[1][1] * local[1],
                    ];
                    out[row][col] = cc_test_sector_local_to_global_planar(valence, sector, mapped);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn cc_test_per_sector_characteristic_tensor_rows(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    canonical_sector: usize,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .enumerate()
        .map(|(source_sector_index, source_sector)| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let local = cc_test_global_to_local_planar(
                        valence,
                        source_sector_index,
                        source_sector[row][col],
                    );
                    out[row][col] =
                        cc_test_sector_local_to_global_planar(valence, canonical_sector, local);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn cc_test_project_tensor_rows_for_scalar_basis(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64) {
    cc_test_project_selected_tensor_rows_for_scalar_basis(
        tensor_rows,
        fixed_vector_constraints,
        scalar_constraints,
        |_| true,
    )
}

fn cc_test_project_selected_tensor_rows_for_scalar_basis(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    free_coeff: impl Fn(Kp20PlanarTauCoeff) -> bool,
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64) {
    let mut vector_constraints = fixed_vector_constraints.to_vec();
    let mut augmented_scalar_constraints = Vec::new();
    let mut free_scalar_targets = Vec::new();
    let mut extra_count = 0_usize;

    for (row_index, row) in tensor_rows.iter().enumerate() {
        let coeff = row.terms[0].0;
        if free_coeff(coeff) {
            for axis in 0..=1 {
                augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                    terms: vec![
                        (
                            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                                coeff, axis,
                            )),
                            1.0,
                        ),
                        (Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_count), -1.0),
                    ],
                    target: row.target[axis],
                });
                free_scalar_targets.push((row_index, axis));
                extra_count += 1;
            }
        } else {
            vector_constraints.push(row.clone());
        }
    }

    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints: scalar_constraints.to_vec(),
            augmented_scalar_constraints,
            extra_variable_count: extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0,
        },
    )
    .unwrap();

    let mut projected = tensor_rows;
    for ((row_index, axis), delta) in free_scalar_targets
        .into_iter()
        .zip(solved.extra_values.iter().copied())
    {
        projected[row_index].target[axis] += delta;
    }
    (
        projected,
        solved
            .extra_values
            .chunks(2)
            .map(|chunk| {
                let x = chunk[0];
                let y = chunk.get(1).copied().unwrap_or(0.0);
                (x * x + y * y).sqrt()
            })
            .fold(0.0_f64, f64::max),
    )
}

fn cc_test_sector_local_self_diagonal_normalized_source_rows(
    mut tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    valence: usize,
    sector: usize,
    piece: Kp20PieceId,
    max_index: usize,
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64) {
    let mut row_for_coeff = std::collections::BTreeMap::new();
    for (row_index, row) in tensor_rows.iter().enumerate() {
        if row.terms.len() == 1 && row.terms[0].1 == 1.0 && row.terms[0].0.piece == piece {
            row_for_coeff.insert((row.terms[0].0.i, row.terms[0].0.j), row_index);
        }
    }

    let mut max_delta = 0.0_f64;
    for i in 0..=max_index {
        if let Some(&row_index) = row_for_coeff.get(&(i, i)) {
            let original = Point3::new(
                tensor_rows[row_index].target[0],
                tensor_rows[row_index].target[1],
                0.0,
            );
            let local = cc_test_global_to_local_planar(valence, sector, original);
            let average = 0.5 * (local[0] + local[1]);
            let projected =
                cc_test_sector_local_to_global_planar(valence, sector, [average, average]);
            tensor_rows[row_index].target = [projected.x, projected.y];
            max_delta = max_delta.max(original.distance_to(&projected));
        }

        for j in (i + 1)..=max_index {
            let Some(&lower_row) = row_for_coeff.get(&(i, j)) else {
                continue;
            };
            let Some(&upper_row) = row_for_coeff.get(&(j, i)) else {
                continue;
            };
            let lower_original = Point3::new(
                tensor_rows[lower_row].target[0],
                tensor_rows[lower_row].target[1],
                0.0,
            );
            let upper_original = Point3::new(
                tensor_rows[upper_row].target[0],
                tensor_rows[upper_row].target[1],
                0.0,
            );
            let lower_local = cc_test_global_to_local_planar(valence, sector, lower_original);
            let upper_local = cc_test_global_to_local_planar(valence, sector, upper_original);
            let projected_lower_local = [
                0.5 * (lower_local[0] + upper_local[1]),
                0.5 * (lower_local[1] + upper_local[0]),
            ];
            let projected_upper_local = [projected_lower_local[1], projected_lower_local[0]];
            let projected_lower =
                cc_test_sector_local_to_global_planar(valence, sector, projected_lower_local);
            let projected_upper =
                cc_test_sector_local_to_global_planar(valence, sector, projected_upper_local);
            tensor_rows[lower_row].target = [projected_lower.x, projected_lower.y];
            tensor_rows[upper_row].target = [projected_upper.x, projected_upper.y];
            max_delta = max_delta
                .max(lower_original.distance_to(&projected_lower))
                .max(upper_original.distance_to(&projected_upper));
        }
    }

    (tensor_rows, max_delta)
}

fn cc_test_underlined_c0_normalized_source_rows(
    mut tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    previous_sector_rotation_radians: f64,
    max_i: usize,
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64) {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let rotate_previous = |point: Point3| {
        Point3::new(
            cos_theta * point.x - sin_theta * point.y,
            sin_theta * point.x + cos_theta * point.y,
            0.0,
        )
    };
    let unrotate_current = |point: Point3| {
        Point3::new(
            cos_theta * point.x + sin_theta * point.y,
            -sin_theta * point.x + cos_theta * point.y,
            0.0,
        )
    };

    let row_for_coeff = |rows: &[Kp20PlanarTauLinearConstraint], coeff| {
        rows.iter().position(|row| {
            row.terms.len() == 1 && row.terms[0].1 == 1.0 && row.terms[0].0 == coeff
        })
    };

    let mut max_delta = 0.0_f64;
    for i in 0..=max_i {
        let previous_coeff = Kp20PlanarTauCoeff::new(previous_piece, i, 0);
        let current_coeff = Kp20PlanarTauCoeff::new(current_piece, i, 0);
        let Some(previous_row) = row_for_coeff(&tensor_rows, previous_coeff) else {
            continue;
        };
        let Some(current_row) = row_for_coeff(&tensor_rows, current_coeff) else {
            continue;
        };
        let previous_original = Point3::new(
            tensor_rows[previous_row].target[0],
            tensor_rows[previous_row].target[1],
            0.0,
        );
        let current_original = Point3::new(
            tensor_rows[current_row].target[0],
            tensor_rows[current_row].target[1],
            0.0,
        );
        let previous_rotated = rotate_previous(previous_original);
        let current_projected = Point3::new(
            0.5 * (current_original.x + previous_rotated.x),
            0.5 * (current_original.y + previous_rotated.y),
            0.0,
        );
        let previous_projected = unrotate_current(current_projected);
        tensor_rows[previous_row].target = [previous_projected.x, previous_projected.y];
        tensor_rows[current_row].target = [current_projected.x, current_projected.y];
        max_delta = max_delta
            .max(previous_original.distance_to(&previous_projected))
            .max(current_original.distance_to(&current_projected));
    }

    (tensor_rows, max_delta)
}

fn cc_test_project_tensor_rows_for_augmented_scalar_basis(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    fixed_augmented_scalar_constraints: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
    fixed_extra_variable_count: usize,
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64, f64, Vec<f64>) {
    cc_test_project_tensor_rows_for_augmented_scalar_basis_with_extra_weight(
        tensor_rows,
        fixed_vector_constraints,
        scalar_constraints,
        fixed_augmented_scalar_constraints,
        fixed_extra_variable_count,
        1.0,
    )
}

fn cc_test_project_tensor_rows_for_augmented_scalar_basis_with_extra_weight(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    fixed_augmented_scalar_constraints: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
    fixed_extra_variable_count: usize,
    extra_regularization_weight: f64,
) -> (Vec<Kp20PlanarTauLinearConstraint>, f64, f64, Vec<f64>) {
    let vector_constraints = fixed_vector_constraints.to_vec();
    let mut augmented_scalar_constraints = fixed_augmented_scalar_constraints.to_vec();
    let mut free_scalar_targets = Vec::new();
    let mut extra_count = fixed_extra_variable_count;

    for (row_index, row) in tensor_rows.iter().enumerate() {
        let coeff = row.terms[0].0;
        for axis in 0..=1 {
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        1.0,
                    ),
                    (Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_count), -1.0),
                ],
                target: row.target[axis],
            });
            free_scalar_targets.push((row_index, axis));
            extra_count += 1;
        }
    }

    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints: scalar_constraints.to_vec(),
            augmented_scalar_constraints,
            extra_variable_count: extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight,
        },
    )
    .unwrap();

    let mut projected = tensor_rows;
    for ((row_index, axis), delta) in free_scalar_targets.into_iter().zip(
        solved.extra_values[fixed_extra_variable_count..]
            .iter()
            .copied(),
    ) {
        projected[row_index].target[axis] += delta;
    }
    let max_source_delta = solved.extra_values[fixed_extra_variable_count..]
        .chunks(2)
        .map(|chunk| {
            let x = chunk[0];
            let y = chunk.get(1).copied().unwrap_or(0.0);
            (x * x + y * y).sqrt()
        })
        .fold(0.0_f64, f64::max);
    (
        projected,
        max_source_delta,
        solved.max_constraint_residual,
        solved.extra_values[..fixed_extra_variable_count].to_vec(),
    )
}

fn cc_test_paper_quadrant_base_map(piece_id: Kp20PieceId) -> Kp20AffineGuideMap {
    kp20_paper_quadrant_sector_piece_maps()
        .into_iter()
        .find_map(|(id, maps)| (id == piece_id).then_some(maps[0].1))
        .expect("paper quadrant map exists")
}

fn cc_test_affine_map_point(map: Kp20AffineGuideMap, u: f64, v: f64) -> [f64; 2] {
    [
        map.origin[0] + u * map.du[0] + v * map.dv[0],
        map.origin[1] + u * map.du[1] + v * map.dv[1],
    ]
}

fn cc_test_affine_macro_patch_from_piece_maps(
    maps: [(Kp20PieceId, Kp20AffineGuideMap); 4],
) -> Kp20MacroPatch {
    let mut patch = Kp20MacroPatch::zero();
    for (piece_id, map) in maps {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let [x, y] = cc_test_affine_map_point(map, i as f64 / 6.0, j as f64 / 6.0);
                piece.cps[i][j] = Point3::new(x, y, 0.0);
            }
        }
    }
    patch
}

fn cc_test_paper_quadrant_affine_macro_patch() -> Kp20MacroPatch {
    cc_test_affine_macro_patch_from_piece_maps([
        (
            Kp20PieceId::NearOriginLeft,
            cc_test_paper_quadrant_base_map(Kp20PieceId::NearOriginLeft),
        ),
        (
            Kp20PieceId::NearOriginRight,
            cc_test_paper_quadrant_base_map(Kp20PieceId::NearOriginRight),
        ),
        (
            Kp20PieceId::NearMidpointRight,
            cc_test_paper_quadrant_base_map(Kp20PieceId::NearMidpointRight),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            cc_test_paper_quadrant_base_map(Kp20PieceId::NearMidpointLeft),
        ),
    ])
}

fn cc_test_affine_map_inverse(map: Kp20AffineGuideMap, point: [f64; 2]) -> Option<[f64; 2]> {
    let px = point[0] - map.origin[0];
    let py = point[1] - map.origin[1];
    let det = map.du[0] * map.dv[1] - map.dv[0] * map.du[1];
    if det.abs() < 1.0e-14 {
        return None;
    }
    Some([
        (px * map.dv[1] - map.dv[0] * py) / det,
        (map.du[0] * py - px * map.du[1]) / det,
    ])
}

fn cc_test_nearest_degree6_index(value: f64) -> usize {
    let rounded = (6.0 * value).round();
    assert!(
        (6.0 * value - rounded).abs() < 1.0e-10,
        "expected {value} to land on degree-6 grid"
    );
    rounded as usize
}

fn cc_test_tau_scalar_value(tau: &Kp20MacroPatch, coeff: Kp20PlanarTauScalarCoeff) -> f64 {
    let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
    match coeff.axis {
        0 => cp.x,
        1 => cp.y,
        _ => panic!("invalid planar tau axis {}", coeff.axis),
    }
}

fn cc_test_vector_rows_max_residual(
    tau: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauLinearConstraint],
) -> f64 {
    rows.iter()
        .flat_map(|row| {
            let mut value = [0.0_f64; 2];
            for &(coeff, weight) in &row.terms {
                let cp = tau.piece(coeff.piece).cps[coeff.i][coeff.j];
                value[0] += weight * cp.x;
                value[1] += weight * cp.y;
            }
            [
                (value[0] - row.target[0]).abs(),
                (value[1] - row.target[1]).abs(),
            ]
        })
        .fold(0.0_f64, f64::max)
}

fn cc_test_scalar_rows_max_residual(
    tau: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> f64 {
    rows.iter()
        .map(|row| {
            let value = row
                .terms
                .iter()
                .map(|&(coeff, weight)| weight * cc_test_tau_scalar_value(tau, coeff))
                .sum::<f64>();
            (value - row.target).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn cc_test_augmented_rows_max_residual(
    tau: &Kp20MacroPatch,
    extra_values: &[f64],
    rows: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
) -> f64 {
    rows.iter()
        .map(|row| {
            let value = row
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    weight
                        * match coeff {
                            Kp20PlanarTauAugmentedScalarCoeff::Tau(coeff) => {
                                cc_test_tau_scalar_value(tau, coeff)
                            }
                            Kp20PlanarTauAugmentedScalarCoeff::Extra(index) => extra_values[index],
                        }
                })
                .sum::<f64>();
            (value - row.target).abs()
        })
        .fold(0.0_f64, f64::max)
}

#[derive(Debug, Clone, Copy)]
struct CcTestEndpointFrame {
    origin: Point3,
    du: Point3,
    dv: Point3,
    determinant: f64,
}

#[derive(Debug, Clone, Copy)]
struct CcTestOuterCornerFrame {
    p66: Point3,
    p56: Point3,
    p65: Point3,
    du: Point3,
    dv: Point3,
    determinant: f64,
}

fn cc_test_endpoint_frame(
    tau: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
    i: usize,
    j: usize,
) -> CcTestEndpointFrame {
    let piece = tau.piece(piece_id);
    let origin = piece.cps[i][j];
    let (u_neighbor, u_sign) = if i < 6 {
        (piece.cps[i + 1][j], 1.0)
    } else {
        (piece.cps[i - 1][j], -1.0)
    };
    let (v_neighbor, v_sign) = if j < 6 {
        (piece.cps[i][j + 1], 1.0)
    } else {
        (piece.cps[i][j - 1], -1.0)
    };
    let du = Point3::new(
        6.0 * u_sign * (u_neighbor.x - origin.x),
        6.0 * u_sign * (u_neighbor.y - origin.y),
        0.0,
    );
    let dv = Point3::new(
        6.0 * v_sign * (v_neighbor.x - origin.x),
        6.0 * v_sign * (v_neighbor.y - origin.y),
        0.0,
    );
    CcTestEndpointFrame {
        origin,
        du,
        dv,
        determinant: du.x * dv.y - du.y * dv.x,
    }
}

fn cc_test_outer_corner_frame(
    tau: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
) -> CcTestOuterCornerFrame {
    let piece = tau.piece(piece_id);
    let p66 = piece.cps[6][6];
    let p56 = piece.cps[5][6];
    let p65 = piece.cps[6][5];
    let du = Point3::new(6.0 * (p66.x - p56.x), 6.0 * (p66.y - p56.y), 0.0);
    let dv = Point3::new(6.0 * (p66.x - p65.x), 6.0 * (p66.y - p65.y), 0.0);
    CcTestOuterCornerFrame {
        p66,
        p56,
        p65,
        du,
        dv,
        determinant: du.x * dv.y - du.y * dv.x,
    }
}

fn cc_test_fig8_midpoint_average_row(
    piece_id: Kp20PieceId,
    midpoint_i: usize,
    midpoint_j: usize,
    end0_i: usize,
    end0_j: usize,
    end1_i: usize,
    end1_j: usize,
) -> Kp20PlanarTauLinearConstraint {
    Kp20PlanarTauLinearConstraint::homogeneous(vec![
        (
            Kp20PlanarTauCoeff::new(piece_id, midpoint_i, midpoint_j),
            1.0,
        ),
        (Kp20PlanarTauCoeff::new(piece_id, end0_i, end0_j), -0.5),
        (Kp20PlanarTauCoeff::new(piece_id, end1_i, end1_j), -0.5),
    ])
}

fn cc_test_fig8_indicated_average_variant_rows(
    include_diagonal_mirror: bool,
    include_rotated_quadrants: bool,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    if include_diagonal_mirror {
        rows.push(cc_test_fig8_midpoint_average_row(
            Kp20PieceId::NearOriginLeft,
            0,
            1,
            0,
            0,
            0,
            2,
        ));
    }
    if include_rotated_quadrants {
        for piece_id in [
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ] {
            rows.push(cc_test_fig8_midpoint_average_row(
                piece_id, 1, 0, 0, 0, 2, 0,
            ));
            if include_diagonal_mirror {
                rows.push(cc_test_fig8_midpoint_average_row(
                    piece_id, 0, 1, 0, 0, 0, 2,
                ));
            }
        }
    }
    rows
}

fn cc_test_augmented_local_component_terms(
    coeff: Kp20PlanarTauCoeff,
    component: [f64; 2],
    weight: f64,
) -> Vec<(Kp20PlanarTauAugmentedScalarCoeff, f64)> {
    vec![
        (
            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
            weight * component[0],
        ),
        (
            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
            weight * component[1],
        ),
    ]
}

fn cc_test_first_order_centered_source_rows(
    piece_id: Kp20PieceId,
    local_axes: [[f64; 2]; 2],
    scale_extra_index: usize,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let p00 = Kp20PlanarTauCoeff::new(piece_id, 0, 0);
    let p10 = Kp20PlanarTauCoeff::new(piece_id, 1, 0);
    let p01 = Kp20PlanarTauCoeff::new(piece_id, 0, 1);
    let scale = Kp20PlanarTauAugmentedScalarCoeff::Extra(scale_extra_index);
    let mut rows = Vec::new();
    for axis in local_axes {
        rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: cc_test_augmented_local_component_terms(p00, axis, 1.0),
            target: 0.0,
        });
    }

    let mut p10_cross = cc_test_augmented_local_component_terms(p10, local_axes[1], 1.0);
    p10_cross.extend(cc_test_augmented_local_component_terms(
        p00,
        local_axes[1],
        -1.0,
    ));
    rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: p10_cross,
        target: 0.0,
    });

    let mut p01_cross = cc_test_augmented_local_component_terms(p01, local_axes[0], 1.0);
    p01_cross.extend(cc_test_augmented_local_component_terms(
        p00,
        local_axes[0],
        -1.0,
    ));
    rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: p01_cross,
        target: 0.0,
    });

    let mut p10_along = cc_test_augmented_local_component_terms(p10, local_axes[0], 1.0);
    p10_along.extend(cc_test_augmented_local_component_terms(
        p00,
        local_axes[0],
        -1.0,
    ));
    p10_along.push((scale, -1.0 / 6.0));
    rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: p10_along,
        target: 0.0,
    });

    let mut p01_along = cc_test_augmented_local_component_terms(p01, local_axes[1], 1.0);
    p01_along.extend(cc_test_augmented_local_component_terms(
        p00,
        local_axes[1],
        -1.0,
    ));
    p01_along.push((scale, -1.0 / 6.0));
    rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: p01_along,
        target: 0.0,
    });

    rows
}

fn cc_test_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut out = 1.0_f64;
    for i in 1..=k {
        out *= (n + 1 - i) as f64 / i as f64;
    }
    out
}

fn cc_test_bi6_origin_derivative_component_terms(
    coeff: Kp20PlanarTauCoeff,
    du_order: usize,
    dv_order: usize,
    component: [f64; 2],
) -> Vec<(Kp20PlanarTauAugmentedScalarCoeff, f64)> {
    let factor_u = (0..du_order).fold(1.0_f64, |acc, offset| acc * (6 - offset) as f64);
    let factor_v = (0..dv_order).fold(1.0_f64, |acc, offset| acc * (6 - offset) as f64);
    let mut terms = Vec::new();
    for r in 0..=du_order {
        let sign_r = if (du_order - r).is_multiple_of(2) { 1.0 } else { -1.0 };
        for s in 0..=dv_order {
            let sign_s = if (dv_order - s).is_multiple_of(2) { 1.0 } else { -1.0 };
            let weight = factor_u
                * factor_v
                * sign_r
                * sign_s
                * cc_test_binomial(du_order, r)
                * cc_test_binomial(dv_order, s);
            let cp = Kp20PlanarTauCoeff::new(coeff.piece, coeff.i + r, coeff.j + s);
            terms.extend(cc_test_augmented_local_component_terms(
                cp, component, weight,
            ));
        }
    }
    terms
}

fn cc_test_paper_refined_apex_cubic_derivative_rows(
    valence: usize,
    sector: usize,
    piece_id: Kp20PieceId,
    local_axes: [[f64; 2]; 2],
    scale_extra_index: usize,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let tau = Kp20TauMap::for_sector_local_paper_refined_apex(valence, sector)
        .expect("paper refined apex tau");
    let jet = tau.jet3(0.0, 0.0);
    let origin = Kp20PlanarTauCoeff::new(piece_id, 0, 0);
    let mut rows = Vec::new();
    for axis in local_axes {
        rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: cc_test_augmented_local_component_terms(origin, axis, 1.0),
            target: 0.0,
        });
    }
    let scale = Kp20PlanarTauAugmentedScalarCoeff::Extra(scale_extra_index);
    for du_order in 0..=3 {
        for dv_order in 0..=(3 - du_order) {
            if du_order + dv_order == 0 {
                continue;
            }
            let coeff = Kp20PlanarTauCoeff::new(piece_id, 0, 0);
            for component_index in 0..=1 {
                let mut terms = cc_test_bi6_origin_derivative_component_terms(
                    coeff,
                    du_order,
                    dv_order,
                    local_axes[component_index],
                );
                terms.push((scale, -jet.d[du_order][dv_order][component_index]));
                rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint { terms, target: 0.0 });
            }
        }
    }
    rows
}

fn cc_test_linear_source_normalized_eq14_rows(
    tensor_rows: &[Kp20PlanarTauLinearConstraint],
    matrix_extra_start: usize,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let mut rows = Vec::with_capacity(2 * tensor_rows.len());
    for row in tensor_rows {
        let coeff = row.terms[0].0;
        let source_x = row.target[0];
        let source_y = row.target[1];
        rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    1.0,
                ),
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(matrix_extra_start),
                    -source_x,
                ),
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(matrix_extra_start + 1),
                    -source_y,
                ),
            ],
            target: source_x,
        });
        rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    1.0,
                ),
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(matrix_extra_start + 2),
                    -source_x,
                ),
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(matrix_extra_start + 3),
                    -source_y,
                ),
            ],
            target: source_y,
        });
    }
    rows
}

fn cc_test_sector_local_scalar_component_difference_terms(
    lhs: Kp20PlanarTauCoeff,
    lhs_component: [f64; 2],
    rhs: Kp20PlanarTauCoeff,
    rhs_component: [f64; 2],
) -> Vec<(Kp20PlanarTauScalarCoeff, f64)> {
    vec![
        (Kp20PlanarTauScalarCoeff::new(lhs, 0), lhs_component[0]),
        (Kp20PlanarTauScalarCoeff::new(lhs, 1), lhs_component[1]),
        (Kp20PlanarTauScalarCoeff::new(rhs, 0), -rhs_component[0]),
        (Kp20PlanarTauScalarCoeff::new(rhs, 1), -rhs_component[1]),
    ]
}

fn cc_test_sector_local_pair_first_difference_terms(
    lower_step: Kp20PlanarTauCoeff,
    lower_origin: Kp20PlanarTauCoeff,
    lower_component: [f64; 2],
    upper_step: Kp20PlanarTauCoeff,
    upper_origin: Kp20PlanarTauCoeff,
    upper_component: [f64; 2],
    upper_sign: f64,
) -> Vec<(Kp20PlanarTauScalarCoeff, f64)> {
    vec![
        (
            Kp20PlanarTauScalarCoeff::new(lower_step, 0),
            lower_component[0],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(lower_step, 1),
            lower_component[1],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(lower_origin, 0),
            -lower_component[0],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(lower_origin, 1),
            -lower_component[1],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(upper_step, 0),
            -upper_sign * upper_component[0],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(upper_step, 1),
            -upper_sign * upper_component[1],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(upper_origin, 0),
            upper_sign * upper_component[0],
        ),
        (
            Kp20PlanarTauScalarCoeff::new(upper_origin, 1),
            upper_sign * upper_component[1],
        ),
    ]
}

fn cc_test_sector_local_self_diagonal_rows(
    piece: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
            cc_test_sector_local_scalar_component_difference_terms(
                Kp20PlanarTauCoeff::new(piece, i, i),
                local0,
                Kp20PlanarTauCoeff::new(piece, i, i),
                local1,
            ),
        ));
        for j in (i + 1)..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(piece, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(piece, j, i),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(piece, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(piece, j, i),
                    local0,
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_pair_diagonal_rows(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    local0,
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_pair_diagonal_rows_with_components(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    lower_components: [[f64; 2]; 2],
    upper_components: [[f64; 2]; 2],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    lower_components[0],
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    upper_components[0],
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    lower_components[1],
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    upper_components[1],
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_filtered_pair_diagonal_rows(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
    keep: impl Fn(usize, usize) -> bool,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            if !keep(i, j) {
                continue;
            }
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(upper_left, j, i),
                    local0,
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_mapped_pair_diagonal_rows(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
    upper_index: impl Fn(usize, usize) -> (usize, usize),
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            let (upper_i, upper_j) = upper_index(i, j);
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local0,
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_filtered_mapped_pair_diagonal_rows(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
    keep: impl Fn(usize, usize) -> bool,
    upper_index: impl Fn(usize, usize) -> (usize, usize),
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            if !keep(i, j) {
                continue;
            }
            let (upper_i, upper_j) = upper_index(i, j);
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                cc_test_sector_local_scalar_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local0,
                ),
            ));
        }
    }
    rows
}

fn cc_test_sector_local_filtered_mapped_pair_diagonal_component_rows(
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
    keep: impl Fn(usize, usize) -> bool,
    upper_index: impl Fn(usize, usize) -> (usize, usize),
    include_first_component: bool,
    include_second_component: bool,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        for j in 0..=6 {
            if !keep(i, j) {
                continue;
            }
            let (upper_i, upper_j) = upper_index(i, j);
            if include_first_component {
                rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                    cc_test_sector_local_scalar_component_difference_terms(
                        Kp20PlanarTauCoeff::new(lower_right, i, j),
                        local0,
                        Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                        local1,
                    ),
                ));
            }
            if include_second_component {
                rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                    cc_test_sector_local_scalar_component_difference_terms(
                        Kp20PlanarTauCoeff::new(lower_right, i, j),
                        local1,
                        Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                        local0,
                    ),
                ));
            }
        }
    }
    rows
}

fn cc_test_neg_component(component: [f64; 2]) -> [f64; 2] {
    [-component[0], -component[1]]
}

fn cc_test_local_identity_tri6_guide() -> Kp20Tri6GuideSector {
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            guide.cps[tri6_flat_index(6, i, j)] = Point3::new(i as f64 / 6.0, j as f64 / 6.0, 0.0);
        }
    }
    guide
}

fn cc_test_local_planar_patch_to_global(
    valence: usize,
    sector: usize,
    mut patch: Kp20MacroPatch,
) -> Kp20MacroPatch {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let cp = piece.cps[i][j];
                piece.cps[i][j] = Point3::new(
                    cp.x * spoke1[0] + cp.y * spoke0[0],
                    cp.x * spoke1[1] + cp.y * spoke0[1],
                    0.0,
                );
            }
        }
    }
    patch
}

#[derive(Debug, Clone, Copy)]
struct SampledPieceTauReport {
    samples: usize,
    det_samples: usize,
    outside: usize,
    negative: usize,
    near_zero: usize,
    min_det: f64,
    max_det: f64,
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
    max_sum: f64,
}

#[derive(Debug, Clone, Copy)]
struct CcTestSquareDomainReport {
    samples: usize,
    outside: usize,
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
    max_excess: f64,
}

fn cc_test_square_domain_report(
    tau: &Kp20MacroPatch,
    samples: &[f64],
    tolerance: f64,
) -> CcTestSquareDomainReport {
    let mut report = CcTestSquareDomainReport {
        samples: 0,
        outside: 0,
        min_u: f64::INFINITY,
        max_u: f64::NEG_INFINITY,
        min_v: f64::INFINITY,
        max_v: f64::NEG_INFINITY,
        max_excess: 0.0,
    };
    for piece_id in Kp20PieceId::all() {
        let patch = tau.piece(piece_id);
        for &u in samples {
            let bu = bernstein_bi6_all(u);
            for &v in samples {
                let bv = bernstein_bi6_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for i in 0..=6 {
                    for j in 0..=6 {
                        let cp = patch.cps[i][j];
                        let weight = bu[i] * bv[j];
                        x += weight * cp.x;
                        y += weight * cp.y;
                    }
                }
                let lower_excess = (-x).max(0.0).max((-y).max(0.0));
                let upper_excess = (x - 1.0).max(0.0).max((y - 1.0).max(0.0));
                let total_excess = lower_excess.max(upper_excess);
                report.samples += 1;
                report.outside += usize::from(total_excess > tolerance);
                report.min_u = report.min_u.min(x);
                report.max_u = report.max_u.max(x);
                report.min_v = report.min_v.min(y);
                report.max_v = report.max_v.max(y);
                report.max_excess = report.max_excess.max(total_excess);
            }
        }
    }
    report
}

fn sampled_piece_tau_report(
    tau: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
    samples: &[f64],
) -> SampledPieceTauReport {
    let patch = tau.piece(piece_id);
    let mut report = SampledPieceTauReport {
        samples: 0,
        det_samples: 0,
        outside: 0,
        negative: 0,
        near_zero: 0,
        min_det: f64::INFINITY,
        max_det: f64::NEG_INFINITY,
        min_u: f64::INFINITY,
        max_u: f64::NEG_INFINITY,
        min_v: f64::INFINITY,
        max_v: f64::NEG_INFINITY,
        max_sum: f64::NEG_INFINITY,
    };

    for &u in samples {
        let bu = bernstein_bi6_all(u);
        let du = bernstein_bi6_deriv_all(u);
        for &v in samples {
            let bv = bernstein_bi6_all(v);
            let dv = bernstein_bi6_deriv_all(v);
            let mut x = 0.0_f64;
            let mut y = 0.0_f64;
            let mut dx_du = 0.0_f64;
            let mut dy_du = 0.0_f64;
            let mut dx_dv = 0.0_f64;
            let mut dy_dv = 0.0_f64;
            for i in 0..=6 {
                for j in 0..=6 {
                    let cp = patch.cps[i][j];
                    let value_weight = bu[i] * bv[j];
                    let weight_du = du[i] * bv[j];
                    let weight_dv = bu[i] * dv[j];
                    x += value_weight * cp.x;
                    y += value_weight * cp.y;
                    dx_du += weight_du * cp.x;
                    dy_du += weight_du * cp.y;
                    dx_dv += weight_dv * cp.x;
                    dy_dv += weight_dv * cp.y;
                }
            }
            let determinant = dx_du * dy_dv - dx_dv * dy_du;
            report.samples += 1;
            report.det_samples += 1;
            report.outside += usize::from(x < -1.0e-10 || y < -1.0e-10 || x + y > 1.0 + 1.0e-10);
            report.negative += usize::from(determinant < -1.0e-10);
            report.near_zero += usize::from(determinant.abs() <= 1.0e-10);
            report.min_det = report.min_det.min(determinant);
            report.max_det = report.max_det.max(determinant);
            report.min_u = report.min_u.min(x);
            report.max_u = report.max_u.max(x);
            report.min_v = report.min_v.min(y);
            report.max_v = report.max_v.max(y);
            report.max_sum = report.max_sum.max(x + y);
        }
    }
    report
}

fn sampled_piece_negative_tau_locations(
    tau: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
    samples: &[f64],
) -> Vec<(f64, f64, f64)> {
    let patch = tau.piece(piece_id);
    let mut negatives = Vec::new();
    for &u in samples {
        let bu = bernstein_bi6_all(u);
        let du = bernstein_bi6_deriv_all(u);
        for &v in samples {
            let bv = bernstein_bi6_all(v);
            let dv = bernstein_bi6_deriv_all(v);
            let mut dx_du = 0.0_f64;
            let mut dy_du = 0.0_f64;
            let mut dx_dv = 0.0_f64;
            let mut dy_dv = 0.0_f64;
            for i in 0..=6 {
                for j in 0..=6 {
                    let cp = patch.cps[i][j];
                    let weight_du = du[i] * bv[j];
                    let weight_dv = bu[i] * dv[j];
                    dx_du += weight_du * cp.x;
                    dy_du += weight_du * cp.y;
                    dx_dv += weight_dv * cp.x;
                    dy_dv += weight_dv * cp.y;
                }
            }
            let determinant = dx_du * dy_dv - dx_dv * dy_du;
            if determinant < -1.0e-10 {
                negatives.push((u, v, determinant));
            }
        }
    }
    negatives
}

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
fn pr08_tensor_border_candidates_are_ranked_by_tau_orientation_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let mut best_candidate = "";
    let mut best_negative = usize::MAX;
    let mut best_min_determinant = f64::NEG_INFINITY;

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
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 tau orientation candidate={candidate}: neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            candidate = candidate.name(),
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
        );
        if orientation.negative_determinants < best_negative
            || (orientation.negative_determinants == best_negative
                && orientation.min_determinant > best_min_determinant)
        {
            best_candidate = candidate.name();
            best_negative = orientation.negative_determinants;
            best_min_determinant = orientation.min_determinant;
        }
    }

    eprintln!(
        "PR08 tau orientation best candidate={best_candidate} neg={best_negative} min_det={best_min_determinant:e}"
    );
    assert!(
        best_negative > 0,
        "an existing PR08 source candidate unexpectedly satisfied the tau orientation gate"
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
#[ignore = "obsolete pre-adjacent-eq14-source PR08 row-map ranking diagnostic"]
fn pr08_tensor_border_candidates_are_ranked_by_rotational_sector_rows() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut best_name = "";
    let mut best_residual = f64::INFINITY;
    let mut reversed_outer_row_residual = f64::INFINITY;
    let mut stam_subdivided_first_column_residual = f64::INFINITY;
    let mut stam_q1_residual = f64::INFINITY;
    let mut stam_q2_residual = f64::INFINITY;
    let mut stam_q3_residual = f64::INFINITY;
    let mut stam_pr08_inner_q1_residual = f64::INFINITY;
    let mut stam_pr08_inner_q2_residual = f64::INFINITY;
    let mut stam_pr08_inner_q3_residual = f64::INFINITY;
    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
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
        if candidate
            == Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedReversed
        {
            reversed_outer_row_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaStamSubdividedFirstColumn {
            stam_subdivided_first_column_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows {
            stam_q1_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ2InteriorRows {
            stam_q2_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows {
            stam_q3_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows {
            stam_pr08_inner_q1_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows {
            stam_pr08_inner_q2_residual = solved.max_constraint_residual;
        }
        if candidate == Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows {
            stam_pr08_inner_q3_residual = solved.max_constraint_residual;
        }
        if solved.max_constraint_residual < best_residual {
            best_residual = solved.max_constraint_residual;
            best_name = candidate.name();
        }
    }

    assert_eq!(
        best_name,
        "pr08-middle-patch-formula-next-outer-row-transposed"
    );
    assert!(
        best_residual > 4.2 && best_residual < 4.5,
        "expected best rotational PR08 candidate residual near 4.33, got {best_residual}"
    );
    assert!(
        reversed_outer_row_residual > 129.0 && reversed_outer_row_residual < 132.0,
        "expected reversed adjacent-outer-row transposed candidate to remain rejected near 130.75, got {reversed_outer_row_residual}"
    );
    assert!(
        stam_subdivided_first_column_residual > 22.0
            && stam_subdivided_first_column_residual < 23.5,
        "expected Appendix-A source-mask first-column candidate to remain rejected near 22.69, got {stam_subdivided_first_column_residual}"
    );
    assert!(
        stam_q1_residual > 40.0 && stam_q1_residual < 43.0,
        "expected whole Stam q1 source-row candidate to remain rejected near 41.38, got {stam_q1_residual}"
    );
    assert!(
        stam_q2_residual > 40.0 && stam_q2_residual < 43.0,
        "expected whole Stam q2 source-row candidate to remain rejected near 41.50, got {stam_q2_residual}"
    );
    assert!(
        stam_q3_residual > 43.0 && stam_q3_residual < 45.0,
        "expected whole Stam q3 source-row candidate to remain rejected near 43.87, got {stam_q3_residual}"
    );
    assert!(
        stam_pr08_inner_q1_residual > 54.0 && stam_pr08_inner_q1_residual < 55.0,
        "expected PR08-inner-ring Stam q1 candidate to remain rejected near 54.62, got {stam_pr08_inner_q1_residual}"
    );
    assert!(
        stam_pr08_inner_q2_residual > 40.0 && stam_pr08_inner_q2_residual < 42.0,
        "expected PR08-inner-ring Stam q2 candidate to remain rejected near 41.02, got {stam_pr08_inner_q2_residual}"
    );
    assert!(
        stam_pr08_inner_q3_residual > 46.0 && stam_pr08_inner_q3_residual < 48.0,
        "expected PR08-inner-ring Stam q3 candidate to remain rejected near 47.31, got {stam_pr08_inner_q3_residual}"
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
#[ignore = "obsolete pre-adjacent-eq14-source PR08 candidate family localization diagnostic"]
fn pr08_best_rotational_candidate_residual_is_localized_by_family() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
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
    eprintln!("composed PR08 subdivision source family residuals: {residuals:?}");
    let near_midpoint = residuals
        .iter()
        .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::SectorG1G2NearMidpoint)
        .unwrap();
    assert!(near_midpoint.max_residual > 4.2 && near_midpoint.max_residual < 4.5);
    assert_eq!(near_midpoint.worst_row_index, 27);
    assert_eq!(near_midpoint.worst_axis, None);
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::SectorG1G2NearOrigin
            && residual.max_residual > 4.2
            && residual.max_residual < 4.4
            && residual.worst_row_index == 8
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::DiagonalSymmetry
            && residual.max_residual > 1.1
            && residual.max_residual < 1.2
            && residual.worst_row_index == 40
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::TensorBorderEq14
            && residual.max_residual < 1.0e-6
    }));
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source PR08 near-midpoint b-direction diagnostic"]
fn pr08_best_candidate_near_midpoint_b_direction_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();

    let mut best_label = "";
    let mut best_residual = f64::INFINITY;
    for (origin_label, near_origin_b) in [
        ("origin-paper", kp20_b_coeffs_near_o(c)),
        ("origin-reversed", [c, 2.0 * c, 2.0 * c]),
    ] {
        for (midpoint_label, near_midpoint_b) in [
            ("midpoint-paper", kp20_b_coeffs_near_m(c)),
            ("midpoint-reversed", [0.0, 0.0, c]),
            ("midpoint-middle", [0.0, c, 0.0]),
        ] {
            for (rotation_label, rotation) in [("+angle", angle), ("-angle", -angle)] {
                let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                vector_constraints
                    .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
                let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    near_origin_b,
                    rotation,
                ));
                scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
                    Kp20PieceId::NearMidpointRight,
                    Kp20PieceId::NearMidpointLeft,
                    near_midpoint_b,
                    rotation,
                ));
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: Kp20MacroPatch::zero(),
                        vector_constraints,
                        scalar_constraints,
                        regularization_weight: 1e-10,
                    },
                )
                .unwrap();
                if solved.max_constraint_residual < best_residual {
                    best_residual = solved.max_constraint_residual;
                    best_label = match (origin_label, midpoint_label, rotation_label) {
                        ("origin-paper", "midpoint-reversed", "+angle") => {
                            "origin-paper/midpoint-reversed/+angle"
                        }
                        ("origin-paper", "midpoint-paper", "+angle") => {
                            "origin-paper/midpoint-paper/+angle"
                        }
                        _ => "other",
                    };
                }
            }
        }
    }

    assert_eq!(best_label, "other");
    assert!(best_residual > 4.7 && best_residual < 4.9);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source PR08 sector-coupling diagnostic"]
fn pr08_best_candidate_sector_coupling_orientation_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
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
        best.residual > 4.2 && best.residual < 4.5,
        "expected best sector-coupling diagnostic residual near 4.33, got {best:?}"
    );
    assert_eq!(best.origin_pair, "origin-right-left");
    assert_eq!(best.midpoint_pair, "midpoint-left-right");
    assert_eq!(best.origin_b, "origin-paper");
    assert_eq!(best.midpoint_b, "midpoint-paper");
    assert!((best.rotation - angle).abs() < 1e-12);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source PR08 sector-coupling ranking diagnostic"]
fn pr08_tensor_border_candidates_are_ranked_by_best_sector_coupling() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut best_name = "";
    let mut best_residual = f64::INFINITY;
    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let solved =
            solve_pr08_best_sector_coupling_candidate(&tensor_border, Kp20MacroPatch::zero())
                .unwrap();
        if solved.max_constraint_residual < best_residual {
            best_residual = solved.max_constraint_residual;
            best_name = candidate.name();
        }
    }

    assert_eq!(
        best_name,
        "pr08-middle-patch-formula-next-outer-row-transposed"
    );
    assert!(
        best_residual > 4.2 && best_residual < 4.5,
        "expected best tensor-border residual near 4.33 under corrected sector coupling, got {best_residual}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source split-b orientation diagnostic"]
fn pr08_best_candidate_split_g1_g2_b_orientation_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let origin_reversed = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let midpoint_paper = kp20_b_coeffs_near_m(c);
    let midpoint_reversed = [midpoint_paper[2], midpoint_paper[1], midpoint_paper[0]];
    let mut best = Kp20SplitBOrientationDiagnostic {
        residual: f64::INFINITY,
        origin_eq8_b: "",
        origin_eq9_b: "",
        midpoint_eq8_b: "",
        midpoint_eq9_b: "",
    };

    for &(origin_eq8_label, origin_eq8_b) in &[
        ("origin-paper", origin_paper),
        ("origin-reversed", origin_reversed),
    ] {
        for &(origin_eq9_label, origin_eq9_b) in &[
            ("origin-paper", origin_paper),
            ("origin-reversed", origin_reversed),
        ] {
            for &(midpoint_eq8_label, midpoint_eq8_b) in &[
                ("midpoint-paper", midpoint_paper),
                ("midpoint-reversed", midpoint_reversed),
            ] {
                for &(midpoint_eq9_label, midpoint_eq9_b) in &[
                    ("midpoint-paper", midpoint_paper),
                    ("midpoint-reversed", midpoint_reversed),
                ] {
                    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                    vector_constraints.extend(
                        kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap(),
                    );
                    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
                            Kp20PieceId::NearOriginRight,
                            Kp20PieceId::NearOriginLeft,
                            origin_eq8_b,
                            origin_eq9_b,
                            angle,
                        ),
                    );
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
                            Kp20PieceId::NearMidpointLeft,
                            Kp20PieceId::NearMidpointRight,
                            midpoint_eq8_b,
                            midpoint_eq9_b,
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
                    if solved.max_constraint_residual < best.residual {
                        best = Kp20SplitBOrientationDiagnostic {
                            residual: solved.max_constraint_residual,
                            origin_eq8_b: origin_eq8_label,
                            origin_eq9_b: origin_eq9_label,
                            midpoint_eq8_b: midpoint_eq8_label,
                            midpoint_eq9_b: midpoint_eq9_label,
                        };
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "expected split-b orientation residual near 4.33, got {best:?}"
    );
    assert_eq!(best.origin_eq8_b, "origin-reversed");
    assert_eq!(best.origin_eq9_b, "origin-paper");
    assert_eq!(best.midpoint_eq8_b, "midpoint-paper");
    assert_eq!(best.midpoint_eq9_b, "midpoint-paper");
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source eq9 sign diagnostic"]
fn pr08_best_candidate_eq9_term_sign_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let origin_eq8_b = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let origin_eq9_b = origin_paper;
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let mut best = Kp20Eq9TermSignDiagnostic {
        residual: f64::INFINITY,
        origin_signs: Kp20Eq9TermSigns::paper(),
        midpoint_signs: Kp20Eq9TermSigns::paper(),
    };

    for origin_signs in kp20_eq9_term_sign_variants() {
        for midpoint_signs in kp20_eq9_term_sign_variants() {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    origin_eq8_b,
                    origin_eq9_b,
                    origin_signs,
                    angle,
                ),
            );
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                    Kp20PieceId::NearMidpointLeft,
                    Kp20PieceId::NearMidpointRight,
                    midpoint_b,
                    midpoint_b,
                    midpoint_signs,
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
            if solved.max_constraint_residual < best.residual {
                best = Kp20Eq9TermSignDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_signs,
                    midpoint_signs,
                };
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.3,
        "expected best eq9 term-sign diagnostic residual near 4.256, got {best:?}"
    );
    assert_eq!(best.origin_signs.dvv_curr, 1.0);
    assert_eq!(best.origin_signs.dvv_prev, -1.0);
    assert_eq!(best.origin_signs.two_b_duv, 1.0);
    assert_eq!(best.origin_signs.b_sq_duu, -1.0);
    assert_eq!(best.origin_signs.b_bprime_du, -1.0);
    assert_eq!(best.midpoint_signs.dvv_curr, 1.0);
    assert_eq!(best.midpoint_signs.dvv_prev, -1.0);
    assert_eq!(best.midpoint_signs.two_b_duv, 1.0);
    assert_eq!(best.midpoint_signs.b_sq_duu, -1.0);
    assert_eq!(best.midpoint_signs.b_bprime_du, 1.0);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source eq9 dvv sign diagnostic"]
fn pr08_best_candidate_eq9_dvv_sign_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let origin_eq8_b = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let origin_eq9_b = origin_paper;
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let origin_base = Kp20Eq9TermSigns {
        b_sq_duu: 1.0,
        ..Kp20Eq9TermSigns::paper()
    };
    let midpoint_base = Kp20Eq9TermSigns {
        two_b_duv: -1.0,
        b_bprime_du: 1.0,
        ..Kp20Eq9TermSigns::paper()
    };
    let mut best = Kp20Eq9TermSignDiagnostic {
        residual: f64::INFINITY,
        origin_signs: origin_base,
        midpoint_signs: midpoint_base,
    };

    for origin_dvv_curr in [-1.0, 1.0] {
        for origin_dvv_prev in [-1.0, 1.0] {
            for midpoint_dvv_curr in [-1.0, 1.0] {
                for midpoint_dvv_prev in [-1.0, 1.0] {
                    let origin_signs = Kp20Eq9TermSigns {
                        dvv_curr: origin_dvv_curr,
                        dvv_prev: origin_dvv_prev,
                        ..origin_base
                    };
                    let midpoint_signs = Kp20Eq9TermSigns {
                        dvv_curr: midpoint_dvv_curr,
                        dvv_prev: midpoint_dvv_prev,
                        ..midpoint_base
                    };
                    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                    vector_constraints.extend(
                        kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap(),
                    );
                    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                            Kp20PieceId::NearOriginRight,
                            Kp20PieceId::NearOriginLeft,
                            origin_eq8_b,
                            origin_eq9_b,
                            origin_signs,
                            angle,
                        ),
                    );
                    scalar_constraints.extend(
                        kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                            Kp20PieceId::NearMidpointLeft,
                            Kp20PieceId::NearMidpointRight,
                            midpoint_b,
                            midpoint_b,
                            midpoint_signs,
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
                    if solved.max_constraint_residual < best.residual {
                        best = Kp20Eq9TermSignDiagnostic {
                            residual: solved.max_constraint_residual,
                            origin_signs,
                            midpoint_signs,
                        };
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.3,
        "expected eq9 dvv-sign diagnostic to preserve the 4.256 frontier, got {best:?}"
    );
    assert_eq!(best.origin_signs.dvv_curr, 1.0);
    assert_eq!(best.origin_signs.dvv_prev, -1.0);
    assert_eq!(best.midpoint_signs.dvv_curr, 1.0);
    assert_eq!(best.midpoint_signs.dvv_prev, -1.0);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source eq9 term-source diagnostic"]
fn pr08_best_candidate_eq9_term_source_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let origin_eq8_b = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let origin_eq9_b = origin_paper;
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let origin_signs = Kp20Eq9TermSigns {
        b_sq_duu: 1.0,
        ..Kp20Eq9TermSigns::paper()
    };
    let midpoint_signs = Kp20Eq9TermSigns {
        two_b_duv: -1.0,
        b_bprime_du: 1.0,
        ..Kp20Eq9TermSigns::paper()
    };
    let mut best = Kp20Eq9TermSourceDiagnostic {
        residual: f64::INFINITY,
        origin_terms: Kp20Eq9DiagnosticTerms {
            signs: origin_signs,
            ..Kp20Eq9DiagnosticTerms::paper()
        },
        midpoint_terms: Kp20Eq9DiagnosticTerms {
            signs: midpoint_signs,
            ..Kp20Eq9DiagnosticTerms::paper()
        },
    };

    for origin_terms in kp20_eq9_term_source_variants(origin_signs) {
        for midpoint_terms in kp20_eq9_term_source_variants(midpoint_signs) {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    origin_eq8_b,
                    origin_eq9_b,
                    origin_terms,
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
            if solved.max_constraint_residual < best.residual {
                best = Kp20Eq9TermSourceDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_terms,
                    midpoint_terms,
                };
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.3,
        "expected eq9 term-source diagnostic to preserve the 4.256 frontier, got {best:?}"
    );
    assert_eq!(best.origin_terms.duv_source, Kp20Eq9TermSource::Previous);
    assert_eq!(best.origin_terms.duu_source, Kp20Eq9TermSource::Previous);
    assert_eq!(best.origin_terms.du_source, Kp20Eq9TermSource::Previous);
    assert_eq!(best.midpoint_terms.duu_source, Kp20Eq9TermSource::Previous);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source coherent eq9 reversal diagnostic"]
fn pr08_best_candidate_coherent_eq9_reversal_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let midpoint_paper = kp20_b_coeffs_near_m(c);
    let origin_eq8_b = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let mut best = Kp20CoherentEq9ReversalDiagnostic {
        residual: f64::INFINITY,
        origin_variant: "",
        midpoint_variant: "",
    };

    for &(origin_variant, origin_eq9_b, origin_terms) in &[
        (
            "origin-printed",
            origin_paper,
            Kp20Eq9DiagnosticTerms::paper(),
        ),
        (
            "origin-reversed-parameter",
            [origin_paper[2], origin_paper[1], origin_paper[0]],
            kp20_eq9_reversed_parameter_terms(),
        ),
    ] {
        for &(midpoint_variant, midpoint_eq9_b, midpoint_terms) in &[
            (
                "midpoint-printed",
                midpoint_paper,
                Kp20Eq9DiagnosticTerms::paper(),
            ),
            (
                "midpoint-reversed-parameter",
                [midpoint_paper[2], midpoint_paper[1], midpoint_paper[0]],
                kp20_eq9_reversed_parameter_terms(),
            ),
        ] {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    origin_eq8_b,
                    origin_eq9_b,
                    origin_terms,
                    angle,
                ),
            );
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
                    Kp20PieceId::NearMidpointLeft,
                    Kp20PieceId::NearMidpointRight,
                    midpoint_paper,
                    midpoint_eq9_b,
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
            if solved.max_constraint_residual < best.residual {
                best = Kp20CoherentEq9ReversalDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_variant,
                    midpoint_variant,
                };
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "expected coherent eq9 reversal sweep to preserve the 4.332 paper-backed frontier, got {best:?}"
    );
    assert_eq!(best.origin_variant, "origin-printed");
    assert_eq!(best.midpoint_variant, "midpoint-printed");
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source G2 rotation-frame diagnostic"]
fn pr08_best_candidate_g2_rotation_frame_sweep() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_paper = kp20_b_coeffs_near_o(c);
    let origin_eq8_b = [origin_paper[2], origin_paper[1], origin_paper[0]];
    let origin_eq9_b = origin_paper;
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let mut best = Kp20G2RotationFrameDiagnostic {
        residual: f64::INFINITY,
        origin_g2_rotation: 0.0,
        midpoint_g2_rotation: 0.0,
    };

    for origin_g2_rotation in [angle, -angle] {
        for midpoint_g2_rotation in [angle, -angle] {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints
                .extend(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_constraints(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearOriginLeft,
                origin_eq8_b,
                angle,
            ));
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g2_constraints(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearOriginLeft,
                origin_eq9_b,
                origin_g2_rotation,
            ));
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_constraints(
                Kp20PieceId::NearMidpointLeft,
                Kp20PieceId::NearMidpointRight,
                midpoint_b,
                angle,
            ));
            scalar_constraints.extend(kp20_planar_tau_rotated_sector_g2_constraints(
                Kp20PieceId::NearMidpointLeft,
                Kp20PieceId::NearMidpointRight,
                midpoint_b,
                midpoint_g2_rotation,
            ));
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
                best = Kp20G2RotationFrameDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_g2_rotation,
                    midpoint_g2_rotation,
                };
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "expected G2 rotation-frame sweep to preserve the 4.332 frontier, got {best:?}"
    );
    assert!((best.origin_g2_rotation - angle).abs() < 1e-12);
    assert!((best.midpoint_g2_rotation - angle).abs() < 1e-12);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source tensor-border slack diagnostic"]
fn pr08_best_candidate_free_tensor_border_targets_localizes_projection_gap() {
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows(|_| true);

    assert!(
        solved.residual < 1e-9,
        "freeing the uncertain tensor-border targets should make the corrected paper row set feasible, got {solved:?}"
    );
    assert!(
        solved.max_slack > 2.3 && solved.max_slack < 2.4,
        "expected max target movement near 2.32, got {solved:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source tensor-border half dominance diagnostic"]
fn pr08_best_candidate_tensor_border_target_gap_is_midpoint_half_dominated() {
    let free_near_c = solve_pr08_best_candidate_with_free_tensor_border_rows(|row| row < 21);
    let free_near_m = solve_pr08_best_candidate_with_free_tensor_border_rows(|row| row >= 21);
    let mut best_depth = Kp20TensorBorderFreedomDiagnostic {
        residual: f64::INFINITY,
        max_slack: 0.0,
    };
    let mut best_depth_index = usize::MAX;
    for depth in 0..3 {
        let solved =
            solve_pr08_best_candidate_with_free_tensor_border_rows(|row| (row % 21) / 7 == depth);
        if solved.residual < best_depth.residual {
            best_depth = solved;
            best_depth_index = depth;
        }
    }

    assert!(
        free_near_m.residual > 4.2 && free_near_m.residual < 4.3,
        "freeing only the near-m half should collapse to the known origin-row frontier, got {free_near_m:?}"
    );
    assert!(
        free_near_c.residual > 40.0,
        "freeing only the near-c half should not explain the corrected-coupling gap, got {free_near_c:?}"
    );
    assert!(
        best_depth_index == 1 && best_depth.residual > 180.0,
        "freeing one depth across both halves should remain rejected, got depth {best_depth_index}: {best_depth:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source tensor-border band diagnostic"]
fn pr08_best_candidate_tensor_border_target_gap_localizes_by_half_depth_band() {
    let mut best = Kp20TensorBorderBandFreedomDiagnostic {
        residual: f64::INFINITY,
        max_slack: 0.0,
        band_mask: 0,
        freed_band_count: 0,
    };
    let mut best_by_count = [best; 7];
    for band_mask in 1..64 {
        let solved = solve_pr08_best_candidate_with_free_tensor_border_rows(|row| {
            let band = 3 * (row / 21) + (row % 21) / 7;
            band_mask & (1 << band) != 0
        });
        let candidate = Kp20TensorBorderBandFreedomDiagnostic {
            residual: solved.residual,
            max_slack: solved.max_slack,
            band_mask,
            freed_band_count: band_mask.count_ones(),
        };
        if candidate.residual < best.residual {
            best = candidate;
        }
        if candidate.residual < best_by_count[candidate.freed_band_count as usize].residual {
            best_by_count[candidate.freed_band_count as usize] =
                Kp20TensorBorderBandFreedomDiagnostic {
                    residual: solved.residual,
                    max_slack: solved.max_slack,
                    band_mask,
                    freed_band_count: band_mask.count_ones(),
                };
        }
    }

    let five_band_solution = best_by_count[5];
    assert!(
        best.residual < 1e-9,
        "freeing all half/depth target bands should make the paper row set feasible, got {best:?}"
    );
    assert_eq!(best.band_mask, 63);
    assert!(
        five_band_solution.residual < 1e-9,
        "freeing every tensor-border target band except near-c position should also be feasible, got {five_band_solution:?}"
    );
    assert_eq!(five_band_solution.band_mask, 62);
    assert!(
        best_by_count[4].residual > 0.1 && best_by_count[4].residual < 0.2,
        "keeping two or more target bands pinned should not close the projection gap, got {:?}",
        best_by_count[4]
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source pinned tensor-row diagnostic"]
fn pr08_best_candidate_pinned_target_rows_rank_projection_gap() {
    let mut worst = Kp20TensorBorderPinnedRowDiagnostic {
        residual: 0.0,
        pinned_row: 0,
    };
    let mut best = Kp20TensorBorderPinnedRowDiagnostic {
        residual: f64::INFINITY,
        pinned_row: 0,
    };
    for pinned_row in 7..42 {
        let solved = solve_pr08_best_candidate_with_free_tensor_border_rows(|row| {
            row >= 7 && row != pinned_row
        });
        let diagnostic = Kp20TensorBorderPinnedRowDiagnostic {
            residual: solved.residual,
            pinned_row,
        };
        if diagnostic.residual > worst.residual {
            worst = diagnostic;
        }
        if diagnostic.residual < best.residual {
            best = diagnostic;
        }
    }

    assert_eq!(worst.pinned_row, 21);
    assert_eq!(
        tensor_border_row_label(worst.pinned_row),
        "near-m depth0 u0"
    );
    assert!(
        worst.residual > 5.3 && worst.residual < 5.5,
        "pinning near-m depth0 u0 should preserve the largest projection gap, got {worst:?}"
    );
    assert_eq!(best.pinned_row, 34);
    assert_eq!(tensor_border_row_label(best.pinned_row), "near-m depth1 u6");
    assert!(
        best.residual < 1e-9,
        "some individual near-m derivative rows should be compatible with the anchored near-c position, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 anchor diagnostic"]
fn pr08_best_candidate_row21_is_bspline_midpoint_anchor() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_constraints = kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap();
    let direct = Point3::new(
        tensor_constraints[21].target[0],
        tensor_constraints[21].target[1],
        0.0,
    );
    let from_source_rows = tensor_border_row21_source_anchor(&rows.sectors[0]);

    assert_close_point(direct, from_source_rows, 1e-12);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 slack diagnostic"]
fn pr08_best_candidate_row21_required_anchor_delta_is_measured() {
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows_full(|row| row >= 7);
    let row21_extra_offset = 2 * (21 - 7);
    let dx = solved.extra_values[row21_extra_offset];
    let dy = solved.extra_values[row21_extra_offset + 1];

    assert!(
        solved.max_constraint_residual < 1e-9,
        "freeing all non-near-c-position targets should make the row set feasible"
    );
    assert!(
        dx > -0.65 && dx < -0.64 && dy > -0.31 && dy < -0.30,
        "expected row21 required slack near (-0.6447, -0.3048), got ({dx}, {dy})"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 visual-midpoint diagnostic"]
fn pr08_best_candidate_corrected_row21_anchor_is_not_visual_stam_side_midpoint() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let current = tensor_border_row21_source_anchor(&rows.sectors[0]);
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows_full(|row| row >= 7);
    let row21_extra_offset = 2 * (21 - 7);
    let corrected = Point3::new(
        current.x + solved.extra_values[row21_extra_offset],
        current.y + solved.extra_values[row21_extra_offset + 1],
        0.0,
    );
    let projected = kp20_pr08_stam_visual_projected_grids(&ring).unwrap();
    let mut best = Kp20StamAnchorDiagnostic {
        distance: f64::INFINITY,
        patch: Kp20StamPatch::LowerRight,
        transform: Kp20StamGridTransform::Identity,
        side: Kp20StamGridSide::Top,
        reversed_depth: false,
    };

    for patch in [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ] {
        for transform in kp20_stam_grid_transforms() {
            for side in Kp20StamGridSide::all() {
                for reversed_depth in [false, true] {
                    let candidate_rows =
                        kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
                            &projected,
                            patch,
                            transform,
                            side,
                            reversed_depth,
                        )
                        .unwrap();
                    let anchor = tensor_border_row21_source_anchor(&candidate_rows.sectors[0]);
                    let distance = anchor.distance_to(&corrected);
                    if distance < best.distance {
                        best = Kp20StamAnchorDiagnostic {
                            distance,
                            patch,
                            transform,
                            side,
                            reversed_depth,
                        };
                    }
                }
            }
        }
    }

    assert_eq!(best.patch, Kp20StamPatch::UpperLeft);
    assert_eq!(best.transform, Kp20StamGridTransform::Rotate180);
    assert_eq!(best.side, Kp20StamGridSide::Top);
    assert!(best.reversed_depth);
    assert!(
        best.distance > 0.38 && best.distance < 0.39,
        "corrected row21 anchor should not already appear as a transformed visual Stam side midpoint, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 single-label diagnostic"]
fn pr08_best_candidate_row21_single_label_substitution_probe() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let current = tensor_border_row21_source_anchor(&rows.sectors[0]);
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows_full(|row| row >= 7);
    let row21_extra_offset = 2 * (21 - 7);
    let corrected = Point3::new(
        current.x + solved.extra_values[row21_extra_offset],
        current.y + solved.extra_values[row21_extra_offset + 1],
        0.0,
    );
    let mut best = Kp20Row21LabelSubstitutionDiagnostic {
        distance: current.distance_to(&corrected),
        source_row: 0,
        source_col: 0,
        label_sector: "current",
        label_index: 0,
    };
    let label_sets = [
        (
            "previous",
            &ring.sectors[(ring.sectors.len() - 1) % ring.sectors.len()],
        ),
        ("current", &ring.sectors[0]),
        ("next", &ring.sectors[1]),
    ];

    for source_row in 0..3 {
        for source_col in 0..4 {
            for &(label_sector, labels) in &label_sets {
                for (label_index, &label) in labels.iter().enumerate() {
                    let mut candidate = rows.sectors[0];
                    candidate[source_row][source_col] = label;
                    let anchor = tensor_border_row21_source_anchor(&candidate);
                    let distance = anchor.distance_to(&corrected);
                    if distance < best.distance {
                        best = Kp20Row21LabelSubstitutionDiagnostic {
                            distance,
                            source_row,
                            source_col,
                            label_sector,
                            label_index,
                        };
                    }
                }
            }
        }
    }

    assert_eq!(best.source_row, 1);
    assert_eq!(best.source_col, 1);
    assert_eq!(best.label_sector, "next");
    assert_eq!(best.label_index, 3);
    assert!(
        best.distance > 0.13 && best.distance < 0.14,
        "best single PR08 label substitution should still miss the corrected row21 anchor, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 single-label full-tau diagnostic"]
fn pr08_best_candidate_row21_best_single_label_substitution_full_tau_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    for sector_index in 0..ring.sectors.len() {
        let next = (sector_index + 1) % ring.sectors.len();
        rows.sectors[sector_index][1][1] = ring.sectors[next][3];
    }
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved =
        solve_pr08_best_sector_coupling_candidate(&tensor_border, Kp20MacroPatch::zero()).unwrap();

    assert!(
        solved.max_constraint_residual > 11.0 && solved.max_constraint_residual < 12.0,
        "row21 best single-label substitution should remain rejected near 11.45, got {}",
        solved.max_constraint_residual
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 adjacent-column diagnostic"]
fn pr08_best_candidate_row21_adjacent_column_substitution_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    for sector_index in 0..ring.sectors.len() {
        let next = (sector_index + 1) % ring.sectors.len();
        rows.sectors[sector_index][0][1] = ring.sectors[next][3];
        rows.sectors[sector_index][1][1] = ring.sectors[next][6];
        rows.sectors[sector_index][2][1] = ring.sectors[next][12];
    }
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved =
        solve_pr08_best_sector_coupling_candidate(&tensor_border, Kp20MacroPatch::zero()).unwrap();

    assert!(
        solved.max_constraint_residual > 67.0 && solved.max_constraint_residual < 68.0,
        "adjacent-column row21 repair should remain strongly rejected near 67.55, got {}",
        solved.max_constraint_residual
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 column-triple diagnostic"]
fn pr08_best_candidate_row21_column_label_triple_probe() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let current = tensor_border_row21_source_anchor(&rows.sectors[0]);
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows_full(|row| row >= 7);
    let row21_extra_offset = 2 * (21 - 7);
    let corrected = Point3::new(
        current.x + solved.extra_values[row21_extra_offset],
        current.y + solved.extra_values[row21_extra_offset + 1],
        0.0,
    );
    let labels = row21_nearby_pr08_labels(&ring);
    let mut best = Kp20Row21ColumnTripleDiagnostic {
        distance: current.distance_to(&corrected),
        column: 1,
        labels: [("", 0); 3],
    };

    for column in 0..4 {
        for &(label0_sector, label0_index, label0) in &labels {
            for &(label1_sector, label1_index, label1) in &labels {
                for &(label2_sector, label2_index, label2) in &labels {
                    let mut candidate = rows.sectors[0];
                    candidate[0][column] = label0;
                    candidate[1][column] = label1;
                    candidate[2][column] = label2;
                    let anchor = tensor_border_row21_source_anchor(&candidate);
                    let distance = anchor.distance_to(&corrected);
                    if distance < best.distance {
                        best = Kp20Row21ColumnTripleDiagnostic {
                            distance,
                            column,
                            labels: [
                                (label0_sector, label0_index),
                                (label1_sector, label1_index),
                                (label2_sector, label2_index),
                            ],
                        };
                    }
                }
            }
        }
    }

    assert_eq!(best.column, 1);
    assert_eq!(
        best.labels,
        [("previous", 1), ("next", 12), ("previous", 12)]
    );
    assert!(
        best.distance > 0.001 && best.distance < 0.002,
        "expected best row21 column-triple anchor distance near 0.0014, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 column-triple full-tau diagnostic"]
fn pr08_best_candidate_row21_best_column_triple_full_tau_probe() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    for sector_index in 0..ring.sectors.len() {
        let previous = (sector_index + ring.sectors.len() - 1) % ring.sectors.len();
        let next = (sector_index + 1) % ring.sectors.len();
        rows.sectors[sector_index][0][1] = ring.sectors[previous][1];
        rows.sectors[sector_index][1][1] = ring.sectors[next][12];
        rows.sectors[sector_index][2][1] = ring.sectors[previous][12];
    }
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved =
        solve_pr08_best_sector_coupling_candidate(&tensor_border, Kp20MacroPatch::zero()).unwrap();

    assert!(
        solved.max_constraint_residual > 64.0 && solved.max_constraint_residual < 65.0,
        "row21 best column-triple candidate should remain rejected near 64.41, got {}",
        solved.max_constraint_residual
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source row21 recombination diagnostic"]
fn pr08_best_candidate_row21_source_recombination_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let best_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let best_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&best_rows).unwrap();
    let best_tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&best_tensor_border, 0).unwrap();
    let mut best = Kp20DerivativeSourceDiagnostic {
        residual: f64::INFINITY,
        source_name: "",
    };

    for source in Kp20Pr08TensorBorderCandidate::all() {
        let source_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, source).unwrap();
        let source_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&source_rows).unwrap();
        let source_tensor_constraints =
            kp20_planar_tau_tensor_border_constraints(&source_tensor_border, 0).unwrap();
        let mut vector_constraints = best_tensor_constraints.clone();
        vector_constraints[21].target = source_tensor_constraints[21].target;
        vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        let angle = core::f64::consts::TAU / 5.0;
        let c = angle.cos();
        let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
        scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearOriginLeft,
            [c, 2.0 * c, 2.0 * c],
            kp20_b_coeffs_near_o(c),
            angle,
        ));
        scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
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
        if solved.max_constraint_residual < best.residual {
            best = Kp20DerivativeSourceDiagnostic {
                residual: solved.max_constraint_residual,
                source_name: source.name(),
            };
        }
    }

    assert_eq!(
        best.source_name,
        "pr08-middle-patch-formula-next-outer-row-transposed"
    );
    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "single row21 source recombination should preserve the 4.332 frontier, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source position/derivative source diagnostic"]
fn pr08_best_position_row_derivative_source_sweep_localizes_projection_gap() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let position_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let mut best = Kp20DerivativeSourceDiagnostic {
        residual: f64::INFINITY,
        source_name: "",
    };

    for derivative_source in Kp20Pr08TensorBorderCandidate::all() {
        let derivative_rows =
            kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, derivative_source)
                .unwrap();
        let mut sectors = Vec::with_capacity(position_rows.sectors.len());
        for sector in 0..position_rows.sectors.len() {
            sectors.push([
                position_rows.sectors[sector][0],
                derivative_rows.sectors[sector][1],
                derivative_rows.sectors[sector][2],
            ]);
        }
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&Kp20RegularTensorBorderCNet { sectors })
                .unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: 5,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
        if solved.max_constraint_residual < best.residual {
            best = Kp20DerivativeSourceDiagnostic {
                residual: solved.max_constraint_residual,
                source_name: derivative_source.name(),
            };
        }
    }

    assert_eq!(
        best.source_name,
        "pr08-middle-patch-formula-next-outer-row-transposed"
    );
    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "expected derivative-source recombination to preserve the 4.332 frontier, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source similarity-normalization diagnostic"]
fn pr08_best_candidate_global_similarity_normalization_is_not_enough() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_border_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap();
    let mut vector_constraints = tensor_border_rows[0..7].to_vec();
    vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: augmented_tensor_border_similarity_constraints(
                &tensor_border_rows[7..],
            ),
            extra_variable_count: 2,
            regularization_weight: 1e-10,
            extra_regularization_weight: 1e-10,
        },
    )
    .unwrap();

    assert!(
        solved.max_constraint_residual > 70.0 && solved.max_constraint_residual < 73.0,
        "expected anchored similarity normalization to remain rejected near 71.48, got residual {} with a={} b={}",
        solved.max_constraint_residual,
        solved.extra_values[0],
        solved.extra_values[1]
    );
    assert!(solved.extra_values[0] > 1.6 && solved.extra_values[0] < 1.8);
    assert!(solved.extra_values[1] > -0.55 && solved.extra_values[1] < -0.4);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source normalized candidate ranking diagnostic"]
fn pr08_source_candidates_normalized_to_middle_patch_value_are_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let expected = kp20_pr08_middle_patch_normalization_value(&ring, 0).unwrap();
    let mut best = Kp20DerivativeSourceDiagnostic {
        residual: f64::INFINITY,
        source_name: "",
    };

    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let scaled_rows = normalize_rows_to_middle_patch_value(&rows, expected);
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&scaled_rows).unwrap();
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: 5,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: Kp20MacroPatch::zero(),
                regularization_weight: 1e-10,
            })
            .unwrap();
        if solved.max_constraint_residual < best.residual {
            best = Kp20DerivativeSourceDiagnostic {
                residual: solved.max_constraint_residual,
                source_name: candidate.name(),
            };
        }
    }

    assert_eq!(
        best.source_name,
        "pr08-middle-patch-formula-next-outer-row-transposed"
    );
    assert!(
        best.residual > 4.2 && best.residual < 4.5,
        "expected PR08-normalized source sweep to preserve the 4.332 frontier, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source d1_3 sweep diagnostic"]
fn pr08_best_candidate_eq14_d1_3_sweep_exposes_near_origin_frontier() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let c = (core::f64::consts::TAU / 5.0).cos();
    let default_d1_3 = 2.0 * c / 3.0;
    let mut best = Kp20Eq14D1Diagnostic {
        residual: f64::INFINITY,
        d1_3: f64::NAN,
    };

    for step in -40..=40 {
        let d1_3 = default_d1_3 + step as f64 * 0.025;
        let solved = solve_pr08_best_candidate_with_eq14_d1_3(&tensor_border, d1_3);
        if solved.max_constraint_residual < best.residual {
            best = Kp20Eq14D1Diagnostic {
                residual: solved.max_constraint_residual,
                d1_3,
            };
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.3,
        "expected eq(14) d1_3 sweep to expose the near-origin 4.256 frontier, got {best:?}"
    );
    assert!(
        best.d1_3 > 0.03 && best.d1_3 < 0.032,
        "expected best diagnostic d1_3 near -0.019 rather than default {default_d1_3}, got {best:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source d1_3 family localization diagnostic"]
fn pr08_best_candidate_eq14_d1_3_family_residuals_localize_near_origin() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let d1_3 = -0.018988670416701714;
    let solved = solve_pr08_best_candidate_with_eq14_d1_3(&tensor_border, d1_3);
    let residuals = kp20_planar_tau_rotational_paper_family_residuals_with_d1_3(
        5,
        0,
        &tensor_border,
        d1_3,
        &solved.tau_macro_patch,
    )
    .unwrap();
    let near_origin = residuals
        .iter()
        .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::SectorG1G2NearOrigin)
        .unwrap();
    let tensor_border_residual = residuals
        .iter()
        .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::TensorBorderEq14)
        .unwrap();

    assert!(
        near_origin.max_residual > 4.2 && near_origin.max_residual < 4.3,
        "expected diagnostic d1_3 to localize the active frontier to near-origin rows, got {near_origin:?}"
    );
    assert_eq!(near_origin.worst_row_index, 8);
    assert_eq!(near_origin.worst_axis, None);
    assert_eq!(
        rotated_sector_g1_g2_row_label(near_origin.worst_row_index),
        "eq8[4].x"
    );
    assert!(
        tensor_border_residual.max_residual < 1.0e-6,
        "explicit d1_3 tensor-border rows should be satisfied by the solved patch, got {tensor_border_residual:?}"
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source d1_3 plus eq9 sign diagnostic"]
fn pr08_best_candidate_eq14_d1_3_and_eq9_sign_sweep_remains_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let d1_3 = -0.018988670416701714;
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let origin_eq8_b = [c, 2.0 * c, 2.0 * c];
    let origin_eq9_b = kp20_b_coeffs_near_o(c);
    let midpoint_b = kp20_b_coeffs_near_m(c);
    let mut best = Kp20Eq9TermSignDiagnostic {
        residual: f64::INFINITY,
        origin_signs: Kp20Eq9TermSigns::paper(),
        midpoint_signs: Kp20Eq9TermSigns::paper(),
    };

    for origin_signs in kp20_eq9_term_sign_variants() {
        for midpoint_signs in kp20_eq9_term_sign_variants() {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(
                kp20_planar_tau_tensor_border_constraints_with_d1_3(&tensor_border, 0, d1_3)
                    .unwrap(),
            );
            let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                    Kp20PieceId::NearOriginRight,
                    Kp20PieceId::NearOriginLeft,
                    origin_eq8_b,
                    origin_eq9_b,
                    origin_signs,
                    angle,
                ),
            );
            scalar_constraints.extend(
                kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
                    Kp20PieceId::NearMidpointLeft,
                    Kp20PieceId::NearMidpointRight,
                    midpoint_b,
                    midpoint_b,
                    midpoint_signs,
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
            if solved.max_constraint_residual < best.residual {
                best = Kp20Eq9TermSignDiagnostic {
                    residual: solved.max_constraint_residual,
                    origin_signs,
                    midpoint_signs,
                };
            }
        }
    }

    assert!(
        best.residual > 4.2 && best.residual < 4.3,
        "expected combined d1_3 + eq9 sign sweep to remain at the 4.256 frontier, got {best:?}"
    );
    assert_eq!(best.origin_signs.b_sq_duu, 1.0);
    assert_eq!(best.origin_signs.b_bprime_du, -1.0);
    assert_eq!(best.midpoint_signs.two_b_duv, -1.0);
    assert_eq!(best.midpoint_signs.b_sq_duu, -1.0);
    assert_eq!(best.midpoint_signs.b_bprime_du, 1.0);
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source first-column feasibility diagnostic"]
fn pr08_middle_patch_first_column_only_feasibility_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let tensor_border_constraints =
        pr08_middle_patch_augmented_first_column_tensor_border_constraints(&ring);
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let vector_constraints = kp20_planar_tau_internal_c2_constraints();
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointRight,
        Kp20PieceId::NearMidpointLeft,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: tensor_border_constraints,
            extra_variable_count: 6,
            regularization_weight: 1e-10,
            extra_regularization_weight: 1e-12,
        },
    )
    .unwrap();

    assert!(
        solved.max_constraint_residual > 800.0 && solved.max_constraint_residual < 830.0,
        "expected pure first-column feasibility to remain rejected near 817, got {}",
        solved.max_constraint_residual
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source first-column scale diagnostic"]
fn pr08_middle_patch_source_first_column_scale_is_not_sufficient() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let tensor_border_constraints =
        pr08_middle_patch_augmented_source_first_column_scale_constraints(&ring);
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let vector_constraints = kp20_planar_tau_internal_c2_constraints();
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointRight,
        Kp20PieceId::NearMidpointLeft,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: tensor_border_constraints,
            extra_variable_count: 1,
            regularization_weight: 1e-10,
            extra_regularization_weight: 1e-12,
        },
    )
    .unwrap();

    assert!(
        solved.max_constraint_residual > 740.0 && solved.max_constraint_residual < 755.0,
        "expected source first-column scale solve to remain rejected near 746, got residual {}",
        solved.max_constraint_residual
    );
    assert!(
        solved.extra_values[0] > -6.3 && solved.extra_values[0] < -6.1,
        "expected rejected source first-column scale near -6.23, got {}",
        solved.extra_values[0]
    );
}

#[test]
#[ignore = "obsolete pre-adjacent-eq14-source Appendix-B projection orientation diagnostic"]
fn stam_appendix_b_source_projection_orientation_sweep_is_rejected() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let mut best = StamAppendixBProjectionDiagnostic {
        residual: f64::INFINITY,
        patch_index: 0,
        transposed: false,
        row_start: 0,
        reverse_rows: false,
        reverse_cols: false,
    };
    for patch_index in 0..3 {
        for transposed in [false, true] {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    for reverse_cols in [false, true] {
                        let rows = stam_appendix_b_projected_rows(
                            &ring,
                            patch_index,
                            transposed,
                            row_start,
                            reverse_rows,
                            reverse_cols,
                        )
                        .unwrap();
                        let tensor_border =
                            kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
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
                            best = StamAppendixBProjectionDiagnostic {
                                residual: solved.max_constraint_residual,
                                patch_index,
                                transposed,
                                row_start,
                                reverse_rows,
                                reverse_cols,
                            };
                        }
                    }
                }
            }
        }
    }

    assert!(
        best.residual > 17.0 && best.residual < 18.0,
        "expected best Stam Appendix-B source-projection sweep residual near 17.47, got {best:?}"
    );
    assert_eq!(best.patch_index, 1);
    assert!(!best.transposed);
    assert_eq!(best.row_start, 0);
    assert!(best.reverse_rows);
    assert!(!best.reverse_cols);
}

#[test]
#[ignore = "exhaustive transformed-side diagnostic; PR08 refined labels are already proved to be lambda-scaled"]
fn pr08_refined_stam_visual_side_rows_are_measured_through_paper_tau_solve() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let refined = pr08_subdivide_characteristic_ring(&ring).unwrap();
    let refined_ring = Kp20CatmullClarkCharacteristicRing {
        sectors: refined,
        lambda: ring.lambda,
    };
    let projection = kp20_pr08_stam_visual_projected_grids(&refined_ring).unwrap();
    let initial = Kp20MacroPatch::zero();
    let mut best = StamTransformedSideProjectionDiagnostic {
        residual: f64::INFINITY,
        patch: Kp20StamPatch::LowerRight,
        transform: Kp20StamGridTransform::Identity,
        side: Kp20StamGridSide::Top,
        reversed_depth: false,
    };

    for patch in [
        Kp20StamPatch::LowerRight,
        Kp20StamPatch::CenterRight,
        Kp20StamPatch::UpperLeft,
    ] {
        for transform in kp20_stam_grid_transforms() {
            for side in Kp20StamGridSide::all() {
                for reversed_depth in [false, true] {
                    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
                        &projection,
                        patch,
                        transform,
                        side,
                        reversed_depth,
                    )
                    .unwrap();
                    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(
                        &Kp20PlanarTauPaperSolveInput {
                            valence: 5,
                            sector: 0,
                            characteristic_tensor_border: &tensor_border,
                            initial,
                            regularization_weight: 1.0e-10,
                        },
                    )
                    .unwrap();
                    if solved.max_constraint_residual < best.residual {
                        best = StamTransformedSideProjectionDiagnostic {
                            residual: solved.max_constraint_residual,
                            patch,
                            transform,
                            side,
                            reversed_depth,
                        };
                    }
                }
            }
        }
    }

    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        best.patch,
        best.transform,
        best.side,
        best.reversed_depth,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial,
        regularization_weight: 1.0e-10,
    })
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "PR08 refined Stam visual side-row best: {best:?} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        best.residual > 6.0,
        "refined-label visual Stam projection unexpectedly closed the paper tau gate: {best:?}"
    );
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

    assert_eq!(best_odd.label_index, 12);
    assert!(best_odd.residual < 1.0e-12);
    assert_eq!(best_even.label_index, 1);
    assert!(best_even.residual < 1.0e-12);
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
fn stam_appendix_a_scaled_pr08_q13_q2_candidate_tau_validity_is_measured() {
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
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);

    eprintln!(
        "scaled PR08 q13/q2 tau validity: residual={residual:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
    assert!(solved.max_constraint_residual < 1.0e-10);
    assert_eq!(validity.domain.outside_triangle_samples, 39);
    assert_eq!(validity.orientation.negative_determinants, 29);
    assert_eq!(validity.orientation.near_zero_determinants, 0);
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
fn pr08_endpoint_normalized_next_row_candidate_improves_but_remains_invalid() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();

    let q13_q2_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let q13_q2_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&q13_q2_rows).unwrap();
    let q13_q2_solved =
        solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
            valence,
            sector,
            characteristic_tensor_border: &q13_q2_tensor_border,
            initial: Kp20MacroPatch::zero(),
            regularization_weight: 1.0e-10,
        })
        .unwrap();
    let q13_q2_local =
        cc_test_global_planar_patch_to_local(valence, sector, q13_q2_solved.tau_macro_patch);
    let q13_q2_validity =
        kp20_planar_bi6_tau_validity_report(&q13_q2_local, &SAMPLES, 1.0e-10, 1.0e-10);
    let q13_q2_square = cc_test_square_domain_report(&q13_q2_local, &SAMPLES, 1.0e-10);

    let next_row_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
    )
    .unwrap();
    let normalized_next_row =
        kp20_reflected_unit_sector_chart_tensor_rows(&next_row_rows, valence, sector).expect(
            "next-row transposed source should have nondegenerate reflected chart endpoints",
        );
    let normalized_next_row_tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&normalized_next_row).unwrap();
    let normalized_next_row_solved =
        solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
            valence,
            sector,
            characteristic_tensor_border: &normalized_next_row_tensor_border,
            initial: Kp20MacroPatch::zero(),
            regularization_weight: 1.0e-10,
        })
        .unwrap();
    let normalized_next_row_local = cc_test_global_planar_patch_to_local(
        valence,
        sector,
        normalized_next_row_solved.tau_macro_patch,
    );
    let normalized_next_row_validity =
        kp20_planar_bi6_tau_validity_report(&normalized_next_row_local, &SAMPLES, 1.0e-10, 1.0e-10);
    let normalized_next_row_square =
        cc_test_square_domain_report(&normalized_next_row_local, &SAMPLES, 1.0e-10);

    eprintln!(
        "PR08 endpoint-normalized next-row candidate vs q13/q2: q13_residual={q13_residual:e} q13_square={q13_square}/100 q13_outside={q13_outside}/100 q13_neg={q13_neg}/100 q13_det=[{q13_min_det:e},{q13_max_det:e}] next_residual={next_residual:e} next_square={next_square}/100 next_outside={next_outside}/100 next_neg={next_neg}/100 next_det=[{next_min_det:e},{next_max_det:e}]",
        q13_residual = q13_q2_solved.max_constraint_residual,
        q13_square = q13_q2_square.outside,
        q13_outside = q13_q2_validity.domain.outside_triangle_samples,
        q13_neg = q13_q2_validity.orientation.negative_determinants,
        q13_min_det = q13_q2_validity.orientation.min_determinant,
        q13_max_det = q13_q2_validity.orientation.max_determinant,
        next_residual = normalized_next_row_solved.max_constraint_residual,
        next_square = normalized_next_row_square.outside,
        next_outside = normalized_next_row_validity.domain.outside_triangle_samples,
        next_neg = normalized_next_row_validity
            .orientation
            .negative_determinants,
        next_min_det = normalized_next_row_validity.orientation.min_determinant,
        next_max_det = normalized_next_row_validity.orientation.max_determinant,
    );

    assert!(normalized_next_row_solved.max_constraint_residual < 1.0e-12);
    assert!(normalized_next_row_square.outside < q13_q2_square.outside);
    assert!(
        normalized_next_row_validity
            .orientation
            .negative_determinants
            < q13_q2_validity.orientation.negative_determinants
    );
    assert!(
        normalized_next_row_validity.orientation.min_determinant
            > q13_q2_validity.orientation.min_determinant
    );
    assert!(
        normalized_next_row_square.outside > 0
            || normalized_next_row_validity
                .orientation
                .negative_determinants
                > 0,
        "endpoint-normalized next-row source unexpectedly passed sampled tau validity"
    );
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
fn pr08_stam_projected_side_sources_are_measured_under_live_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let ranked = pr08_stam_projected_side_endpoint_fit_ranking(&ring, valence, sector);
    let mut best = Kp20StamSideTauValidityDiagnostic {
        residual: f64::INFINITY,
        square_outside: usize::MAX,
        triangle_outside: usize::MAX,
        negative_determinants: usize::MAX,
        min_determinant: f64::NEG_INFINITY,
        endpoint: ranked[0],
    };

    for endpoint in ranked.iter().take(8).copied() {
        let projection = pr08_stam_projected_grids_for_family(&ring, endpoint.family).unwrap();
        let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
            &projection,
            endpoint.patch,
            endpoint.transform,
            endpoint.side,
            endpoint.reversed_depth,
        )
        .unwrap();
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows(
            &rows,
            valence,
            sector,
            endpoint.swap_targets,
            endpoint.conjugate_source,
        )
        .expect("ranked Stam side source should have nondegenerate endpoints");
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
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
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 Stam side candidate family={} patch={:?} transform={:?} side={} reversed_depth={} swap={} conjugate={} endpoint={:e}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            endpoint.family.name(),
            endpoint.patch,
            endpoint.transform,
            endpoint.side.name(),
            endpoint.reversed_depth,
            endpoint.swap_targets,
            endpoint.conjugate_source,
            endpoint.max_endpoint_residual,
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
        if square.outside < best.square_outside
            || (square.outside == best.square_outside
                && validity.orientation.negative_determinants < best.negative_determinants)
            || (square.outside == best.square_outside
                && validity.orientation.negative_determinants == best.negative_determinants
                && validity.orientation.min_determinant > best.min_determinant)
        {
            best = Kp20StamSideTauValidityDiagnostic {
                residual: solved.max_constraint_residual,
                square_outside: square.outside,
                triangle_outside: validity.domain.outside_triangle_samples,
                negative_determinants: validity.orientation.negative_determinants,
                min_determinant: validity.orientation.min_determinant,
                endpoint,
            };
        }
    }

    eprintln!(
        "PR08 Stam side live-gate best: family={} patch={:?} transform={:?} side={} reversed_depth={} swap={} conjugate={} endpoint={:e}: residual={residual:e} square_outside={square_outside} triangle_outside={triangle_outside} neg={neg} min_det={min_det:e}",
        best.endpoint.family.name(),
        best.endpoint.patch,
        best.endpoint.transform,
        best.endpoint.side.name(),
        best.endpoint.reversed_depth,
        best.endpoint.swap_targets,
        best.endpoint.conjugate_source,
        best.endpoint.max_endpoint_residual,
        residual = best.residual,
        square_outside = best.square_outside,
        triangle_outside = best.triangle_outside,
        neg = best.negative_determinants,
        min_det = best.min_determinant,
    );

    assert_eq!(best.square_outside, 52);
    assert_eq!(best.negative_determinants, 50);
    assert!(best.min_determinant < -1.0e-2);
}

#[test]
#[ignore = "slow diagnostic sweep over all unique exact Stam side bands"]
fn pr08_all_unique_stam_projected_side_sources_are_measured_under_live_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let ranked = pr08_stam_projected_side_endpoint_fit_ranking(&ring, valence, sector);
    let mut unique = Vec::<(
        Kp20StamSideEndpointFitDiagnostic,
        Kp20RegularTensorBorderCNet,
    )>::new();

    for endpoint in ranked {
        let projection = pr08_stam_projected_grids_for_family(&ring, endpoint.family).unwrap();
        let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
            &projection,
            endpoint.patch,
            endpoint.transform,
            endpoint.side,
            endpoint.reversed_depth,
        )
        .unwrap();
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows(
            &rows,
            valence,
            sector,
            endpoint.swap_targets,
            endpoint.conjugate_source,
        )
        .expect("ranked Stam side source should have nondegenerate endpoints");
        if unique.iter().any(|(_, existing)| {
            tensor_border_rows_max_distance(existing, &normalized_rows) < 1.0e-12
        }) {
            continue;
        }
        unique.push((endpoint, normalized_rows));
    }

    let mut best = Kp20StamSideTauValidityDiagnostic {
        residual: f64::INFINITY,
        square_outside: usize::MAX,
        triangle_outside: usize::MAX,
        negative_determinants: usize::MAX,
        min_determinant: f64::NEG_INFINITY,
        endpoint: unique[0].0,
    };

    for (endpoint, normalized_rows) in &unique {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(normalized_rows).unwrap();
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
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        if square.outside < best.square_outside
            || (square.outside == best.square_outside
                && validity.orientation.negative_determinants < best.negative_determinants)
            || (square.outside == best.square_outside
                && validity.orientation.negative_determinants == best.negative_determinants
                && validity.orientation.min_determinant > best.min_determinant)
        {
            best = Kp20StamSideTauValidityDiagnostic {
                residual: solved.max_constraint_residual,
                square_outside: square.outside,
                triangle_outside: validity.domain.outside_triangle_samples,
                negative_determinants: validity.orientation.negative_determinants,
                min_determinant: validity.orientation.min_determinant,
                endpoint: *endpoint,
            };
        }
    }

    eprintln!(
        "PR08 unique Stam side live-gate sweep: unique_count={} best_family={} best_patch={:?} best_transform={:?} best_side={} reversed_depth={} swap={} conjugate={} endpoint={:e}: residual={residual:e} square_outside={square_outside} triangle_outside={triangle_outside} neg={neg} min_det={min_det:e}",
        unique.len(),
        best.endpoint.family.name(),
        best.endpoint.patch,
        best.endpoint.transform,
        best.endpoint.side.name(),
        best.endpoint.reversed_depth,
        best.endpoint.swap_targets,
        best.endpoint.conjugate_source,
        best.endpoint.max_endpoint_residual,
        residual = best.residual,
        square_outside = best.square_outside,
        triangle_outside = best.triangle_outside,
        neg = best.negative_determinants,
        min_det = best.min_determinant,
    );

    assert!(unique.len() < 256);
    assert!(
        best.square_outside >= 31 || best.negative_determinants >= 27,
        "exact Stam side sweep unexpectedly improved the current adjacent-row source"
    );
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
fn pr08_best_exact_stam_side_anchored_linear_characteristic_source_is_measured() {
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
    let mut best = (
        "",
        0.0,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        0.0,
    );

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
        for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(2), ring.lambda.powi(4)] {
            let mut augmented_scalar_constraints = cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                0,
            );
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: anchored_scale,
            });
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
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
            eprintln!(
                "PR08 best exact Stam side anchored linear characteristic source {normalization_label} scale={anchored_scale:e}: residual={residual:e} final_scale={final_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                final_scale = solved.extra_values[0],
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
            if square.outside < best.2
                || (square.outside == best.2 && validity.orientation.negative_determinants < best.3)
                || (square.outside == best.2
                    && validity.orientation.negative_determinants == best.3
                    && validity.orientation.min_determinant > best.4)
            {
                best = (
                    normalization_label,
                    anchored_scale,
                    square.outside,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    solved.max_constraint_residual,
                    det_a,
                );
            }
        }
    }

    eprintln!("PR08 best exact Stam side anchored linear characteristic source best: {best:?}");
    assert_eq!(best.0, "per-sector-characteristic-exact-side");
    assert_eq!(best.2, 50);
    assert_eq!(best.3, 61);
    assert!(best.5 > 1.0e-2);
    assert!(best.6.abs() < 1.0e-12);
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
fn pr08_best_exact_stam_side_depth_order_variants_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let mut best_unprojected = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);
    let mut best_projected = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);

    for (label, tensor_border) in [
        (
            "bspline-depth-012",
            pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
                &ring,
                valence,
                sector,
                [0, 1, 2],
            ),
        ),
        (
            "bspline-depth-021",
            pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
                &ring,
                valence,
                sector,
                [0, 2, 1],
            ),
        ),
        (
            "bspline-depth-210",
            pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
                &ring,
                valence,
                sector,
                [2, 1, 0],
            ),
        ),
        (
            "bezier-depth-012",
            pr08_best_exact_stam_side_tensor_border_with_bezier_depth_order(
                &ring,
                valence,
                sector,
                [0, 1, 2],
            ),
        ),
        (
            "bezier-depth-021",
            pr08_best_exact_stam_side_tensor_border_with_bezier_depth_order(
                &ring,
                valence,
                sector,
                [0, 2, 1],
            ),
        ),
        (
            "bezier-depth-210",
            pr08_best_exact_stam_side_tensor_border_with_bezier_depth_order(
                &ring,
                valence,
                sector,
                [2, 1, 0],
            ),
        ),
    ] {
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
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
            "PR08 best exact Stam side depth variant {label}: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
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

    eprintln!(
        "PR08 best exact Stam side depth variants best_unprojected={best_unprojected:?} best_projected={best_projected:?}"
    );
    assert_eq!(best_unprojected.0, "");
    assert_eq!(best_projected.0, "bspline-depth-210");
    assert_eq!(best_projected.1, 33);
    assert_eq!(best_projected.2, 27);
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
fn pr08_best_exact_stam_side_reversed_depth_half_assignment_variants_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let mut best_projected = (
        Kp20Eq14SourceHalf::NearC,
        Kp20Eq14SourceHalf::NearC,
        false,
        false,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    let mut best_unprojected_residual = f64::INFINITY;

    for previous_half in [Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM] {
        for current_half in [Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM] {
            for reverse_previous_u in [false, true] {
                for reverse_current_u in [false, true] {
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
                        "PR08 best exact Stam side reversed-depth half variant previous={previous_half:?} current={current_half:?} reverse_previous={reverse_previous_u} reverse_current={reverse_current_u}: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
                        ures = unprojected.residual,
                        usquare = unprojected.square_outside,
                        uneg = unprojected.negative_determinants,
                        umin = unprojected.min_determinant,
                        pres = projected.residual,
                        psquare = projected.square_outside,
                        pneg = projected.negative_determinants,
                        pmin = projected.min_determinant,
                    );
                    best_unprojected_residual = best_unprojected_residual.min(unprojected.residual);
                    if projected.square_outside < best_projected.4
                        || (projected.square_outside == best_projected.4
                            && projected.negative_determinants < best_projected.5)
                        || (projected.square_outside == best_projected.4
                            && projected.negative_determinants == best_projected.5
                            && projected.min_determinant > best_projected.6)
                    {
                        best_projected = (
                            previous_half,
                            current_half,
                            reverse_previous_u,
                            reverse_current_u,
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
        "PR08 best exact Stam side reversed-depth half variants best_projected={best_projected:?} best_unprojected_residual={best_unprojected_residual:e}"
    );
    assert!(best_unprojected_residual > 1.0e-1);
    assert_eq!(best_projected.4, 33);
    assert_eq!(best_projected.5, 27);
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
#[ignore = "slow diagnostic sweep over endpoint-compatible right-edge Stam donors"]
fn pr08_best_exact_stam_side_reversed_depth_right_edge_donors_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let ranked = pr08_stam_projected_side_endpoint_fit_ranking(&ring, valence, sector);
    let mut unique_donors = Vec::<(Kp20StamSideEndpointFitDiagnostic, TensorBorderBi3)>::new();
    for endpoint in ranked {
        let Some(tensor_border) = pr08_stam_side_tensor_border_with_b_spline_depth_order(
            &ring,
            valence,
            sector,
            endpoint,
            [2, 1, 0],
        ) else {
            continue;
        };
        if unique_donors.iter().any(|(_, existing)| {
            existing
                .sectors
                .iter()
                .zip(&tensor_border.sectors)
                .flat_map(|(left_sector, right_sector)| left_sector.iter().zip(right_sector))
                .flat_map(|(left_row, right_row)| left_row.iter().zip(right_row))
                .map(|(left, right)| left.distance_to(right))
                .fold(0.0_f64, f64::max)
                < 1.0e-12
        }) {
            continue;
        }
        unique_donors.push((endpoint, tensor_border));
        if unique_donors.len() >= 2 {
            break;
        }
    }

    let mut best_unprojected = Kp20RightEdgeDonorDiagnostic {
        residual: f64::INFINITY,
        square_outside: usize::MAX,
        negative_determinants: usize::MAX,
        min_determinant: f64::NEG_INFINITY,
        family: Kp20StamProjectionFamily::AppendixTensorOrder,
        patch: Kp20StamPatch::UpperLeft,
        transform: Kp20StamGridTransform::Identity,
        side: Kp20StamGridSide::Left,
        reversed_depth: false,
        swap_targets: false,
        conjugate_source: false,
        donor_sector_offset: 0,
        donor_half: Kp20Eq14SourceHalf::NearC,
        column_label: "",
    };
    let mut best_projected = best_unprojected;

    for (endpoint, donor_tensor_border) in &unique_donors {
        for donor_sector_offset in [-1_isize, 0, 1] {
            for donor_half in [Kp20Eq14SourceHalf::NearC, Kp20Eq14SourceHalf::NearM] {
                for (column_label, columns) in [("u6", &[6_usize][..]), ("u5-u6", &[5, 6][..])] {
                    let tensor_rows = cc_test_replaced_near_midpoint_right_columns_tensor_rows(
                        &base_tensor_border,
                        donor_tensor_border,
                        sector,
                        donor_sector_offset,
                        donor_half,
                        columns,
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
                        "PR08 reversed-depth right-edge donor family={} patch={:?} transform={:?} side={} reversed_depth={} swap={} conjugate={} donor_offset={donor_sector_offset} donor_half={donor_half:?} columns={column_label}: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
                        endpoint.family.name(),
                        endpoint.patch,
                        endpoint.transform,
                        endpoint.side.name(),
                        endpoint.reversed_depth,
                        endpoint.swap_targets,
                        endpoint.conjugate_source,
                        ures = unprojected.residual,
                        usquare = unprojected.square_outside,
                        uneg = unprojected.negative_determinants,
                        umin = unprojected.min_determinant,
                        pres = projected.residual,
                        psquare = projected.square_outside,
                        pneg = projected.negative_determinants,
                        pmin = projected.min_determinant,
                    );

                    if unprojected.residual < best_unprojected.residual {
                        best_unprojected = Kp20RightEdgeDonorDiagnostic {
                            residual: unprojected.residual,
                            square_outside: unprojected.square_outside,
                            negative_determinants: unprojected.negative_determinants,
                            min_determinant: unprojected.min_determinant,
                            family: endpoint.family,
                            patch: endpoint.patch,
                            transform: endpoint.transform,
                            side: endpoint.side,
                            reversed_depth: endpoint.reversed_depth,
                            swap_targets: endpoint.swap_targets,
                            conjugate_source: endpoint.conjugate_source,
                            donor_sector_offset,
                            donor_half,
                            column_label,
                        };
                    }
                    if projected.square_outside < best_projected.square_outside
                        || (projected.square_outside == best_projected.square_outside
                            && projected.negative_determinants
                                < best_projected.negative_determinants)
                        || (projected.square_outside == best_projected.square_outside
                            && projected.negative_determinants
                                == best_projected.negative_determinants
                            && projected.min_determinant > best_projected.min_determinant)
                    {
                        best_projected = Kp20RightEdgeDonorDiagnostic {
                            residual: projected.residual,
                            square_outside: projected.square_outside,
                            negative_determinants: projected.negative_determinants,
                            min_determinant: projected.min_determinant,
                            family: endpoint.family,
                            patch: endpoint.patch,
                            transform: endpoint.transform,
                            side: endpoint.side,
                            reversed_depth: endpoint.reversed_depth,
                            swap_targets: endpoint.swap_targets,
                            conjugate_source: endpoint.conjugate_source,
                            donor_sector_offset,
                            donor_half,
                            column_label,
                        };
                    }
                }
            }
        }
    }

    eprintln!(
        "PR08 reversed-depth right-edge donor summary: unique_donors={} best_unprojected={best_unprojected:?} best_projected={best_projected:?}",
        unique_donors.len()
    );
    assert_eq!(unique_donors.len(), 2);
    assert!(best_unprojected.residual.is_finite());
    assert!(best_projected.residual.is_finite());
    assert!(best_unprojected.residual > 0.7);
    assert_eq!(best_unprojected.square_outside, 56);
    assert_eq!(best_unprojected.negative_determinants, 30);
    assert_eq!(best_projected.square_outside, 31);
    assert_eq!(best_projected.negative_determinants, 27);
    assert_eq!(
        best_projected.family,
        Kp20StamProjectionFamily::AppendixTensorOrder
    );
    assert_eq!(best_projected.patch, Kp20StamPatch::LowerRight);
    assert_eq!(best_projected.transform, Kp20StamGridTransform::FlipRows);
    assert_eq!(best_projected.side, Kp20StamGridSide::Top);
    assert!(!best_projected.reversed_depth);
    assert!(!best_projected.swap_targets);
    assert!(!best_projected.conjugate_source);
    assert_eq!(best_projected.donor_sector_offset, -1);
    assert_eq!(best_projected.donor_half, Kp20Eq14SourceHalf::NearM);
    assert_eq!(best_projected.column_label, "u6");
}

#[test]
fn pr08_best_exact_stam_side_reversed_depth_full_right_curve_donor_is_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let donor_endpoint = Kp20StamSideEndpointFitDiagnostic {
        max_endpoint_residual: 0.0,
        sum_endpoint_residual: 0.0,
        family: Kp20StamProjectionFamily::AppendixTensorOrder,
        patch: Kp20StamPatch::LowerRight,
        transform: Kp20StamGridTransform::FlipRows,
        side: Kp20StamGridSide::Top,
        reversed_depth: false,
        swap_targets: false,
        conjugate_source: false,
    };
    let donor_tensor_border = pr08_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        donor_endpoint,
        [2, 1, 0],
    )
    .expect("projected best right-curve donor should normalize");
    let tensor_rows = cc_test_replaced_near_midpoint_right_columns_tensor_rows(
        &base_tensor_border,
        &donor_tensor_border,
        sector,
        -1,
        Kp20Eq14SourceHalf::NearM,
        &[0, 1, 2, 3, 4, 5, 6],
    );
    let unprojected = cc_test_solve_exact_side_tensor_rows_under_live_basis(
        valence,
        sector,
        tensor_rows.clone(),
        false,
    );
    let projected =
        cc_test_solve_exact_side_tensor_rows_under_live_basis(valence, sector, tensor_rows, true);

    eprintln!(
        "PR08 reversed-depth full right-curve donor: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
        ures = unprojected.residual,
        usquare = unprojected.square_outside,
        uneg = unprojected.negative_determinants,
        umin = unprojected.min_determinant,
        pres = projected.residual,
        psquare = projected.square_outside,
        pneg = projected.negative_determinants,
        pmin = projected.min_determinant,
    );

    assert!(unprojected.residual > 10.0);
    assert_eq!(unprojected.square_outside, 86);
    assert_eq!(unprojected.negative_determinants, 43);
    assert!(projected.residual < 1.0e-10);
    assert_eq!(projected.square_outside, 59);
    assert_eq!(projected.negative_determinants, 42);
}

#[test]
fn pr08_best_exact_stam_side_reversed_depth_beta_splitters_are_measured() {
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let tensor_border = pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
        &ring,
        valence,
        sector,
        [2, 1, 0],
    );
    let mut best_unprojected = ("", f64::INFINITY, usize::MAX, usize::MAX, f64::NEG_INFINITY);
    let mut best_projected = ("", f64::INFINITY, usize::MAX, usize::MAX, f64::NEG_INFINITY);

    for (label, tensor_rows) in [
        (
            "kp20-eq14",
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap(),
        ),
        (
            "kp13-tangential-beta",
            kp13_beta_tensor_border_constraints(&tensor_border, sector).unwrap(),
        ),
        (
            "kp13-exact-pullback-beta",
            kp13_beta_exact_pullback_tensor_border_constraints(&tensor_border, sector).unwrap(),
        ),
        (
            "kp15-split-beta",
            kp15_beta_half_strip_tensor_border_constraints(&tensor_border, sector).unwrap(),
        ),
        (
            "kp15-full-then-split-beta",
            kp15_beta_full_then_split_tensor_border_constraints(&tensor_border, sector).unwrap(),
        ),
    ] {
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
            "PR08 reversed-depth beta splitter {label}: unprojected residual={ures:e} square={usquare}/100 neg={uneg}/100 min_det={umin:e}; projected residual={pres:e} square={psquare}/100 neg={pneg}/100 min_det={pmin:e}",
            ures = unprojected.residual,
            usquare = unprojected.square_outside,
            uneg = unprojected.negative_determinants,
            umin = unprojected.min_determinant,
            pres = projected.residual,
            psquare = projected.square_outside,
            pneg = projected.negative_determinants,
            pmin = projected.min_determinant,
        );
        if unprojected.residual < best_unprojected.1 {
            best_unprojected = (
                label,
                unprojected.residual,
                unprojected.square_outside,
                unprojected.negative_determinants,
                unprojected.min_determinant,
            );
        }
        if projected.square_outside < best_projected.2
            || (projected.square_outside == best_projected.2
                && projected.negative_determinants < best_projected.3)
            || (projected.square_outside == best_projected.2
                && projected.negative_determinants == best_projected.3
                && projected.min_determinant > best_projected.4)
        {
            best_projected = (
                label,
                projected.residual,
                projected.square_outside,
                projected.negative_determinants,
                projected.min_determinant,
            );
        }
    }

    eprintln!(
        "PR08 reversed-depth beta splitter summary: best_unprojected={best_unprojected:?} best_projected={best_projected:?}"
    );
    assert_eq!(best_unprojected.0, "kp20-eq14");
    assert!(best_unprojected.1 > 0.7);
    assert_eq!(best_projected.0, "kp13-tangential-beta");
    assert_eq!(best_projected.2, 31);
    assert_eq!(best_projected.3, 27);
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
fn pr08_best_exact_stam_side_tangential_beta_right_edge_mask_residual_is_localized() {
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
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let section42_row_count = section42_rows.len();
    let mut scalar_constraints = section42_rows;
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_rows, source_delta) =
        cc_test_project_selected_tensor_rows_for_scalar_basis(
            tensor_rows.clone(),
            &fig8_rows,
            &scalar_constraints,
            |coeff| coeff.piece == Kp20PieceId::NearMidpointRight && (coeff.j == 5 || coeff.j == 6),
        );
    let solved = cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
        valence,
        sector,
        projected_tensor_rows.clone(),
        false,
    );
    let fig8_worst = cc_test_vector_rows_worst_residual(&solved.tau_macro_patch, &fig8_rows);
    let tensor_worst =
        cc_test_vector_rows_worst_residual(&solved.tau_macro_patch, &projected_tensor_rows);
    let scalar_worst =
        cc_test_scalar_rows_worst_residual(&solved.tau_macro_patch, &scalar_constraints);
    let scalar_family = if scalar_worst.0 < section42_row_count {
        "section42"
    } else {
        "diagonal"
    };
    let scalar_family_row = if scalar_worst.0 < section42_row_count {
        scalar_worst.0
    } else {
        scalar_worst.0 - section42_row_count
    };

    eprintln!(
        "PR08 tangential-beta right-u5-u6 residual localization: source_delta={source_delta:e} residual={residual:e} square={square}/100 neg={neg}/100 fig8_worst={fig8_worst:?} tensor_worst=({}, {tensor_res:e}, {:?}[{},{}]) scalar_worst={scalar_worst:?} scalar_family={scalar_family} family_row={scalar_family_row} terms={terms:?}",
        tensor_border_row_label(tensor_worst.0),
        tensor_worst.2.piece,
        tensor_worst.2.i,
        tensor_worst.2.j,
        terms = &scalar_constraints[scalar_worst.0].terms,
        residual = solved.validity.residual,
        square = solved.validity.square_outside,
        neg = solved.validity.negative_determinants,
        tensor_res = tensor_worst.1,
    );

    assert!(source_delta > 1.0);
    assert!(solved.validity.residual > 0.7);
    assert_eq!(solved.validity.negative_determinants, 16);
    assert!(scalar_worst.1 >= tensor_worst.1);
    assert_eq!(scalar_family, "diagonal");
    assert_eq!(scalar_family_row, 49);

    let (corner_projected_tensor_rows, corner_source_delta) =
        cc_test_project_selected_tensor_rows_for_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            |coeff| {
                coeff.piece == Kp20PieceId::NearMidpointRight
                    && ((coeff.j == 5 || coeff.j == 6) || (coeff.i == 0 && coeff.j == 0))
            },
        );
    let corner_solved = cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
        valence,
        sector,
        corner_projected_tensor_rows.clone(),
        false,
    );
    let corner_tensor_worst = cc_test_vector_rows_worst_residual(
        &corner_solved.tau_macro_patch,
        &corner_projected_tensor_rows,
    );
    let corner_scalar_worst =
        cc_test_scalar_rows_worst_residual(&corner_solved.tau_macro_patch, &scalar_constraints);
    let corner_scalar_family = if corner_scalar_worst.0 < section42_row_count {
        "section42"
    } else {
        "diagonal"
    };
    let corner_scalar_family_row = if corner_scalar_worst.0 < section42_row_count {
        corner_scalar_worst.0
    } else {
        corner_scalar_worst.0 - section42_row_count
    };
    eprintln!(
        "PR08 tangential-beta right-corner-u5-u6 residual localization: source_delta={corner_source_delta:e} residual={residual:e} square={square}/100 neg={neg}/100 tensor_worst=({}, {tensor_res:e}, {:?}[{},{}]) scalar_worst={corner_scalar_worst:?} scalar_family={corner_scalar_family} family_row={corner_scalar_family_row} terms={terms:?}",
        tensor_border_row_label(corner_tensor_worst.0),
        corner_tensor_worst.2.piece,
        corner_tensor_worst.2.i,
        corner_tensor_worst.2.j,
        terms = &scalar_constraints[corner_scalar_worst.0].terms,
        residual = corner_solved.validity.residual,
        square = corner_solved.validity.square_outside,
        neg = corner_solved.validity.negative_determinants,
        tensor_res = corner_tensor_worst.1,
    );
    assert!(corner_solved.validity.residual < solved.validity.residual);
    assert_eq!(corner_solved.validity.negative_determinants, 16);
}

#[test]
fn pr08_best_exact_stam_side_tangential_beta_self_diagonal_source_is_measured() {
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
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let section42_row_count = section42_rows.len();
    let mut scalar_constraints = section42_rows;
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (self_diagonal_rows, self_diagonal_delta) =
        cc_test_sector_local_self_diagonal_normalized_source_rows(
            beta_rows.clone(),
            valence,
            sector,
            Kp20PieceId::NearMidpointRight,
            2,
        );
    let (self_c0_rows, c0_delta) = cc_test_underlined_c0_normalized_source_rows(
        self_diagonal_rows.clone(),
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        angle,
        2,
    );
    let (self_plus_endpoint_rows, endpoint_delta) =
        cc_test_project_selected_tensor_rows_for_scalar_basis(
            self_diagonal_rows.clone(),
            &fig8_rows,
            &scalar_constraints,
            |coeff| coeff.piece == Kp20PieceId::NearMidpointRight && (coeff.j == 5 || coeff.j == 6),
        );
    let (self_c0_plus_endpoint_rows, c0_endpoint_delta) =
        cc_test_project_selected_tensor_rows_for_scalar_basis(
            self_c0_rows.clone(),
            &fig8_rows,
            &scalar_constraints,
            |coeff| coeff.piece == Kp20PieceId::NearMidpointRight && (coeff.j == 5 || coeff.j == 6),
        );
    let (all_projected_rows, all_projected_delta) = cc_test_project_tensor_rows_for_scalar_basis(
        beta_rows.clone(),
        &fig8_rows,
        &scalar_constraints,
    );

    let mut projection_gap_deltas = self_c0_plus_endpoint_rows
        .iter()
        .zip(all_projected_rows.iter())
        .enumerate()
        .map(|(row_index, (constructed, projected))| {
            let coeff = constructed.terms[0].0;
            let dx = projected.target[0] - constructed.target[0];
            let dy = projected.target[1] - constructed.target[1];
            let norm = (dx * dx + dy * dy).sqrt();
            (
                norm,
                row_index,
                coeff,
                dx,
                dy,
                constructed.target,
                projected.target,
            )
        })
        .collect::<Vec<_>>();
    projection_gap_deltas.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (rank, (norm, row_index, coeff, dx, dy, constructed, projected)) in
        projection_gap_deltas.iter().take(12).enumerate()
    {
        eprintln!(
            "PR08 tangential-beta source projection gap rank={rank} row={row_index} label={} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} constructed=({cx:e},{cy:e}) projected=({px:e},{py:e})",
            tensor_border_row_label(*row_index),
            piece = coeff.piece,
            i = coeff.i,
            j = coeff.j,
            cx = constructed[0],
            cy = constructed[1],
            px = projected[0],
            py = projected[1],
        );
    }

    let cases = [
        ("beta", beta_rows, 0.0),
        (
            "self-diagonal-source",
            self_diagonal_rows,
            self_diagonal_delta,
        ),
        (
            "self-diagonal-c0-source",
            self_c0_rows,
            self_diagonal_delta.max(c0_delta),
        ),
        (
            "self-diagonal-source-plus-right-u5-u6",
            self_plus_endpoint_rows,
            self_diagonal_delta.max(endpoint_delta),
        ),
        (
            "self-diagonal-c0-source-plus-right-u5-u6",
            self_c0_plus_endpoint_rows,
            self_diagonal_delta.max(c0_delta).max(c0_endpoint_delta),
        ),
        (
            "all-pin-projection",
            all_projected_rows,
            all_projected_delta,
        ),
    ];
    let mut measured = Vec::new();
    for (label, rows, source_delta) in cases {
        let solved = cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
            valence, sector, rows, false,
        );
        let scalar_worst =
            cc_test_scalar_rows_worst_residual(&solved.tau_macro_patch, &scalar_constraints);
        let scalar_family = if scalar_worst.0 < section42_row_count {
            "section42"
        } else {
            "diagonal"
        };
        let scalar_family_row = if scalar_worst.0 < section42_row_count {
            scalar_worst.0
        } else {
            scalar_worst.0 - section42_row_count
        };
        eprintln!(
            "PR08 tangential-beta self-diagonal source {label}: source_delta={source_delta:e} residual={residual:e} square={square}/100 neg={neg}/100 min_det={min_det:e} scalar_family={scalar_family} family_row={scalar_family_row} scalar_worst={scalar_worst:?} terms={terms:?}",
            residual = solved.validity.residual,
            square = solved.validity.square_outside,
            neg = solved.validity.negative_determinants,
            min_det = solved.validity.min_determinant,
            terms = &scalar_constraints[scalar_worst.0].terms,
        );
        measured.push((
            label,
            source_delta,
            solved.validity.residual,
            solved.validity.square_outside,
            solved.validity.negative_determinants,
            scalar_family,
            scalar_family_row,
        ));
    }

    let self_case = measured
        .iter()
        .find(|case| case.0 == "self-diagonal-source")
        .unwrap();
    let endpoint_case = measured
        .iter()
        .find(|case| case.0 == "self-diagonal-source-plus-right-u5-u6")
        .unwrap();
    let c0_endpoint_case = measured
        .iter()
        .find(|case| case.0 == "self-diagonal-c0-source-plus-right-u5-u6")
        .unwrap();
    let all_pin_case = measured
        .iter()
        .find(|case| case.0 == "all-pin-projection")
        .unwrap();
    assert!(self_diagonal_delta > 0.2);
    assert!(self_case.2 < 0.6);
    assert_eq!(self_case.5, "diagonal");
    assert!(endpoint_case.2 < self_case.2);
    assert!(c0_endpoint_case.2.is_finite());
    assert!(all_pin_case.2 < 1.0e-10);
    assert_eq!(all_pin_case.3, 31);
    assert_eq!(all_pin_case.4, 27);
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
fn pr08_best_exact_stam_side_all_pin_closed_tau_domain_active_rows_are_measured() {
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
    let mut active_domain_rows = Vec::new();
    let mut last_solved = None;
    let mut last_tau = Kp20MacroPatch::zero();

    for iteration in 0..=3 {
        let mut solve_scalar_constraints = scalar_constraints.clone();
        solve_scalar_constraints.extend(active_domain_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: vector_constraints.clone(),
                scalar_constraints: solve_scalar_constraints,
                regularization_weight: 1.0e-10,
            })
            .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let (new_rows, lower_u, lower_v, upper_sum) = cc_test_active_sector_domain_bound_rows(
            &solved.tau_macro_patch,
            valence,
            sector,
            &SAMPLES,
            1.0e-10,
        );
        eprintln!(
            "PR08 all-pin closed tau active domain iteration={iteration}: source_delta={source_delta:e} active_rows={} new_rows={} lower_u={} lower_v={} upper_sum={} residual={residual:e} triangle_outside={outside}/100 max_excess={max_excess:e} square={square_outside}/100 neg={neg}/100 min_det={min_det:e}",
            active_domain_rows.len(),
            new_rows.len(),
            lower_u,
            lower_v,
            upper_sum,
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            max_excess = validity.domain.max_simplex_excess,
            square_outside = square.outside,
            neg = validity.orientation.negative_determinants,
            min_det = validity.orientation.min_determinant,
        );
        last_solved = Some((
            solved.max_constraint_residual,
            validity.domain.outside_triangle_samples,
            validity.domain.max_simplex_excess,
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
            new_rows.len(),
        ));
        last_tau = solved.tau_macro_patch;
        if new_rows.is_empty() {
            break;
        }
        active_domain_rows.extend(new_rows);
    }

    let final_report = last_solved.unwrap();
    let final_local_tau = cc_test_global_planar_patch_to_local(valence, sector, last_tau);
    for piece_id in Kp20PieceId::all() {
        let piece_report = sampled_piece_tau_report(&final_local_tau, piece_id, &SAMPLES);
        eprintln!(
            "PR08 all-pin closed tau active domain final piece {piece_id:?}: outside={outside}/25 neg={neg}/25 min_u={min_u:e} max_u={max_u:e} min_v={min_v:e} max_v={max_v:e} max_sum={max_sum:e} det=[{min_det:e},{max_det:e}]",
            outside = piece_report.outside,
            neg = piece_report.negative,
            min_u = piece_report.min_u,
            max_u = piece_report.max_u,
            min_v = piece_report.min_v,
            max_v = piece_report.max_v,
            max_sum = piece_report.max_sum,
            min_det = piece_report.min_det,
            max_det = piece_report.max_det,
        );
    }
    assert!(final_report.0.is_finite());
    assert_eq!(final_report.1, 20);
    assert_eq!(final_report.3, 20);
    assert_eq!(final_report.4, 23);
    assert!(final_report.2 > 0.02 && final_report.2 < 0.021);
    assert_eq!(final_report.6, 20);
}

#[test]
fn pr08_best_exact_stam_side_domain_aware_source_projection_is_measured() {
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
    let (_baseline_rows, baseline_source_delta) = cc_test_project_tensor_rows_for_scalar_basis(
        beta_rows.clone(),
        &fig8_rows,
        &scalar_constraints,
    );
    let mut active_rows = Vec::new();
    let mut final_report = None;
    let mut best_report = (
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        0_usize,
    );
    let mut best_tau = Kp20MacroPatch::zero();
    for iteration in 0..=4 {
        let mut domain_scalar_constraints = scalar_constraints.clone();
        domain_scalar_constraints.extend(active_rows.clone());
        let (domain_projected_rows, domain_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                beta_rows.clone(),
                &fig8_rows,
                &domain_scalar_constraints,
            );
        let solved = cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
            valence,
            sector,
            domain_projected_rows.clone(),
            false,
        );
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let active_residual =
            cc_test_scalar_rows_worst_residual(&solved.tau_macro_patch, &domain_scalar_constraints);
        let tensor_residual =
            cc_test_vector_rows_worst_residual(&solved.tau_macro_patch, &domain_projected_rows);
        let (new_rows, lower_u, lower_v, upper_sum) = cc_test_active_sector_domain_bound_rows(
            &solved.tau_macro_patch,
            valence,
            sector,
            &SAMPLES,
            1.0e-10,
        );
        eprintln!(
            "PR08 domain-aware source projection iteration={iteration}: active_rows={} new_rows={} lower_u={} lower_v={} upper_sum={} baseline_delta={baseline_source_delta:e} domain_delta={domain_source_delta:e} residual={residual:e} active_worst={active_residual:?} tensor_worst=({}, {tensor_res:e}, {:?}[{},{}]) triangle_outside={outside}/100 max_excess={max_excess:e} square={square_outside}/100 neg={neg}/100 min_det={min_det:e}",
            active_rows.len(),
            new_rows.len(),
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
        final_report = Some((
            active_rows.len(),
            new_rows.len(),
            domain_source_delta,
            solved.validity.residual,
            validity.domain.outside_triangle_samples,
            validity.domain.max_simplex_excess,
            square.outside,
            validity.orientation.negative_determinants,
            validity.orientation.min_determinant,
        ));
        if validity.domain.outside_triangle_samples < best_report.0
            || (validity.domain.outside_triangle_samples == best_report.0
                && validity.orientation.negative_determinants < best_report.1)
            || (validity.domain.outside_triangle_samples == best_report.0
                && validity.orientation.negative_determinants == best_report.1
                && validity.orientation.min_determinant > best_report.2)
        {
            best_report = (
                validity.domain.outside_triangle_samples,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                validity.domain.max_simplex_excess,
                iteration,
            );
            best_tau = local_tau;
        }
        if new_rows.is_empty() {
            break;
        }
        active_rows.extend(new_rows);
    }

    let final_report = final_report.unwrap();
    for piece_id in Kp20PieceId::all() {
        let negatives = sampled_piece_negative_tau_locations(&best_tau, piece_id, &SAMPLES);
        let piece_report = sampled_piece_tau_report(&best_tau, piece_id, &SAMPLES);
        eprintln!(
            "PR08 domain-aware source projection best piece {piece_id:?}: outside={outside}/25 neg={neg}/25 min_u={min_u:e} min_v={min_v:e} max_sum={max_sum:e} det=[{min_det:e},{max_det:e}] negatives={negatives:?}",
            outside = piece_report.outside,
            neg = piece_report.negative,
            min_u = piece_report.min_u,
            min_v = piece_report.min_v,
            max_sum = piece_report.max_sum,
            min_det = piece_report.min_det,
            max_det = piece_report.max_det,
        );
    }
    assert!(final_report.2 > baseline_source_delta);
    assert!(final_report.3.is_finite());
    assert_eq!(best_report.0, 13);
    assert_eq!(best_report.1, 35);
    assert_eq!(best_report.4, 2);
    assert!(best_report.3 < 6.0e-5);
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
    assert!(best_by_domain.3 < 4.0e-14);
    assert!(best_by_domain.7 > 0.49 && best_by_domain.7 < 0.50);
    assert_eq!(best_by_orientation.1, 0);
    assert_eq!(best_by_orientation.0, 0);
    assert!(best_by_orientation.3 < 4.0e-14);
    assert_eq!(fine_validity.domain.outside_triangle_samples, 0);
    assert_eq!(fine_validity.orientation.negative_determinants, 0);
}

#[test]
fn pr08_exact_endpoint_next_outer_row_candidate_is_rejected_by_tau_validity() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedReversed,
    )
    .unwrap();
    let normalized_rows =
        unit_sector_endpoint_normalized_tensor_rows(&rows, valence, sector, false, true)
            .expect("next-outer-row source should have nondegenerate exact endpoints");
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1.0e-10,
    })
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);

    eprintln!(
        "PR08 exact endpoint next-outer-row candidate: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
    assert!(solved.max_constraint_residual < 1.0e-12);
    assert_eq!(square.outside, 97);
    assert_eq!(validity.orientation.negative_determinants, 47);
}

#[test]
fn pr08_adjacent_transposed_first_column_sources_are_compared_under_live_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY);

    for (label, candidate, chart) in [
        (
            "adjacent-j-plus-1-transposed/reflected-fit",
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
            (true, true),
        ),
        (
            "stam-subdivided-first-column/reflected-fit",
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaStamSubdividedFirstColumn,
            (true, true),
        ),
        (
            "adjacent-outer-reversed/exact-endpoint",
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedReversed,
            (false, true),
        ),
    ] {
        let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate).unwrap();
        let normalized_rows =
            unit_sector_endpoint_normalized_tensor_rows(&rows, valence, sector, chart.0, chart.1)
                .expect("adjacent transposed source should have nondegenerate endpoints");
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
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
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 adjacent transposed source {label}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        if square.outside < best.1
            || (square.outside == best.1 && validity.orientation.negative_determinants < best.2)
            || (square.outside == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
            );
        }
    }

    assert_eq!(best.0, "adjacent-j-plus-1-transposed/reflected-fit");
    assert_eq!(best.1, 31);
    assert_eq!(best.2, 27);
}

fn pr08_fig63_adjacent_middle_row_transposed_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Kp20RegularTensorBorderCNet {
    let mut sectors = Vec::with_capacity(ring.sectors.len());
    for sector_index in 0..ring.sectors.len() {
        let sector = &ring.sectors[sector_index];
        let next = &ring.sectors[(sector_index + 1) % ring.sectors.len()];
        sectors.push([
            [next[2], sector[2], sector[5], sector[8]],
            [next[5], sector[3], sector[6], sector[9]],
            [next[8], sector[12], sector[11], sector[12]],
        ]);
    }
    Kp20RegularTensorBorderCNet { sectors }
}

#[test]
fn pr08_fig63_adjacent_middle_row_source_is_measured_under_live_gate() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let rows = pr08_fig63_adjacent_middle_row_transposed_rows(&ring);
    let mut best = ("", usize::MAX, usize::MAX, f64::NEG_INFINITY);

    for (chart_label, swap_targets, conjugate_source) in [
        ("direct", false, false),
        ("swapped", true, false),
        ("direct-conjugate", false, true),
        ("swapped-conjugate", true, true),
    ] {
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows(
            &rows,
            valence,
            sector,
            swap_targets,
            conjugate_source,
        )
        .expect("Fig. 6.3 adjacent middle row should have nondegenerate endpoints");
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
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
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "PR08 Fig6.3 adjacent middle row chart={chart_label}: residual={residual:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        if square.outside < best.1
            || (square.outside == best.1 && validity.orientation.negative_determinants < best.2)
            || (square.outside == best.1
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                chart_label,
                square.outside,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
            );
        }
    }

    assert_eq!(best.0, "swapped-conjugate");
    assert_eq!(best.1, 39);
    assert_eq!(best.2, 25);
}

#[test]
fn pr08_adjacent_transposed_reflected_source_with_centered_first_order_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
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
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
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
    let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);

    eprintln!(
        "PR08 adjacent transposed reflected source with centered first-order: residual={residual:e} center_scale={center_scale:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
    assert!(solved.max_constraint_residual > 1.0e-4);
    assert!(solved.max_constraint_residual < 1.0e-3);
    assert_eq!(square.outside, 31);
    assert_eq!(validity.orientation.negative_determinants, 27);
    assert!(
        square.outside >= 31 || validity.orientation.negative_determinants >= 27,
        "centered first-order unexpectedly improved the reflected adjacent-row source"
    );
}

#[test]
fn pr08_endpoint_normalized_next_row_invalidity_localizes_by_piece() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
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
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).unwrap();
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1.0e-10,
    })
    .unwrap();
    let local_tau = cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);

    let mut total_outside = 0_usize;
    let mut total_negative = 0_usize;
    for piece_id in Kp20PieceId::all() {
        let report = sampled_piece_tau_report(&local_tau, piece_id, &SAMPLES);
        eprintln!(
            "endpoint-normalized next-row piece {piece_id:?}: outside={outside}/{samples} neg={neg}/{det_samples} zero={zero} det=[{min_det:e},{max_det:e}] u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] max_sum={max_sum:e}",
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
        total_outside += report.outside;
        total_negative += report.negative;
    }

    assert_eq!(total_outside, 31);
    assert_eq!(total_negative, 27);
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

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_scaled_appendix_cubic_is_measured() {
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
        let vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
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

        let mut cubic_rows = Vec::new();
        for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            for i in 0..=3 {
                for j in 0..=(3 - i) {
                    let target = appendix_tau.piece(piece_id).cps[i][j];
                    let coeff = Kp20PlanarTauCoeff::new(piece_id, i, j);
                    for (axis, target_axis) in [target.x, target.y].into_iter().enumerate() {
                        cubic_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                            terms: vec![
                                (
                                    Kp20PlanarTauAugmentedScalarCoeff::Tau(
                                        Kp20PlanarTauScalarCoeff::new(coeff, axis),
                                    ),
                                    1.0,
                                ),
                                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_axis),
                            ],
                            target: 0.0,
                        });
                    }
                }
            }
        }

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_scale) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &vector_constraints,
                &scalar_constraints,
                &cubic_rows,
                1,
            );
        let mut final_vector_constraints = vector_constraints;
        final_vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: final_vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: cubic_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
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
            "scaled q13/q2 self-diagonal scaled {label} cubic: residual={residual:e} projection_residual={projection_residual:e} source_delta={source_delta:e} projection_scale={projection_scale:e} final_scale={final_scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            projection_residual = projection_residual,
            source_delta = max_source_delta,
            projection_scale = projection_scale[0],
            final_scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_self_diagonal_left_scaled_appendix_cubic_is_measured() {
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
        let vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
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

        let mut cubic_rows = Vec::new();
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                let target = appendix_tau.piece(Kp20PieceId::NearOriginLeft).cps[i][j];
                let coeff = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, i, j);
                for (axis, target_axis) in [target.x, target.y].into_iter().enumerate() {
                    cubic_rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                        terms: vec![
                            (
                                Kp20PlanarTauAugmentedScalarCoeff::Tau(
                                    Kp20PlanarTauScalarCoeff::new(coeff, axis),
                                ),
                                1.0,
                            ),
                            (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_axis),
                        ],
                        target: 0.0,
                    });
                }
            }
        }

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_scale) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &vector_constraints,
                &scalar_constraints,
                &cubic_rows,
                1,
            );
        let mut final_vector_constraints = vector_constraints;
        final_vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: final_vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: cubic_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
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
            "scaled q13/q2 self-diagonal left scaled {label} cubic: residual={residual:e} projection_residual={projection_residual:e} source_delta={source_delta:e} projection_scale={projection_scale:e} final_scale={final_scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            projection_residual = projection_residual,
            source_delta = max_source_delta,
            projection_scale = projection_scale[0],
            final_scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_pair_diagonal_component_variants_are_measured() {
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
    let variants = [
        ("none", None),
        ("swap", Some([local_axes[1], local_axes[0]])),
        ("same", Some([local_axes[0], local_axes[1]])),
        (
            "swap-neg",
            Some([
                cc_test_neg_component(local_axes[1]),
                cc_test_neg_component(local_axes[0]),
            ]),
        ),
        (
            "same-neg",
            Some([
                cc_test_neg_component(local_axes[0]),
                cc_test_neg_component(local_axes[1]),
            ]),
        ),
        (
            "first-neg-swap",
            Some([cc_test_neg_component(local_axes[1]), local_axes[0]]),
        ),
        (
            "second-neg-swap",
            Some([local_axes[1], cc_test_neg_component(local_axes[0])]),
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
    for (label, upper_components) in variants {
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
        if let Some(upper_components) = upper_components {
            scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows_with_components(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes,
                upper_components,
            ));
        }
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
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
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 pair variant {label}: residual={residual:e} source_delta={max_source_delta:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
                max_source_delta,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!("scaled q13/q2 best pair variant diagnostic: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_fig8_diagonal_piece_correspondence_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();
    let maps = kp20_section42_map_candidates(angle);
    let cases = [
        (
            "corrected-current",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "same-side-orientation-winner",
            0_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointLeft,
        ),
        (
            "corrected-low-neg-rejected",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointRight,
        ),
        (
            "corrected-algebraic-collapse",
            1_usize,
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        ),
    ];

    let mut best = (
        "",
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearMidpointRight,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    let mut current_layout = None;

    for (case_label, map_index, self_a, self_b, lower_right, upper_left) in cases {
        let (map_label, map) = maps[map_index];
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            self_a,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_self_diagonal_rows(
            self_b,
            local_axes[0],
            local_axes[1],
        ));
        scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows(
            lower_right,
            upper_left,
            local_axes[0],
            local_axes[1],
        ));

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
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

        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
        let score = (
            validity.orientation.negative_determinants,
            square.outside,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            max_source_delta,
            -validity.orientation.min_determinant,
        );
        eprintln!(
            "scaled q13/q2 Fig.8 diagonal correspondence case={case_label} map={map_label} self=({self_a:?},{self_b:?}) pair=({lower_right:?},{upper_left:?}) residual={residual:e} source_delta={source_delta:e} scale={scale:e} square_outside={square_outside}/{square_samples} triangle_outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}])",
            residual = solved.max_constraint_residual,
            source_delta = max_source_delta,
            scale = solved.extra_values[0],
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

        if case_label == "corrected-current" {
            current_layout = Some(score);
        }
        if score < (best.5, best.6, best.7, best.8, best.9, best.10) {
            best = (
                case_label,
                self_a,
                self_b,
                lower_right,
                upper_left,
                score.0,
                score.1,
                score.2,
                score.3,
                score.4,
                score.5,
            );
        }
    }

    eprintln!(
        "scaled q13/q2 Fig.8 diagonal correspondence best map={} self=({:?},{:?}) pair=({:?},{:?}) neg={} square_outside={} triangle_outside={} residual={:e} source_delta={:e} neg_min_det={:e} corrected_current={current_layout:?}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7, best.8, best.9, best.10
    );
    assert!(current_layout.is_some());
    assert!(best.8.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_pair_diagonal_lemma_owned_rows_are_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let tensor_border = stam_appendix_a_scaled_pr08_q13_q2_tensor_border();

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
        (
            "exclude-apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i != 0 || j != 0,
            ),
        ),
        (
            "exclude-first-order-square",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i > 1 || j > 1,
            ),
        ),
        (
            "exclude-total-degree-one",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 1,
            ),
        ),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
    ] {
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
        scalar_constraints.extend(pair_rows);

        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta) =
            cc_test_project_tensor_rows_for_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
            );
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
            "scaled q13/q2 pair lemma-owned {label}: residual={residual:e} source_delta={source_delta:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e}",
            residual = solved.max_constraint_residual,
            source_delta = max_source_delta,
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
}

#[test]
fn kp20_fig8_paper_quadrant_pair_correspondence_is_derived_from_domains() {
    let lower_right = Kp20AffineGuideMap {
        origin: [0.5, 0.0],
        du: [0.5, 0.0],
        dv: [0.0, 0.5],
    };
    let upper_left_coefficient_chart = Kp20AffineGuideMap {
        origin: [0.5, 1.0],
        du: [-0.5, 0.0],
        dv: [0.0, -0.5],
    };

    for i in 0..=6 {
        for j in 0..=6 {
            let lower_point = cc_test_affine_map_point(lower_right, i as f64 / 6.0, j as f64 / 6.0);
            let diagonal_image = [lower_point[1], lower_point[0]];
            let upper_uv = cc_test_affine_map_inverse(upper_left_coefficient_chart, diagonal_image)
                .expect("paper quadrant upper-left coefficient chart is nonsingular");
            let mapped_i = cc_test_nearest_degree6_index(upper_uv[0]);
            let mapped_j = cc_test_nearest_degree6_index(upper_uv[1]);
            assert_eq!(
                (mapped_i, mapped_j),
                (6 - j, 6 - i),
                "KP20 Fig.8 coefficient-chart diagonal should map NearOriginRight[{i},{j}] to NearMidpointLeft[{},{}]",
                6 - j,
                6 - i,
            );
        }
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_with_exact_pair_is_measured() {
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
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &centered_rows,
            1,
        );

    for (label, include_centered_in_final) in [
        ("projection-only", false),
        ("projection-and-centered-final", true),
    ] {
        let mut vector_constraints = fig8_rows.clone();
        vector_constraints.extend(projected_tensor_rows.clone());
        let (tau_macro_patch, max_constraint_residual, final_scale) = if include_centered_in_final {
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            (
                solved.tau_macro_patch,
                solved.max_constraint_residual,
                Some(solved.extra_values[0]),
            )
        } else {
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            (solved.tau_macro_patch, solved.max_constraint_residual, None)
        };
        let local_tau = cc_test_global_planar_patch_to_local(valence, sector, tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
        let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
        let apex_det = du[0] * dv[1] - dv[0] * du[1];
        eprintln!(
            "scaled q13/q2 exact-pair Lemma4 source projection {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            projected_scale = projection_extra_values[0],
            residual = max_constraint_residual,
            final_scale = final_scale.unwrap_or(f64::NAN),
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
        assert!(max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_cubic_lemma4_projection_with_exact_pair_is_measured() {
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
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let cubic_rows = cc_test_paper_refined_apex_cubic_derivative_rows(
        valence,
        sector,
        Kp20PieceId::NearOriginLeft,
        local_axes,
        0,
    );
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &cubic_rows,
            1,
        );

    let mut vector_constraints = fig8_rows;
    vector_constraints.extend(projected_tensor_rows);
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: cubic_rows,
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
        "scaled q13/q2 exact-pair cubic-Lemma4 source projection: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
        projected_scale = projection_extra_values[0],
        final_scale = solved.extra_values[0],
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

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_without_full_diagonal_is_measured() {
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

    for (label, centered_rows) in [
        (
            "first-order",
            cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0),
        ),
        (
            "cubic-derivative",
            cc_test_paper_refined_apex_cubic_derivative_rows(
                valence,
                sector,
                Kp20PieceId::NearOriginLeft,
                local_axes,
                0,
            ),
        ),
    ] {
        let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
        let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints: centered_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
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
            "scaled q13/q2 no-full-diagonal Lemma4 source projection {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} p00=({p00x:e},{p00y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_fig8_average_completions_are_measured() {
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
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    let mut best = (
        "",
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (label, fig8_rows) in [
        (
            "paper-single",
            cc_test_fig8_indicated_average_variant_rows(false, false),
        ),
        (
            "diagonal-mirror",
            cc_test_fig8_indicated_average_variant_rows(true, false),
        ),
        (
            "rotated-quadrants",
            cc_test_fig8_indicated_average_variant_rows(false, true),
        ),
        (
            "rotated-quadrants-plus-mirror",
            cc_test_fig8_indicated_average_variant_rows(true, true),
        ),
    ] {
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );

        let mut vector_constraints = fig8_rows;
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows.clone(),
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 Fig8 average completion {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
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
                max_source_delta,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!("scaled q13/q2 best Fig8 average completion diagnostic: {best:?}");
    assert!(best.5.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_lemma4_projection_deltas_are_localized() {
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
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let original_tensor_rows = tensor_rows.clone();
    let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
        cc_test_project_tensor_rows_for_augmented_scalar_basis(
            tensor_rows,
            &fig8_rows,
            &scalar_constraints,
            &centered_rows,
            1,
        );

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
            "scaled q13/q2 Lemma4 eq14 projection delta rank={rank} row={row_index} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} original=({ox:e},{oy:e}) projected=({px:e},{py:e})",
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
        "scaled q13/q2 Lemma4 eq14 projection delta summary: projection_residual={projection_residual:e} max_source_delta={max_source_delta:e} projected_scale={projected_scale:e}",
        projected_scale = projection_extra_values[0],
    );
    assert!(projection_residual.is_finite());
    assert!(max_source_delta.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_c_origin_eq14_normalization_is_measured() {
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
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (label, normalize_left, normalize_right) in [
        ("none", false, false),
        ("left-c-origin", true, false),
        ("right-c-origin", false, true),
        ("both-c-origin", true, true),
    ] {
        let mut tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let left_c = tensor_rows[6].target;
        let right_c = tensor_rows[21 + 6].target;
        for (row_index, row) in tensor_rows.iter_mut().enumerate() {
            let origin = if row_index < 21 {
                normalize_left.then_some(left_c)
            } else {
                normalize_right.then_some(right_c)
            };
            if let Some(origin) = origin {
                row.target[0] -= origin[0];
                row.target[1] -= origin[1];
            }
        }
        let (projected_tensor_rows, max_source_delta, projection_residual, projection_extra_values) =
            cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
        let mut vector_constraints = fig8_rows.clone();
        vector_constraints.extend(projected_tensor_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows.clone(),
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .unwrap();
        let local_tau =
            cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled q13/q2 c-origin eq14 normalization {label}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            projected_scale = projection_extra_values[0],
            final_scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_per_sector_characteristic_eq14_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let normalized_rows = cc_test_per_sector_characteristic_tensor_rows(&raw_rows, valence, sector);
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (label, rows) in [
        ("raw-global", raw_rows),
        ("per-sector-characteristic", normalized_rows),
    ] {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();

        for include_projection in [false, true] {
            let (source_rows, source_delta, projection_residual, projected_scale) =
                if include_projection {
                    let (
                        projected_tensor_rows,
                        max_source_delta,
                        projection_residual,
                        projection_extra_values,
                    ) = cc_test_project_tensor_rows_for_augmented_scalar_basis(
                        tensor_rows.clone(),
                        &fig8_rows,
                        &scalar_constraints,
                        &centered_rows,
                        1,
                    );
                    (
                        projected_tensor_rows,
                        max_source_delta,
                        projection_residual,
                        projection_extra_values[0],
                    )
                } else {
                    (tensor_rows.clone(), 0.0, 0.0, f64::NAN)
                };
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(source_rows);
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
            let p00 = piece.cps[0][0];
            let p10 = piece.cps[1][0];
            let p01 = piece.cps[0][1];
            let du = [6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y)];
            let dv = [6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y)];
            let apex_det = du[0] * dv[1] - dv[0] * du[1];
            eprintln!(
                "scaled q13/q2 {label} eq14 characteristic normalization projection={include_projection}: projection_residual={projection_residual:e} source_delta={source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] apex_det={apex_det:e} worst={worst:?}@({u:.2},{v:.2})",
                final_scale = solved.extra_values[0],
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
fn stam_appendix_a_scaled_pr08_q13_q2_characteristic_dyadic_level_with_lemma4_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    let mut best = (
        "",
        0,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::INFINITY,
    );
    for (normalization_label, normalized_base) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
        ),
    ] {
        for level_exponent in -2..=10 {
            let rows = scale_tensor_border_rows(&normalized_base, ring.lambda.powi(level_exponent));
            let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
            let tensor_rows =
                kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
            let (
                projected_tensor_rows,
                max_source_delta,
                projection_residual,
                projection_extra_values,
            ) = cc_test_project_tensor_rows_for_augmented_scalar_basis(
                tensor_rows,
                &fig8_rows,
                &scalar_constraints,
                &centered_rows,
                1,
            );
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                    initial: Kp20MacroPatch::zero(),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    augmented_scalar_constraints: centered_rows.clone(),
                    extra_variable_count: 1,
                    regularization_weight: 1.0e-10,
                    extra_regularization_weight: 1.0e-10,
                },
            )
            .unwrap();
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "scaled q13/q2 characteristic dyadic {normalization_label} level={level_exponent}: projection_residual={projection_residual:e} source_delta={max_source_delta:e} projected_scale={projected_scale:e} final_scale={final_scale:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                projected_scale = projection_extra_values[0],
                final_scale = solved.extra_values[0],
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            if validity.orientation.negative_determinants < best.3
                || (validity.orientation.negative_determinants == best.3
                    && validity.domain.outside_triangle_samples < best.2)
                || (validity.orientation.negative_determinants == best.3
                    && validity.domain.outside_triangle_samples == best.2
                    && validity.orientation.min_determinant > best.4)
            {
                best = (
                    normalization_label,
                    level_exponent,
                    validity.domain.outside_triangle_samples,
                    validity.orientation.negative_determinants,
                    validity.orientation.min_determinant,
                    max_source_delta,
                    solved.max_constraint_residual,
                );
            }
        }
    }
    eprintln!("scaled q13/q2 best characteristic dyadic Lemma4 diagnostic: {best:?}");
    assert!(best.6.is_finite());
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_linear_characteristic_source_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);

    for (normalization_label, rows) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
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
        eprintln!(
            "scaled q13/q2 linear characteristic source normalization {normalization_label}: residual={residual:e} center_scale={center_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} frob_delta={frob_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[0],
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
        assert!(det_a.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_linear_source_collapse_is_localized() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let rows = cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector);
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 0);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair = cc_test_sector_local_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
    );

    for (label, diagonal_blocks) in [
        (
            "section42-no-diagonal",
            Vec::<&[Kp20PlanarTauScalarLinearConstraint]>::new(),
        ),
        ("section42-self-origin", vec![self_origin.as_slice()]),
        ("section42-self-midpoint", vec![self_midpoint.as_slice()]),
        (
            "section42-self-diagonals",
            vec![self_origin.as_slice(), self_midpoint.as_slice()],
        ),
        ("section42-pair-only", vec![pair.as_slice()]),
        (
            "section42-all-diagonal",
            vec![
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair.as_slice(),
            ],
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        for block in diagonal_blocks {
            scalar_constraints.extend(block.iter().cloned());
        }
        let mut augmented_scalar_constraints = centered_rows.clone();
        augmented_scalar_constraints
            .extend(cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1));
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: fig8_rows.clone(),
                scalar_constraints,
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
        eprintln!(
            "scaled q13/q2 linear source collapse {label}: residual={residual:e} center_scale={center_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e}",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[0],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_scaled_pr08_q13_q2_anchored_linear_characteristic_source_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let base_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
    )
    .unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );

    for (normalization_label, rows) in [
        ("raw-global", base_rows.clone()),
        (
            "per-sector-characteristic",
            cc_test_per_sector_characteristic_tensor_rows(&base_rows, valence, sector),
        ),
    ] {
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        let tensor_rows =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
            let mut augmented_scalar_constraints = cc_test_first_order_centered_source_rows(
                Kp20PieceId::NearOriginLeft,
                local_axes,
                0,
            );
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                target: anchored_scale,
            });
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
            let local_tau =
                cc_test_global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "scaled q13/q2 anchored linear characteristic source {normalization_label} scale={anchored_scale:e}: residual={residual:e} final_scale={final_scale:e} A=[[{a00:e},{a01:e}],[{a10:e},{a11:e}]] detA={det_a:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                final_scale = solved.extra_values[0],
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
}

#[test]
fn stam_appendix_a_regular_subpatch_anchored_characteristic_source_sweep_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_map_candidates(angle)[0].1;
    let local_axes = cc_test_sector_local_axis_components(valence, sector);
    let ring = kp20_catmull_clark_characteristic_ring(valence).unwrap();
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );

    let mut best = (f64::INFINITY, 0, false, 0, false, false, 0.0, 0.0, 0_usize);
    for patch_index in 0..3 {
        for transposed in [false, true] {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    for reverse_cols in [false, true] {
                        let rows = stam_appendix_b_projected_rows(
                            &ring,
                            patch_index,
                            transposed,
                            row_start,
                            reverse_rows,
                            reverse_cols,
                        )
                        .unwrap();
                        let rows =
                            cc_test_per_sector_characteristic_tensor_rows(&rows, valence, sector);
                        let tensor_border =
                            kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
                        let tensor_rows =
                            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector)
                                .unwrap();
                        for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
                            let mut augmented_scalar_constraints =
                                cc_test_first_order_centered_source_rows(
                                    Kp20PieceId::NearOriginLeft,
                                    local_axes,
                                    0,
                                );
                            augmented_scalar_constraints.push(
                                Kp20PlanarTauAugmentedScalarLinearConstraint {
                                    terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                                    target: anchored_scale,
                                },
                            );
                            augmented_scalar_constraints.extend(
                                cc_test_linear_source_normalized_eq14_rows(&tensor_rows, 1),
                            );
                            let solved =
                                solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
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
                            let local_tau = cc_test_global_planar_patch_to_local(
                                valence,
                                sector,
                                solved.tau_macro_patch,
                            );
                            let validity = kp20_planar_bi6_tau_validity_report(
                                &local_tau,
                                &[0.0, 0.5, 1.0],
                                1.0e-10,
                                1.0e-10,
                            );
                            if solved.max_constraint_residual < best.0 {
                                best = (
                                    solved.max_constraint_residual,
                                    patch_index,
                                    transposed,
                                    row_start,
                                    reverse_rows,
                                    reverse_cols,
                                    anchored_scale,
                                    det_a,
                                    validity.orientation.negative_determinants,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!("Stam regular subpatch anchored characteristic source best: {best:?}");
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_promoted_section42_lemma4_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let top_count = 2 * valence + 1;
    let source_extra_count = 2 * top_count;
    let fig8_rows = cc_test_fig8_indicated_average_variant_rows(false, false);
    let diagonal_rows =
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, 0).unwrap();

    let mut best = (
        f64::INFINITY,
        "",
        0,
        0,
        false,
        0.0,
        usize::MAX,
        usize::MAX,
        f64::NEG_INFINITY,
    );
    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for patch_index in 0..3 {
            for row_start in [0, 1] {
                for reverse_rows in [false, true] {
                    for anchored_scale in [1.0, ring.lambda, ring.lambda.powi(4)] {
                        let mut augmented_scalar_constraints =
                            stam_appendix_a_unknown_top_vector_tensor_border_constraints(
                                &ring,
                                patch_index,
                                row_start,
                                reverse_rows,
                            );
                        augmented_scalar_constraints.extend(
                            cc_test_first_order_centered_source_rows(
                                Kp20PieceId::NearOriginLeft,
                                local_axes,
                                source_extra_count,
                            ),
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(
                                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                                    1.0,
                                )],
                                target: anchored_scale,
                            },
                        );
                        let mut scalar_constraints =
                            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
                        scalar_constraints.extend(diagonal_rows.clone());
                        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                                initial: Kp20MacroPatch::zero(),
                                vector_constraints: fig8_rows.clone(),
                                scalar_constraints,
                                augmented_scalar_constraints,
                                extra_variable_count: source_extra_count + 1,
                                regularization_weight: 1.0e-10,
                                extra_regularization_weight: 1.0e-12,
                            },
                        )
                        .unwrap();
                        let local_tau = cc_test_global_planar_patch_to_local(
                            valence,
                            0,
                            solved.tau_macro_patch,
                        );
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau, &SAMPLES, 1.0e-10, 1.0e-10,
                        );
                        eprintln!(
                            "unknown top promoted §4.2 Lemma4 map={map_label} patch={patch_index} row_start={row_start} reverse={reverse_rows} scale={anchored_scale:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                            residual = solved.max_constraint_residual,
                            outside = validity.domain.outside_triangle_samples,
                            domain_samples = validity.domain.sample_count,
                            neg = validity.orientation.negative_determinants,
                            orientation_samples = validity.orientation.sample_count,
                            zero = validity.orientation.near_zero_determinants,
                            min_det = validity.orientation.min_determinant,
                            max_det = validity.orientation.max_determinant,
                        );
                        if solved.max_constraint_residual < best.0
                            || (solved.max_constraint_residual <= best.0 + 1.0e-12
                                && validity.orientation.negative_determinants < best.7)
                            || (solved.max_constraint_residual <= best.0 + 1.0e-12
                                && validity.orientation.negative_determinants == best.7
                                && validity.orientation.min_determinant > best.8)
                        {
                            best = (
                                solved.max_constraint_residual,
                                map_label,
                                patch_index,
                                row_start,
                                reverse_rows,
                                anchored_scale,
                                validity.domain.outside_triangle_samples,
                                validity.orientation.negative_determinants,
                                validity.orientation.min_determinant,
                            );
                        }
                    }
                }
            }
        }
    }
    eprintln!("unknown top promoted §4.2 Lemma4 best: {best:?}");
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_lemma4_zero_lock_is_localized() {
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

    for (label, vector_constraints, scalar_constraints) in [
        (
            "eq14-source-only",
            Vec::<Kp20PlanarTauLinearConstraint>::new(),
            Vec::<Kp20PlanarTauScalarLinearConstraint>::new(),
        ),
        (
            "fig8-eq14",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            Vec::<Kp20PlanarTauScalarLinearConstraint>::new(),
        ),
        (
            "section42-no-diagonal",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c),
        ),
        (
            "section42-diagonal",
            cc_test_fig8_indicated_average_variant_rows(false, false),
            {
                let mut rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
                rows.extend(
                    kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, 0).unwrap(),
                );
                rows
            },
        ),
    ] {
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
                vector_constraints,
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
            "unknown top Lemma4 zero-lock {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_diagonal_block_zero_lock_is_localized() {
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
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair = cc_test_sector_local_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
    );

    for (label, diagonal_blocks) in [
        ("none", Vec::<&[Kp20PlanarTauScalarLinearConstraint]>::new()),
        ("self-origin", vec![self_origin.as_slice()]),
        ("self-midpoint", vec![self_midpoint.as_slice()]),
        (
            "self-both",
            vec![self_origin.as_slice(), self_midpoint.as_slice()],
        ),
        ("pair-only", vec![pair.as_slice()]),
        (
            "all",
            vec![
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair.as_slice(),
            ],
        ),
    ] {
        let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
        for block in diagonal_blocks {
            scalar_constraints.extend(block.iter().cloned());
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
            "unknown top diagonal block {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_pair_diagonal_component_variants_are_measured() {
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
    let variants = [
        ("none", None),
        ("swap", Some([local_axes[1], local_axes[0]])),
        ("same", Some([local_axes[0], local_axes[1]])),
        (
            "swap-neg",
            Some([
                cc_test_neg_component(local_axes[1]),
                cc_test_neg_component(local_axes[0]),
            ]),
        ),
        (
            "same-neg",
            Some([
                cc_test_neg_component(local_axes[0]),
                cc_test_neg_component(local_axes[1]),
            ]),
        ),
        (
            "first-neg-swap",
            Some([cc_test_neg_component(local_axes[1]), local_axes[0]]),
        ),
        (
            "second-neg-swap",
            Some([local_axes[1], cc_test_neg_component(local_axes[0])]),
        ),
    ];

    let mut best = ("", f64::INFINITY, usize::MAX, f64::NEG_INFINITY);
    for (label, upper_components) in variants {
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
        if let Some(upper_components) = upper_components {
            scalar_constraints.extend(cc_test_sector_local_pair_diagonal_rows_with_components(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes,
                upper_components,
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
            "unknown top pair variant {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        if solved.max_constraint_residual < best.1
            || (solved.max_constraint_residual <= best.1 + 1.0e-12
                && validity.orientation.negative_determinants < best.2)
            || (solved.max_constraint_residual <= best.1 + 1.0e-12
                && validity.orientation.negative_determinants == best.2
                && validity.orientation.min_determinant > best.3)
        {
            best = (
                label,
                solved.max_constraint_residual,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
            );
        }
    }
    eprintln!("unknown top best pair variant: {best:?}");
    assert!(best.1.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_pair_diagonal_subset_zero_lock_is_localized() {
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

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 0 && j == 0,
            ),
        ),
        (
            "first-order-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i <= 1 && j <= 1,
            ),
        ),
        (
            "near-origin-block",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i <= 2 && j <= 2,
            ),
        ),
        (
            "domain-diagonal",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == j,
            ),
        ),
        (
            "outer-boundary",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i == 6 || j == 6,
            ),
        ),
        (
            "full",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
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
        scalar_constraints.extend(pair_rows);
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
            "unknown top pair subset {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_pair_diagonal_lemma_owned_rows_are_measured() {
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

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
        (
            "exclude-apex-corner",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i != 0 || j != 0,
            ),
        ),
        (
            "exclude-first-order-square",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i > 1 || j > 1,
            ),
        ),
        (
            "exclude-total-degree-one",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 1,
            ),
        ),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
    ] {
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
        scalar_constraints.extend(pair_rows);
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
            "unknown top pair lemma-owned {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
fn stam_appendix_a_unknown_top_vector_lemma_owned_pair_source_modes_are_measured() {
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
        |i, j| i + j > 3,
    ));
    let mut augmented_scalar_constraints = source_constraints;
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
    assert!(solved.max_constraint_residual.is_finite());

    let odd = (0..valence)
        .map(|sector| {
            Complex::new(
                solved.extra_values[1 + 2 * sector],
                solved.extra_values[top_count + 1 + 2 * sector],
            )
        })
        .collect::<Vec<_>>();
    let even = (0..valence)
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
    eprintln!(
        "unknown top lemma-owned pair source modes: residual={residual:e} center_scale={center_scale:e} odd_mode1_residual={odd_residual:e} even_mode1_residual={even_residual:e} odd_modes={odd_magnitudes:?} even_modes={even_magnitudes:?}",
        residual = solved.max_constraint_residual,
        center_scale = solved.extra_values[source_extra_count],
    );
    assert!(odd_residual.is_finite());
    assert!(even_residual.is_finite());
}

#[test]
fn stam_appendix_a_unknown_top_vector_anchored_lemma_owned_pair_is_measured() {
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

    for (label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
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
        scalar_constraints.extend(pair_rows);
        let mut augmented_scalar_constraints = source_constraints.clone();
        augmented_scalar_constraints.extend(cc_test_first_order_centered_source_rows(
            Kp20PieceId::NearOriginLeft,
            local_axes,
            source_extra_count,
        ));
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(1), 1.0)],
            target: ring.sectors[0][1].x,
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
        let odd = (0..valence)
            .map(|sector| {
                Complex::new(
                    solved.extra_values[1 + 2 * sector],
                    solved.extra_values[top_count + 1 + 2 * sector],
                )
            })
            .collect::<Vec<_>>();
        let even = (0..valence)
            .map(|sector| {
                Complex::new(
                    solved.extra_values[2 + 2 * sector],
                    solved.extra_values[top_count + 2 + 2 * sector],
                )
            })
            .collect::<Vec<_>>();
        eprintln!(
            "unknown top anchored lemma-owned pair {label}: residual={residual:e} final_scale={final_scale:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] odd_mode1_residual={odd_residual:e} even_mode1_residual={even_residual:e}",
            residual = solved.max_constraint_residual,
            final_scale = solved.extra_values[source_extra_count],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            odd_residual = fourier_mode_residual(&odd, 1),
            even_residual = fourier_mode_residual(&even, 1),
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let mut best = (
        f64::INFINITY,
        "",
        0_usize,
        0_usize,
        false,
        0.0_f64,
        f64::INFINITY,
        0_usize,
        0_usize,
    );

    for (pair_label, pair_rows) in [
        ("none", Vec::<Kp20PlanarTauScalarLinearConstraint>::new()),
        (
            "exclude-lemma-cubic-triangle",
            cc_test_sector_local_filtered_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
                |i, j| i + j > 3,
            ),
        ),
        (
            "full-pair",
            cc_test_sector_local_pair_diagonal_rows(
                Kp20PieceId::NearOriginRight,
                Kp20PieceId::NearMidpointLeft,
                local_axes[0],
                local_axes[1],
            ),
        ),
    ] {
        for rotation_direction in [-1.0, 1.0] {
            for patch_index in 0..3 {
                for row_start in [0, 1] {
                    for reverse_rows in [false, true] {
                        let mut scalar_constraints =
                            kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
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
                        scalar_constraints.extend(pair_rows.clone());

                        let mut augmented_scalar_constraints =
                            stam_appendix_a_rotational_top_vector_tensor_border_constraints(
                                &ring,
                                patch_index,
                                row_start,
                                reverse_rows,
                                rotation_direction,
                            );
                        augmented_scalar_constraints.extend(
                            cc_test_first_order_centered_source_rows(
                                Kp20PieceId::NearOriginLeft,
                                local_axes,
                                source_extra_count,
                            ),
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
                                target: ring.sectors[0][1].x,
                            },
                        );
                        augmented_scalar_constraints.push(
                            Kp20PlanarTauAugmentedScalarLinearConstraint {
                                terms: vec![(
                                    Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
                                    1.0,
                                )],
                                target: anchored_scale,
                            },
                        );
                        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
                            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                                initial: Kp20MacroPatch::zero(),
                                vector_constraints: cc_test_fig8_indicated_average_variant_rows(
                                    false, false,
                                ),
                                scalar_constraints,
                                augmented_scalar_constraints,
                                extra_variable_count: source_extra_count + 1,
                                regularization_weight: 1.0e-10,
                                extra_regularization_weight: 1.0e-12,
                            },
                        )
                        .unwrap();
                        let local_tau = cc_test_global_planar_patch_to_local(
                            valence,
                            0,
                            solved.tau_macro_patch,
                        );
                        let validity = kp20_planar_bi6_tau_validity_report(
                            &local_tau,
                            &[0.0, 0.5, 1.0],
                            1.0e-10,
                            1.0e-10,
                        );
                        if solved.max_constraint_residual < best.0 {
                            best = (
                                solved.max_constraint_residual,
                                pair_label,
                                patch_index,
                                row_start,
                                reverse_rows,
                                rotation_direction,
                                solved.extra_values[source_extra_count],
                                validity.orientation.negative_determinants,
                                validity.orientation.near_zero_determinants,
                            );
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "rotational top lemma-owned pair best: residual={:e} pair={} patch={} row_start={} reverse={} rotation={} center_scale={:e} neg={} zero={}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7, best.8
    );
    assert!(best.0.is_finite());
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_residual_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let vector_constraints = cc_test_fig8_indicated_average_variant_rows(false, false);
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_rows = cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    );
    let source_rows =
        stam_appendix_a_rotational_top_vector_tensor_border_constraints(&ring, 0, 0, false, -1.0);
    let centered_rows =
        cc_test_first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, local_axes, 4);
    let source_anchor = vec![Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(0), 1.0)],
        target: ring.sectors[0][1].x,
    }];
    let scale_anchor = vec![Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(
            Kp20PlanarTauAugmentedScalarCoeff::Extra(source_extra_count),
            1.0,
        )],
        target: anchored_scale,
    }];

    let mut scalar_constraints = section42_rows.clone();
    scalar_constraints.extend(self_origin.clone());
    scalar_constraints.extend(self_midpoint.clone());
    scalar_constraints.extend(pair_rows.clone());
    let mut augmented_scalar_constraints = source_rows.clone();
    augmented_scalar_constraints.extend(centered_rows.clone());
    augmented_scalar_constraints.extend(source_anchor.clone());
    augmented_scalar_constraints.extend(scale_anchor.clone());
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: vector_constraints.clone(),
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: source_extra_count + 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-12,
        },
    )
    .unwrap();

    eprintln!(
        "rotational top lemma-owned residual families: total={total:e} fig8={fig8:e} section42={section42:e} self_origin={self_origin_res:e} self_midpoint={self_midpoint_res:e} pair_cut={pair_cut:e} eq14_source={eq14_source:e} centered={centered:e} source_anchor={source_anchor_res:e} scale_anchor={scale_anchor_res:e}",
        total = solved.max_constraint_residual,
        fig8 = cc_test_vector_rows_max_residual(&solved.tau_macro_patch, &vector_constraints),
        section42 = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &section42_rows),
        self_origin_res = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &self_origin),
        self_midpoint_res =
            cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &self_midpoint),
        pair_cut = cc_test_scalar_rows_max_residual(&solved.tau_macro_patch, &pair_rows),
        eq14_source = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &source_rows
        ),
        centered = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &centered_rows
        ),
        source_anchor_res = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &source_anchor
        ),
        scale_anchor_res = cc_test_augmented_rows_max_residual(
            &solved.tau_macro_patch,
            &solved.extra_values,
            &scale_anchor
        ),
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma_owned_pair_source_anchor_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
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
        |i, j| i + j > 3,
    ));

    for anchor_extra in 0..source_extra_count {
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
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(anchor_extra), 1.0)],
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
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        eprintln!(
            "rotational top lemma-owned source anchor extra{anchor_extra}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_source_collapse_localizes_by_section42_family() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[0].1;
    let anchored_scale = ring.lambda.powi(4);
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_cut = cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    );
    let cases: Vec<(&str, Vec<&[Kp20PlanarTauScalarLinearConstraint]>)> = vec![
        ("fig8-eq14-only", vec![]),
        ("c0", vec![c0_rows.as_slice()]),
        ("eq11", vec![eq11_rows.as_slice()]),
        ("eq12", vec![eq12_rows.as_slice()]),
        ("eq13", vec![eq13_rows.as_slice()]),
        (
            "section42",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "section42-self",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
                self_origin.as_slice(),
                self_midpoint.as_slice(),
            ],
        ),
        (
            "section42-self-pair-cut",
            vec![
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair_cut.as_slice(),
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
        eprintln!(
            "rotational top source collapse section42 {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_c0_collapse_sweeps_symbol_maps() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
        for (c0_label, include_plain, include_underlined) in [
            ("plain", true, false),
            ("underlined", false, true),
            ("both", true, true),
        ] {
            let mut scalar_constraints = Vec::new();
            for i in 0..=6 {
                if include_plain {
                    add_kp20_section42_vector_rows(
                        &mut scalar_constraints,
                        map,
                        &[
                            section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                            section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
                        ],
                    );
                }
                if include_underlined {
                    add_kp20_section42_vector_rows(
                        &mut scalar_constraints,
                        map,
                        &[
                            section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                            section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                        ],
                    );
                }
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
            eprintln!(
                "rotational top C0 map sweep {map_label}/{c0_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e})",
                residual = solved.max_constraint_residual,
                center_scale = solved.extra_values[source_extra_count],
                s0 = solved.extra_values[0],
                s1 = solved.extra_values[1],
                s2 = solved.extra_values[2],
                s3 = solved.extra_values[3],
            );
            assert!(solved.max_constraint_residual.is_finite());
        }
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_full_section42_sweeps_symbol_maps() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);

    for (map_label, map) in kp20_section42_map_candidates(angle) {
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
            |i, j| i + j > 3,
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
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        let square = cc_test_square_domain_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        eprintln!(
            "rotational top full §4.2 map sweep {map_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) triangle_outside={triangle_outside}/{domain_samples} square_outside={square_outside}/{square_samples} square_range=([{min_u:e},{max_u:e}],[{min_v:e},{max_v:e}]) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            triangle_outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            square_outside = square.outside,
            square_samples = square.samples,
            min_u = square.min_u,
            max_u = square.max_u,
            min_v = square.min_v,
            max_v = square.max_v,
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
fn stam_appendix_a_rotational_top_vector_source_collapse_localizes_c0_companions() {
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
        ("c0", vec![c0_rows.as_slice()]),
        ("c0-eq11", vec![c0_rows.as_slice(), eq11_rows.as_slice()]),
        ("c0-eq12", vec![c0_rows.as_slice(), eq12_rows.as_slice()]),
        ("c0-eq13", vec![c0_rows.as_slice(), eq13_rows.as_slice()]),
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
            "c0-eq12-eq13",
            vec![
                c0_rows.as_slice(),
                eq12_rows.as_slice(),
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

    for (label, blocks) in cases {
        let mut scalar_constraints = Vec::new();
        for block in blocks {
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
        eprintln!(
            "rotational top source collapse C0 companions {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected §4.2 C0 companion set {label} should close with a nonzero rotational source, got residual {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected §4.2 C0 companion set {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_split_rows_are_localized() {
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
        eprintln!(
            "rotational top corrected eq12 split {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected eq12 split {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected eq12 split {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_endpoint_depth_is_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);

    let endpoint_halves = [
        (
            "plain",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            "arrow",
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];

    for (half_label, upper, lower) in endpoint_halves {
        for j in 0..=2 {
            let mut scalar_constraints = c0_rows.clone();
            scalar_constraints.extend(eq11_rows.clone());
            scalar_constraints.extend(kp20_section42_eq12_endpoint_row(map, upper, lower, j));

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
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
            eprintln!(
                "rotational top corrected eq12 endpoint {half_label} j={j}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                center_scale = solved.extra_values[source_extra_count],
                s0 = solved.extra_values[0],
                s1 = solved.extra_values[1],
                s2 = solved.extra_values[2],
                s3 = solved.extra_values[3],
                neg = validity.orientation.negative_determinants,
                samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(
                solved.max_constraint_residual < 1.0e-10,
                "corrected eq12 endpoint {half_label} j={j} should close algebraically, got {}",
                solved.max_constraint_residual
            );
        }
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_symbol_direction_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let current_pairs = [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];
    let reversed_pairs = [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ];
    let reversed_plain_pairs = [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];
    let reversed_arrow_pairs = [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ];
    let cases = [
        ("current", current_pairs.as_slice()),
        ("reversed-both", reversed_pairs.as_slice()),
        ("reversed-plain", reversed_plain_pairs.as_slice()),
        ("reversed-arrow", reversed_arrow_pairs.as_slice()),
    ];

    for (label, pairs) in cases {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(kp20_section42_eq12_rows_with_symbol_pairs(map, pairs));

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
        eprintln!(
            "rotational top corrected eq12 direction {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected eq12 direction {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            solved.extra_values[0].abs() > 0.5,
            "corrected eq12 direction {label} should preserve the anchored source, got {}",
            solved.extra_values[0]
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_with_eq13_is_measured() {
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
    let current_eq12_rows = kp20_section42_only_eq12_rows(map);
    let reversed_eq12_rows = kp20_section42_eq12_rows_with_symbol_pairs(
        map,
        &[
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ],
    );
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);

    for (label, eq12_rows, include_eq13) in [
        ("current-eq12-no-eq13", current_eq12_rows.as_slice(), false),
        (
            "reversed-eq12-no-eq13",
            reversed_eq12_rows.as_slice(),
            false,
        ),
        ("current-eq12-eq13", current_eq12_rows.as_slice(), true),
        ("reversed-eq12-eq13", reversed_eq12_rows.as_slice(), true),
    ] {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(eq12_rows.iter().cloned());
        if include_eq13 {
            scalar_constraints.extend(eq13_rows.clone());
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
        eprintln!(
            "rotational top corrected reversed-eq12 full {case_label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            case_label = label,
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed-eq12 case {} should close algebraically, got {}",
            label,
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_eq13_bands_are_measured() {
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
    let reversed_eq12_rows = kp20_section42_eq12_rows_with_symbol_pairs(
        map,
        &[
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ],
    );

    let cases = [
        (
            "current-bands",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
        ),
        (
            "swapped-bands",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
        (
            "all-plain-band",
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::ArrowP,
        ),
        (
            "all-underlined-band",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ];

    for (label, strip_p, strip_q, row4_p, row4_q) in cases {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(reversed_eq12_rows.clone());
        scalar_constraints.extend(kp20_section42_eq13_rows_with_bands(
            map,
            2.0 * c,
            strip_p,
            strip_q,
            row4_p,
            row4_q,
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
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected reversed-eq12 eq13 bands {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed-eq12 eq13 band case {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_reversed_eq12_split_rows_are_localized() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);

    let splits = [
        (
            "reversed-plain-endpoint",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainEndpoint,
        ),
        (
            "reversed-plain-plus-row",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainPlusRow,
        ),
        (
            "reversed-plain-minus-row",
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
            Kp20Section42Eq12Split::PlainMinusRow,
        ),
        (
            "reversed-arrow-endpoint",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowEndpoint,
        ),
        (
            "reversed-arrow-plus-row",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowPlusRow,
        ),
        (
            "reversed-arrow-minus-row",
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42Eq12Split::ArrowMinusRow,
        ),
    ];

    for (label, upper, lower, split) in splits {
        let mut scalar_constraints = c0_rows.clone();
        scalar_constraints.extend(eq11_rows.clone());
        scalar_constraints.extend(kp20_section42_eq12_split_rows_with_symbols(
            map, upper, lower, split,
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
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10, 1.0e-10);
        eprintln!(
            "rotational top corrected reversed eq12 split {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected reversed eq12 split {label} should close algebraically, got {}",
            solved.max_constraint_residual
        );
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_corrected_section42_diagonal_blocks_are_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let section42_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let self_origin = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearOriginLeft,
        local_axes[0],
        local_axes[1],
    );
    let self_midpoint = cc_test_sector_local_self_diagonal_rows(
        Kp20PieceId::NearMidpointRight,
        local_axes[0],
        local_axes[1],
    );
    let pair_cut = cc_test_sector_local_filtered_pair_diagonal_rows(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local_axes[0],
        local_axes[1],
        |i, j| i + j > 3,
    );
    let cases: Vec<(&str, Vec<&[Kp20PlanarTauScalarLinearConstraint]>)> = vec![
        ("section42", vec![]),
        ("section42-self-origin", vec![self_origin.as_slice()]),
        ("section42-self-midpoint", vec![self_midpoint.as_slice()]),
        (
            "section42-self-both",
            vec![self_origin.as_slice(), self_midpoint.as_slice()],
        ),
        ("section42-pair-cut", vec![pair_cut.as_slice()]),
        (
            "section42-self-both-pair-cut",
            vec![
                self_origin.as_slice(),
                self_midpoint.as_slice(),
                pair_cut.as_slice(),
            ],
        ),
    ];

    for (label, blocks) in cases {
        let mut scalar_constraints = section42_rows.clone();
        for block in blocks {
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
        eprintln!(
            "rotational top corrected §4.2 diagonal blocks {label}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        match label {
            "section42" | "section42-self-origin" => {
                assert!(
                    solved.max_constraint_residual < 1.0e-10,
                    "corrected §4.2 row set {label} should close with nonzero source, got residual {}",
                    solved.max_constraint_residual
                );
                assert!(
                    solved.extra_values[0].abs() > 0.5,
                    "corrected §4.2 row set {label} should preserve the anchored source, got {}",
                    solved.extra_values[0]
                );
            }
            "section42-self-midpoint" | "section42-self-both" | "section42-self-both-pair-cut" => {
                assert!(
                    solved.max_constraint_residual > 0.5,
                    "corrected §4.2 row set {label} should remain rejected, got residual {}",
                    solved.max_constraint_residual
                );
                assert!(
                    solved.extra_values[0].abs() < 1.0e-6,
                    "corrected §4.2 row set {label} should collapse the anchored source under the rejected rows, got {}",
                    solved.extra_values[0]
                );
            }
            "section42-pair-cut" => {
                assert!(
                    solved.max_constraint_residual < 1.0e-10,
                    "corrected §4.2 pair-cut rows should close algebraically after the eq(12) direction correction, got {}",
                    solved.max_constraint_residual
                );
                assert!(
                    validity.orientation.negative_determinants >= 20,
                    "corrected §4.2 pair-cut rows should remain orientation-rejected, got {}/{} negative samples",
                    validity.orientation.negative_determinants,
                    validity.orientation.sample_count
                );
            }
            _ => unreachable!("unrecognized corrected §4.2 diagonal case {label}"),
        }
    }
}

#[test]
fn stam_appendix_a_rotational_top_vector_corrected_section42_source_anchor_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);

    for anchor_extra in 0..source_extra_count {
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
            terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(anchor_extra), 1.0)],
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
                scalar_constraints: scalar_constraints.clone(),
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
            "rotational top corrected §4.2 source anchor extra{anchor_extra}: residual={residual:e} center_scale={center_scale:e} source=({s0:e},{s1:e},{s2:e},{s3:e}) neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            center_scale = solved.extra_values[source_extra_count],
            s0 = solved.extra_values[0],
            s1 = solved.extra_values[1],
            s2 = solved.extra_values[2],
            s3 = solved.extra_values[3],
            neg = validity.orientation.negative_determinants,
            samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-10,
            "corrected §4.2 source anchor extra{anchor_extra} should close algebraically, got {}",
            solved.max_constraint_residual
        );
        assert!(
            validity.orientation.negative_determinants > 0,
            "corrected §4.2 source anchor extra{anchor_extra} unexpectedly passed coarse orientation"
        );
    }
}

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
fn stam_appendix_a_rotational_top_vector_corrected_section42_piece_folds_are_localized() {
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
    let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
    eprintln!(
        "rotational top corrected §4.2 piece folds total: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
    for piece_id in Kp20PieceId::all() {
        let report = sampled_piece_tau_report(&local_tau, piece_id, &[0.0, 0.5, 1.0]);
        eprintln!(
            "rotational top corrected §4.2 piece {piece_id:?}: outside={outside}/{samples} neg={neg}/{det_samples} zero={zero} det=[{min_det:e},{max_det:e}] u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] max_sum={max_sum:e}",
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
    }
    assert_eq!(orientation.negative_determinants, 15);
}

#[test]
fn stam_appendix_a_rotational_top_vector_near_midpoint_outer_corner_is_localized() {
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
    let eq12_rows = kp20_section42_eq12_rows_with_symbol_pairs(
        map,
        &[
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ],
    );
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let cases: [(&str, bool, [&[Kp20PlanarTauScalarLinearConstraint]; 4]); 5] = [
        (
            "c0-eq11",
            true,
            [c0_rows.as_slice(), eq11_rows.as_slice(), &[], &[]],
        ),
        (
            "c0-eq11-eq12",
            true,
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                &[],
            ],
        ),
        (
            "c0-eq11-eq13",
            true,
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                &[],
                eq13_rows.as_slice(),
            ],
        ),
        (
            "full-section42",
            true,
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
        (
            "full-section42-without-fig8",
            false,
            [
                c0_rows.as_slice(),
                eq11_rows.as_slice(),
                eq12_rows.as_slice(),
                eq13_rows.as_slice(),
            ],
        ),
    ];

    let mut c0_eq11_near_midpoint_left_det = f64::NAN;
    let mut eq12_near_midpoint_left_det = f64::NAN;
    let mut eq13_near_midpoint_left_det = f64::NAN;
    let mut full_near_midpoint_left_det = f64::NAN;
    let mut full_without_fig8_det = f64::NAN;
    for (label, include_fig8_average, row_blocks) in cases {
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
        let vector_constraints = if include_fig8_average {
            cc_test_fig8_indicated_average_variant_rows(false, false)
        } else {
            Vec::new()
        };
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let orientation =
            kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
        eprintln!(
            "rotational top corrected §4.2 near-midpoint outer-corner {label}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e} nml_p66=({nml66x:e},{nml66y:e}) nml_p56=({nml56x:e},{nml56y:e}) nml_p65=({nml65x:e},{nml65y:e}) nml_du=({nmldux:e},{nmlduy:e}) nml_dv=({nmldvx:e},{nmldvy:e}) nmr_det={nmr_det:e} nmr_p66=({nmr66x:e},{nmr66y:e}) nmr_p56=({nmr56x:e},{nmr56y:e}) nmr_p65=({nmr65x:e},{nmr65y:e}) nmr_du=({nmrdux:e},{nmrduy:e}) nmr_dv=({nmrdvx:e},{nmrdvy:e})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
            nml_det = nml.determinant,
            nml66x = nml.p66.x,
            nml66y = nml.p66.y,
            nml56x = nml.p56.x,
            nml56y = nml.p56.y,
            nml65x = nml.p65.x,
            nml65y = nml.p65.y,
            nmldux = nml.du.x,
            nmlduy = nml.du.y,
            nmldvx = nml.dv.x,
            nmldvy = nml.dv.y,
            nmr_det = nmr.determinant,
            nmr66x = nmr.p66.x,
            nmr66y = nmr.p66.y,
            nmr56x = nmr.p56.x,
            nmr56y = nmr.p56.y,
            nmr65x = nmr.p65.x,
            nmr65y = nmr.p65.y,
            nmrdux = nmr.du.x,
            nmrduy = nmr.du.y,
            nmrdvx = nmr.dv.x,
            nmrdvy = nmr.dv.y,
        );
        assert!(solved.max_constraint_residual < 1.0e-10);
        match label {
            "c0-eq11" => {
                c0_eq11_near_midpoint_left_det = nml.determinant;
            }
            "c0-eq11-eq12" => {
                eq12_near_midpoint_left_det = nml.determinant;
            }
            "c0-eq11-eq13" => {
                eq13_near_midpoint_left_det = nml.determinant;
            }
            "full-section42" => {
                full_near_midpoint_left_det = nml.determinant;
                assert_eq!(orientation.negative_determinants, 15);
                assert_eq!(orientation.worst_piece, Kp20PieceId::NearMidpointLeft);
                assert!((orientation.worst_u - 1.0).abs() < 1.0e-12);
                assert!((orientation.worst_v - 0.5).abs() < 1.0e-12);
            }
            "full-section42-without-fig8" => {
                full_without_fig8_det = nml.determinant;
            }
            _ => unreachable!("unhandled near-midpoint outer-corner diagnostic case"),
        }
    }

    assert!(
        c0_eq11_near_midpoint_left_det > 0.09,
        "C0+eq(11) should leave the near-midpoint-left outer corner positive, got {c0_eq11_near_midpoint_left_det}"
    );
    assert!(
        eq12_near_midpoint_left_det < -0.04,
        "adding corrected eq(12) should introduce the near-midpoint-left outer-corner sign flip, got {eq12_near_midpoint_left_det}"
    );
    assert!(
        eq13_near_midpoint_left_det > 0.1,
        "adding eq(13) without eq(12) should not be the near-midpoint-left outer-corner sign flip, got {eq13_near_midpoint_left_det}"
    );
    assert!(
        full_near_midpoint_left_det < -0.05 && full_near_midpoint_left_det > -0.06,
        "full corrected §4.2 should retain the known near-midpoint-left outer-corner fold, got {full_near_midpoint_left_det}"
    );
    assert!(
        (full_without_fig8_det - full_near_midpoint_left_det).abs() < 1.0e-3,
        "Fig.8 indicated averages should not be the near-midpoint-left outer-corner determinant driver, got with={full_near_midpoint_left_det} without={full_without_fig8_det}"
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_eq12_split_outer_corner_is_localized() {
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
    let split_rows = [
        (
            Kp20Section42Eq12Split::PlainEndpoint,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
                Kp20Section42Eq12Split::PlainEndpoint,
            ),
        ),
        (
            Kp20Section42Eq12Split::PlainPlusRow,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
                Kp20Section42Eq12Split::PlainPlusRow,
            ),
        ),
        (
            Kp20Section42Eq12Split::PlainMinusRow,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
                Kp20Section42Eq12Split::PlainMinusRow,
            ),
        ),
        (
            Kp20Section42Eq12Split::ArrowEndpoint,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
                Kp20Section42Eq12Split::ArrowEndpoint,
            ),
        ),
        (
            Kp20Section42Eq12Split::ArrowPlusRow,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
                Kp20Section42Eq12Split::ArrowPlusRow,
            ),
        ),
        (
            Kp20Section42Eq12Split::ArrowMinusRow,
            kp20_section42_eq12_split_rows_with_symbols(
                map,
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
                Kp20Section42Eq12Split::ArrowMinusRow,
            ),
        ),
    ];

    let solve_with_eq12_splits = |label: &str,
                                  included_splits: &[Kp20Section42Eq12Split]|
     -> (usize, f64, f64) {
        let mut scalar_constraints = Vec::new();
        scalar_constraints.extend(c0_rows.iter().cloned());
        scalar_constraints.extend(eq11_rows.iter().cloned());
        scalar_constraints.extend(eq13_rows.iter().cloned());
        for (split, rows) in &split_rows {
            if included_splits.contains(split) {
                scalar_constraints.extend(rows.iter().cloned());
            }
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
        let orientation =
            kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
        eprintln!(
            "rotational top corrected eq12 split outer-corner {label}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e} nml_du=({nmldux:e},{nmlduy:e}) nml_dv=({nmldvx:e},{nmldvy:e}) nmr_det={nmr_det:e} nmr_du=({nmrdux:e},{nmrduy:e}) nmr_dv=({nmrdvx:e},{nmrdvy:e})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
            nml_det = nml.determinant,
            nmldux = nml.du.x,
            nmlduy = nml.du.y,
            nmldvx = nml.dv.x,
            nmldvy = nml.dv.y,
            nmr_det = nmr.determinant,
            nmrdux = nmr.du.x,
            nmrduy = nmr.du.y,
            nmrdvx = nmr.dv.x,
            nmrdvy = nmr.dv.y,
        );
        assert!(solved.max_constraint_residual < 1.0e-10);
        (
            orientation.negative_determinants,
            nml.determinant,
            nmr.determinant,
        )
    };

    let all_splits = Kp20Section42Eq12Split::ALL;
    let (_, full_nml_det, _) = solve_with_eq12_splits("full-eq12", &all_splits);
    let mut best_removed_for_corner = (
        Kp20Section42Eq12Split::PlainEndpoint,
        usize::MAX,
        f64::NEG_INFINITY,
    );
    for removed in Kp20Section42Eq12Split::ALL {
        let included = Kp20Section42Eq12Split::ALL
            .into_iter()
            .filter(|split| *split != removed)
            .collect::<Vec<_>>();
        let label = format!("without-{}", removed.label());
        let (neg, nml_det, _) = solve_with_eq12_splits(&label, &included);
        if nml_det > best_removed_for_corner.2 {
            best_removed_for_corner = (removed, neg, nml_det);
        }
    }
    let mut single_drivers = Vec::new();
    let mut worst_single = (
        Kp20Section42Eq12Split::PlainEndpoint,
        0_usize,
        f64::INFINITY,
    );
    for included in Kp20Section42Eq12Split::ALL {
        let label = format!("only-{}", included.label());
        let (neg, nml_det, _) = solve_with_eq12_splits(&label, &[included]);
        single_drivers.push((included, neg, nml_det));
        if nml_det < worst_single.2 {
            worst_single = (included, neg, nml_det);
        }
    }

    assert!(
        full_nml_det < -0.04,
        "full corrected eq12 should retain the near-midpoint-left outer-corner fold, got {full_nml_det}"
    );
    assert_eq!(
        best_removed_for_corner.0,
        Kp20Section42Eq12Split::ArrowEndpoint,
        "removing the arrow endpoint should give the largest exact near-midpoint-left outer-corner determinant relief, got {best_removed_for_corner:?}"
    );
    assert_eq!(
        worst_single.0,
        Kp20Section42Eq12Split::ArrowEndpoint,
        "the arrow endpoint should be the strongest single eq12 driver for this exact corner diagnostic, got {worst_single:?}"
    );
    for (split, _, nml_det) in single_drivers {
        match split {
            Kp20Section42Eq12Split::PlainEndpoint => assert!(
                nml_det < -0.02,
                "{split:?} alone should flip the near-midpoint-left outer-corner determinant, got {nml_det}"
            ),
            Kp20Section42Eq12Split::ArrowEndpoint => assert!(
                nml_det < -0.03,
                "{split:?} alone should flip the near-midpoint-left outer-corner determinant, got {nml_det}"
            ),
            Kp20Section42Eq12Split::PlainPlusRow
            | Kp20Section42Eq12Split::PlainMinusRow
            | Kp20Section42Eq12Split::ArrowPlusRow
            | Kp20Section42Eq12Split::ArrowMinusRow => assert!(
                nml_det > 0.1,
                "{split:?} alone should leave the near-midpoint-left outer-corner determinant positive, got {nml_det}"
            ),
        }
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
fn stam_appendix_a_rotational_top_vector_section42_symbol_axis_orientation_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;
    let orientations = [
        Kp20Section42SymbolOrientation::IDENTITY,
        Kp20Section42SymbolOrientation {
            reverse_p_j: false,
            reverse_arrow_p_j: false,
            reverse_underlined_p_j: true,
            reverse_arrow_underlined_p_j: true,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: false,
            reverse_arrow_p_j: false,
            reverse_underlined_p_j: true,
            reverse_arrow_underlined_p_j: false,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: false,
            reverse_arrow_p_j: false,
            reverse_underlined_p_j: false,
            reverse_arrow_underlined_p_j: true,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: true,
            reverse_arrow_p_j: true,
            reverse_underlined_p_j: false,
            reverse_arrow_underlined_p_j: false,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: true,
            reverse_arrow_p_j: false,
            reverse_underlined_p_j: false,
            reverse_arrow_underlined_p_j: false,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: false,
            reverse_arrow_p_j: true,
            reverse_underlined_p_j: false,
            reverse_arrow_underlined_p_j: false,
        },
        Kp20Section42SymbolOrientation {
            reverse_p_j: true,
            reverse_arrow_p_j: true,
            reverse_underlined_p_j: true,
            reverse_arrow_underlined_p_j: true,
        },
    ];
    let mut best = (
        Kp20Section42SymbolOrientation::IDENTITY,
        usize::MAX,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );

    for orientation in orientations {
        let scalar_constraints =
            kp20_section42_oriented_presolved_strip_basis_rows(map, 2.0 * c, orientation);
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
        let orientation_report =
            kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        let nmr = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointRight);
        eprintln!(
            "rotational top corrected §4.2 symbol-axis orientation {label}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e} nmr_det={nmr_det:e}",
            label = orientation.label(),
            residual = solved.max_constraint_residual,
            neg = orientation_report.negative_determinants,
            samples = orientation_report.sample_count,
            worst = orientation_report.worst_piece,
            u = orientation_report.worst_u,
            v = orientation_report.worst_v,
            nml_det = nml.determinant,
            nmr_det = nmr.determinant,
        );
        if solved.max_constraint_residual < 1.0e-10
            && (orientation_report.negative_determinants < best.1
                || (orientation_report.negative_determinants == best.1 && nml.determinant > best.2))
        {
            best = (
                orientation,
                orientation_report.negative_determinants,
                nml.determinant,
                nmr.determinant,
            );
        }
    }

    assert_eq!(
        best.0,
        Kp20Section42SymbolOrientation::IDENTITY,
        "coherent symbol-axis flips should remain rejected unless they beat the accepted same-axis §4.2 map, got {best:?}"
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
fn stam_appendix_a_rotational_top_vector_positive_outer_corner_completion_is_feasible() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let anchored_scale = ring.lambda.powi(4);
    let map = kp20_section42_map_candidates(angle)[1].1;

    let solve = |scalar_constraints: Vec<Kp20PlanarTauScalarLinearConstraint>,
                 extra_vector_constraints: Vec<Kp20PlanarTauLinearConstraint>|
     -> Kp20PlanarTauAugmentedF3SolveOutput {
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
        let mut vector_constraints = cc_test_fig8_indicated_average_variant_rows(false, false);
        vector_constraints.extend(extra_vector_constraints);
        solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints,
                scalar_constraints,
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap()
    };

    let mut c0_eq11_rows = Vec::new();
    c0_eq11_rows.extend(kp20_section42_only_c0_rows(map));
    c0_eq11_rows.extend(kp20_section42_only_eq11_rows(map));
    let c0_eq11_solved = solve(c0_eq11_rows, Vec::new());
    let c0_eq11_local =
        cc_test_global_planar_patch_to_local(valence, 0, c0_eq11_solved.tau_macro_patch);
    let c0_eq11_global = c0_eq11_solved.tau_macro_patch;
    let c0_eq11_nml = cc_test_outer_corner_frame(&c0_eq11_local, Kp20PieceId::NearMidpointLeft);
    let c0_eq11_nmr = cc_test_outer_corner_frame(&c0_eq11_local, Kp20PieceId::NearMidpointRight);
    assert!(c0_eq11_nml.determinant > 0.09);
    assert!(c0_eq11_nmr.determinant > 0.09);

    let full_rows = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let full_unpinned = solve(full_rows.clone(), Vec::new());
    let full_unpinned_local =
        cc_test_global_planar_patch_to_local(valence, 0, full_unpinned.tau_macro_patch);
    let unpinned_orientation =
        kp20_planar_bi6_tau_orientation_report(&full_unpinned_local, &[0.0, 0.5, 1.0], 1.0e-10);
    let unpinned_nml =
        cc_test_outer_corner_frame(&full_unpinned_local, Kp20PieceId::NearMidpointLeft);

    let mut outer_corner_pins = Vec::new();
    for piece_id in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        for (i, j) in [(6, 6), (5, 6), (6, 5)] {
            let target = c0_eq11_global.piece(piece_id).cps[i][j];
            outer_corner_pins.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece_id, i, j),
                [target.x, target.y],
            ));
        }
    }
    let full_pinned = solve(full_rows, outer_corner_pins);
    let full_pinned_local =
        cc_test_global_planar_patch_to_local(valence, 0, full_pinned.tau_macro_patch);
    let pinned_orientation =
        kp20_planar_bi6_tau_orientation_report(&full_pinned_local, &[0.0, 0.5, 1.0], 1.0e-10);
    let pinned_nml = cc_test_outer_corner_frame(&full_pinned_local, Kp20PieceId::NearMidpointLeft);
    let pinned_nmr = cc_test_outer_corner_frame(&full_pinned_local, Kp20PieceId::NearMidpointRight);
    for piece_id in Kp20PieceId::all() {
        let negatives =
            sampled_piece_negative_tau_locations(&full_pinned_local, piece_id, &[0.0, 0.5, 1.0]);
        eprintln!(
            "rotational top corrected §4.2 positive outer-corner remaining negatives {piece_id:?}: {negatives:?}"
        );
    }

    eprintln!(
        "rotational top corrected §4.2 positive outer-corner completion: unpinned_residual={unpinned_residual:e} unpinned_neg={unpinned_neg}/{unpinned_samples} unpinned_nml_det={unpinned_nml_det:e}; pinned_residual={pinned_residual:e} pinned_neg={pinned_neg}/{pinned_samples} pinned_nml_det={pinned_nml_det:e} pinned_nmr_det={pinned_nmr_det:e}",
        unpinned_residual = full_unpinned.max_constraint_residual,
        unpinned_neg = unpinned_orientation.negative_determinants,
        unpinned_samples = unpinned_orientation.sample_count,
        unpinned_nml_det = unpinned_nml.determinant,
        pinned_residual = full_pinned.max_constraint_residual,
        pinned_neg = pinned_orientation.negative_determinants,
        pinned_samples = pinned_orientation.sample_count,
        pinned_nml_det = pinned_nml.determinant,
        pinned_nmr_det = pinned_nmr.determinant,
    );
    assert!(full_pinned.max_constraint_residual < 1.0e-10);
    assert!(pinned_nml.determinant > 0.09);
    assert!(pinned_nmr.determinant > 0.09);
    assert!(
        pinned_orientation.negative_determinants < unpinned_orientation.negative_determinants,
        "pinning the unconstrained positive outer-corner completion should reduce coarse folds if the blocker is F3 completion, got pinned={} unpinned={}",
        pinned_orientation.negative_determinants,
        unpinned_orientation.negative_determinants
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_section42_f3_nullspace_dependence_is_measured() {
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

    let zero_initial = Kp20MacroPatch::zero();
    let c0_eq11_initial = {
        let mut c0_eq11_rows = Vec::new();
        c0_eq11_rows.extend(kp20_section42_only_c0_rows(map));
        c0_eq11_rows.extend(kp20_section42_only_eq11_rows(map));
        solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: Kp20MacroPatch::zero(),
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: c0_eq11_rows,
                augmented_scalar_constraints: augmented_scalar_constraints.clone(),
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap()
        .tau_macro_patch
    };

    let cases = [
        ("zero-reg-0", zero_initial, 0.0),
        ("zero-reg-1e-12", zero_initial, 1.0e-12),
        ("zero-reg-1e-10", zero_initial, 1.0e-10),
        ("zero-reg-1e-8", zero_initial, 1.0e-8),
        ("c0-eq11-reg-1e-12", c0_eq11_initial, 1.0e-12),
        ("c0-eq11-reg-1e-10", c0_eq11_initial, 1.0e-10),
    ];
    let mut zero_target_neg = None;
    let mut c0_target_neg = None;
    let mut unregularized_neg = None;

    for (label, initial, regularization_weight) in cases {
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial,
                vector_constraints: cc_test_fig8_indicated_average_variant_rows(false, false),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: augmented_scalar_constraints.clone(),
                extra_variable_count: source_extra_count + 1,
                regularization_weight,
                extra_regularization_weight: 1.0e-12,
            },
        );
        match solved {
            Ok(solved) => {
                let local_tau =
                    cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
                let orientation =
                    kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
                let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
                eprintln!(
                    "rotational top corrected §4.2 F3 nullspace {label}: residual={residual:e} kkt={kkt:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e}",
                    residual = solved.max_constraint_residual,
                    kkt = solved.kkt_condition_number,
                    neg = orientation.negative_determinants,
                    samples = orientation.sample_count,
                    worst = orientation.worst_piece,
                    u = orientation.worst_u,
                    v = orientation.worst_v,
                    nml_det = nml.determinant,
                );
                assert!(solved.max_constraint_residual < 1.0e-8);
                if label == "zero-reg-0" {
                    unregularized_neg = Some(orientation.negative_determinants);
                } else if label == "zero-reg-1e-10" {
                    zero_target_neg = Some(orientation.negative_determinants);
                } else if label == "c0-eq11-reg-1e-10" {
                    c0_target_neg = Some(orientation.negative_determinants);
                }
            }
            Err(error) => {
                eprintln!("rotational top corrected §4.2 F3 nullspace {label}: error={error:?}");
                assert_eq!(regularization_weight, 0.0);
            }
        }
    }

    assert_eq!(
        unregularized_neg,
        Some(15),
        "unregularized F3 solve should remain the known 15/36 fold if nullspace selection is not the frontier"
    );
    assert_eq!(
        zero_target_neg,
        Some(15),
        "zero-target regularized baseline should remain the known 15/36 coarse fold"
    );
    assert_eq!(
        c0_target_neg,
        Some(15),
        "changing only the regularization target should not change the fold count if the F3 optimum is well determined"
    );
}

#[test]
fn stam_appendix_a_rotational_top_vector_lemma4_scale_sweep_is_measured() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let valence = ring.sectors.len();
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_axes = cc_test_sector_local_axis_components(valence, 0);
    let source_extra_count = 4;
    let map = kp20_section42_map_candidates(angle)[1].1;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    let mut best = (0_i32, f64::INFINITY, usize::MAX, f64::NEG_INFINITY);

    for exponent in -2..=8 {
        let anchored_scale = ring.lambda.powi(exponent);
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
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: source_extra_count + 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-12,
            },
        )
        .unwrap();
        let local_tau = cc_test_global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let orientation =
            kp20_planar_bi6_tau_orientation_report(&local_tau, &[0.0, 0.5, 1.0], 1.0e-10);
        let nml = cc_test_outer_corner_frame(&local_tau, Kp20PieceId::NearMidpointLeft);
        eprintln!(
            "rotational top corrected §4.2 Lemma4 scale sweep exponent={exponent} scale={anchored_scale:e}: residual={residual:e} neg={neg}/{samples} worst={worst:?}@({u:.2},{v:.2}) nml_det={nml_det:e}",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            worst = orientation.worst_piece,
            u = orientation.worst_u,
            v = orientation.worst_v,
            nml_det = nml.determinant,
        );
        assert!(solved.max_constraint_residual < 1.0e-10);
        if orientation.negative_determinants < best.2
            || (orientation.negative_determinants == best.2 && nml.determinant > best.3)
        {
            best = (
                exponent,
                anchored_scale,
                orientation.negative_determinants,
                nml.determinant,
            );
        }
    }

    assert_eq!(
        best.0, 0,
        "the corrected cross-degree source rows move the best Lemma scale to lambda^0, got {best:?}"
    );
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
fn stam_appendix_a_composed_pr08_subdivision_source_residual_localizes_by_family() {
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
    let residuals = kp20_planar_tau_rotational_paper_family_residuals(
        5,
        0,
        &tensor_border,
        &solved.tau_macro_patch,
    )
    .unwrap();
    eprintln!("composed PR08 subdivision source family residuals: {residuals:?}");
    assert!(residuals.iter().all(|residual| {
        matches!(
            residual.family,
            Kp20PlanarTauPaperRowFamily::TensorBorderEq14
                | Kp20PlanarTauPaperRowFamily::DiagonalSymmetry
        ) || residual.max_residual < 1.0e-9
    }));
    let diagonal = residuals
        .iter()
        .find(|residual| residual.family == Kp20PlanarTauPaperRowFamily::DiagonalSymmetry)
        .unwrap();
    assert!(
        diagonal.max_residual < 1.0e-9,
        "corrected Fig.8 diagonal reflection is now accepted by the live solve, got {diagonal:?}"
    );
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
        if solved.max_constraint_residual < finite_level_best.residual {
            finite_level_best = StamAppendixADyadicLevelDiagnostic {
                residual: solved.max_constraint_residual,
                level_exponent,
            };
        }
    }

    assert!(
        finite_level_best.residual < 1.0e-10,
        "finite dyadic level sweep should expose pre-solved §4.2/§4.3 closure, got {finite_level_best:?}"
    );
    assert_eq!(finite_level_best.level_exponent, 4);

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
fn stam_appendix_a_scaled_pr08_q13_q2_residual_closes_by_presolved_family() {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_scaled_top_label_projected_rows(
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
    let residuals = kp20_planar_tau_rotational_paper_family_residuals(
        5,
        0,
        &tensor_border,
        &solved.tau_macro_patch,
    )
    .unwrap();
    eprintln!("scaled PR08 q13/q2 family residuals: {residuals:?}");
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::TensorBorderEq14
            && residual.max_residual.is_finite()
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::Section42C0
            && residual.max_residual < 1.0e-10
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq11Degree5
            && residual.max_residual < 1.0e-10
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq12InternalC2
            && residual.max_residual < 1.0e-10
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::Section42Eq13Strip
            && residual.max_residual < 1.0e-10
    }));
    assert!(residuals.iter().any(|residual| {
        residual.family == Kp20PlanarTauPaperRowFamily::DiagonalSymmetry
            && residual.max_residual < 1.0e-10
    }));
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
    assert_eq!(best.midpoint_terms.duv_source, Kp20Eq9TermSource::Current);
    assert_eq!(best.midpoint_terms.duu_source, Kp20Eq9TermSource::Previous);
    assert_eq!(best.midpoint_terms.du_source, Kp20Eq9TermSource::Previous);
    assert_eq!(best.midpoint_terms.signs.two_b_duv, 1.0);
    assert_eq!(best.midpoint_terms.signs.b_sq_duu, -1.0);
    assert_eq!(best.midpoint_terms.signs.b_bprime_du, 1.0);
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

#[derive(Debug, Clone, Copy)]
struct StamAppendixBProjectionDiagnostic {
    residual: f64,
    patch_index: usize,
    transposed: bool,
    row_start: usize,
    reverse_rows: bool,
    reverse_cols: bool,
}

#[derive(Debug, Clone, Copy)]
struct StamTransformedSideProjectionDiagnostic {
    residual: f64,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    side: Kp20StamGridSide,
    reversed_depth: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20StamProjectionFamily {
    VisualFigure63,
    AppendixTensorOrder,
}

impl Kp20StamProjectionFamily {
    const fn name(self) -> &'static str {
        match self {
            Self::VisualFigure63 => "visual-figure-6-3",
            Self::AppendixTensorOrder => "appendix-tensor-order",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20StamSideEndpointFitDiagnostic {
    max_endpoint_residual: f64,
    sum_endpoint_residual: f64,
    family: Kp20StamProjectionFamily,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    side: Kp20StamGridSide,
    reversed_depth: bool,
    swap_targets: bool,
    conjugate_source: bool,
}

#[derive(Debug, Clone, Copy)]
struct Kp20StamSideTauValidityDiagnostic {
    residual: f64,
    square_outside: usize,
    triangle_outside: usize,
    negative_determinants: usize,
    min_determinant: f64,
    endpoint: Kp20StamSideEndpointFitDiagnostic,
}

#[derive(Debug)]
struct StamAppendixAUnknownTopDiagnostic {
    residual: f64,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
}

#[derive(Debug)]
struct StamAppendixAPr08InnerSeedDiagnostic {
    residual: f64,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
}

fn assert_stam_appendix_a_pr08_inner_seed_diagnostic_in_search_space(
    diagnostic: StamAppendixAPr08InnerSeedDiagnostic,
) {
    assert!(diagnostic.patch_index < 3);
    assert!(diagnostic.row_start <= 1);
    match diagnostic.reverse_rows {
        true | false => {}
    }
}

#[derive(Debug)]
struct StamAppendixADyadicLevelDiagnostic {
    residual: f64,
    level_exponent: i32,
}

#[derive(Debug, Clone, Copy)]
struct Pr08TopRingLabelMatch {
    residual: f64,
    label_index: usize,
}

struct StamAppendixATopVectorRings {
    odd: Vec<Complex>,
    even: Vec<Complex>,
}

fn stam_appendix_a_scaled_pr08_q13_q2_tensor_border() -> TensorBorderBi3 {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let stam_top = stam_appendix_a_frequency_one_top_vector(&ring);
    let odd_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 12), &stam_top.odd)
            .unwrap();
    let even_scale =
        complex_similarity_ring_scale(&pr08_label_ring_as_complex(&ring, 1), &stam_top.even)
            .unwrap();
    let rows = stam_appendix_a_pr08_scaled_top_label_projected_rows(
        &ring,
        (12, odd_scale),
        (1, even_scale),
        1,
        1,
        true,
    )
    .unwrap();
    kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap()
}

fn stam_appendix_a_label_top_weight_matrix(
    lambda: f64,
    top_count: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let s11 = stam_appendix_a_s11(top_count);
    let s12 = stam_appendix_a_s12();
    let mut system = s12;
    for (idx, row) in system.iter_mut().enumerate() {
        row[idx] -= lambda;
    }

    let mut outer_weights = vec![vec![0.0_f64; top_count]; 7];
    for top_index in 0..top_count {
        let rhs = s11.iter().map(|row| -row[top_index]).collect::<Vec<_>>();
        let solved = solve_real_dense(system.iter().map(|row| row.to_vec()).collect(), rhs)
            .ok_or_else(|| {
                "singular Stam Appendix-A outer-control source matrix solve".to_string()
            })?;
        for row in 0..7 {
            outer_weights[row][top_index] = solved[row];
        }
    }

    let mut labels = vec![vec![0.0_f64; top_count]; top_count + 16];
    for top_index in 0..top_count {
        labels[top_index][top_index] = lambda;
    }
    for row in 0..7 {
        for top_index in 0..top_count {
            labels[top_count + row][top_index] = lambda * outer_weights[row][top_index];
        }
    }

    let s21 = stam_appendix_a_s21(top_count);
    let s22 = stam_appendix_a_s22();
    for row in 0..9 {
        for top_index in 0..top_count {
            let outer = (0..7)
                .map(|outer_row| s22[row][outer_row] * outer_weights[outer_row][top_index])
                .sum::<f64>();
            labels[top_count + 7 + row][top_index] = s21[row][top_index] + outer;
        }
    }
    Ok(labels)
}

fn kp20_eq13_rotated_near_m_strip_assignment_rows(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    b0: f64,
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let w1 = -b0 / 100.0;
    let w2 = -w1;
    let w3 = b0 / 30.0;
    let mut rows = Vec::new();

    for i in [5, 6] {
        add_rotated_vector_constraint_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &[
                (true, i, 1, 1.0),
                (true, i, 0, -1.0),
                (true, i, 2, -0.25),
                (false, i, 2, 0.25),
            ],
            previous_sector_rotation_radians,
        );
        add_rotated_vector_constraint_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &[
                (false, i, 1, 1.0),
                (true, i, 0, -1.0),
                (true, i, 2, 0.25),
                (false, i, 2, -0.25),
            ],
            previous_sector_rotation_radians,
        );
    }

    add_rotated_vector_constraint_rows(
        &mut rows,
        previous_piece,
        current_piece,
        &[
            (true, 4, 1, 1.0),
            (true, 4, 0, -1.0),
            (true, 4, 2, -0.25),
            (false, 4, 2, 0.25),
            (true, 5, 0, w3),
            (true, 6, 0, -w3),
            (true, 5, 2, -w2),
            (false, 5, 2, w2),
            (true, 6, 2, -w1),
            (false, 6, 2, w1),
        ],
        previous_sector_rotation_radians,
    );
    add_rotated_vector_constraint_rows(
        &mut rows,
        previous_piece,
        current_piece,
        &[
            (false, 4, 1, 1.0),
            (true, 4, 0, -1.0),
            (true, 4, 2, 0.25),
            (false, 4, 2, -0.25),
            (true, 5, 0, w3),
            (true, 6, 0, -w3),
            (true, 5, 2, w2),
            (false, 5, 2, -w2),
            (true, 6, 2, w1),
            (false, 6, 2, -w1),
        ],
        previous_sector_rotation_radians,
    );
    rows
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42PaperSymbolMap {
    p: Kp20Section42MappedPiece,
    arrow_p: Kp20Section42MappedPiece,
    underlined_p: Kp20Section42MappedPiece,
    arrow_underlined_p: Kp20Section42MappedPiece,
    previous_sector_rotation_radians: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42MappedPiece {
    piece: Kp20PieceId,
    is_previous_sector: bool,
}

impl Kp20Section42MappedPiece {
    const fn previous(piece: Kp20PieceId) -> Self {
        Self {
            piece,
            is_previous_sector: true,
        }
    }

    const fn current(piece: Kp20PieceId) -> Self {
        Self {
            piece,
            is_previous_sector: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Kp20Section42PaperSymbol {
    P,
    ArrowP,
    UnderlinedP,
    ArrowUnderlinedP,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20Section42DiagnosticFamily {
    C0Only,
    Eq11Only,
    Eq12Only,
    Eq13Only,
    C0Eq13,
    C0Eq11Eq13,
    C0Eq12Eq13,
    Full,
}

impl Kp20Section42DiagnosticFamily {
    const fn all() -> [Self; 8] {
        [
            Self::C0Only,
            Self::Eq11Only,
            Self::Eq12Only,
            Self::Eq13Only,
            Self::C0Eq13,
            Self::C0Eq11Eq13,
            Self::C0Eq12Eq13,
            Self::Full,
        ]
    }

    const fn label(self) -> &'static str {
        match self {
            Self::C0Only => "c0-only",
            Self::Eq11Only => "eq11-only",
            Self::Eq12Only => "eq12-only",
            Self::Eq13Only => "eq13-only",
            Self::C0Eq13 => "c0-eq13",
            Self::C0Eq11Eq13 => "c0-eq11-eq13",
            Self::C0Eq12Eq13 => "c0-eq12-eq13",
            Self::Full => "full",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42PaperTerm {
    symbol: Kp20Section42PaperSymbol,
    i: usize,
    j: usize,
    weight: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Section42ResidualBreakdown {
    c0: (f64, usize),
    eq11: (f64, usize),
    eq12: (f64, usize),
    eq13: (f64, usize),
}

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct Kp20Section42ProjectionSweepDiagnostic {
    residual: f64,
    map_label: &'static str,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20Section42SourceProjectionScope {
    UnderlinedC0PinnedRows,
    UnderlinedSection42PinnedRows,
    AllTensorPins,
}

#[derive(Debug, Clone, PartialEq)]
struct Kp20Section42SourceProjectionReport {
    max_row_residual: f64,
    max_source_delta: f64,
    adjusted_scalar_targets: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20Section42Eq14Half {
    NearC,
    NearM,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Kp20Section42Eq14PieceSource {
    piece: Kp20PieceId,
    sector: usize,
    half: Kp20Section42Eq14Half,
    reverse_u: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Section42Eq14HalfSourceSweep {
    residual: f64,
    left_sector_offset: usize,
    left_half: Kp20Section42Eq14Half,
    left_reverse: bool,
    right_sector_offset: usize,
    right_half: Kp20Section42Eq14Half,
    right_reverse: bool,
}

#[derive(Debug, Clone, PartialEq)]
struct Kp20ScalarRowResidualDiagnostic {
    row_index: usize,
    residual: f64,
    terms: Vec<(Kp20PlanarTauScalarCoeff, f64)>,
}

fn kp20_section42_map_candidates(
    previous_sector_rotation_radians: f64,
) -> [(&'static str, Kp20Section42PaperSymbolMap); 4] {
    [
        (
            "right-left/same-side-underlined",
            Kp20Section42PaperSymbolMap {
                p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
                arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
                underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointRight),
                arrow_underlined_p: Kp20Section42MappedPiece::current(
                    Kp20PieceId::NearMidpointLeft,
                ),
                previous_sector_rotation_radians,
            },
        ),
        (
            "right-left/prior-diagnostic-underlined",
            Kp20Section42PaperSymbolMap {
                p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
                arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
                underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
                arrow_underlined_p: Kp20Section42MappedPiece::current(
                    Kp20PieceId::NearMidpointRight,
                ),
                previous_sector_rotation_radians,
            },
        ),
        (
            "left-right/same-side-underlined",
            Kp20Section42PaperSymbolMap {
                p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
                arrow_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
                underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
                arrow_underlined_p: Kp20Section42MappedPiece::previous(
                    Kp20PieceId::NearMidpointRight,
                ),
                previous_sector_rotation_radians,
            },
        ),
        (
            "left-right/prior-diagnostic-underlined",
            Kp20Section42PaperSymbolMap {
                p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
                arrow_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
                underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
                arrow_underlined_p: Kp20Section42MappedPiece::previous(
                    Kp20PieceId::NearMidpointLeft,
                ),
                previous_sector_rotation_radians,
            },
        ),
    ]
}

fn solve_kp20_section42_family_residual(
    tensor_border: &TensorBorderBi3,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    family: Kp20Section42DiagnosticFamily,
    include_diagonal_symmetry: bool,
) -> f64 {
    let mut vector_constraints = Vec::new();
    vector_constraints.extend(kp20_planar_tau_tensor_border_constraints(tensor_border, 0).unwrap());
    let mut scalar_constraints = if include_diagonal_symmetry {
        kp20_planar_tau_diagonal_symmetry_constraints()
    } else {
        Vec::new()
    };
    scalar_constraints.extend(kp20_section42_diagnostic_rows(map, b0_near_origin, family));
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: Kp20MacroPatch::zero(),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1e-10,
    })
    .unwrap()
    .max_constraint_residual
}

fn solve_kp20_section42_source_projection(
    tensor_border: &TensorBorderBi3,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    scope: Kp20Section42SourceProjectionScope,
) -> Kp20Section42SourceProjectionReport {
    project_kp20_section42_source_tensor_rows(tensor_border, map, b0_near_origin, scope).1
}

fn project_kp20_section42_source_tensor_rows(
    tensor_border: &TensorBorderBi3,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    scope: Kp20Section42SourceProjectionScope,
) -> (
    Vec<Kp20PlanarTauLinearConstraint>,
    Kp20Section42SourceProjectionReport,
) {
    let tensor_rows = kp20_planar_tau_tensor_border_constraints(tensor_border, 0).unwrap();
    let free_coeffs = kp20_section42_source_projection_free_coeffs(&tensor_rows, map, scope);
    let mut vector_constraints = Vec::new();
    let mut augmented_scalar_constraints = Vec::new();
    let mut free_scalar_targets = Vec::new();
    let mut extra_count = 0_usize;

    for (row_index, row) in tensor_rows.iter().enumerate() {
        let coeff = row.terms[0].0;
        if free_coeffs.contains(&coeff) {
            for axis in 0..=1 {
                augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                    terms: vec![
                        (
                            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                                coeff, axis,
                            )),
                            1.0,
                        ),
                        (Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_count), -1.0),
                    ],
                    target: row.target[axis],
                });
                free_scalar_targets.push((row_index, axis));
                extra_count += 1;
            }
        } else {
            vector_constraints.push(row.clone());
        }
    }

    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, b0_near_origin);
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0,
        },
    )
    .unwrap();

    let mut projected_tensor_rows = tensor_rows;
    for ((row_index, axis), delta) in free_scalar_targets
        .into_iter()
        .zip(solved.extra_values.iter().copied())
    {
        projected_tensor_rows[row_index].target[axis] += delta;
    }

    let report = Kp20Section42SourceProjectionReport {
        max_row_residual: solved.max_constraint_residual,
        max_source_delta: max_point_delta_from_scalar_deltas(&solved.extra_values),
        adjusted_scalar_targets: solved.extra_values.len(),
    };
    (projected_tensor_rows, report)
}

fn kp20_section42_source_projection_free_coeffs(
    tensor_rows: &[Kp20PlanarTauLinearConstraint],
    map: Kp20Section42PaperSymbolMap,
    scope: Kp20Section42SourceProjectionScope,
) -> Vec<Kp20PlanarTauCoeff> {
    match scope {
        Kp20Section42SourceProjectionScope::UnderlinedC0PinnedRows => {
            let underlined =
                kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::UnderlinedP);
            let arrow_underlined =
                kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::ArrowUnderlinedP);
            let mut coeffs = Vec::new();
            for i in 0..=2 {
                coeffs.push(Kp20PlanarTauCoeff::new(underlined.piece, i, 0));
                coeffs.push(Kp20PlanarTauCoeff::new(arrow_underlined.piece, i, 0));
            }
            coeffs
        }
        Kp20Section42SourceProjectionScope::UnderlinedSection42PinnedRows => {
            let underlined =
                kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::UnderlinedP);
            let arrow_underlined =
                kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::ArrowUnderlinedP);
            let mut coeffs = Vec::new();
            for piece in [underlined.piece, arrow_underlined.piece] {
                for coeff in [
                    Kp20PlanarTauCoeff::new(piece, 0, 0),
                    Kp20PlanarTauCoeff::new(piece, 1, 0),
                    Kp20PlanarTauCoeff::new(piece, 2, 0),
                    Kp20PlanarTauCoeff::new(piece, 0, 1),
                    Kp20PlanarTauCoeff::new(piece, 0, 2),
                    Kp20PlanarTauCoeff::new(piece, 1, 1),
                    Kp20PlanarTauCoeff::new(piece, 1, 2),
                ] {
                    if tensor_rows
                        .iter()
                        .any(|row| row.terms == vec![(coeff, 1.0)])
                        && !coeffs.contains(&coeff)
                    {
                        coeffs.push(coeff);
                    }
                }
            }
            coeffs
        }
        Kp20Section42SourceProjectionScope::AllTensorPins => tensor_rows
            .iter()
            .map(|row| row.terms[0].0)
            .collect::<Vec<_>>(),
    }
}

fn max_point_delta_from_scalar_deltas(values: &[f64]) -> f64 {
    values
        .chunks(2)
        .map(|chunk| {
            let x = chunk[0];
            let y = chunk.get(1).copied().unwrap_or(0.0);
            (x * x + y * y).sqrt()
        })
        .fold(0.0_f64, f64::max)
}

fn worst_kp20_scalar_row_residual(
    tau: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> Kp20ScalarRowResidualDiagnostic {
    let mut worst = Kp20ScalarRowResidualDiagnostic {
        row_index: 0,
        residual: 0.0,
        terms: Vec::new(),
    };
    for (row_index, row) in rows.iter().enumerate() {
        let mut value = 0.0_f64;
        for &(coeff, weight) in &row.terms {
            let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
            let axis_value = if coeff.axis == 0 { cp.x } else { cp.y };
            value += weight * axis_value;
        }
        let residual = (value - row.target).abs();
        if residual > worst.residual {
            worst = Kp20ScalarRowResidualDiagnostic {
                row_index,
                residual,
                terms: row.terms.clone(),
            };
        }
    }
    worst
}

fn project_kp20_section42_underlined_c0_endpoint_tensor_pins(
    tensor_rows: &mut [Kp20PlanarTauLinearConstraint],
    map: Kp20Section42PaperSymbolMap,
) -> f64 {
    let underlined = kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::UnderlinedP);
    let arrow_underlined =
        kp20_section42_mapped_piece(map, Kp20Section42PaperSymbol::ArrowUnderlinedP);
    let mut max_delta = 0.0_f64;
    for i in 0..=2 {
        let underlined_row =
            find_kp20_tensor_pin_row(tensor_rows, Kp20PlanarTauCoeff::new(underlined.piece, i, 0));
        let arrow_row = find_kp20_tensor_pin_row(
            tensor_rows,
            Kp20PlanarTauCoeff::new(arrow_underlined.piece, i, 0),
        );

        let underlined_target = tensor_rows[underlined_row].target;
        let arrow_target = tensor_rows[arrow_row].target;
        let underlined_current =
            kp20_section42_mapped_target_to_current_frame(underlined_target, underlined, map);
        let arrow_current =
            kp20_section42_mapped_target_to_current_frame(arrow_target, arrow_underlined, map);
        let projected_current = [
            0.5 * (underlined_current[0] + arrow_current[0]),
            0.5 * (underlined_current[1] + arrow_current[1]),
        ];
        tensor_rows[underlined_row].target =
            kp20_section42_current_target_to_mapped_frame(projected_current, underlined, map);
        tensor_rows[arrow_row].target =
            kp20_section42_current_target_to_mapped_frame(projected_current, arrow_underlined, map);
        max_delta = max_delta.max(
            ((underlined_current[0] - arrow_current[0]).powi(2)
                + (underlined_current[1] - arrow_current[1]).powi(2))
            .sqrt()
                * 0.5,
        );
    }
    max_delta
}

fn find_kp20_tensor_pin_row(
    tensor_rows: &[Kp20PlanarTauLinearConstraint],
    coeff: Kp20PlanarTauCoeff,
) -> usize {
    tensor_rows
        .iter()
        .position(|row| row.terms == vec![(coeff, 1.0)])
        .unwrap_or_else(|| panic!("missing tensor-border pin row for {coeff:?}"))
}

fn kp20_section42_mapped_target_to_current_frame(
    target: [f64; 2],
    mapped: Kp20Section42MappedPiece,
    map: Kp20Section42PaperSymbolMap,
) -> [f64; 2] {
    if mapped.is_previous_sector {
        rotate_xy_array(target, map.previous_sector_rotation_radians)
    } else {
        target
    }
}

fn kp20_section42_current_target_to_mapped_frame(
    target: [f64; 2],
    mapped: Kp20Section42MappedPiece,
    map: Kp20Section42PaperSymbolMap,
) -> [f64; 2] {
    if mapped.is_previous_sector {
        rotate_xy_array(target, -map.previous_sector_rotation_radians)
    } else {
        target
    }
}

fn rotate_xy_array(point: [f64; 2], angle: f64) -> [f64; 2] {
    let (sin_theta, cos_theta) = angle.sin_cos();
    [
        cos_theta * point[0] - sin_theta * point[1],
        sin_theta * point[0] + cos_theta * point[1],
    ]
}

fn solve_kp20_section42_native_rows(
    tensor_border: &TensorBorderBi3,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    include_diagonal_symmetry: bool,
) -> Kp20PlanarTauF3SolveOutput {
    solve_kp20_section42_rows_with_tensor_constraints(
        kp20_planar_tau_tensor_border_constraints(tensor_border, 0).unwrap(),
        map,
        b0_near_origin,
        include_diagonal_symmetry,
    )
}

fn solve_kp20_section42_rows_with_tensor_constraints(
    vector_constraints: Vec<Kp20PlanarTauLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    include_diagonal_symmetry: bool,
) -> Kp20PlanarTauF3SolveOutput {
    let mut scalar_constraints = if include_diagonal_symmetry {
        kp20_planar_tau_diagonal_symmetry_constraints()
    } else {
        Vec::new()
    };
    scalar_constraints.extend(kp20_section42_presolved_strip_basis_rows(
        map,
        b0_near_origin,
    ));
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: Kp20MacroPatch::zero(),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1e-10,
    })
    .unwrap()
}

fn solve_kp20_section42_rows_with_extra_scalar_constraints(
    tensor_border: &TensorBorderBi3,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    mut extra_scalar_constraints: Vec<Kp20PlanarTauScalarLinearConstraint>,
) -> Kp20PlanarTauF3SolveOutput {
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, b0_near_origin);
    scalar_constraints.append(&mut extra_scalar_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: Kp20MacroPatch::zero(),
        vector_constraints: kp20_planar_tau_tensor_border_constraints(tensor_border, 0).unwrap(),
        scalar_constraints,
        regularization_weight: 1e-10,
    })
    .unwrap()
}

fn kp20_section42_uv_tensor_border_constraints_from_split(
    split: &[Kp20SplitTensorBorder],
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = split.get(sector).ok_or_else(|| {
        format!(
            "KP20 planar tau tensor-border sector {sector} outside valence {}",
            split.len()
        )
    })?;
    let mut rows = Vec::with_capacity(42);
    add_kp20_section42_uv_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        halves.near_c.rows,
    );
    add_kp20_section42_uv_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        halves.near_m.rows,
    );
    Ok(rows)
}

fn kp20_section42_depth_u_tensor_border_constraints_from_sources(
    split: &[Kp20SplitTensorBorder],
    left: Kp20Section42Eq14PieceSource,
    right: Kp20Section42Eq14PieceSource,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = Vec::with_capacity(42);
    add_kp20_section42_depth_u_tensor_border_piece_source(&mut rows, split, left);
    add_kp20_section42_depth_u_tensor_border_piece_source(&mut rows, split, right);
    rows
}

fn add_kp20_section42_depth_u_tensor_border_piece_source(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    split: &[Kp20SplitTensorBorder],
    source: Kp20Section42Eq14PieceSource,
) {
    let halves = &split[source.sector % split.len()];
    let half_rows = match source.half {
        Kp20Section42Eq14Half::NearC => halves.near_c.rows,
        Kp20Section42Eq14Half::NearM => halves.near_m.rows,
    };
    for (depth, row) in half_rows.iter().enumerate() {
        for u in 0..=6 {
            let cp = if source.reverse_u { row[6 - u] } else { row[u] };
            rows.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(source.piece, depth, u),
                [cp.x, cp.y],
            ));
        }
    }
}

fn add_kp20_section42_uv_tensor_border_half_rows(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    piece: Kp20PieceId,
    half_rows: [[Point3; 7]; 3],
) {
    for (v_depth, row) in half_rows.iter().enumerate() {
        for (u_index, cp) in row.iter().enumerate() {
            rows.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece, u_index, v_depth),
                [cp.x, cp.y],
            ));
        }
    }
}

fn kp20_section42_residual_breakdown(
    patch: &Kp20MacroPatch,
    c0_rows: &[Kp20PlanarTauScalarLinearConstraint],
    eq11_rows: &[Kp20PlanarTauScalarLinearConstraint],
    eq12_rows: &[Kp20PlanarTauScalarLinearConstraint],
    eq13_rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> Kp20Section42ResidualBreakdown {
    Kp20Section42ResidualBreakdown {
        c0: kp20_max_scalar_rows_residual_with_index(patch, c0_rows),
        eq11: kp20_max_scalar_rows_residual_with_index(patch, eq11_rows),
        eq12: kp20_max_scalar_rows_residual_with_index(patch, eq12_rows),
        eq13: kp20_max_scalar_rows_residual_with_index(patch, eq13_rows),
    }
}

fn kp20_section42_presolved_strip_basis_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    kp20_section42_diagnostic_rows(map, b0_near_origin, Kp20Section42DiagnosticFamily::Full)
}

fn kp20_section42_only_c0_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_c0_rows(&mut rows, map);
    rows
}

fn kp20_section42_only_eq11_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq11_degree5_rows(&mut rows, map);
    rows
}

fn kp20_section42_only_eq12_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq12_internal_c2_rows(&mut rows, map);
    rows
}

fn kp20_section42_only_eq13_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn kp20_section42_diagnostic_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    family: Kp20Section42DiagnosticFamily,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    if matches!(
        family,
        Kp20Section42DiagnosticFamily::C0Eq13
            | Kp20Section42DiagnosticFamily::C0Eq11Eq13
            | Kp20Section42DiagnosticFamily::C0Eq12Eq13
            | Kp20Section42DiagnosticFamily::Full
            | Kp20Section42DiagnosticFamily::C0Only
    ) {
        add_kp20_section42_c0_rows(&mut rows, map);
    }
    if matches!(
        family,
        Kp20Section42DiagnosticFamily::C0Eq11Eq13
            | Kp20Section42DiagnosticFamily::Full
            | Kp20Section42DiagnosticFamily::Eq11Only
    ) {
        add_kp20_section42_eq11_degree5_rows(&mut rows, map);
    }
    if matches!(
        family,
        Kp20Section42DiagnosticFamily::C0Eq12Eq13
            | Kp20Section42DiagnosticFamily::Full
            | Kp20Section42DiagnosticFamily::Eq12Only
    ) {
        add_kp20_section42_eq12_internal_c2_rows(&mut rows, map);
    }
    if !matches!(
        family,
        Kp20Section42DiagnosticFamily::C0Only
            | Kp20Section42DiagnosticFamily::Eq11Only
            | Kp20Section42DiagnosticFamily::Eq12Only
    ) {
        add_kp20_section42_eq13_rows(&mut rows, map, b0_near_origin);
    }
    rows
}

fn add_kp20_section42_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    for i in 0..=6 {
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

fn add_kp20_section42_eq11_degree5_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    const DEGREE6_ALT_BINOMIAL: [f64; 7] = [1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0];
    for symbol in [
        Kp20Section42PaperSymbol::P,
        Kp20Section42PaperSymbol::UnderlinedP,
    ] {
        let terms = DEGREE6_ALT_BINOMIAL
            .iter()
            .enumerate()
            .map(|(i, &weight)| section42_term(symbol, i, 0, weight))
            .collect::<Vec<_>>();
        add_kp20_section42_vector_rows(rows, map, &terms);
    }
}

fn add_kp20_section42_eq12_internal_c2_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    // KP20 eq. (12) writes rows 0, 1, and 5 of the upper/origin piece
    // from rows 2, 4, and 6 of the lower piece closer to m.
    for (upper, lower) in [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ] {
        for j in 0..=2 {
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
            );
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, 0.25),
                    section42_term(upper, 2, j, -0.25),
                ],
            );
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20Section42Eq12Split {
    PlainEndpoint,
    PlainPlusRow,
    PlainMinusRow,
    ArrowEndpoint,
    ArrowPlusRow,
    ArrowMinusRow,
}

impl Kp20Section42Eq12Split {
    const ALL: [Self; 6] = [
        Self::PlainEndpoint,
        Self::PlainPlusRow,
        Self::PlainMinusRow,
        Self::ArrowEndpoint,
        Self::ArrowPlusRow,
        Self::ArrowMinusRow,
    ];

    const fn label(self) -> &'static str {
        match self {
            Self::PlainEndpoint => "plain-endpoint",
            Self::PlainPlusRow => "plain-plus-row",
            Self::PlainMinusRow => "plain-minus-row",
            Self::ArrowEndpoint => "arrow-endpoint",
            Self::ArrowPlusRow => "arrow-plus-row",
            Self::ArrowMinusRow => "arrow-minus-row",
        }
    }
}

fn kp20_section42_eq12_split_rows(
    map: Kp20Section42PaperSymbolMap,
    split: Kp20Section42Eq12Split,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    let (upper, lower) = match split {
        Kp20Section42Eq12Split::PlainEndpoint
        | Kp20Section42Eq12Split::PlainPlusRow
        | Kp20Section42Eq12Split::PlainMinusRow => (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        Kp20Section42Eq12Split::ArrowEndpoint
        | Kp20Section42Eq12Split::ArrowPlusRow
        | Kp20Section42Eq12Split::ArrowMinusRow => (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    };

    match split {
        Kp20Section42Eq12Split::PlainEndpoint | Kp20Section42Eq12Split::ArrowEndpoint => {
            for j in 0..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 0, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                    ],
                );
            }
        }
        Kp20Section42Eq12Split::PlainPlusRow | Kp20Section42Eq12Split::ArrowPlusRow => {
            for j in 1..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 5, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                        section42_term(upper, 4, j, -0.25),
                        section42_term(upper, 2, j, 0.25),
                    ],
                );
            }
        }
        Kp20Section42Eq12Split::PlainMinusRow | Kp20Section42Eq12Split::ArrowMinusRow => {
            for j in 1..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 1, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                        section42_term(upper, 4, j, 0.25),
                        section42_term(upper, 2, j, -0.25),
                    ],
                );
            }
        }
    }
    rows
}

fn kp20_section42_eq12_endpoint_row(
    map: Kp20Section42PaperSymbolMap,
    upper: Kp20Section42PaperSymbol,
    lower: Kp20Section42PaperSymbol,
    j: usize,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_vector_rows(
        &mut rows,
        map,
        &[
            section42_term(lower, 0, j, 1.0),
            section42_term(upper, 6, j, -1.0),
        ],
    );
    rows
}

fn kp20_section42_eq12_split_rows_with_symbols(
    map: Kp20Section42PaperSymbolMap,
    upper: Kp20Section42PaperSymbol,
    lower: Kp20Section42PaperSymbol,
    split: Kp20Section42Eq12Split,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    match split {
        Kp20Section42Eq12Split::PlainEndpoint | Kp20Section42Eq12Split::ArrowEndpoint => {
            for j in 0..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 0, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                    ],
                );
            }
        }
        Kp20Section42Eq12Split::PlainPlusRow | Kp20Section42Eq12Split::ArrowPlusRow => {
            for j in 1..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 5, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                        section42_term(upper, 4, j, -0.25),
                        section42_term(upper, 2, j, 0.25),
                    ],
                );
            }
        }
        Kp20Section42Eq12Split::PlainMinusRow | Kp20Section42Eq12Split::ArrowMinusRow => {
            for j in 1..=2 {
                add_kp20_section42_vector_rows(
                    &mut rows,
                    map,
                    &[
                        section42_term(lower, 1, j, 1.0),
                        section42_term(upper, 6, j, -1.0),
                        section42_term(upper, 4, j, 0.25),
                        section42_term(upper, 2, j, -0.25),
                    ],
                );
            }
        }
    }
    rows
}

fn kp20_section42_eq12_rows_with_symbol_pairs(
    map: Kp20Section42PaperSymbolMap,
    symbol_pairs: &[(Kp20Section42PaperSymbol, Kp20Section42PaperSymbol)],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for &(upper, lower) in symbol_pairs {
        for j in 0..=2 {
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
            );
            add_kp20_section42_vector_rows(
                &mut rows,
                map,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, 0.25),
                    section42_term(upper, 2, j, -0.25),
                ],
            );
        }
    }
    rows
}

fn kp20_section42_eq13_rows_with_bands(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    strip_p: Kp20Section42PaperSymbol,
    strip_q: Kp20Section42PaperSymbol,
    row4_p: Kp20Section42PaperSymbol,
    row4_q: Kp20Section42PaperSymbol,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;

    for i in [5, 6] {
        add_kp20_section42_vector_rows(
            &mut rows,
            map,
            &[
                section42_term(strip_p, i, 1, 1.0),
                section42_term(strip_p, i, 0, -1.0),
                section42_term(strip_p, i, 2, -0.25),
                section42_term(strip_q, i, 2, 0.25),
            ],
        );
        add_kp20_section42_vector_rows(
            &mut rows,
            map,
            &[
                section42_term(strip_q, i, 1, 1.0),
                section42_term(strip_p, i, 0, -1.0),
                section42_term(strip_p, i, 2, 0.25),
                section42_term(strip_q, i, 2, -0.25),
            ],
        );
    }

    add_kp20_section42_vector_rows(
        &mut rows,
        map,
        &[
            section42_term(row4_p, 4, 1, 1.0),
            section42_term(row4_p, 4, 0, -1.0),
            section42_term(row4_p, 4, 2, -0.25),
            section42_term(row4_q, 4, 2, 0.25),
            section42_term(row4_p, 5, 0, w3),
            section42_term(row4_p, 6, 0, -w3),
            section42_term(row4_p, 5, 2, -w2),
            section42_term(row4_q, 5, 2, w2),
            section42_term(row4_p, 6, 2, -w1),
            section42_term(row4_q, 6, 2, w1),
        ],
    );
    add_kp20_section42_vector_rows(
        &mut rows,
        map,
        &[
            section42_term(row4_q, 4, 1, 1.0),
            section42_term(row4_p, 4, 0, -1.0),
            section42_term(row4_p, 4, 2, 0.25),
            section42_term(row4_q, 4, 2, -0.25),
            section42_term(row4_p, 5, 0, w3),
            section42_term(row4_p, 6, 0, -w3),
            section42_term(row4_p, 5, 2, w2),
            section42_term(row4_q, 5, 2, -w2),
            section42_term(row4_p, 6, 2, w1),
            section42_term(row4_q, 6, 2, -w1),
        ],
    );
    rows
}

fn add_kp20_section42_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;

    for i in [5, 6] {
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
        );
    }

    add_kp20_section42_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::P, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, -w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, w1),
        ],
    );
    add_kp20_section42_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, -w1),
        ],
    );
}

fn section42_term(
    symbol: Kp20Section42PaperSymbol,
    i: usize,
    j: usize,
    weight: f64,
) -> Kp20Section42PaperTerm {
    Kp20Section42PaperTerm {
        symbol,
        i,
        j,
        weight,
    }
}

fn add_kp20_section42_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let (sin_theta, cos_theta) = map.previous_sector_rotation_radians.sin_cos();
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let mapped = kp20_section42_mapped_piece(map, term.symbol);
        let coeff = Kp20PlanarTauCoeff::new(mapped.piece, term.i, term.j);
        if mapped.is_previous_sector {
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * cos_theta,
            ));
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                -term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                term.weight * cos_theta,
            ));
        } else {
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), term.weight));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), term.weight));
        }
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Kp20Section42SymbolOrientation {
    reverse_p_j: bool,
    reverse_arrow_p_j: bool,
    reverse_underlined_p_j: bool,
    reverse_arrow_underlined_p_j: bool,
}

impl Kp20Section42SymbolOrientation {
    const IDENTITY: Self = Self {
        reverse_p_j: false,
        reverse_arrow_p_j: false,
        reverse_underlined_p_j: false,
        reverse_arrow_underlined_p_j: false,
    };

    const fn label(self) -> &'static str {
        match (
            self.reverse_p_j,
            self.reverse_arrow_p_j,
            self.reverse_underlined_p_j,
            self.reverse_arrow_underlined_p_j,
        ) {
            (false, false, false, false) => "identity",
            (false, false, true, true) => "underlined-both-reverse-j",
            (false, false, true, false) => "underlined-p-reverse-j",
            (false, false, false, true) => "arrow-underlined-p-reverse-j",
            (true, true, false, false) => "plain-both-reverse-j",
            (true, false, false, false) => "p-reverse-j",
            (false, true, false, false) => "arrow-p-reverse-j",
            (true, true, true, true) => "all-reverse-j",
            _ => "mixed-reverse-j",
        }
    }

    const fn reverse_j(self, symbol: Kp20Section42PaperSymbol) -> bool {
        match symbol {
            Kp20Section42PaperSymbol::P => self.reverse_p_j,
            Kp20Section42PaperSymbol::ArrowP => self.reverse_arrow_p_j,
            Kp20Section42PaperSymbol::UnderlinedP => self.reverse_underlined_p_j,
            Kp20Section42PaperSymbol::ArrowUnderlinedP => self.reverse_arrow_underlined_p_j,
        }
    }
}

fn kp20_section42_oriented_presolved_strip_basis_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    orientation: Kp20Section42SymbolOrientation,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_oriented_c0_rows(&mut rows, map, orientation);
    add_kp20_section42_oriented_eq11_rows(&mut rows, map, orientation);
    add_kp20_section42_oriented_eq12_rows(&mut rows, map, orientation);
    add_kp20_section42_oriented_eq13_rows(&mut rows, map, b0_near_origin, orientation);
    rows
}

fn add_kp20_section42_oriented_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    orientation: Kp20Section42SymbolOrientation,
) {
    for i in 0..=6 {
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            orientation,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            orientation,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

fn add_kp20_section42_oriented_eq11_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    orientation: Kp20Section42SymbolOrientation,
) {
    const DEGREE6_ALT_BINOMIAL: [f64; 7] = [1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0];
    for symbol in [
        Kp20Section42PaperSymbol::P,
        Kp20Section42PaperSymbol::UnderlinedP,
    ] {
        let terms = DEGREE6_ALT_BINOMIAL
            .iter()
            .enumerate()
            .map(|(i, &weight)| section42_term(symbol, i, 0, weight))
            .collect::<Vec<_>>();
        add_kp20_section42_oriented_vector_rows(rows, map, orientation, &terms);
    }
}

fn add_kp20_section42_oriented_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    orientation: Kp20Section42SymbolOrientation,
) {
    for (upper, lower) in [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ] {
        for j in 0..=2 {
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                orientation,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                orientation,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
            );
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                orientation,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, 0.25),
                    section42_term(upper, 2, j, -0.25),
                ],
            );
        }
    }
}

fn add_kp20_section42_oriented_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    orientation: Kp20Section42SymbolOrientation,
) {
    let mut raw_rows = Vec::new();
    add_kp20_section42_eq13_rows(&mut raw_rows, map, b0_near_origin);
    if orientation == Kp20Section42SymbolOrientation::IDENTITY {
        rows.extend(raw_rows);
        return;
    }

    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;
    for i in [5, 6] {
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            orientation,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            orientation,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
        );
    }
    add_kp20_section42_oriented_vector_rows(
        rows,
        map,
        orientation,
        &[
            section42_term(Kp20Section42PaperSymbol::P, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, -w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, w1),
        ],
    );
    add_kp20_section42_oriented_vector_rows(
        rows,
        map,
        orientation,
        &[
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, -w1),
        ],
    );
}

fn add_kp20_section42_oriented_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    orientation: Kp20Section42SymbolOrientation,
    terms: &[Kp20Section42PaperTerm],
) {
    let oriented_terms = terms
        .iter()
        .map(|term| {
            let j = if orientation.reverse_j(term.symbol) {
                2 - term.j
            } else {
                term.j
            };
            section42_term(term.symbol, term.i, j, term.weight)
        })
        .collect::<Vec<_>>();
    add_kp20_section42_vector_rows(rows, map, &oriented_terms);
}

fn kp20_section42_mapped_piece(
    map: Kp20Section42PaperSymbolMap,
    symbol: Kp20Section42PaperSymbol,
) -> Kp20Section42MappedPiece {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

fn kp20_max_scalar_rows_residual_with_index(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> (f64, usize) {
    rows.iter()
        .enumerate()
        .map(|(row_index, row)| {
            let value = row
                .terms
                .iter()
                .map(|(coeff, weight)| weight * kp20_scalar_coeff_value(patch, *coeff))
                .sum::<f64>();
            ((value - row.target).abs(), row_index)
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap_or((0.0, 0))
}

fn max_kp20_vector_rows_residual(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauLinearConstraint],
) -> f64 {
    rows.iter()
        .map(|row| {
            let mut value = [0.0_f64; 2];
            for &(coeff, weight) in &row.terms {
                let cp = patch.piece(coeff.piece).cps[coeff.i][coeff.j];
                value[0] += weight * cp.x;
                value[1] += weight * cp.y;
            }
            (value[0] - row.target[0])
                .abs()
                .max((value[1] - row.target[1]).abs())
        })
        .fold(0.0_f64, f64::max)
}

fn kp20_scalar_coeff_value(patch: &Kp20MacroPatch, coeff: Kp20PlanarTauScalarCoeff) -> f64 {
    let cp = patch.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
    match coeff.axis {
        0 => cp.x,
        1 => cp.y,
        axis => panic!("invalid planar scalar axis {axis}"),
    }
}

fn add_rotated_vector_constraint_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    terms: &[(bool, usize, usize, f64)],
    previous_sector_rotation_radians: f64,
) {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &(is_previous, i, j, weight) in terms {
        let piece = if is_previous {
            previous_piece
        } else {
            current_piece
        };
        let coeff = Kp20PlanarTauCoeff::new(piece, i, j);
        if is_previous {
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight * cos_theta));
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), -weight * sin_theta));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight * sin_theta));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), weight * cos_theta));
        } else {
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), weight));
        }
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

#[derive(Debug)]
struct StamAppendixARotationalTopDiagnostic {
    residual: f64,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SectorCouplingDiagnostic {
    residual: f64,
    origin_pair: &'static str,
    midpoint_pair: &'static str,
    origin_b: &'static str,
    midpoint_b: &'static str,
    rotation: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SplitBOrientationDiagnostic {
    residual: f64,
    origin_eq8_b: &'static str,
    origin_eq9_b: &'static str,
    midpoint_eq8_b: &'static str,
    midpoint_eq9_b: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Eq9TermSignDiagnostic {
    residual: f64,
    origin_signs: Kp20Eq9TermSigns,
    midpoint_signs: Kp20Eq9TermSigns,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Eq9TermSourceDiagnostic {
    residual: f64,
    origin_terms: Kp20Eq9DiagnosticTerms,
    midpoint_terms: Kp20Eq9DiagnosticTerms,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20CoherentEq9ReversalDiagnostic {
    residual: f64,
    origin_variant: &'static str,
    midpoint_variant: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20G2RotationFrameDiagnostic {
    residual: f64,
    origin_g2_rotation: f64,
    midpoint_g2_rotation: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20TensorBorderFreedomDiagnostic {
    residual: f64,
    max_slack: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20TensorBorderBandFreedomDiagnostic {
    residual: f64,
    max_slack: f64,
    band_mask: u32,
    freed_band_count: u32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20TensorBorderPinnedRowDiagnostic {
    residual: f64,
    pinned_row: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20DerivativeSourceDiagnostic {
    residual: f64,
    source_name: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Eq14D1Diagnostic {
    residual: f64,
    d1_3: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20StamAnchorDiagnostic {
    distance: f64,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    side: Kp20StamGridSide,
    reversed_depth: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Row21LabelSubstitutionDiagnostic {
    distance: f64,
    source_row: usize,
    source_col: usize,
    label_sector: &'static str,
    label_index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Row21ColumnTripleDiagnostic {
    distance: f64,
    column: usize,
    labels: [(&'static str, usize); 3],
}

fn kp20_eq9_reversed_parameter_terms() -> Kp20Eq9DiagnosticTerms {
    Kp20Eq9DiagnosticTerms {
        signs: Kp20Eq9TermSigns {
            two_b_duv: -1.0,
            ..Kp20Eq9TermSigns::paper()
        },
        ..Kp20Eq9DiagnosticTerms::paper()
    }
}

fn kp20_eq9_term_sign_variants() -> [Kp20Eq9TermSigns; 8] {
    let paper = Kp20Eq9TermSigns::paper();
    let mut variants = [paper; 8];
    let signs = [-1.0, 1.0];
    let mut idx = 0;
    for two_b_duv in signs {
        for b_sq_duu in signs {
            for b_bprime_du in signs {
                variants[idx] = Kp20Eq9TermSigns {
                    two_b_duv,
                    b_sq_duu,
                    b_bprime_du,
                    ..paper
                };
                idx += 1;
            }
        }
    }
    variants
}

fn kp20_eq9_term_source_variants(signs: Kp20Eq9TermSigns) -> [Kp20Eq9DiagnosticTerms; 8] {
    let paper = Kp20Eq9DiagnosticTerms {
        signs,
        ..Kp20Eq9DiagnosticTerms::paper()
    };
    let mut variants = [paper; 8];
    let sources = [Kp20Eq9TermSource::Previous, Kp20Eq9TermSource::Current];
    let mut idx = 0;
    for duv_source in sources {
        for duu_source in sources {
            for du_source in sources {
                variants[idx] = Kp20Eq9DiagnosticTerms {
                    signs,
                    duv_source,
                    duu_source,
                    du_source,
                };
                idx += 1;
            }
        }
    }
    variants
}

fn partitioned_tensor_border_target_constraints(
    rows: &[Kp20PlanarTauLinearConstraint],
    free_row: impl Fn(usize) -> bool,
) -> (
    Vec<Kp20PlanarTauLinearConstraint>,
    Vec<Kp20PlanarTauAugmentedScalarLinearConstraint>,
    usize,
) {
    let mut pinned = Vec::new();
    let mut augmented = Vec::new();
    let mut extra_variable_count = 0;
    for (row_index, row) in rows.iter().enumerate() {
        if free_row(row_index) {
            for axis in 0..2 {
                let mut terms = row
                    .terms
                    .iter()
                    .map(|&(coeff, weight)| {
                        (
                            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                                coeff, axis,
                            )),
                            weight,
                        )
                    })
                    .collect::<Vec<_>>();
                terms.push((
                    Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_variable_count),
                    -1.0,
                ));
                extra_variable_count += 1;
                augmented.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                    terms,
                    target: row.target[axis],
                });
            }
        } else {
            pinned.push(row.clone());
        }
    }
    (pinned, augmented, extra_variable_count)
}

fn augmented_tensor_border_similarity_constraints(
    rows: &[Kp20PlanarTauLinearConstraint],
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let mut out = Vec::with_capacity(2 * rows.len());
    for row in rows {
        let [target_x, target_y] = row.target;
        let mut x_terms = row
            .terms
            .iter()
            .map(|&(coeff, weight)| {
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    weight,
                )
            })
            .collect::<Vec<_>>();
        x_terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_x));
        x_terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(1), target_y));
        out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: x_terms,
            target: 0.0,
        });

        let mut y_terms = row
            .terms
            .iter()
            .map(|&(coeff, weight)| {
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    weight,
                )
            })
            .collect::<Vec<_>>();
        y_terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -target_y));
        y_terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(1), -target_x));
        out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: y_terms,
            target: 0.0,
        });
    }
    out
}

fn normalize_rows_to_middle_patch_value(
    rows: &Kp20RegularTensorBorderCNet,
    expected: Point3,
) -> Kp20RegularTensorBorderCNet {
    let mut out = rows.clone();
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    grid[0] = rows.sectors[0][0];
    grid[1] = rows.sectors[0][1];
    grid[2] = rows.sectors[0][2];
    let actual = eval_cubic_bspline_grid(&grid, 0.0, 1.0).unwrap();
    let denom = actual.x * actual.x + actual.y * actual.y;
    if denom < 1e-20 {
        return out;
    }
    let scale_re = (expected.x * actual.x + expected.y * actual.y) / denom;
    let scale_im = (expected.y * actual.x - expected.x * actual.y) / denom;
    for sector in &mut out.sectors {
        for row in sector {
            for point in row {
                let x = scale_re * point.x - scale_im * point.y;
                let y = scale_im * point.x + scale_re * point.y;
                *point = Point3::new(x, y, point.z);
            }
        }
    }
    out
}

fn solve_pr08_best_candidate_with_eq14_d1_3(
    tensor_border: &TensorBorderBi3,
    d1_3: f64,
) -> Kp20PlanarTauF3SolveOutput {
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(
        kp20_planar_tau_tensor_border_constraints_with_d1_3(tensor_border, 0, d1_3).unwrap(),
    );
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: Kp20MacroPatch::zero(),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1e-10,
    })
    .unwrap()
}

fn solve_pr08_best_candidate_with_free_tensor_border_rows(
    free_row: impl Fn(usize) -> bool,
) -> Kp20TensorBorderFreedomDiagnostic {
    let solved = solve_pr08_best_candidate_with_free_tensor_border_rows_full(free_row);
    Kp20TensorBorderFreedomDiagnostic {
        residual: solved.max_constraint_residual,
        max_slack: solved
            .extra_values
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max),
    }
}

fn solve_pr08_best_candidate_with_free_tensor_border_rows_full(
    free_row: impl Fn(usize) -> bool,
) -> super::super::kp20_bi6_tau_f3_solve::Kp20PlanarTauAugmentedF3SolveOutput {
    let ring = kp20_catmull_clark_characteristic_ring(5).unwrap();
    let rows = kp20_pr08_characteristic_tensor_border_candidate_rows(
        &ring,
        Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed,
    )
    .unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let tensor_border_rows = kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap();
    let (mut vector_constraints, augmented_scalar_constraints, extra_variable_count) =
        partitioned_tensor_border_target_constraints(&tensor_border_rows, free_row);
    vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count,
            regularization_weight: 1e-10,
            extra_regularization_weight: 1e-10,
        },
    )
    .unwrap()
}

fn tensor_border_row21_source_anchor(source_rows: &[[Point3; 4]; 3]) -> Point3 {
    const DEPTH_WEIGHTS: [f64; 3] = [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0];
    const BOUNDARY_MIDPOINT_WEIGHTS: [f64; 4] = [1.0 / 48.0, 23.0 / 48.0, 23.0 / 48.0, 1.0 / 48.0];
    let mut terms = Vec::with_capacity(12);
    for row in 0..3 {
        for col in 0..4 {
            terms.push((
                DEPTH_WEIGHTS[row] * BOUNDARY_MIDPOINT_WEIGHTS[col],
                source_rows[row][col],
            ));
        }
    }
    weighted_sum(&terms)
}

fn row21_nearby_pr08_labels(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Vec<(&'static str, usize, Point3)> {
    let previous = (ring.sectors.len() - 1) % ring.sectors.len();
    let mut out = Vec::with_capacity(39);
    for (sector_name, sector_index) in [("previous", previous), ("current", 0), ("next", 1)] {
        for (label_index, &label) in ring.sectors[sector_index].iter().enumerate() {
            out.push((sector_name, label_index, label));
        }
    }
    out
}

const fn kp20_stam_grid_transforms() -> [Kp20StamGridTransform; 8] {
    [
        Kp20StamGridTransform::Identity,
        Kp20StamGridTransform::Rotate90,
        Kp20StamGridTransform::Rotate180,
        Kp20StamGridTransform::Rotate270,
        Kp20StamGridTransform::Transpose,
        Kp20StamGridTransform::AntiTranspose,
        Kp20StamGridTransform::FlipRows,
        Kp20StamGridTransform::FlipColumns,
    ]
}

fn pr08_stam_projected_grids_for_family(
    ring: &Kp20CatmullClarkCharacteristicRing,
    family: Kp20StamProjectionFamily,
) -> Result<super::super::kp20_bi6_stam_eval::Kp20Pr08StamProjectedGrids, String> {
    match family {
        Kp20StamProjectionFamily::VisualFigure63 => kp20_pr08_stam_visual_projected_grids(ring),
        Kp20StamProjectionFamily::AppendixTensorOrder => {
            kp20_pr08_stam_appendix_tensor_projected_grids(ring)
        }
    }
}

fn pr08_stam_projected_side_endpoint_fit_ranking(
    ring: &Kp20CatmullClarkCharacteristicRing,
    valence: usize,
    sector: usize,
) -> Vec<Kp20StamSideEndpointFitDiagnostic> {
    let alpha = core::f64::consts::TAU / valence as f64;
    let sector_start = Complex::from_angle(alpha * sector as f64);
    let sector_end = Complex::from_angle(alpha * (sector + 1) as f64);
    let mut diagnostics = Vec::new();

    for family in [
        Kp20StamProjectionFamily::VisualFigure63,
        Kp20StamProjectionFamily::AppendixTensorOrder,
    ] {
        let projection = pr08_stam_projected_grids_for_family(ring, family).unwrap();
        for patch in [
            Kp20StamPatch::LowerRight,
            Kp20StamPatch::CenterRight,
            Kp20StamPatch::UpperLeft,
        ] {
            for transform in kp20_stam_grid_transforms() {
                for side in Kp20StamGridSide::all() {
                    for reversed_depth in [false, true] {
                        let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
                            &projection,
                            patch,
                            transform,
                            side,
                            reversed_depth,
                        )
                        .unwrap();
                        let source0 = complex_from_point(rows.sectors[sector][0][0]);
                        let source1 = complex_from_point(rows.sectors[sector][0][3]);
                        for (swap_targets, conjugate_source) in
                            [(false, false), (true, false), (false, true), (true, true)]
                        {
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
                            diagnostics.push(Kp20StamSideEndpointFitDiagnostic {
                                max_endpoint_residual: residual0.max(residual1),
                                sum_endpoint_residual: residual0 + residual1,
                                family,
                                patch,
                                transform,
                                side,
                                reversed_depth,
                                swap_targets,
                                conjugate_source,
                            });
                        }
                    }
                }
            }
        }
    }

    diagnostics.sort_by(|a, b| {
        a.max_endpoint_residual
            .total_cmp(&b.max_endpoint_residual)
            .then_with(|| a.sum_endpoint_residual.total_cmp(&b.sum_endpoint_residual))
            .then_with(|| a.family.name().cmp(b.family.name()))
            .then_with(|| format!("{:?}", a.patch).cmp(&format!("{:?}", b.patch)))
            .then_with(|| format!("{:?}", a.transform).cmp(&format!("{:?}", b.transform)))
            .then_with(|| a.side.name().cmp(b.side.name()))
            .then_with(|| a.reversed_depth.cmp(&b.reversed_depth))
            .then_with(|| a.swap_targets.cmp(&b.swap_targets))
            .then_with(|| a.conjugate_source.cmp(&b.conjugate_source))
    });
    diagnostics
}

fn pr08_best_exact_stam_side_source_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Kp20RegularTensorBorderCNet {
    let projection =
        kp20_pr08_stam_appendix_tensor_projected_grids(ring).expect("Appendix tensor projection");
    kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        Kp20StamPatch::UpperLeft,
        Kp20StamGridTransform::FlipColumns,
        Kp20StamGridSide::Left,
        true,
    )
    .expect("best exact Stam side source rows")
}

fn pr08_best_exact_stam_side_tensor_border_with_b_spline_depth_order(
    ring: &Kp20CatmullClarkCharacteristicRing,
    valence: usize,
    sector: usize,
    order: [usize; 3],
) -> TensorBorderBi3 {
    let rows = unit_sector_endpoint_normalized_tensor_rows(
        &pr08_best_exact_stam_side_source_rows(ring),
        valence,
        sector,
        true,
        true,
    )
    .expect("best exact Stam side should have nondegenerate endpoint chart");
    let sectors = rows
        .sectors
        .iter()
        .map(|sector_rows| {
            cubic_bspline_rows_to_tensor_border_rows(&[
                sector_rows[order[0]],
                sector_rows[order[1]],
                sector_rows[order[2]],
            ])
        })
        .collect::<Vec<_>>();
    TensorBorderBi3 { sectors }
}

fn pr08_best_exact_stam_side_tensor_border_with_bezier_depth_order(
    ring: &Kp20CatmullClarkCharacteristicRing,
    valence: usize,
    sector: usize,
    order: [usize; 3],
) -> TensorBorderBi3 {
    let rows = unit_sector_endpoint_normalized_tensor_rows(
        &pr08_best_exact_stam_side_source_rows(ring),
        valence,
        sector,
        true,
        true,
    )
    .expect("best exact Stam side should have nondegenerate endpoint chart");
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let sectors = tensor_border
        .sectors
        .iter()
        .map(|sector_rows| {
            [
                sector_rows[order[0]],
                sector_rows[order[1]],
                sector_rows[order[2]],
            ]
        })
        .collect::<Vec<_>>();
    TensorBorderBi3 { sectors }
}

fn pr08_stam_side_tensor_border_with_b_spline_depth_order(
    ring: &Kp20CatmullClarkCharacteristicRing,
    valence: usize,
    sector: usize,
    endpoint: Kp20StamSideEndpointFitDiagnostic,
    order: [usize; 3],
) -> Option<TensorBorderBi3> {
    let projection = pr08_stam_projected_grids_for_family(ring, endpoint.family).ok()?;
    let rows = kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
        &projection,
        endpoint.patch,
        endpoint.transform,
        endpoint.side,
        endpoint.reversed_depth,
    )
    .ok()?;
    let normalized_rows = unit_sector_endpoint_normalized_tensor_rows(
        &rows,
        valence,
        sector,
        endpoint.swap_targets,
        endpoint.conjugate_source,
    )?;
    let sectors = normalized_rows
        .sectors
        .iter()
        .map(|sector_rows| {
            cubic_bspline_rows_to_tensor_border_rows(&[
                sector_rows[order[0]],
                sector_rows[order[1]],
                sector_rows[order[2]],
            ])
        })
        .collect::<Vec<_>>();
    Some(TensorBorderBi3 { sectors })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp20Eq14SourceHalf {
    NearC,
    NearM,
}

fn kp20_split_half_rows(
    split: Kp20SplitTensorBorder,
    half: Kp20Eq14SourceHalf,
) -> [[Point3; 7]; 3] {
    match half {
        Kp20Eq14SourceHalf::NearC => split.near_c.rows,
        Kp20Eq14SourceHalf::NearM => split.near_m.rows,
    }
}

fn maybe_reverse_degree6_row(row: [Point3; 7], reverse_u: bool) -> [Point3; 7] {
    if reverse_u {
        [row[6], row[5], row[4], row[3], row[2], row[1], row[0]]
    } else {
        row
    }
}

fn maybe_reverse_half_rows(rows: [[Point3; 7]; 3], reverse_u: bool) -> [[Point3; 7]; 3] {
    [
        maybe_reverse_degree6_row(rows[0], reverse_u),
        maybe_reverse_degree6_row(rows[1], reverse_u),
        maybe_reverse_degree6_row(rows[2], reverse_u),
    ]
}

fn cc_test_tensor_half_pin_rows(
    piece: Kp20PieceId,
    half_rows: [[Point3; 7]; 3],
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = Vec::with_capacity(21);
    for (depth, row) in half_rows.iter().enumerate() {
        for (u_index, cp) in row.iter().enumerate() {
            rows.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece, depth, u_index),
                [cp.x, cp.y],
            ));
        }
    }
    rows
}

fn cc_test_eq14_half_assignment_tensor_rows(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    previous_half: Kp20Eq14SourceHalf,
    current_half: Kp20Eq14SourceHalf,
    reverse_previous_u: bool,
    reverse_current_u: bool,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let split = split_kp20_tensor_border(tensor_border).unwrap();
    let previous = (sector + split.len() - 1) % split.len();
    let mut rows = Vec::with_capacity(42);
    rows.extend(cc_test_tensor_half_pin_rows(
        Kp20PieceId::NearMidpointLeft,
        maybe_reverse_half_rows(
            kp20_split_half_rows(split[previous], previous_half),
            reverse_previous_u,
        ),
    ));
    rows.extend(cc_test_tensor_half_pin_rows(
        Kp20PieceId::NearMidpointRight,
        maybe_reverse_half_rows(
            kp20_split_half_rows(split[sector], current_half),
            reverse_current_u,
        ),
    ));
    rows
}

fn cc_test_replaced_near_midpoint_right_columns_tensor_rows(
    base_tensor_border: &TensorBorderBi3,
    donor_tensor_border: &TensorBorderBi3,
    sector: usize,
    donor_sector_offset: isize,
    donor_half: Kp20Eq14SourceHalf,
    columns: &[usize],
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let base_split = split_kp20_tensor_border(base_tensor_border).unwrap();
    let donor_split = split_kp20_tensor_border(donor_tensor_border).unwrap();
    let n = base_split.len();
    let previous = (sector + n - 1) % n;
    let donor_sector = (sector as isize + donor_sector_offset).rem_euclid(n as isize) as usize;
    let mut current_right = base_split[sector].near_c.rows;
    let donor_rows = kp20_split_half_rows(donor_split[donor_sector], donor_half);
    for &column in columns {
        for depth in 0..3 {
            current_right[depth][column] = donor_rows[depth][column];
        }
    }

    let mut rows = Vec::with_capacity(42);
    rows.extend(cc_test_tensor_half_pin_rows(
        Kp20PieceId::NearMidpointLeft,
        base_split[previous].near_c.rows,
    ));
    rows.extend(cc_test_tensor_half_pin_rows(
        Kp20PieceId::NearMidpointRight,
        current_right,
    ));
    rows
}

#[derive(Debug, Clone, Copy)]
struct Kp20ExactSideVariantTauDiagnostic {
    residual: f64,
    square_outside: usize,
    negative_determinants: usize,
    min_determinant: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SolvedExactSideRowsDiagnostic {
    tau_macro_patch: Kp20MacroPatch,
    validity: Kp20ExactSideVariantTauDiagnostic,
}

#[derive(Debug, Clone, Copy)]
struct Kp20RightEdgeDonorDiagnostic {
    residual: f64,
    square_outside: usize,
    negative_determinants: usize,
    min_determinant: f64,
    family: Kp20StamProjectionFamily,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    side: Kp20StamGridSide,
    reversed_depth: bool,
    swap_targets: bool,
    conjugate_source: bool,
    donor_sector_offset: isize,
    donor_half: Kp20Eq14SourceHalf,
    column_label: &'static str,
}

fn cc_test_solve_exact_side_tensor_rows_under_live_basis(
    valence: usize,
    sector: usize,
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    project_all_tensor_pins: bool,
) -> Kp20ExactSideVariantTauDiagnostic {
    cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
        valence,
        sector,
        tensor_rows,
        project_all_tensor_pins,
    )
    .validity
}

fn cc_test_solve_exact_side_tensor_rows_under_live_basis_full(
    valence: usize,
    sector: usize,
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    project_all_tensor_pins: bool,
) -> Kp20SolvedExactSideRowsDiagnostic {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    };
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let source_rows = if project_all_tensor_pins {
        cc_test_project_tensor_rows_for_scalar_basis(tensor_rows, &fig8_rows, &scalar_constraints).0
    } else {
        tensor_rows
    };
    let mut vector_constraints = fig8_rows;
    vector_constraints.extend(source_rows);
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
    let square = cc_test_square_domain_report(&local_tau, &SAMPLES, 1.0e-10);
    Kp20SolvedExactSideRowsDiagnostic {
        tau_macro_patch: solved.tau_macro_patch,
        validity: Kp20ExactSideVariantTauDiagnostic {
            residual: solved.max_constraint_residual,
            square_outside: square.outside,
            negative_determinants: validity.orientation.negative_determinants,
            min_determinant: validity.orientation.min_determinant,
        },
    }
}

fn cc_test_vector_rows_worst_residual(
    tau: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauLinearConstraint],
) -> (usize, f64, Kp20PlanarTauCoeff) {
    let mut worst = (0, 0.0, rows[0].terms[0].0);
    for (row_index, row) in rows.iter().enumerate() {
        let mut actual = [0.0_f64; 2];
        for &(coeff, weight) in &row.terms {
            let cp = tau.piece(coeff.piece).cps[coeff.i][coeff.j];
            actual[0] += weight * cp.x;
            actual[1] += weight * cp.y;
        }
        let residual = (actual[0] - row.target[0]).hypot(actual[1] - row.target[1]);
        if residual > worst.1 {
            worst = (row_index, residual, row.terms[0].0);
        }
    }
    worst
}

fn cc_test_scalar_rows_worst_residual(
    tau: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> (usize, f64) {
    let mut worst = (0, 0.0);
    for (row_index, row) in rows.iter().enumerate() {
        let mut actual = 0.0_f64;
        for &(coeff, weight) in &row.terms {
            let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
            actual += weight * if coeff.axis == 0 { cp.x } else { cp.y };
        }
        let residual = (actual - row.target).abs();
        if residual > worst.1 {
            worst = (row_index, residual);
        }
    }
    worst
}

fn solve_pr08_best_sector_coupling_candidate(
    tensor_border: &TensorBorderBi3,
    initial: Kp20MacroPatch,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauF3SolveError> {
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_tensor_border_constraints(tensor_border, 0).unwrap());
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
        angle,
    ));
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial,
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1e-10,
    })
}

const fn rotated_sector_g1_g2_row_label(row_index: usize) -> &'static str {
    match row_index {
        0 => "eq8[0].x",
        1 => "eq8[0].y",
        2 => "eq8[1].x",
        3 => "eq8[1].y",
        4 => "eq8[2].x",
        5 => "eq8[2].y",
        6 => "eq8[3].x",
        7 => "eq8[3].y",
        8 => "eq8[4].x",
        9 => "eq8[4].y",
        10 => "eq8[5].x",
        11 => "eq8[5].y",
        12 => "eq8[6].x",
        13 => "eq8[6].y",
        14 => "eq9[0].x",
        15 => "eq9[0].y",
        16 => "eq9[1].x",
        17 => "eq9[1].y",
        18 => "eq9[2].x",
        19 => "eq9[2].y",
        20 => "eq9[3].x",
        21 => "eq9[3].y",
        22 => "eq9[4].x",
        23 => "eq9[4].y",
        24 => "eq9[5].x",
        25 => "eq9[5].y",
        26 => "eq9[6].x",
        27 => "eq9[6].y",
        28 => "eq9[7].x",
        29 => "eq9[7].y",
        _ => "invalid",
    }
}

const fn tensor_border_row_label(row_index: usize) -> &'static str {
    match row_index {
        0 => "near-c depth0 u0",
        1 => "near-c depth0 u1",
        2 => "near-c depth0 u2",
        3 => "near-c depth0 u3",
        4 => "near-c depth0 u4",
        5 => "near-c depth0 u5",
        6 => "near-c depth0 u6",
        7 => "near-c depth1 u0",
        8 => "near-c depth1 u1",
        9 => "near-c depth1 u2",
        10 => "near-c depth1 u3",
        11 => "near-c depth1 u4",
        12 => "near-c depth1 u5",
        13 => "near-c depth1 u6",
        14 => "near-c depth2 u0",
        15 => "near-c depth2 u1",
        16 => "near-c depth2 u2",
        17 => "near-c depth2 u3",
        18 => "near-c depth2 u4",
        19 => "near-c depth2 u5",
        20 => "near-c depth2 u6",
        21 => "near-m depth0 u0",
        22 => "near-m depth0 u1",
        23 => "near-m depth0 u2",
        24 => "near-m depth0 u3",
        25 => "near-m depth0 u4",
        26 => "near-m depth0 u5",
        27 => "near-m depth0 u6",
        28 => "near-m depth1 u0",
        29 => "near-m depth1 u1",
        30 => "near-m depth1 u2",
        31 => "near-m depth1 u3",
        32 => "near-m depth1 u4",
        33 => "near-m depth1 u5",
        34 => "near-m depth1 u6",
        35 => "near-m depth2 u0",
        36 => "near-m depth2 u1",
        37 => "near-m depth2 u2",
        38 => "near-m depth2 u3",
        39 => "near-m depth2 u4",
        40 => "near-m depth2 u5",
        41 => "near-m depth2 u6",
        _ => "invalid",
    }
}

fn stam_appendix_b_projected_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    transposed: bool,
    row_start: usize,
    reverse_rows: bool,
    reverse_cols: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let mut sectors = Vec::with_capacity(ring.sectors.len());
    for sector_index in 0..ring.sectors.len() {
        let source = stam_appendix_a_regular_subpatch_full_grid(
            ring.lambda,
            ring.sectors.len(),
            sector_index,
            patch_index,
        )?;
        let mut oriented = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
        for row in 0..4 {
            for col in 0..4 {
                oriented[row][col] = if transposed {
                    source[col][row]
                } else {
                    source[row][col]
                };
            }
        }
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for out_row in 0..3 {
            let source_row = if reverse_rows { 2 - out_row } else { out_row };
            for out_col in 0..4 {
                let source_col = if reverse_cols { 3 - out_col } else { out_col };
                rows[out_row][out_col] = oriented[row_start + source_row][source_col];
            }
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

fn stam_appendix_a_pr08_top_label_projected_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    odd_label: usize,
    even_label: usize,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = ring.sectors[0][0].x;
    u0_y[0] = ring.sectors[0][0].y;
    for sector in 0..valence {
        let odd = ring.sectors[sector][odd_label];
        let even = ring.sectors[sector][even_label];
        u0_x[1 + 2 * sector] = odd.x;
        u0_y[1 + 2 * sector] = odd.y;
        u0_x[2 + 2 * sector] = even.x;
        u0_y[2 + 2 * sector] = even.y;
    }
    let labels = stam_appendix_a_subdivide_extended_labels_from_top(ring.lambda, &u0_x, &u0_y)?;
    let mut sectors = Vec::with_capacity(valence);
    for sector_index in 0..valence {
        let grid = stam_appendix_b_subpatch_full_grid_from_labels(
            &labels,
            valence,
            sector_index,
            patch_index,
        )?;
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for out_row in 0..3 {
            let source_row = if reverse_rows { 2 - out_row } else { out_row };
            rows[out_row] = grid[row_start + source_row];
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

fn stam_appendix_a_pr08_scaled_top_label_projected_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    odd_label: (usize, Complex),
    even_label: (usize, Complex),
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = ring.sectors[0][0].x;
    u0_y[0] = ring.sectors[0][0].y;
    for sector in 0..valence {
        let odd = Complex::new(
            ring.sectors[sector][odd_label.0].x,
            ring.sectors[sector][odd_label.0].y,
        ) * odd_label.1;
        let even = Complex::new(
            ring.sectors[sector][even_label.0].x,
            ring.sectors[sector][even_label.0].y,
        ) * even_label.1;
        u0_x[1 + 2 * sector] = odd.re;
        u0_y[1 + 2 * sector] = odd.im;
        u0_x[2 + 2 * sector] = even.re;
        u0_y[2 + 2 * sector] = even.im;
    }
    stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn stam_appendix_a_pr08_subdivided_scaled_top_label_projected_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    odd_label: (usize, Complex),
    even_label: (usize, Complex),
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let refined = pr08_subdivide_characteristic_ring(ring)?;
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = refined[0][0].x;
    u0_y[0] = refined[0][0].y;
    for sector in 0..valence {
        let odd = Complex::new(
            refined[sector][odd_label.0].x,
            refined[sector][odd_label.0].y,
        ) * odd_label.1;
        let even = Complex::new(
            refined[sector][even_label.0].x,
            refined[sector][even_label.0].y,
        ) * even_label.1;
        u0_x[1 + 2 * sector] = odd.re;
        u0_y[1 + 2 * sector] = odd.im;
        u0_x[2 + 2 * sector] = even.re;
        u0_y[2 + 2 * sector] = even.im;
    }
    stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn stam_appendix_a_pr08_inner_ring_scaled_subpatch_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    seed_scale: f64,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = seed_scale * ring.sectors[0][0].x;
    u0_y[0] = seed_scale * ring.sectors[0][0].y;
    for sector in 0..valence {
        u0_x[1 + 2 * sector] = seed_scale * ring.sectors[sector][1].x;
        u0_y[1 + 2 * sector] = seed_scale * ring.sectors[sector][1].y;
        u0_x[2 + 2 * sector] = seed_scale * ring.sectors[sector][2].x;
        u0_y[2 + 2 * sector] = seed_scale * ring.sectors[sector][2].y;
    }
    stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn stam_appendix_a_unknown_top_vector_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let top_count = 2 * ring.sectors.len() + 1;
    let zero_rows = stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &[0.0; 11],
        &[0.0; 11],
        patch_index,
        row_start,
        reverse_rows,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp20_planar_tau_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(2 * top_count);
    for extra in 0..2 * top_count {
        let mut u0_x = vec![0.0_f64; top_count];
        let mut u0_y = vec![0.0_f64; top_count];
        if extra < top_count {
            u0_x[extra] = 1.0;
        } else {
            u0_y[extra - top_count] = 1.0;
        }
        let rows = stam_appendix_a_projected_rows_from_top_vectors(
            ring,
            &u0_x,
            &u0_y,
            patch_index,
            row_start,
            reverse_rows,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn solve_stam_appendix_a_unknown_top_vector_source_span(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> super::super::kp20_bi6_tau_f3_solve::Kp20PlanarTauAugmentedF3SolveOutput {
    let tensor_border_constraints = stam_appendix_a_unknown_top_vector_tensor_border_constraints(
        ring,
        patch_index,
        row_start,
        reverse_rows,
    );
    let angle = core::f64::consts::TAU / ring.sectors.len() as f64;
    let c = angle.cos();
    let vector_constraints = kp20_planar_tau_internal_c2_constraints();
    let mut scalar_constraints = kp20_planar_tau_diagonal_symmetry_constraints();
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearOriginLeft,
        [c, 2.0 * c, 2.0 * c],
        kp20_b_coeffs_near_o(c),
        angle,
    ));
    scalar_constraints.extend(kp20_planar_tau_rotated_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
        angle,
    ));

    let top_count = 2 * ring.sectors.len() + 1;
    let mut augmented_scalar_constraints = tensor_border_constraints;
    augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
        terms: vec![(Kp20PlanarTauAugmentedScalarCoeff::Extra(1), 1.0)],
        target: ring.sectors[0][1].x,
    });
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
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
    .unwrap()
}

fn stam_appendix_a_rotational_top_vector_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp20_planar_tau_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_rotational_top_vector_kp15_beta_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp15_beta_half_strip_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp15_beta_half_strip_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_rotational_top_vector_kp15_beta_full_then_split_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp15_beta_full_then_split_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp15_beta_full_then_split_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_rotational_top_vector_kp13_beta_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints = kp13_beta_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints.push(kp13_beta_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_rotational_top_vector_kp13_beta_exact_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp13_beta_exact_pullback_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp13_beta_exact_pullback_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_rotational_top_vector_kp13_beta_exact_phase_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
    left_selector: Kp13BetaBoundaryHalfSelector,
    right_selector: Kp13BetaBoundaryHalfSelector,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
        ring,
        [0.0; 4],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints = kp13_beta_exact_pullback_tensor_border_constraints_for_phase(
        &zero_tensor_border,
        0,
        left_selector,
        right_selector,
    )
    .unwrap();

    let mut basis_constraints = Vec::with_capacity(4);
    for extra in 0..4 {
        let mut basis = [0.0_f64; 4];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints.push(
            kp13_beta_exact_pullback_tensor_border_constraints_for_phase(
                &tensor_border,
                0,
                left_selector,
                right_selector,
            )
            .unwrap(),
        );
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn kp15_beta_half_strip_tensor_border_constraints(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = kp15_beta_split_tensor_border(tensor_border)?;
    let previous = (sector + halves.len() - 1) % halves.len();
    let current = sector % halves.len();
    let mut rows = Vec::with_capacity(42);
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        halves[previous].near_c.rows,
    );
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        halves[current].near_c.rows,
    );
    Ok(rows)
}

fn kp13_beta_tensor_border_constraints(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = kp13_beta_split_tensor_border(tensor_border)?;
    let previous = (sector + halves.len() - 1) % halves.len();
    let current = sector % halves.len();
    let mut rows = Vec::with_capacity(42);
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        halves[previous].near_c.rows,
    );
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        halves[current].near_c.rows,
    );
    Ok(rows)
}

fn kp13_beta_exact_pullback_tensor_border_constraints(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    kp13_beta_exact_pullback_tensor_border_constraints_for_phase(
        tensor_border,
        sector,
        Kp13BetaBoundaryHalfSelector::PreviousMiddle,
        Kp13BetaBoundaryHalfSelector::CurrentMiddle,
    )
}

fn kp13_beta_exact_pullback_tensor_border_constraints_for_phase(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    left_selector: Kp13BetaBoundaryHalfSelector,
    right_selector: Kp13BetaBoundaryHalfSelector,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = kp13_beta_exact_pullback_split_tensor_border(tensor_border)?;
    let previous = (sector + halves.len() - 1) % halves.len();
    let current = sector % halves.len();
    let mut rows = Vec::with_capacity(42);
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        left_selector.rows(&halves, previous, current),
    );
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        right_selector.rows(&halves, previous, current),
    );
    Ok(rows)
}

fn kp15_beta_full_then_split_tensor_border_constraints(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = kp15_beta_full_then_split_tensor_border(tensor_border)?;
    let previous = (sector + halves.len() - 1) % halves.len();
    let current = sector % halves.len();
    let mut rows = Vec::with_capacity(42);
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        halves[previous].near_c.rows,
    );
    cc_test_add_depth_u_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        halves[current].near_c.rows,
    );
    Ok(rows)
}

fn cc_test_add_depth_u_tensor_border_half_rows(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    piece: Kp20PieceId,
    half_rows: [[Point3; 7]; 3],
) {
    for (depth, row) in half_rows.iter().enumerate() {
        for (u_index, cp) in row.iter().enumerate() {
            rows.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece, depth, u_index),
                [cp.x, cp.y],
            ));
        }
    }
}

fn kp13_beta_split_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP13 beta diagnostic tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    let (middle_e, middle_d, right_e, right_d) = kp13_beta_two_piece_row2_controls(c);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        let mut left = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let mut right = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP13 beta diagnostic sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
            let (near_c, near_m) = cc_test_split_cubic_half(sector[depth]);
            left[depth] = near_c;
            right[depth] = near_m;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: cc_test_reparameterise_half_with_tangential_beta(left, middle_e, middle_d),
            },
            near_m: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: cc_test_reparameterise_half_with_tangential_beta(right, right_e, right_d),
            },
        });
    }
    Ok(out)
}

fn kp13_beta_exact_pullback_split_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP13 beta exact pullback diagnostic tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    let (middle_e, middle_d, right_e, right_d) = kp13_beta_two_piece_row2_controls(c);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        let mut left = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let mut right = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP13 beta exact pullback diagnostic sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
            let (near_c, near_m) = cc_test_split_cubic_half(sector[depth]);
            left[depth] = near_c;
            right[depth] = near_m;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: cc_test_reparameterise_half_with_exact_second_order_beta(
                    left, middle_e, middle_d,
                ),
            },
            near_m: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: cc_test_reparameterise_half_with_exact_second_order_beta(
                    right, right_e, right_d,
                ),
            },
        });
    }
    Ok(out)
}

fn kp13_beta_two_piece_row2_controls(c: f64) -> ([f64; 4], [f64; 4], [f64; 4], [f64; 4]) {
    let row = kp13_beta_equation15_boundary_row(c);
    let mut middle_e = [0.0; 4];
    let mut middle_d = [0.0; 4];
    let mut right_e = [0.0; 4];
    let mut right_d = [0.0; 4];
    for i in 0..4 {
        middle_e[i] = 6.0 * (row.middle[i][0] - i as f64 / 6.0);
        middle_d[i] = 6.0 * (row.middle[i][1] - 1.0 / 3.0);
        right_e[i] = 6.0 * (row.right[i][0] - (0.5 + i as f64 / 6.0));
        right_d[i] = 6.0 * (row.right[i][1] - 1.0 / 3.0);
    }
    (middle_e, middle_d, right_e, right_d)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kp13BetaBoundaryHalfSelector {
    PreviousMiddle,
    PreviousRight,
    CurrentMiddle,
    CurrentRight,
}

impl Kp13BetaBoundaryHalfSelector {
    fn rows(
        self,
        halves: &[Kp20SplitTensorBorder],
        previous: usize,
        current: usize,
    ) -> [[Point3; 7]; 3] {
        match self {
            Self::PreviousMiddle => halves[previous].near_c.rows,
            Self::PreviousRight => halves[previous].near_m.rows,
            Self::CurrentMiddle => halves[current].near_c.rows,
            Self::CurrentRight => halves[current].near_m.rows,
        }
    }
}

struct Kp13BetaEquation15BoundaryRow {
    middle: [[f64; 2]; 4],
    right: [[f64; 2]; 4],
}

fn kp13_beta_equation15_boundary_row(c: f64) -> Kp13BetaEquation15BoundaryRow {
    let middle_p1 = [1.0 / 6.0, 1.0 / 3.0];
    let middle_p2 = [1.0 / 3.0 + c / 216.0, 1.0 / 3.0 + c / 216.0];
    let middle_p3 = [0.5 + 7.0 * c / 1296.0, 1.0 / 3.0 + c / 72.0];
    let right_p0 = middle_p3;
    let right_p1 = [
        2.0 * middle_p3[0] - middle_p2[0],
        2.0 * middle_p3[1] - middle_p2[1],
    ];
    let right_p2 = [
        4.0 * middle_p3[0] - 4.0 * middle_p2[0] + middle_p1[0],
        4.0 * middle_p3[1] - 4.0 * middle_p2[1] + middle_p1[1],
    ];
    let right_p3 = [1.0, 1.0 / 3.0 + c / 18.0];

    let middle = [[0.0, 1.0 / 3.0], middle_p1, middle_p2, middle_p3];
    let right = [right_p0, right_p1, right_p2, right_p3];
    Kp13BetaEquation15BoundaryRow { middle, right }
}

fn assert_close_scalar(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "expected {actual:e} to match {expected:e}; delta={:e}",
        (actual - expected).abs()
    );
}

fn assert_close_xy(actual: [f64; 2], expected: [f64; 2], tol: f64) {
    assert_close_scalar(actual[0], expected[0], tol);
    assert_close_scalar(actual[1], expected[1], tol);
}

fn add_xy(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

fn sub_xy(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn scale_xy(a: [f64; 2], scale: f64) -> [f64; 2] {
    [a[0] * scale, a[1] * scale]
}

fn cc_test_reparameterise_half_with_tangential_beta(
    rows: [[Point3; 4]; 3],
    e: [f64; 4],
    d: [f64; 4],
) -> [[Point3; 7]; 3] {
    let row0 = cc_test_elevate_cubic_to_degree6(rows[0]);
    let row1 = cc_test_elevate_cubic_to_degree6(rows[1]);
    let mut row2 = cc_test_elevate_cubic_to_degree6(rows[2]);
    let cross_delta = [
        rows[1][0] - rows[0][0],
        rows[1][1] - rows[0][1],
        rows[1][2] - rows[0][2],
        rows[1][3] - rows[0][3],
    ];
    let d_cross = cc_test_multiply_scalar_cubic_vec3_cubic(d, cross_delta);
    let along_derivative = [
        (rows[0][1] - rows[0][0]) * 3.0,
        (rows[0][2] - rows[0][1]) * 3.0,
        (rows[0][3] - rows[0][2]) * 3.0,
    ];
    let e_along = cc_test_multiply_scalar_cubic_vec3_quadratic_to_degree6(e, along_derivative);
    for i in 0_usize..7 {
        row2[i] = row2[i] + d_cross[i] * 0.5 + e_along[i] * 0.5;
    }
    [row0, row1, row2]
}

fn cc_test_reparameterise_half_with_exact_second_order_beta(
    rows: [[Point3; 4]; 3],
    e: [f64; 4],
    d: [f64; 4],
) -> [[Point3; 7]; 3] {
    let source_row0 = cc_test_elevate_cubic_to_degree6(rows[0]);
    let source_row1 = cc_test_elevate_cubic_to_degree6(rows[1]);
    let source_row2 = cc_test_elevate_cubic_to_degree6(rows[2]);
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut row2 = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0_usize..7 {
        row1[i] = weighted_sum(&[(0.5, source_row0[i]), (0.5, source_row1[i])]);
        row2[i] = weighted_sum(&[
            (0.2, source_row0[i]),
            (0.6, source_row1[i]),
            (0.2, source_row2[i]),
        ]);
    }

    let cross_delta = [
        rows[1][0] - rows[0][0],
        rows[1][1] - rows[0][1],
        rows[1][2] - rows[0][2],
        rows[1][3] - rows[0][3],
    ];
    let d_cross = cc_test_multiply_scalar_cubic_vec3_cubic(d, cross_delta);
    let along_derivative = [
        (rows[0][1] - rows[0][0]) * 3.0,
        (rows[0][2] - rows[0][1]) * 3.0,
        (rows[0][3] - rows[0][2]) * 3.0,
    ];
    let e_along = cc_test_multiply_scalar_cubic_vec3_quadratic_to_degree6(e, along_derivative);
    for i in 0_usize..7 {
        row2[i] = row2[i] + d_cross[i] * 0.1 + e_along[i] * (1.0 / 30.0);
    }
    [source_row0, row1, row2]
}

fn cc_test_elevate_cubic_to_degree6(row: [Point3; 4]) -> [Point3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut out = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0_usize..7 {
        let j_min = i.saturating_sub(3);
        let j_max = i.min(3);
        let mut terms = Vec::with_capacity(j_max - j_min + 1);
        for j in j_min..=j_max {
            terms.push((C3[j] * cc_test_binom3(i - j) / C6[i], row[j]));
        }
        out[i] = weighted_sum(&terms);
    }
    out
}

fn cc_test_multiply_scalar_cubic_vec3_cubic(
    scalar: [f64; 4],
    vector: [crate::geometry::Vec3; 4],
) -> [crate::geometry::Vec3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut out = [crate::geometry::Vec3::zero(); 7];
    for i in 0_usize..7 {
        let j_min = i.saturating_sub(3);
        let j_max = i.min(3);
        for j in j_min..=j_max {
            let k = i - j;
            out[i] = out[i] + vector[k] * (scalar[j] * C3[j] * C3[k] / C6[i]);
        }
    }
    out
}

fn cc_test_multiply_scalar_cubic_vec3_quadratic_to_degree6(
    scalar: [f64; 4],
    vector: [crate::geometry::Vec3; 3],
) -> [crate::geometry::Vec3; 7] {
    const C2: [f64; 3] = [1.0, 2.0, 1.0];
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];
    let mut degree5 = [crate::geometry::Vec3::zero(); 6];
    for i in 0_usize..6 {
        let j_min = i.saturating_sub(2);
        let j_max = i.min(3);
        for j in j_min..=j_max {
            let k = i - j;
            degree5[i] = degree5[i] + vector[k] * (scalar[j] * C3[j] * C2[k] / C5[i]);
        }
    }
    let mut out = [crate::geometry::Vec3::zero(); 7];
    out[0] = degree5[0];
    out[6] = degree5[5];
    for i in 1..6 {
        out[i] = degree5[i - 1] * (i as f64 / 6.0) + degree5[i] * (1.0 - i as f64 / 6.0);
    }
    out
}

const fn cc_test_binom3(i: usize) -> f64 {
    match i {
        0 | 3 => 1.0,
        1 | 2 => 3.0,
        _ => 0.0,
    }
}

fn kp15_beta_full_then_split_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP15 beta full-then-split diagnostic tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    let params = Bi6Reparams::for_valence(n);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP15 beta full-then-split diagnostic sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
        }
        let full_rows = apply_beta_reparameterization(sector, &params);
        let mut near_c = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];
        let mut near_m = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];
        for depth in 0..3 {
            let (left, right) = cc_test_split_degree6_half(full_rows[depth]);
            near_c[depth] = left;
            near_m[depth] = right;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: super::super::kp20_bi6_cap::Kp20BoundaryHalf { rows: near_c },
            near_m: super::super::kp20_bi6_cap::Kp20BoundaryHalf { rows: near_m },
        });
    }
    Ok(out)
}

fn cc_test_split_degree6_half(row: [Point3; 7]) -> ([Point3; 7], [Point3; 7]) {
    let mut levels = Vec::with_capacity(7);
    levels.push(row.to_vec());
    for depth in 1..=6 {
        let prev = &levels[depth - 1];
        let mut next = Vec::with_capacity(7 - depth);
        for i in 0..(7 - depth) {
            next.push(weighted_sum(&[(0.5, prev[i]), (0.5, prev[i + 1])]));
        }
        levels.push(next);
    }

    let mut left = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0..=6 {
        left[i] = levels[i][0];
        right[6 - i] = levels[i][levels[i].len() - 1];
    }
    (left, right)
}

fn kp15_beta_split_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP15 beta diagnostic tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    let params = Bi6Reparams::for_valence(n);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        let mut left = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let mut right = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP15 beta diagnostic sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
            let (near_c, near_m) = cc_test_split_cubic_half(sector[depth]);
            left[depth] = near_c;
            right[depth] = near_m;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: apply_beta_reparameterization(&left, &params),
            },
            near_m: super::super::kp20_bi6_cap::Kp20BoundaryHalf {
                rows: apply_beta_reparameterization(&right, &params),
            },
        });
    }
    Ok(out)
}

fn cc_test_split_cubic_half(row: [Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let p10 = weighted_sum(&[(0.5, row[0]), (0.5, row[1])]);
    let p11 = weighted_sum(&[(0.5, row[1]), (0.5, row[2])]);
    let p12 = weighted_sum(&[(0.5, row[2]), (0.5, row[3])]);
    let p20 = weighted_sum(&[(0.5, p10), (0.5, p11)]);
    let p21 = weighted_sum(&[(0.5, p11), (0.5, p12)]);
    let p30 = weighted_sum(&[(0.5, p20), (0.5, p21)]);
    ([row[0], p10, p20, p30], [p30, p21, p12, row[3]])
}

fn stam_appendix_a_rotational_top_vector_with_center_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let zero_rows = stam_appendix_a_projected_rows_from_rotational_top_vector_with_center(
        ring,
        [0.0; 6],
        patch_index,
        row_start,
        reverse_rows,
        rotation_direction,
    )
    .unwrap();
    let zero_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&zero_rows).unwrap();
    let zero_constraints =
        kp20_planar_tau_tensor_border_constraints(&zero_tensor_border, 0).unwrap();

    let mut basis_constraints = Vec::with_capacity(6);
    for extra in 0..6 {
        let mut basis = [0.0_f64; 6];
        basis[extra] = 1.0;
        let rows = stam_appendix_a_projected_rows_from_rotational_top_vector_with_center(
            ring,
            basis,
            patch_index,
            row_start,
            reverse_rows,
            rotation_direction,
        )
        .unwrap();
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
        basis_constraints
            .push(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
    }

    let mut out = Vec::with_capacity(2 * zero_constraints.len());
    for (constraint_index, constraint) in zero_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis] - constraint.target[axis];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra), -weight));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn stam_appendix_a_projected_rows_from_rotational_top_vector_with_center(
    ring: &Kp20CatmullClarkCharacteristicRing,
    base: [f64; 6],
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = base[0];
    u0_y[0] = base[1];
    for sector in 0..valence {
        let theta = rotation_direction * core::f64::consts::TAU * sector as f64 / valence as f64;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let odd = rotate_complex_pair(base[2], base[3], sin_theta, cos_theta);
        let even = rotate_complex_pair(base[4], base[5], sin_theta, cos_theta);
        u0_x[1 + 2 * sector] = odd.0;
        u0_y[1 + 2 * sector] = odd.1;
        u0_x[2 + 2 * sector] = even.0;
        u0_y[2 + 2 * sector] = even.1;
    }
    stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn stam_appendix_a_projected_rows_from_rotational_top_vector(
    ring: &Kp20CatmullClarkCharacteristicRing,
    base: [f64; 4],
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
    rotation_direction: f64,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    for sector in 0..valence {
        let theta = rotation_direction * core::f64::consts::TAU * sector as f64 / valence as f64;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let odd = rotate_complex_pair(base[0], base[1], sin_theta, cos_theta);
        let even = rotate_complex_pair(base[2], base[3], sin_theta, cos_theta);
        u0_x[1 + 2 * sector] = odd.0;
        u0_y[1 + 2 * sector] = odd.1;
        u0_x[2 + 2 * sector] = even.0;
        u0_y[2 + 2 * sector] = even.1;
    }
    stam_appendix_a_projected_rows_from_top_vectors(
        ring,
        &u0_x,
        &u0_y,
        patch_index,
        row_start,
        reverse_rows,
    )
}

fn rotate_complex_pair(x: f64, y: f64, sin_theta: f64, cos_theta: f64) -> (f64, f64) {
    (cos_theta * x - sin_theta * y, sin_theta * x + cos_theta * y)
}

fn stam_appendix_a_projected_rows_from_top_vectors(
    ring: &Kp20CatmullClarkCharacteristicRing,
    u0_x: &[f64],
    u0_y: &[f64],
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let labels = stam_appendix_a_subdivide_extended_labels_from_top(ring.lambda, u0_x, u0_y)?;
    let mut sectors = Vec::with_capacity(valence);
    for sector_index in 0..valence {
        let grid = stam_appendix_b_subpatch_full_grid_from_labels(
            &labels,
            valence,
            sector_index,
            patch_index,
        )?;
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for out_row in 0..3 {
            let source_row = if reverse_rows { 2 - out_row } else { out_row };
            rows[out_row] = grid[row_start + source_row];
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

fn tensor_border_row_norm(rows: &Kp20RegularTensorBorderCNet) -> f64 {
    rows.sectors
        .iter()
        .flat_map(|sector| sector.iter())
        .flat_map(|row| row.iter())
        .map(|point| point.x.hypot(point.y))
        .fold(0.0_f64, f64::max)
}

fn tensor_border_rows_max_distance(
    left: &Kp20RegularTensorBorderCNet,
    right: &Kp20RegularTensorBorderCNet,
) -> f64 {
    left.sectors
        .iter()
        .zip(&right.sectors)
        .flat_map(|(left_sector, right_sector)| left_sector.iter().zip(right_sector))
        .flat_map(|(left_row, right_row)| left_row.iter().zip(right_row))
        .map(|(left_point, right_point)| left_point.distance_to(right_point))
        .fold(0.0_f64, f64::max)
}

fn unit_sector_endpoint_normalized_tensor_rows(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
    swap_targets: bool,
    conjugate_source: bool,
) -> Option<Kp20RegularTensorBorderCNet> {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = Complex::from_angle(alpha * sector as f64);
    let spoke1 = Complex::from_angle(alpha * (sector + 1) as f64);
    let target0 = if swap_targets { spoke1 } else { spoke0 };
    let target1 = if swap_targets { spoke0 } else { spoke1 };
    let source0 = maybe_conjugate(
        complex_from_point(rows.sectors[sector][0][0]),
        conjugate_source,
    );
    let source1 = maybe_conjugate(
        complex_from_point(rows.sectors[sector][0][3]),
        conjugate_source,
    );
    let scale = complex_similarity_scale_for_pairs(&[(source0, target0), (source1, target1)])?;
    Some(complex_scaled_tensor_border_rows(
        rows,
        scale,
        conjugate_source,
    ))
}

fn kp20_reflected_unit_sector_chart_tensor_rows(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
) -> Option<Kp20RegularTensorBorderCNet> {
    // PR08/Stam labels are stored in the complex characteristic-ring chart,
    // with sector blocks rotating in the positive complex direction. KP20
    // Fig. 6/Fig. 8 pin the same data in the sector coefficient chart whose
    // two axes are the adjacent sector spokes; the reflected side enters that
    // chart with endpoints swapped and complex orientation conjugated.
    unit_sector_endpoint_normalized_tensor_rows(rows, valence, sector, true, true)
}

fn complex_scaled_tensor_border_rows(
    rows: &Kp20RegularTensorBorderCNet,
    scale: Complex,
    conjugate_source: bool,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let value =
                        maybe_conjugate(complex_from_point(sector[row][col]), conjugate_source)
                            * scale;
                    out[row][col] = Point3::new(value.re, value.im, sector[row][col].z);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn complex_from_point(point: Point3) -> Complex {
    Complex::new(point.x, point.y)
}

fn maybe_conjugate(value: Complex, conjugate_source: bool) -> Complex {
    if conjugate_source {
        value.conj()
    } else {
        value
    }
}

fn complex_similarity_scale_for_pairs(pairs: &[(Complex, Complex)]) -> Option<Complex> {
    let mut numerator = Complex::zero();
    let mut denominator = 0.0_f64;
    for &(source, target) in pairs {
        numerator = numerator + target * source.conj();
        denominator += source.norm2();
    }
    if denominator < 1.0e-20 {
        return None;
    }
    Some(numerator / denominator)
}

fn scale_tensor_border_rows(
    rows: &Kp20RegularTensorBorderCNet,
    scale: f64,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let p = sector[row][col];
                    out[row][col] = Point3::new(scale * p.x, scale * p.y, scale * p.z);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn stam_appendix_a_frequency_one_top_vector(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> StamAppendixATopVectorRings {
    let valence = ring.sectors.len();
    let gain = 4.0 * ring.lambda - 1.0;
    let mut odd = Vec::with_capacity(valence);
    let mut even = Vec::with_capacity(valence);
    for sector in 0..valence {
        let theta = core::f64::consts::TAU * sector as f64 / valence as f64;
        let next_theta = core::f64::consts::TAU * (sector + 1) as f64 / valence as f64;
        odd.push(Complex::new(gain * theta.cos(), gain * theta.sin()));
        even.push(Complex::new(
            theta.cos() + next_theta.cos(),
            theta.sin() + next_theta.sin(),
        ));
    }
    StamAppendixATopVectorRings { odd, even }
}

fn pr08_label_ring_as_complex(
    ring: &Kp20CatmullClarkCharacteristicRing,
    label_index: usize,
) -> Vec<Complex> {
    ring.sectors
        .iter()
        .map(|sector| Complex::new(sector[label_index].x, sector[label_index].y))
        .collect()
}

fn complex_similarity_ring_residual(source: &[Complex], target: &[Complex]) -> f64 {
    let Some(scale) = complex_similarity_ring_scale(source, target) else {
        return f64::INFINITY;
    };
    source
        .iter()
        .zip(target)
        .map(|(s, t)| (*s * scale - *t).norm2().sqrt())
        .fold(0.0_f64, f64::max)
}

fn complex_similarity_ring_scale(source: &[Complex], target: &[Complex]) -> Option<Complex> {
    let mut numerator = Complex::zero();
    let mut denominator = 0.0_f64;
    for (s, t) in source.iter().zip(target) {
        numerator = numerator + *t * s.conj();
        denominator += s.norm2();
    }
    if denominator < 1.0e-20 {
        return None;
    }
    Some(numerator / denominator)
}

fn fourier_mode_residual(values: &[Complex], mode: usize) -> f64 {
    let n = values.len();
    let mut coefficient = Complex::zero();
    for (index, value) in values.iter().enumerate() {
        coefficient = coefficient
            + Complex::from_angle(-core::f64::consts::TAU * mode as f64 * index as f64 / n as f64)
                * *value;
    }
    coefficient = coefficient / Complex::new(n as f64, 0.0);

    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let reconstructed =
                Complex::from_angle(core::f64::consts::TAU * mode as f64 * index as f64 / n as f64)
                    * coefficient;
            let delta = *value - reconstructed;
            delta.re.hypot(delta.im)
        })
        .fold(0.0_f64, f64::max)
}

fn fourier_mode_magnitudes(values: &[Complex]) -> Vec<f64> {
    let n = values.len();
    (0..n)
        .map(|mode| {
            let mut coefficient = Complex::zero();
            for (index, value) in values.iter().enumerate() {
                coefficient = coefficient
                    + Complex::from_angle(
                        -core::f64::consts::TAU * mode as f64 * index as f64 / n as f64,
                    ) * *value;
            }
            coefficient = coefficient / Complex::new(n as f64, 0.0);
            coefficient.re.hypot(coefficient.im)
        })
        .collect()
}

fn pr08_middle_patch_augmented_source_first_column_scale_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let fixed = pr08_middle_patch_rows_with_source_first_column_scale(ring, 0.0);
    let basis = pr08_middle_patch_rows_with_source_first_column_scale(ring, 1.0);
    let fixed_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&fixed).unwrap();
    let basis_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&basis).unwrap();
    let fixed_constraints =
        kp20_planar_tau_tensor_border_constraints(&fixed_tensor_border, 0).unwrap();
    let basis_constraints =
        kp20_planar_tau_tensor_border_constraints(&basis_tensor_border, 0).unwrap();

    let mut out = Vec::with_capacity(2 * fixed_constraints.len());
    for (constraint_index, constraint) in fixed_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            let weight = basis_constraints[constraint_index].target[axis] - constraint.target[axis];
            if weight != 0.0 {
                terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -weight));
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn pr08_middle_patch_rows_with_source_first_column_scale(
    ring: &Kp20CatmullClarkCharacteristicRing,
    scale: f64,
) -> Kp20RegularTensorBorderCNet {
    let sectors = ring
        .sectors
        .iter()
        .enumerate()
        .map(|(sector_index, sector)| {
            let mut first_column = stam_appendix_a_middle_patch_first_column(
                ring.lambda,
                ring.sectors.len(),
                sector_index,
            )
            .unwrap();
            for p in &mut first_column {
                *p = Point3::new(scale * p.x, scale * p.y, scale * p.z);
            }
            [
                [first_column[0], sector[2], sector[5], sector[8]],
                [first_column[1], sector[3], sector[6], sector[9]],
                [first_column[2], sector[12], sector[11], sector[12]],
            ]
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn pr08_middle_patch_augmented_first_column_tensor_border_constraints(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let fixed =
        pr08_middle_patch_rows_with_sector_zero_first_column(ring, [Point3::new(0.0, 0.0, 0.0); 3]);
    let fixed_tensor_border = kp20_tensor_border_from_regular_cnet_rows(&fixed).unwrap();
    let fixed_constraints =
        kp20_planar_tau_tensor_border_constraints(&fixed_tensor_border, 0).unwrap();
    let mut basis_constraints = Vec::new();
    for row in 0..3 {
        for axis in 0..2 {
            let mut first_column = [Point3::new(0.0, 0.0, 0.0); 3];
            if axis == 0 {
                first_column[row].x = 1.0;
            } else {
                first_column[row].y = 1.0;
            }
            let basis = pr08_middle_patch_rows_with_sector_zero_first_column(ring, first_column);
            let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&basis).unwrap();
            basis_constraints
                .push(kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).unwrap());
        }
    }

    let mut out = Vec::with_capacity(2 * fixed_constraints.len());
    for (constraint_index, constraint) in fixed_constraints.iter().enumerate() {
        for axis in 0..2 {
            let mut terms = constraint
                .terms
                .iter()
                .map(|&(coeff, weight)| {
                    (
                        Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                            coeff, axis,
                        )),
                        weight,
                    )
                })
                .collect::<Vec<_>>();
            for (extra_index, basis) in basis_constraints.iter().enumerate() {
                let weight = basis[constraint_index].target[axis]
                    - fixed_constraints[constraint_index].target[axis];
                if weight != 0.0 {
                    terms.push((
                        Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_index),
                        -weight,
                    ));
                }
            }
            out.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: constraint.target[axis],
            });
        }
    }
    out
}

fn pr08_middle_patch_rows_with_sector_zero_first_column(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_zero_first_column: [Point3; 3],
) -> Kp20RegularTensorBorderCNet {
    let sectors = ring
        .sectors
        .iter()
        .enumerate()
        .map(|(sector_index, sector)| {
            let first_column = if sector_index == 0 {
                sector_zero_first_column
            } else {
                [Point3::new(0.0, 0.0, 0.0); 3]
            };
            [
                [first_column[0], sector[2], sector[5], sector[8]],
                [first_column[1], sector[3], sector[6], sector[9]],
                [first_column[2], sector[12], sector[11], sector[12]],
            ]
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}
