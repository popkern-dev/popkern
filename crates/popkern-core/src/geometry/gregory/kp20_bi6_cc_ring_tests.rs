#![allow(unused_imports)]
// Paper-notation test module: loop indices, helper arities, and tuple shapes
// mirror KP15/KP20/PR08 subscripts and equation forms; keep the index form.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

pub(crate) use super::super::kp_polar_boundary_input::TensorBorderBi3;
#[allow(unused_imports)]
use super::super::kp_polar_shared::weighted_sum;
pub(crate) use super::super::kp20_bi6_cap::{
    Kp20MacroPatch, Kp20PieceId, Kp20SplitTensorBorder, split_kp20_tensor_border,
};
pub(crate) use super::super::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, cubic_bspline_rows_to_tensor_border_rows,
    kp20_tensor_border_from_regular_cnet_rows,
};
pub(crate) use super::super::kp20_bi6_guide_sample::Kp20AffineGuideMap;
pub(crate) use super::super::kp20_bi6_local_constraint::{
    Kp20Eq9DiagnosticTerms, Kp20Eq9TermSigns, Kp20Eq9TermSource,
};
pub(crate) use super::super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
pub(crate) use super::super::kp20_bi6_piece_maps::{
    kp20_paper_quadrant_sector_piece_maps, kp20_planar_bi6_tau_orientation_report,
    kp20_planar_bi6_tau_validity_report, sample_kp20_macro_patch_from_tau_family_maps,
};
pub(crate) use super::super::kp20_bi6_stam_eval::{
    Kp20StamGridSide, Kp20StamGridTransform, Kp20StamPatch, eval_cubic_bspline_grid,
    kp20_pr08_middle_patch_normalization_value, kp20_pr08_stam_appendix_tensor_projected_grids,
    kp20_pr08_stam_projected_transformed_tensor_border_side_rows,
    kp20_pr08_stam_visual_projected_grids,
};
pub(crate) use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedScalarCoeff,
    Kp20PlanarTauAugmentedScalarLinearConstraint, Kp20PlanarTauCoeff,
    Kp20PlanarTauF3ScalarSolveInput, Kp20PlanarTauF3SolveOutput, Kp20PlanarTauLinearConstraint,
    Kp20PlanarTauScalarCoeff, Kp20PlanarTauScalarLinearConstraint,
    kp20_planar_tau_diagonal_symmetry_constraints,
    kp20_planar_tau_fig8_indicated_averaging_constraints, kp20_planar_tau_internal_c2_constraints,
    kp20_planar_tau_rotated_sector_g1_g2_constraints,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms,
    kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b,
    kp20_planar_tau_tensor_border_constraints,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
    solve_kp20_planar_tau_f3_with_scalar_constraints,
};
pub(crate) use super::super::kp20_bi6_tau_solve_assemble::{
    Kp20PlanarTauPaperRowFamily, Kp20PlanarTauPaperSolveInput,
    kp20_planar_tau_rotational_paper_family_residuals, kp20_planar_tau_rotational_paper_row_set,
    kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs,
    kp20_planar_tau_rotational_section42_source_projection_report,
    kp20_planar_tau_sector_local_diagonal_symmetry_constraints,
    solve_kp20_planar_tau_paper_sector_f3, solve_kp20_planar_tau_rotational_paper_sector_f3,
};
pub(crate) use super::*;
pub(crate) use crate::geometry::gregory::{
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
        let sign_r = if (du_order - r).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        for s in 0..=dv_order {
            let sign_s = if (dv_order - s).is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
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

#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Eq9TermSourceDiagnostic {
    residual: f64,
    origin_terms: Kp20Eq9DiagnosticTerms,
    midpoint_terms: Kp20Eq9DiagnosticTerms,
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

#[path = "kp20_bi6_cc_ring_tests_part1.rs"]
mod part1;
#[path = "kp20_bi6_cc_ring_tests_part2.rs"]
mod part2;
#[path = "kp20_bi6_cc_ring_tests_part3.rs"]
mod part3;
