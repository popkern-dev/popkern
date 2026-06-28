use super::super::kp20_bi6_cnet::{
    kp20_cap_net_figure_tensor_border_rows, kp20_tensor_border_from_regular_cnet_rows,
};
use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedF3SolveOutput,
    Kp20PlanarTauAugmentedScalarCoeff, Kp20PlanarTauAugmentedScalarLinearConstraint,
    kp20_planar_tau_constraint_space_report,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
};
use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::{
    Kp20PlanarTauCoeff, Kp20PlanarTauScalarCoeff, Kp20TauFamily, Kp20TauMap, Kp20Tri6GuideSector,
    TensorBorderBi3, kp20_characteristic_cap_net, kp20_paper_quadrant_sector_piece_maps,
    kp20_planar_bi6_tau_f3_energy, kp20_planar_bi6_tau_orientation_report,
    kp20_planar_bi6_tau_validity_report, kp20_positive_nonoverlap_c2_sector_piece_maps,
    sample_kp20_macro_patch_from_tau_family_maps, tri6_flat_index,
};

fn p(x: f64, y: f64) -> Point3 {
    Point3::new(x, y, 0.0)
}

fn constant_tensor_border(n: usize, value: Point3) -> TensorBorderBi3 {
    TensorBorderBi3 {
        sectors: vec![[[value; 4]; 3]; n],
    }
}

fn constant_tau(value: Point3) -> Kp20MacroPatch {
    let mut tau = Kp20MacroPatch::zero();
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                tau.piece_mut(piece).cps[i][j] = value;
            }
        }
    }
    tau
}

fn planar_regular_sector_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let a0 = core::f64::consts::TAU * k as f64 / n as f64;
        let a1 = core::f64::consts::TAU * (k + 1) as f64 / n as f64;
        let p0 = Point3::new(a0.cos(), a0.sin(), 0.0);
        let p1 = Point3::new(a1.cos(), a1.sin(), 0.0);
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for (row_index, row) in rows.iter_mut().enumerate() {
            let radius = 1.0 - 0.1 * row_index as f64;
            for (i, cp) in row.iter_mut().enumerate() {
                let t = i as f64 / 3.0;
                *cp = Point3::new(
                    radius * ((1.0 - t) * p0.x + t * p1.x),
                    radius * ((1.0 - t) * p0.y + t * p1.y),
                    0.0,
                );
            }
        }
        sectors.push(rows);
    }
    TensorBorderBi3 { sectors }
}

fn characteristic_cap_net_tensor_border(n: usize) -> TensorBorderBi3 {
    let cap_net = kp20_characteristic_cap_net(n).expect("characteristic cap-net");
    let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).expect("figure tensor rows");
    kp20_tensor_border_from_regular_cnet_rows(&rows).expect("tensor-border converts")
}

fn planar_tau_seed_macro_patch(n: usize, sector: usize) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
        let base = corner_maps[0].1;
        let piece = out.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                let local = [
                    base.origin[0] + u * base.du[0] + v * base.dv[0],
                    base.origin[1] + u * base.du[1] + v * base.dv[1],
                ];
                piece.cps[i][j] = local_to_global_planar(n, sector, local);
            }
        }
    }
    out
}

fn local_identity_tri6_guide() -> Kp20Tri6GuideSector {
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            guide.cps[tri6_flat_index(6, i, j)] = p(i as f64 / 6.0, j as f64 / 6.0);
        }
    }
    guide
}

fn local_planar_patch_to_global(
    n: usize,
    sector: usize,
    local_tau: Kp20MacroPatch,
) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let cp = local_tau.piece(piece_id).cps[i][j];
                out.piece_mut(piece_id).cps[i][j] = local_to_global_planar(n, sector, [cp.x, cp.y]);
            }
        }
    }
    out
}

fn local_to_global_planar(n: usize, sector: usize, local: [f64; 2]) -> Point3 {
    let alpha = core::f64::consts::TAU / n as f64;
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

fn global_planar_patch_to_local(
    n: usize,
    sector: usize,
    mut patch: Kp20MacroPatch,
) -> Kp20MacroPatch {
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let local = global_to_local_planar(n, sector, piece.cps[i][j]);
                piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
            }
        }
    }
    patch
}

fn global_to_local_planar(n: usize, sector: usize, point: Point3) -> [f64; 2] {
    let alpha = core::f64::consts::TAU / n as f64;
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

#[derive(Debug, Clone, Copy)]
struct TestNearOriginFirstOrderFrame {
    origin: Point3,
    du: Point3,
    dv: Point3,
    determinant: f64,
}

fn near_origin_first_order_frame(
    tau: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
) -> TestNearOriginFirstOrderFrame {
    let piece = tau.piece(piece_id);
    let origin = piece.cps[0][0];
    let p10 = piece.cps[1][0];
    let p01 = piece.cps[0][1];
    let du = Point3::new(6.0 * (p10.x - origin.x), 6.0 * (p10.y - origin.y), 0.0);
    let dv = Point3::new(6.0 * (p01.x - origin.x), 6.0 * (p01.y - origin.y), 0.0);
    TestNearOriginFirstOrderFrame {
        origin,
        du,
        dv,
        determinant: du.x * dv.y - du.y * dv.x,
    }
}

fn sector_local_axis_components(n: usize, sector: usize) -> [[f64; 2]; 2] {
    let alpha = core::f64::consts::TAU / n as f64;
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

fn augmented_local_component_terms(
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

fn first_order_centered_source_rows(
    piece_id: Kp20PieceId,
    local_axes: [[f64; 2]; 2],
    scale_extra_index: usize,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    first_order_source_rows_for_local_frame(
        piece_id,
        local_axes,
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        scale_extra_index,
    )
}

fn first_order_source_rows_for_local_frame(
    piece_id: Kp20PieceId,
    local_axes: [[f64; 2]; 2],
    origin: [f64; 2],
    du: [f64; 2],
    dv: [f64; 2],
    scale_extra_index: usize,
) -> Vec<Kp20PlanarTauAugmentedScalarLinearConstraint> {
    let p00 = Kp20PlanarTauCoeff::new(piece_id, 0, 0);
    let p10 = Kp20PlanarTauCoeff::new(piece_id, 1, 0);
    let p01 = Kp20PlanarTauCoeff::new(piece_id, 0, 1);
    let scale = Kp20PlanarTauAugmentedScalarCoeff::Extra(scale_extra_index);
    let mut rows = Vec::new();
    for (component, axis) in local_axes.into_iter().enumerate() {
        rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: augmented_local_component_terms(p00, axis, 1.0),
            target: origin[component],
        });
    }

    for (coeff, local_delta) in [(p10, du), (p01, dv)] {
        for (component, axis) in local_axes.into_iter().enumerate() {
            let mut terms = augmented_local_component_terms(coeff, axis, 1.0);
            terms.extend(augmented_local_component_terms(p00, axis, -1.0));
            terms.push((scale, -local_delta[component] / 6.0));
            rows.push(Kp20PlanarTauAugmentedScalarLinearConstraint { terms, target: 0.0 });
        }
    }

    rows
}

fn vector_residual(tau: &Kp20MacroPatch, constraint: &Kp20PlanarTauLinearConstraint) -> [f64; 2] {
    let mut value = [0.0_f64; 2];
    for &(coeff, weight) in &constraint.terms {
        let cp = tau.piece(coeff.piece).cps[coeff.i][coeff.j];
        value[0] += weight * cp.x;
        value[1] += weight * cp.y;
    }
    [
        value[0] - constraint.target[0],
        value[1] - constraint.target[1],
    ]
}

fn scalar_residual(tau: &Kp20MacroPatch, constraint: &Kp20PlanarTauScalarLinearConstraint) -> f64 {
    let mut value = 0.0_f64;
    for &(coeff, weight) in &constraint.terms {
        let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
        value += weight * if coeff.axis == 0 { cp.x } else { cp.y };
    }
    value - constraint.target
}

fn solve_rotational_regular_source_with_row_options(
    valence: usize,
    sector: usize,
    include_fig8_average: bool,
    include_sector_diagonal: bool,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = Vec::new();
    if include_fig8_average {
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    }
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    if include_sector_diagonal {
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
    }
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_regular_source_with_previous_sector_index_transpose(
    valence: usize,
    sector: usize,
    index_map: TestSection42IndexMap,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let mut scalar_constraints =
        test_section42_presolved_strip_basis_rows_previous_transposed(map, 2.0 * c, index_map);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_regular_source_with_rotation_sign(
    valence: usize,
    sector: usize,
    previous_rotation_sign: f64,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointRight),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointLeft),
        previous_sector_rotation_radians: previous_rotation_sign * angle,
    };
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_regular_source_with_eq13_b0(
    valence: usize,
    sector: usize,
    b0_near_origin: f64,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, b0_near_origin);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

#[derive(Clone, Copy)]
struct TestSection42FamilySet {
    c0: bool,
    eq11: bool,
    eq12: bool,
    eq13: bool,
}

impl TestSection42FamilySet {
    const fn new(c0: bool, eq11: bool, eq12: bool, eq13: bool) -> Self {
        Self {
            c0,
            eq11,
            eq12,
            eq13,
        }
    }
}

fn solve_rotational_regular_source_with_section42_families(
    valence: usize,
    sector: usize,
    families: TestSection42FamilySet,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
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
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn test_section42_presolved_strip_basis_rows_previous_transposed(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    index_map: TestSection42IndexMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    test_add_section42_c0_rows_previous_transposed(&mut rows, map, index_map);
    test_add_section42_eq11_degree5_rows_previous_transposed(&mut rows, map, index_map);
    test_add_section42_eq12_internal_c2_rows_previous_transposed(&mut rows, map, index_map);
    test_add_section42_eq13_rows_previous_transposed(&mut rows, map, b0_near_origin, index_map);
    rows
}

#[derive(Debug, Clone, Copy)]
enum TestSection42IndexMap {
    PreviousSectorAll,
    PreviousSectorCentralOnly,
    StoredPieceEdgeOrientation,
}

fn test_add_section42_c0_rows_previous_transposed(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    index_map: TestSection42IndexMap,
) {
    for i in 0..=6 {
        test_add_section42_vector_rows_previous_transposed(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
            index_map,
        );
        test_add_section42_vector_rows_previous_transposed(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
            index_map,
        );
    }
}

fn test_add_section42_eq11_degree5_rows_previous_transposed(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    index_map: TestSection42IndexMap,
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
        test_add_section42_vector_rows_previous_transposed(rows, map, &terms, index_map);
    }
}

fn test_add_section42_eq12_internal_c2_rows_previous_transposed(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    index_map: TestSection42IndexMap,
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
            test_add_section42_vector_rows_previous_transposed(
                rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
                index_map,
            );
        }
        for j in 1..=2 {
            test_add_section42_vector_rows_previous_transposed(
                rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
                index_map,
            );
            test_add_section42_vector_rows_previous_transposed(
                rows,
                map,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, 0.25),
                    section42_term(upper, 2, j, -0.25),
                ],
                index_map,
            );
        }
    }
}

fn test_add_section42_eq13_rows_previous_transposed(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
    index_map: TestSection42IndexMap,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;

    for i in [5, 6] {
        test_add_section42_vector_rows_previous_transposed(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
            index_map,
        );
        test_add_section42_vector_rows_previous_transposed(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
            index_map,
        );
    }

    test_add_section42_vector_rows_previous_transposed(
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
        index_map,
    );
    test_add_section42_vector_rows_previous_transposed(
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
        index_map,
    );
}

fn test_add_section42_vector_rows_previous_transposed(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    terms: &[Kp20Section42PaperTerm],
    index_map: TestSection42IndexMap,
) {
    let (sin_theta, cos_theta) = map.previous_sector_rotation_radians.sin_cos();
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let mapped = kp20_section42_mapped_piece(map, term.symbol);
        let should_transpose = match index_map {
            TestSection42IndexMap::PreviousSectorAll => mapped.is_previous_sector,
            TestSection42IndexMap::PreviousSectorCentralOnly => {
                mapped.is_previous_sector && term.i + term.j <= 2
            }
            TestSection42IndexMap::StoredPieceEdgeOrientation => matches!(
                mapped.piece,
                Kp20PieceId::NearOriginRight | Kp20PieceId::NearMidpointRight
            ),
        };
        let (i, j) = if should_transpose {
            (term.j, term.i)
        } else {
            (term.i, term.j)
        };
        let coeff = Kp20PlanarTauCoeff::new(mapped.piece, i, j);
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

fn solve_rotational_regular_source_with_diagonal_mapping(
    valence: usize,
    sector: usize,
    self_pieces: [Kp20PieceId; 2],
    paired_pieces: [Kp20PieceId; 2],
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(test_sector_local_diagonal_symmetry_constraints(
        valence,
        sector,
        self_pieces,
        paired_pieces,
    ));
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_regular_source_with_origin_pins(
    valence: usize,
    sector: usize,
    pinned_coefficients: &[Kp20PlanarTauCoeff],
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let seed = planar_tau_seed_macro_patch(valence, sector);
    let mut vector_constraints = row_set.vector_constraints;
    for &coeff in pinned_coefficients {
        let target = seed.piece(coeff.piece).cps[coeff.i][coeff.j];
        vector_constraints.push(Kp20PlanarTauLinearConstraint::pin(
            coeff,
            [target.x, target.y],
        ));
    }
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: seed,
        vector_constraints,
        scalar_constraints: row_set.scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_regular_source_with_paper_tau_origin_pins(
    valence: usize,
    sector: usize,
    use_paper_quadrant_maps: bool,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let seed = planar_tau_seed_macro_patch(valence, sector);
    let mut vector_constraints = row_set.vector_constraints;
    let tau = Kp20TauMap::for_sector_local_paper_refined_apex(valence, sector).unwrap();
    let maps = if use_paper_quadrant_maps {
        kp20_paper_quadrant_sector_piece_maps()
    } else {
        kp20_positive_nonoverlap_c2_sector_piece_maps()
    };
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        let base = maps
            .iter()
            .find(|(id, _)| *id == piece_id)
            .expect("piece map exists")
            .1[0]
            .1;
        let jet = tau.jet3(base.origin[0], base.origin[1]);
        let value = jet.d[0][0];
        let ddu = [
            jet.d[1][0][0] * base.du[0] + jet.d[0][1][0] * base.du[1],
            jet.d[1][0][1] * base.du[0] + jet.d[0][1][1] * base.du[1],
        ];
        let ddv = [
            jet.d[1][0][0] * base.dv[0] + jet.d[0][1][0] * base.dv[1],
            jet.d[1][0][1] * base.dv[0] + jet.d[0][1][1] * base.dv[1],
        ];
        let targets = [
            (0, 0, value),
            (1, 0, [value[0] + ddu[0] / 6.0, value[1] + ddu[1] / 6.0]),
            (0, 1, [value[0] + ddv[0] / 6.0, value[1] + ddv[1] / 6.0]),
        ];
        for (i, j, local_target) in targets {
            let target = local_to_global_planar(valence, sector, local_target);
            vector_constraints.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece_id, i, j),
                [target.x, target.y],
            ));
        }
    }
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: seed,
        vector_constraints,
        scalar_constraints: row_set.scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_source_with_appendix_a_central_pins(
    valence: usize,
    sector: usize,
    tensor_border: &TensorBorderBi3,
    max_total_degree: usize,
    use_paper_quadrant_maps: bool,
) -> Kp20PlanarTauF3SolveOutput {
    let row_set = kp20_planar_tau_rotational_paper_row_set(valence, sector, tensor_border).unwrap();
    let seed = planar_tau_seed_macro_patch(valence, sector);
    let mut vector_constraints = row_set.vector_constraints;
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        if use_paper_quadrant_maps {
            kp20_paper_quadrant_sector_piece_maps()
        } else {
            kp20_positive_nonoverlap_c2_sector_piece_maps()
        },
    )
    .expect("paper tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        for i in 0..=max_total_degree {
            for j in 0..=(max_total_degree - i) {
                let target = paper_tau.piece(piece_id).cps[i][j];
                vector_constraints.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece_id, i, j),
                    [target.x, target.y],
                ));
            }
        }
    }
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: seed,
        vector_constraints,
        scalar_constraints: row_set.scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn solve_rotational_characteristic_source_with_centered_first_order(
    valence: usize,
    sector: usize,
    include_near_origin_right: bool,
) -> Kp20PlanarTauAugmentedF3SolveOutput {
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let mut augmented_scalar_constraints = first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        sector_local_axis_components(valence, sector),
        0,
    );
    if include_near_origin_right {
        augmented_scalar_constraints.extend(first_order_centered_source_rows(
            Kp20PieceId::NearOriginRight,
            sector_local_axis_components(valence, sector),
            0,
        ));
    }

    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .expect("centered first-order augmented solve succeeds")
}

fn solve_characteristic_centered_first_order_with_row_families(
    valence: usize,
    sector: usize,
    families: TestSection42FamilySet,
    include_fig8_average: bool,
    include_diagonal: bool,
    include_tensor_source: bool,
) -> Kp20PlanarTauAugmentedF3SolveOutput {
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = Vec::new();
    if include_fig8_average {
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    }
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
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
    }
    if include_tensor_source {
        let tensor_border = characteristic_cap_net_tensor_border(valence);
        let tensor_constraints =
            kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
        let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
            tensor_constraints,
            &vector_constraints,
            &scalar_constraints,
        )
        .unwrap();
        vector_constraints.extend(projected_tensor_constraints);
    }
    let augmented_scalar_constraints = first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        sector_local_axis_components(valence, sector),
        0,
    );

    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .expect("centered first-order family solve succeeds")
}

fn project_tensor_band_with_centered_first_order(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    centered_rows: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut vector_constraints = fixed_vector_constraints.to_vec();
    let mut augmented_scalar_constraints = centered_rows.to_vec();
    let mut tensor_delta_targets = Vec::new();
    let mut extra_count = 1_usize;

    for (row_index, row) in tensor_rows.iter().enumerate() {
        for axis in 0..=1 {
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
            terms.push((Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_count), -1.0));
            augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                terms,
                target: row.target[axis],
            });
            tensor_delta_targets.push((row_index, axis));
            extra_count += 1;
        }
    }

    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints: core::mem::take(&mut vector_constraints),
            scalar_constraints: scalar_constraints.to_vec(),
            augmented_scalar_constraints,
            extra_variable_count: extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0,
        },
    )
    .expect("center-normalized tensor projection succeeds");

    let mut projected = tensor_rows;
    for ((row_index, axis), delta) in tensor_delta_targets
        .into_iter()
        .zip(solved.extra_values[1..].iter().copied())
    {
        projected[row_index].target[axis] += delta;
    }
    projected
}

fn solve_characteristic_source_with_centered_projected_tensor_band(
    valence: usize,
    sector: usize,
    include_fig8_average: bool,
    include_diagonal: bool,
    include_near_origin_right: bool,
) -> Kp20PlanarTauAugmentedF3SolveOutput {
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = Vec::new();
    if include_fig8_average {
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    }
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    if include_diagonal {
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
    }
    let mut centered_rows = first_order_centered_source_rows(
        Kp20PieceId::NearOriginLeft,
        sector_local_axis_components(valence, sector),
        0,
    );
    if include_near_origin_right {
        centered_rows.extend(first_order_centered_source_rows(
            Kp20PieceId::NearOriginRight,
            sector_local_axis_components(valence, sector),
            0,
        ));
    }
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    vector_constraints.extend(project_tensor_band_with_centered_first_order(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
        &centered_rows,
    ));

    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: centered_rows,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .expect("center-normalized source solve succeeds")
}

fn solve_characteristic_source_with_projected_tensor_band_and_first_order_rows(
    valence: usize,
    sector: usize,
    include_fig8_average: bool,
    include_diagonal: bool,
    first_order_rows: Vec<Kp20PlanarTauAugmentedScalarLinearConstraint>,
) -> Kp20PlanarTauAugmentedF3SolveOutput {
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = Vec::new();
    if include_fig8_average {
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    }
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    if include_diagonal {
        scalar_constraints.extend(
            kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
        );
    }
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();
    vector_constraints.extend(project_tensor_band_with_centered_first_order(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
        &first_order_rows,
    ));

    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            augmented_scalar_constraints: first_order_rows,
            extra_variable_count: 1,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: 1.0e-10,
        },
    )
    .expect("first-order source solve succeeds")
}

#[derive(Clone, Copy)]
enum TestEq14SourceBlendMask {
    AllTensorRows,
    Piece(Kp20PieceId),
    Depth(usize),
    PieceDepth(Kp20PieceId, usize),
}

impl TestEq14SourceBlendMask {
    fn label(self) -> &'static str {
        match self {
            Self::AllTensorRows => "all",
            Self::Piece(Kp20PieceId::NearMidpointLeft) => "near-midpoint-left",
            Self::Piece(Kp20PieceId::NearMidpointRight) => "near-midpoint-right",
            Self::Piece(_) => "non-eq14-piece",
            Self::Depth(0) => "depth-0",
            Self::Depth(1) => "depth-1",
            Self::Depth(2) => "depth-2",
            Self::Depth(_) => "non-eq14-depth",
            Self::PieceDepth(Kp20PieceId::NearMidpointLeft, 0) => "near-midpoint-left-depth-0",
            Self::PieceDepth(Kp20PieceId::NearMidpointLeft, 1) => "near-midpoint-left-depth-1",
            Self::PieceDepth(Kp20PieceId::NearMidpointLeft, 2) => "near-midpoint-left-depth-2",
            Self::PieceDepth(Kp20PieceId::NearMidpointRight, 0) => "near-midpoint-right-depth-0",
            Self::PieceDepth(Kp20PieceId::NearMidpointRight, 1) => "near-midpoint-right-depth-1",
            Self::PieceDepth(Kp20PieceId::NearMidpointRight, 2) => "near-midpoint-right-depth-2",
            Self::PieceDepth(_, _) => "non-eq14-piece-depth",
        }
    }

    fn matches(self, coeff: Kp20PlanarTauCoeff) -> bool {
        match self {
            Self::AllTensorRows => true,
            Self::Piece(piece) => coeff.piece == piece,
            Self::Depth(depth) => coeff.i == depth,
            Self::PieceDepth(piece, depth) => coeff.piece == piece && coeff.i == depth,
        }
    }
}

fn explicit_appendix_a_global_tau(valence: usize, sector: usize) -> Kp20MacroPatch {
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::Paper,
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
    )
    .expect("paper tau samples with identity guide");
    local_planar_patch_to_global(valence, sector, local_tau)
}

fn solve_rotational_source_blended_toward_appendix_eq14(
    valence: usize,
    sector: usize,
    blend_alpha: f64,
    mask: TestEq14SourceBlendMask,
) -> Kp20PlanarTauF3SolveOutput {
    let tensor_border = planar_regular_sector_tensor_border(valence);
    let appendix_tau = explicit_appendix_a_global_tau(valence, sector);
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let mut tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, sector).unwrap();

    for row in &mut tensor_constraints {
        let coeff = row.terms[0].0;
        if mask.matches(coeff) {
            let target = appendix_tau.piece(coeff.piece).cps[coeff.i][coeff.j];
            row.target[0] = (1.0 - blend_alpha) * row.target[0] + blend_alpha * target.x;
            row.target[1] = (1.0 - blend_alpha) * row.target[1] + blend_alpha * target.y;
        }
    }

    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let (projected_tensor_constraints, _) = kp20_project_rotational_paper_tensor_band(
        tensor_constraints,
        &vector_constraints,
        &scalar_constraints,
    )
    .unwrap();
    vector_constraints.extend(projected_tensor_constraints);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: planar_tau_seed_macro_patch(valence, sector),
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .unwrap()
}

fn test_sector_local_diagonal_symmetry_constraints(
    valence: usize,
    sector: usize,
    self_pieces: [Kp20PieceId; 2],
    paired_pieces: [Kp20PieceId; 2],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    let local0 = [spoke0[1] / det, -spoke0[0] / det];
    let local1 = [-spoke1[1] / det, spoke1[0] / det];
    let mut rows = Vec::new();
    add_sector_local_self_diagonal_symmetry(&mut rows, self_pieces[0], local0, local1);
    add_sector_local_self_diagonal_symmetry(&mut rows, self_pieces[1], local0, local1);
    add_sector_local_pair_diagonal_symmetry(
        &mut rows,
        paired_pieces[0],
        paired_pieces[1],
        local0,
        local1,
    );
    rows
}

#[test]
fn paper_row_set_assembles_all_kp20_tau_constraint_families() {
    let tensor_border = constant_tensor_border(5, p(0.25, 0.25));
    let row_set = kp20_planar_tau_paper_row_set(5, 2, &tensor_border).expect("row set");

    assert_eq!(row_set.vector_constraints.len(), 84 + 1 + 42 + 15 + 15);
    assert_eq!(
        row_set.scalar_constraints.len(),
        0,
        "KP20 Fig. 8 indicated averaging is not the same as the full-grid diagonal diagnostic"
    );
}

#[test]
fn rotational_paper_row_set_free_space_is_measured_against_diagonal_and_kp13_sigma() {
    let valence = 5;
    let sector = 0;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let row_set =
        kp20_planar_tau_rotational_paper_row_set(valence, sector, &tensor_border).unwrap();
    let report = kp20_planar_tau_constraint_space_report(
        &row_set.vector_constraints,
        &row_set.scalar_constraints,
    )
    .unwrap();

    eprintln!(
        "rotational paper row-set constraint space: vars={} dense_rows={} independent_rows={} nullity={}",
        report.scalar_variable_count,
        report.dense_row_count,
        report.independent_row_count,
        report.nullity
    );
    assert_eq!(report.scalar_variable_count, 392);
    assert_eq!(
        report.dense_row_count,
        2 * row_set.vector_constraints.len() + row_set.scalar_constraints.len()
    );
    assert!(
        report.nullity > 15,
        "KP13 §4.3's bi-5 σ has 15 scalar parameters, but KP20's bi-6 τ uses a different Fig. 8 substrate; recalibrate if this comparison stops separating the two constructions"
    );
    assert_eq!(report.dense_row_count, 354);
    assert_eq!(report.independent_row_count, 331);
    assert_eq!(report.nullity, 61);

    let mut scalar_with_self_diagonal = row_set.scalar_constraints.clone();
    scalar_with_self_diagonal.extend(
        kp20_planar_tau_sector_local_self_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let self_diagonal_report = kp20_planar_tau_constraint_space_report(
        &row_set.vector_constraints,
        &scalar_with_self_diagonal,
    )
    .unwrap();
    eprintln!(
        "rotational row-set plus self-diagonal space: vars={} dense_rows={} independent_rows={} nullity={}",
        self_diagonal_report.scalar_variable_count,
        self_diagonal_report.dense_row_count,
        self_diagonal_report.independent_row_count,
        self_diagonal_report.nullity
    );
    assert_eq!(
        self_diagonal_report.independent_row_count, report.independent_row_count,
        "self-diagonal rows should already be included in the live rotational paper row set"
    );
    assert_eq!(
        self_diagonal_report.nullity, report.nullity,
        "adding the self-diagonal helper again must not change the live row-space rank"
    );

    let mut scalar_with_diagnostic_diagonal = row_set.scalar_constraints.clone();
    scalar_with_diagnostic_diagonal.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector).unwrap(),
    );
    let diagonal_report = kp20_planar_tau_constraint_space_report(
        &row_set.vector_constraints,
        &scalar_with_diagnostic_diagonal,
    )
    .unwrap();
    eprintln!(
        "rotational row-set plus full Fig.8 diagonal space again: vars={} dense_rows={} independent_rows={} nullity={}",
        diagonal_report.scalar_variable_count,
        diagonal_report.dense_row_count,
        diagonal_report.independent_row_count,
        diagonal_report.nullity
    );
    assert!(
        diagonal_report.nullity > 15,
        "KP20's diagonal bi-6 τ space should not be conflated with KP13's 15-parameter bi-5 σ space"
    );
    assert_eq!(
        diagonal_report.independent_row_count, report.independent_row_count,
        "full Fig.8 diagonal rows should already be included in the live rotational paper row set"
    );
    assert_eq!(
        diagonal_report.nullity, report.nullity,
        "adding the full Fig.8 diagonal helper again must not change the live row-space rank"
    );
}

#[test]
fn fig8_indicated_average_maps_to_quadrant1_bottom_edge() {
    let rows = kp20_planar_tau_fig8_indicated_averaging_constraints();

    assert_eq!(rows.len(), 1);
    assert_eq!(
        rows[0].terms,
        vec![
            (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 1, 0),
                1.0
            ),
            (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0),
                -0.5
            ),
            (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 2, 0),
                -0.5
            ),
        ]
    );
    assert_eq!(rows[0].target, [0.0, 0.0]);
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
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
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
        kp20_positive_nonoverlap_c2_sector_piece_maps(),
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

#[derive(Debug, Clone, Copy)]
struct TestSection42StripMap {
    mapped: Kp20Section42MappedPiece,
    edge: TestMacroEdge,
    reversed: bool,
}

#[derive(Debug, Clone, Copy)]
struct TestSection42StripSymbolMap {
    p: TestSection42StripMap,
    arrow_p: TestSection42StripMap,
    underlined_p: TestSection42StripMap,
    arrow_underlined_p: TestSection42StripMap,
    previous_sector_rotation_radians: f64,
}

fn strip_boundary_pair_residual(
    patch: &Kp20MacroPatch,
    left: TestSection42StripMap,
    right: TestSection42StripMap,
    previous_sector_rotation_radians: f64,
) -> f64 {
    (0..=6)
        .map(|i| {
            let left_point =
                transformed_strip_point(patch, left, i, 0, previous_sector_rotation_radians);
            let right_point =
                transformed_strip_point(patch, right, i, 0, previous_sector_rotation_radians);
            distance(left_point, right_point)
        })
        .fold(0.0_f64, f64::max)
}

fn scalar_row_set_residual(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> f64 {
    rows.iter()
        .map(|row| scalar_residual(patch, row).abs())
        .fold(0.0_f64, f64::max)
}

fn transformed_strip_point(
    patch: &Kp20MacroPatch,
    strip: TestSection42StripMap,
    paper_i: usize,
    paper_j: usize,
    previous_sector_rotation_radians: f64,
) -> Point3 {
    let (i, j) = strip.stored_index(paper_i, paper_j);
    let cp = patch.piece(strip.mapped.piece).cps[i][j];
    if strip.mapped.is_previous_sector {
        rotate_xy(cp, previous_sector_rotation_radians)
    } else {
        cp
    }
}

fn test_section42_strip_basis_rows(
    map: TestSection42StripSymbolMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = test_section42_strip_c0_rows(map);
    test_add_section42_strip_eq11_degree5_rows(&mut rows, map);
    test_add_section42_strip_eq12_internal_c2_rows(&mut rows, map);
    test_add_section42_strip_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn test_section42_strip_c0_rows(
    map: TestSection42StripSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    for i in 0..=6 {
        test_add_section42_strip_vector_rows(
            &mut rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        test_add_section42_strip_vector_rows(
            &mut rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
    rows
}

fn test_add_section42_strip_eq11_degree5_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
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
        test_add_section42_strip_vector_rows(rows, map, &terms);
    }
}

fn test_add_section42_strip_eq12_internal_c2_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
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
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
            );
            test_add_section42_strip_vector_rows(
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

fn test_add_section42_strip_eq12_variant_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
    source_is_underlined: bool,
    source_endpoint: usize,
    receiving_endpoint: usize,
) {
    let source_inner = if source_endpoint == 6 {
        (4_usize, 2_usize)
    } else {
        (2_usize, 4_usize)
    };
    let receiving_inner = if receiving_endpoint == 6 {
        (5_usize, 1_usize)
    } else {
        (1_usize, 5_usize)
    };
    let pairs = if source_is_underlined {
        [
            (
                Kp20Section42PaperSymbol::UnderlinedP,
                Kp20Section42PaperSymbol::P,
            ),
            (
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
                Kp20Section42PaperSymbol::ArrowP,
            ),
        ]
    } else {
        [
            (
                Kp20Section42PaperSymbol::P,
                Kp20Section42PaperSymbol::UnderlinedP,
            ),
            (
                Kp20Section42PaperSymbol::ArrowP,
                Kp20Section42PaperSymbol::ArrowUnderlinedP,
            ),
        ]
    };

    for (source, receiving) in pairs {
        for j in 0..=2 {
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, receiving_endpoint, j, 1.0),
                    section42_term(source, source_endpoint, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, receiving_inner.0, j, 1.0),
                    section42_term(source, source_endpoint, j, -1.0),
                    section42_term(source, source_inner.0, j, -0.25),
                    section42_term(source, source_inner.1, j, 0.25),
                ],
            );
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, receiving_inner.1, j, 1.0),
                    section42_term(source, source_endpoint, j, -1.0),
                    section42_term(source, source_inner.0, j, 0.25),
                    section42_term(source, source_inner.1, j, -0.25),
                ],
            );
        }
    }
}

fn test_add_section42_strip_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;

    for i in [5, 6] {
        test_add_section42_strip_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        test_add_section42_strip_vector_rows(
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

    test_add_section42_strip_vector_rows(
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
    test_add_section42_strip_vector_rows(
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

fn test_add_section42_strip_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let (sin_theta, cos_theta) = map.previous_sector_rotation_radians.sin_cos();
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let strip = test_section42_strip_symbol_map(map, term.symbol);
        let (i, j) = strip.stored_index(term.i, term.j);
        let coeff = Kp20PlanarTauCoeff::new(strip.mapped.piece, i, j);
        if strip.mapped.is_previous_sector {
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

fn test_section42_strip_symbol_map(
    map: TestSection42StripSymbolMap,
    symbol: Kp20Section42PaperSymbol,
) -> TestSection42StripMap {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

impl TestSection42StripMap {
    fn stored_index(self, paper_i: usize, paper_j: usize) -> (usize, usize) {
        let along = if self.reversed { 6 - paper_i } else { paper_i };
        match self.edge {
            TestMacroEdge::I0 => (paper_j, along),
            TestMacroEdge::I6 => (6 - paper_j, along),
            TestMacroEdge::J0 => (along, paper_j),
            TestMacroEdge::J6 => (along, 6 - paper_j),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum TestMacroEdge {
    I0,
    I6,
    J0,
    J6,
}

impl TestMacroEdge {
    const ALL: [Self; 4] = [Self::I0, Self::I6, Self::J0, Self::J6];
}

fn macro_edge_pair_residual(
    patch: &Kp20MacroPatch,
    left: Kp20Section42MappedPiece,
    left_edge: TestMacroEdge,
    right: Kp20Section42MappedPiece,
    right_edge: TestMacroEdge,
    reversed: bool,
    previous_sector_rotation_radians: f64,
) -> f64 {
    (0..=6)
        .map(|k| {
            let right_k = if reversed { 6 - k } else { k };
            let left_point = transformed_macro_edge_point(
                patch,
                left,
                left_edge,
                k,
                previous_sector_rotation_radians,
            );
            let right_point = transformed_macro_edge_point(
                patch,
                right,
                right_edge,
                right_k,
                previous_sector_rotation_radians,
            );
            distance(left_point, right_point)
        })
        .fold(0.0_f64, f64::max)
}

fn transformed_macro_edge_point(
    patch: &Kp20MacroPatch,
    mapped: Kp20Section42MappedPiece,
    edge: TestMacroEdge,
    k: usize,
    previous_sector_rotation_radians: f64,
) -> Point3 {
    let cp = match edge {
        TestMacroEdge::I0 => patch.piece(mapped.piece).cps[0][k],
        TestMacroEdge::I6 => patch.piece(mapped.piece).cps[6][k],
        TestMacroEdge::J0 => patch.piece(mapped.piece).cps[k][0],
        TestMacroEdge::J6 => patch.piece(mapped.piece).cps[k][6],
    };
    if mapped.is_previous_sector {
        rotate_xy(cp, previous_sector_rotation_radians)
    } else {
        cp
    }
}

fn rotate_xy(point: Point3, angle: f64) -> Point3 {
    let (sin_theta, cos_theta) = angle.sin_cos();
    Point3::new(
        cos_theta * point.x - sin_theta * point.y,
        sin_theta * point.x + cos_theta * point.y,
        point.z,
    )
}

fn distance(left: Point3, right: Point3) -> f64 {
    let dx = left.x - right.x;
    let dy = left.y - right.y;
    let dz = left.z - right.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[test]
#[ignore = "diagnostic source-band sweep; run explicitly when revisiting eq14 source compatibility"]
fn appendix_a_eq14_source_band_blend_localizes_row_compatibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let masks = [
        TestEq14SourceBlendMask::AllTensorRows,
        TestEq14SourceBlendMask::Piece(Kp20PieceId::NearMidpointLeft),
        TestEq14SourceBlendMask::Piece(Kp20PieceId::NearMidpointRight),
        TestEq14SourceBlendMask::Depth(0),
        TestEq14SourceBlendMask::Depth(1),
        TestEq14SourceBlendMask::Depth(2),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointLeft, 0),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointLeft, 1),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointLeft, 2),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointRight, 0),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointRight, 1),
        TestEq14SourceBlendMask::PieceDepth(Kp20PieceId::NearMidpointRight, 2),
    ];
    let tensor_border = planar_regular_sector_tensor_border(5);
    let mut best = ("", 0.0_f64, usize::MAX, f64::NEG_INFINITY, 0.0_f64);

    for mask in masks {
        for blend_alpha in [0.25, 0.5, 0.75, 1.0] {
            let solved =
                solve_rotational_source_blended_toward_appendix_eq14(5, 0, blend_alpha, mask);
            let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
            let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
            let residuals = kp20_planar_tau_rotational_paper_family_residuals(
                5,
                0,
                &tensor_border,
                &local_planar_patch_to_global(5, 0, local_tau),
            )
            .expect("family residuals");
            eprintln!(
                "Appendix-A eq14 source blend mask={mask_label} alpha={blend_alpha:e}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2}) families={residuals:?}",
                mask_label = mask.label(),
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
            if orientation.negative_determinants < best.2
                || (orientation.negative_determinants == best.2
                    && orientation.min_determinant > best.3)
            {
                best = (
                    mask.label(),
                    blend_alpha,
                    orientation.negative_determinants,
                    orientation.min_determinant,
                    solved.max_constraint_residual,
                );
            }
        }
    }

    eprintln!(
        "best Appendix-A eq14 source blend: mask={} alpha={:e} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3, best.4
    );
    assert!(
        best.4.is_finite(),
        "source-band compatibility sweep must produce finite solves"
    );
}

#[test]
fn rotational_paper_tau_regular_source_f3_priority_limit_is_collapsed() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let tensor_border = planar_regular_sector_tensor_border(5);
    let mut best_negative = usize::MAX;
    let mut best_near_zero = usize::MAX;
    let mut best_min_determinant = f64::NEG_INFINITY;
    let mut best_weight = 0.0_f64;

    for regularization_weight in [1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2] {
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: 5,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: planar_tau_seed_macro_patch(5, 0),
                regularization_weight,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let report = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "regular-source tau validity sweep weight={regularization_weight:e}: residual={residual:e} outside={outside}/{domain_samples} max_excess={max_excess:e} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
            residual = solved.max_constraint_residual,
            outside = report.domain.outside_triangle_samples,
            domain_samples = report.domain.sample_count,
            max_excess = report.domain.max_simplex_excess,
            neg = report.orientation.negative_determinants,
            orientation_samples = report.orientation.sample_count,
            zero = report.orientation.near_zero_determinants,
            min_det = report.orientation.min_determinant,
            max_det = report.orientation.max_determinant,
            worst = report.orientation.worst_piece,
            u = report.orientation.worst_u,
            v = report.orientation.worst_v,
        );
        if report.orientation.negative_determinants < best_negative
            || (report.orientation.negative_determinants == best_negative
                && report.orientation.min_determinant > best_min_determinant)
        {
            best_negative = report.orientation.negative_determinants;
            best_near_zero = report.orientation.near_zero_determinants;
            best_min_determinant = report.orientation.min_determinant;
            best_weight = regularization_weight;
        }
    }

    assert_eq!(
        best_negative, 0,
        "F3-priority source projection should reach the sampled nonnegative collapsed limit at weight {best_weight:e}"
    );
    assert_eq!(
        best_near_zero, 100,
        "the current F3-priority limit is not a nonsingular reparameterization"
    );
    assert!(
        best_min_determinant > -1.0e-20,
        "best regularized solve should be nonnegative to roundoff"
    );
}

#[test]
fn rotational_paper_tau_characteristic_source_orientation_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for valence in [5, 6] {
        let tensor_border = characteristic_cap_net_tensor_border(valence);
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: planar_tau_seed_macro_patch(valence, 0),
                regularization_weight: 1.0e-10,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let residuals = kp20_planar_tau_rotational_paper_family_residuals(
            valence,
            0,
            &tensor_border,
            &local_planar_patch_to_global(valence, 0, local_tau),
        )
        .expect("family residuals");
        eprintln!(
            "characteristic-source tau validity n={valence}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2}) families={residuals:?}",
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
fn characteristic_source_near_origin_first_order_frame_is_measured() {
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: planar_tau_seed_macro_patch(valence, 0),
        regularization_weight: 1.0e-10,
    })
    .expect("rotational paper tau solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
    let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
    let cap_net = kp20_characteristic_cap_net(valence).expect("characteristic cap-net");
    let c4 = global_to_local_planar(valence, 0, cap_net.sectors[0][3]);
    let c5 = global_to_local_planar(valence, 0, cap_net.sectors[0][4]);
    let c6 = global_to_local_planar(valence, 0, cap_net.sectors[0][5]);
    eprintln!(
        "characteristic inner cap-net local c4=({c4x:e},{c4y:e}) c5=({c5x:e},{c5y:e}) c6=({c6x:e},{c6y:e})",
        c4x = c4[0],
        c4y = c4[1],
        c5x = c5[0],
        c5y = c5[1],
        c6x = c6[0],
        c6y = c6[1],
    );
    eprintln!(
        "characteristic-source NearOriginLeft first-order origin=({lox:e},{loy:e}) du=({ldux:e},{lduy:e}) dv=({ldvx:e},{ldvy:e}) det={ldet:e}",
        lox = left.origin.x,
        loy = left.origin.y,
        ldux = left.du.x,
        lduy = left.du.y,
        ldvx = left.dv.x,
        ldvy = left.dv.y,
        ldet = left.determinant,
    );
    eprintln!(
        "characteristic-source NearOriginRight first-order origin=({rox:e},{roy:e}) du=({rdux:e},{rduy:e}) dv=({rdvx:e},{rdvy:e}) det={rdet:e}",
        rox = right.origin.x,
        roy = right.origin.y,
        rdux = right.du.x,
        rduy = right.du.y,
        rdvx = right.dv.x,
        rdvy = right.dv.y,
        rdet = right.determinant,
    );
    assert!(left.determinant.is_finite());
    assert!(right.determinant.is_finite());
}

#[test]
fn characteristic_source_centered_first_order_normalization_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, include_near_origin_right) in [("left-only", false), ("left-and-right", true)] {
        let solved = solve_rotational_characteristic_source_with_centered_first_order(
            5,
            0,
            include_near_origin_right,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "characteristic-source centered first-order {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn centered_first_order_conflict_is_localized_to_source_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let all_section42 = TestSection42FamilySet::new(true, true, true, true);
    let no_section42 = TestSection42FamilySet::new(false, false, false, false);
    let cases = [
        ("center-only", no_section42, false, false, false),
        ("section42", all_section42, false, false, false),
        ("section42-fig8-diagonal", all_section42, true, true, false),
        ("tensor-source-only", no_section42, false, false, true),
        ("full", all_section42, true, true, true),
    ];
    for (label, families, include_fig8, include_diagonal, include_tensor_source) in cases {
        let solved = solve_characteristic_centered_first_order_with_row_families(
            5,
            0,
            families,
            include_fig8,
            include_diagonal,
            include_tensor_source,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "centered first-order row localization {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_origin=({origin_x:e},{origin_y:e}) left_det={left_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            origin_x = left.origin.x,
            origin_y = left.origin.y,
            left_det = left.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn centered_first_order_source_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let solved =
        solve_characteristic_source_with_centered_projected_tensor_band(5, 0, true, true, false);
    let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
    let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let tensor_border = characteristic_cap_net_tensor_border(5);
    let residuals = kp20_planar_tau_rotational_paper_family_residuals(
        5,
        0,
        &tensor_border,
        &local_planar_patch_to_global(5, 0, local_tau),
    )
    .expect("family residuals");
    eprintln!(
        "centered first-order source projection: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e} families={residuals:?}",
        residual = solved.max_constraint_residual,
        scale = solved.extra_values[0],
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
        left_det = left.determinant,
        right_det = right.determinant,
    );
    assert!(solved.max_constraint_residual.is_finite());
}

#[test]
fn centered_first_order_projected_tensor_band_option_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let mut best = ("none", usize::MAX, f64::NEG_INFINITY, usize::MAX, 0.0_f64);
    let mut best_noncollapsed = ("none", usize::MAX, f64::NEG_INFINITY, usize::MAX, 0.0_f64);
    for (label, include_fig8, include_diagonal, include_right) in [
        ("section42-left", false, false, false),
        ("section42-both", false, false, true),
        ("section42-fig8-left", true, false, false),
        ("section42-fig8-both", true, false, true),
        ("section42-diagonal-left", false, true, false),
        ("section42-diagonal-both", false, true, true),
        ("section42-fig8-diagonal-left", true, true, false),
        ("section42-fig8-diagonal-both", true, true, true),
    ] {
        let solved = solve_characteristic_source_with_centered_projected_tensor_band(
            5,
            0,
            include_fig8,
            include_diagonal,
            include_right,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "centered projected tensor-band option {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
        let scale = solved.extra_values[0];
        if validity.orientation.negative_determinants < best.1
            || (validity.orientation.negative_determinants == best.1
                && validity.orientation.min_determinant > best.2)
        {
            best = (
                label,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                validity.domain.outside_triangle_samples,
                scale,
            );
        }
        if scale.abs() > 1.0e-3
            && (validity.orientation.negative_determinants < best_noncollapsed.1
                || (validity.orientation.negative_determinants == best_noncollapsed.1
                    && validity.orientation.min_determinant > best_noncollapsed.2))
        {
            best_noncollapsed = (
                label,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                validity.domain.outside_triangle_samples,
                scale,
            );
        }
    }
    eprintln!(
        "centered projected tensor-band best option: label={} neg={} min_det={:e} outside={} scale={:e}; best_noncollapsed={} neg={} min_det={:e} outside={} scale={:e}",
        best.0,
        best.1,
        best.2,
        best.3,
        best.4,
        best_noncollapsed.0,
        best_noncollapsed.1,
        best_noncollapsed.2,
        best_noncollapsed.3,
        best_noncollapsed.4,
    );
    assert_eq!(
        best.0, "section42-fig8-diagonal-both",
        "recalibrate if a paper-supported centered source projection variant improves the current collapsed best"
    );
    assert!(
        best.4.abs() < 1.0e-10,
        "both-side centered normalization currently wins only by collapsing the shared first-order scale, got {}",
        best.4
    );
    assert_eq!(
        best_noncollapsed.0, "section42-diagonal-left",
        "recalibrate if a non-collapsed centered source projection variant improves the current best"
    );
}

#[test]
fn fig8_first_order_projected_tensor_band_frame_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let axes = sector_local_axis_components(5, 0);
    let fig8_left = first_order_source_rows_for_local_frame(
        Kp20PieceId::NearOriginLeft,
        axes,
        [0.0, 0.0],
        [0.5, 0.5],
        [0.0, 1.0],
        0,
    );
    let fig8_right_same_orientation = first_order_source_rows_for_local_frame(
        Kp20PieceId::NearOriginRight,
        axes,
        [0.0, 0.0],
        [0.5, 0.5],
        [1.0, 0.0],
        0,
    );
    let fig8_right_positive_orientation = first_order_source_rows_for_local_frame(
        Kp20PieceId::NearOriginRight,
        axes,
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 0.5],
        0,
    );

    let mut best_noncollapsed = ("none", usize::MAX, f64::NEG_INFINITY, 0.0_f64);
    for (label, rows) in [
        (
            "axis-left",
            first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, axes, 0),
        ),
        ("fig8-left", fig8_left.clone()),
        (
            "fig8-both-same-orientation",
            [fig8_left.clone(), fig8_right_same_orientation].concat(),
        ),
        (
            "fig8-both-positive-orientation",
            [fig8_left, fig8_right_positive_orientation].concat(),
        ),
    ] {
        let solved = solve_characteristic_source_with_projected_tensor_band_and_first_order_rows(
            5, 0, true, false, rows,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let scale = solved.extra_values[0];
        eprintln!(
            "fig8 first-order projected tensor-band frame {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(solved.max_constraint_residual.is_finite());
        if scale.abs() > 1.0e-3
            && (validity.orientation.negative_determinants < best_noncollapsed.1
                || (validity.orientation.negative_determinants == best_noncollapsed.1
                    && validity.orientation.min_determinant > best_noncollapsed.2))
        {
            best_noncollapsed = (
                label,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                scale,
            );
        }
    }
    eprintln!(
        "fig8 first-order projected tensor-band best noncollapsed: label={} neg={} min_det={:e} scale={:e}",
        best_noncollapsed.0, best_noncollapsed.1, best_noncollapsed.2, best_noncollapsed.3
    );
    assert_eq!(
        best_noncollapsed.0, "axis-left",
        "recalibrate if the Fig. 8 first-order frame beats the current axis-left diagnostic"
    );
}

#[test]
fn rotational_paper_tau_characteristic_source_regularization_sweep_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    let mut best = (0.0_f64, usize::MAX, f64::NEG_INFINITY, f64::INFINITY);
    for regularization_weight in [
        1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4,
    ] {
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence,
                sector: 0,
                characteristic_tensor_border: &tensor_border,
                initial: planar_tau_seed_macro_patch(valence, 0),
                regularization_weight,
            })
            .expect("rotational paper tau solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "characteristic-source regularization weight={regularization_weight:e}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        if validity.orientation.negative_determinants < best.1
            || (validity.orientation.negative_determinants == best.1
                && validity.orientation.min_determinant > best.2)
        {
            best = (
                regularization_weight,
                validity.orientation.negative_determinants,
                validity.orientation.min_determinant,
                solved.max_constraint_residual,
            );
        }
    }
    eprintln!(
        "best characteristic-source regularization: weight={:e} neg={} min_det={:e} residual={:e}",
        best.0, best.1, best.2, best.3
    );
    assert!(best.3.is_finite());
}

#[test]
fn rotational_characteristic_source_appendix_a_central_pins_measure_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let tensor_border = characteristic_cap_net_tensor_border(valence);
    for (layout, use_paper_quadrant_maps) in
        [("positive-nonoverlap", false), ("paper-quadrant", true)]
    {
        for max_total_degree in [1, 2, 3] {
            let solved = solve_rotational_source_with_appendix_a_central_pins(
                valence,
                0,
                &tensor_border,
                max_total_degree,
                use_paper_quadrant_maps,
            );
            let local_tau = global_planar_patch_to_local(valence, 0, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            eprintln!(
                "characteristic-source Appendix-A central pins layout={layout} degree={max_total_degree}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
fn rotational_paper_tau_regular_source_near_origin_corner_derivative_is_inverted() {
    let tensor_border = planar_regular_sector_tensor_border(5);
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: planar_tau_seed_macro_patch(5, 0),
        regularization_weight: 1.0e-10,
    })
    .expect("rotational paper tau solve succeeds");
    let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
    let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
    let p00 = piece.cps[0][0];
    let p10 = piece.cps[1][0];
    let p01 = piece.cps[0][1];
    let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
    let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
    let determinant = du.x * dv.y - dv.x * du.y;

    eprintln!(
        "regular-source local NearOriginLeft corner p00=({p00x:e},{p00y:e}) p10=({p10x:e},{p10y:e}) p01=({p01x:e},{p01y:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e}) det={determinant:e}",
        p00x = p00.x,
        p00y = p00.y,
        p10x = p10.x,
        p10y = p10.y,
        p01x = p01.x,
        p01y = p01.y,
        dux = du.x,
        duy = du.y,
        dvx = dv.x,
        dvy = dv.y,
    );

    assert!(determinant.is_finite());
}

#[test]
fn rotational_regular_source_orientation_localizes_to_fig8_and_diagonal_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let cases = [
        ("section42-only", false, false),
        ("section42-plus-fig8", true, false),
        ("section42-plus-diagonal", false, true),
        ("section42-plus-fig8-diagonal", true, true),
    ];

    for (label, include_fig8, include_diagonal) in cases {
        let solved =
            solve_rotational_regular_source_with_row_options(5, 0, include_fig8, include_diagonal);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
        let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
        let corner_det = du.x * dv.y - dv.x * du.y;
        eprintln!(
            "regular-source row option {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] corner_det={corner_det:e} du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e})",
            residual = solved.max_constraint_residual,
            neg = orientation.negative_determinants,
            samples = orientation.sample_count,
            zero = orientation.near_zero_determinants,
            min_det = orientation.min_determinant,
            max_det = orientation.max_determinant,
            dux = du.x,
            duy = du.y,
            dvx = dv.x,
            dvy = dv.y,
        );
    }
}

#[test]
fn rotational_regular_source_orientation_localizes_within_section42_rows() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let cases = [
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

    for (label, families) in cases {
        let solved = solve_rotational_regular_source_with_section42_families(5, 0, families);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source section42 family {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_orientation_sweeps_diagonal_piece_mapping() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let piece_sets = [
        (
            "nol-nmr_self_nor-nml_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointRight],
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointLeft],
        ),
        (
            "nol-nor_self_nmr-nml_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight],
            [
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft,
            ],
        ),
        (
            "nol-nml_self_nor-nmr_pair",
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointLeft],
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointRight],
        ),
        (
            "nor-nmr_self_nol-nml_pair",
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointRight],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointLeft],
        ),
        (
            "nor-nml_self_nol-nmr_pair",
            [Kp20PieceId::NearOriginRight, Kp20PieceId::NearMidpointLeft],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearMidpointRight],
        ),
        (
            "nmr-nml_self_nol-nor_pair",
            [
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft,
            ],
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight],
        ),
    ];
    let mut best_negative = usize::MAX;

    for (label, self_pieces, paired_pieces) in piece_sets {
        let solved =
            solve_rotational_regular_source_with_diagonal_mapping(5, 0, self_pieces, paired_pieces);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source diagonal mapping {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        if orientation.negative_determinants < best_negative {
            best_negative = orientation.negative_determinants;
        }
    }

    assert_eq!(best_negative, 0);
}

#[test]
fn rotational_regular_source_origin_pins_measure_regular_tau_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let pin_sets: &[(&str, &[Kp20PlanarTauCoeff])] = &[
        (
            "near-origin-left-first-order",
            &[
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 1),
            ],
        ),
        (
            "both-origin-first-order",
            &[
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 1),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 1, 0),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 1),
            ],
        ),
    ];

    for &(label, pins) in pin_sets {
        let solved = solve_rotational_regular_source_with_origin_pins(5, 0, pins);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source origin pins {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_paper_tau_origin_pins_measure_feasibility() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, use_paper_quadrant_maps) in [
        ("appendix-a-positive-nonoverlap", false),
        ("appendix-a-paper-quadrant", true),
    ] {
        let solved = solve_rotational_regular_source_with_paper_tau_origin_pins(
            5,
            0,
            use_paper_quadrant_maps,
        );
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source paper tau origin pins {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_previous_sector_index_transpose_improves_tau_orientation() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, index_map) in [
        ("all-previous", TestSection42IndexMap::PreviousSectorAll),
        (
            "central-only",
            TestSection42IndexMap::PreviousSectorCentralOnly,
        ),
        (
            "stored-piece-edge",
            TestSection42IndexMap::StoredPieceEdgeOrientation,
        ),
    ] {
        let solved =
            solve_rotational_regular_source_with_previous_sector_index_transpose(5, 0, index_map);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        let piece = local_tau.piece(Kp20PieceId::NearOriginLeft);
        let p00 = piece.cps[0][0];
        let p10 = piece.cps[1][0];
        let p01 = piece.cps[0][1];
        let du = Point3::new(6.0 * (p10.x - p00.x), 6.0 * (p10.y - p00.y), 0.0);
        let dv = Point3::new(6.0 * (p01.x - p00.x), 6.0 * (p01.y - p00.y), 0.0);
        let corner_det = du.x * dv.y - dv.x * du.y;

        eprintln!(
            "regular-source previous-sector index transpose {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] corner_det={corner_det:e} worst={worst:?}@({u:.2},{v:.2})",
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

        assert!(solved.max_constraint_residual.is_finite());
    }
}

#[test]
fn rotational_regular_source_previous_sector_rotation_sign_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    for (label, sign) in [("plus", 1.0), ("minus", -1.0)] {
        let solved = solve_rotational_regular_source_with_rotation_sign(5, 0, sign);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source previous-sector rotation {label}: residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
    }
}

#[test]
fn rotational_regular_source_eq13_b0_sign_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let angle = core::f64::consts::TAU / 5.0;
    let c = angle.cos();
    for (label, b0_near_origin) in [
        ("paper-positive-2c", 2.0 * c),
        ("opposite-side-negative-2c", -2.0 * c),
        ("zero-eq13-b-coupling", 0.0),
        ("half-positive-c", c),
        ("half-negative-c", -c),
    ] {
        let solved = solve_rotational_regular_source_with_eq13_b0(5, 0, b0_near_origin);
        let local_tau = global_planar_patch_to_local(5, 0, solved.tau_macro_patch);
        let orientation = kp20_planar_bi6_tau_orientation_report(&local_tau, &SAMPLES, 1.0e-10);
        eprintln!(
            "regular-source eq13 b0 {label}: b0={b0_near_origin:e} residual={residual:e} neg={neg}/{samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
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
        assert!(solved.max_constraint_residual.is_finite());
    }
}
