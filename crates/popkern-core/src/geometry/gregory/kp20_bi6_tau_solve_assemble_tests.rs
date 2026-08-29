#![allow(unused_imports)]
// Paper-notation test module: loop indices, helper arities, and tuple shapes
// mirror KP15/KP20/PR08 subscripts and equation forms; keep the index form.
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::type_complexity
)]

pub(crate) use super::super::kp20_bi6_assemble::{
    Kp20Bi6MacroCapInput, Kp20Bi6MacroCapOutput,
    assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source,
};
pub(crate) use super::super::kp20_bi6_cap::{
    Kp20BoundaryHalf, Kp20SplitTensorBorder, split_kp20_tensor_border,
    split_kp20_tensor_border_with_d1_3,
};
pub(crate) use super::super::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, kp20_cap_net_figure_tensor_border_rows,
    kp20_tensor_border_from_regular_cnet_rows,
};
pub(crate) use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedF3SolveOutput,
    Kp20PlanarTauAugmentedScalarCoeff, Kp20PlanarTauAugmentedScalarLinearConstraint,
    kp20_planar_tau_constraint_space_report,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
};
pub(crate) use super::*;
pub(crate) use crate::geometry::Point3;
pub(crate) use crate::geometry::gregory::kp20_bi6_cc_ring::{
    Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
    kp20_pr08_characteristic_tensor_border_candidate_rows,
};
pub(crate) use crate::geometry::gregory::kp20_bi6_guide_build::{
    Kp20GuideBuildInput, build_kp20_guide_ring,
};
pub(crate) use crate::geometry::gregory::kp20_bi6_local_constraint::{
    Kp20Eq9DiagnosticTerms, Kp20Eq9TermSigns, Kp20Eq9TermSource, kp20_eq8_g1_residual,
    kp20_eq8_jacobian_dense, kp20_eq9_g2_residual, kp20_eq9_jacobian_dense,
    kp20_eq9_jacobian_dense_with_diagnostic_terms, measure_kp20_assembly_constraints,
};
pub(crate) use crate::geometry::gregory::kp20_bi6_math::{
    kp20_b_coeffs_near_m, kp20_b_coeffs_near_o,
};
pub(crate) use crate::geometry::gregory::{
    Bi6Reparams, GregoryPatch, Kp20Bi6Patch, Kp20PlanarTauCoeff, Kp20PlanarTauScalarCoeff,
    Kp20TauFamily, Kp20TauMap, Kp20Tri6GuideSector, RectangularGregoryBiNonic, TensorBorderBi3,
    apply_beta_reparameterization, elevate_bi6_grid_to_bi9, kp20_characteristic_cap_net,
    kp20_paper_quadrant_sector_piece_maps, kp20_planar_bi6_tau_f3_energy,
    kp20_planar_bi6_tau_orientation_report, kp20_planar_bi6_tau_validity_report,
    kp20_positive_nonoverlap_c2_sector_piece_maps,
    kp20_pr08_stam_inner_q1_characteristic_tensor_border,
    kp20_pr08_stam_q1_characteristic_tensor_border,
    sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch,
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

fn eq13_scalar_row_label(row_index: usize) -> &'static str {
    match row_index {
        0 => "underlined-p[5] x",
        1 => "underlined-p[5] y",
        2 => "arrow-underlined-p[5] x",
        3 => "arrow-underlined-p[5] y",
        4 => "underlined-p[6] x",
        5 => "underlined-p[6] y",
        6 => "arrow-underlined-p[6] x",
        7 => "arrow-underlined-p[6] y",
        8 => "p[4] x",
        9 => "p[4] y",
        10 => "arrow-p[4] x",
        11 => "arrow-p[4] y",
        _ => "unknown eq13 row",
    }
}

fn coefficient_is_pinned(pinned: &[Kp20PlanarTauCoeff], coeff: Kp20PlanarTauCoeff) -> bool {
    pinned.contains(&coeff)
}

fn scalar_row_coverage(
    row: &Kp20PlanarTauScalarLinearConstraint,
    pinned: &[Kp20PlanarTauCoeff],
) -> (usize, usize, f64) {
    let mut covered = 0_usize;
    let mut missing_weight = 0.0_f64;
    for &(scalar_coeff, weight) in &row.terms {
        if coefficient_is_pinned(pinned, scalar_coeff.coeff) {
            covered += 1;
        } else {
            missing_weight += weight.abs();
        }
    }
    (covered, row.terms.len(), missing_weight)
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
    let w3 = (0.5 * b0_near_origin) / 30.0;

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

fn fig6c_active_set_source_projection(
    valence: usize,
    sector: usize,
    samples: &[f64],
) -> Option<TestSquareActiveSetProjectionReport> {
    fig6c_active_set_source_projection_with_map(valence, sector, samples).map(|(report, _)| report)
}

fn fig6c_active_set_source_projection_with_map(
    valence: usize,
    sector: usize,
    samples: &[f64],
) -> Option<(
    TestSquareActiveSetProjectionReport,
    TestSection42StripSymbolMap,
)> {
    let angle = core::f64::consts::TAU / valence as f64;
    let piece_maps = kp20_paper_quadrant_sector_piece_maps();
    let refined_local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        piece_maps,
    )
    .expect("refined paper tau samples with identity guide");
    let refined_global_tau = local_planar_patch_to_global(valence, sector, refined_local_tau);
    let strips = appendix_a_candidate_section42_strips();
    let refined_exact_pairs = exact_strip_boundary_pairs(&refined_global_tau, &strips, angle);
    let (report, _, _) = best_fig6c_source_row_basis_projected_graph(
        &refined_global_tau,
        &strips,
        &refined_exact_pairs,
        angle,
    )?;
    let (boundary_projected, _, _) =
        fig6c_boundary_row_basis_projected_source_patch(&refined_global_tau, report.map);
    let active_set = fig6c_square_domain_active_set_projection(
        &boundary_projected,
        report.map,
        valence,
        sector,
        samples,
    );
    Some((active_set, report.map))
}

#[derive(Debug, Clone)]
struct TestComposedCapProjectionReport {
    patch: Kp20MacroPatch,
    alpha: f64,
    valid_alpha: f64,
    valid_cap_residual: f64,
    valid_tau_full: f64,
    max_delta: f64,
    linear_residual: f64,
}

fn fig6c_composed_cap_residual_projection(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
) -> TestComposedCapProjectionReport {
    fig6c_composed_cap_residual_projection_with_weights(
        patch, map, guide, valence, sector, 100.0, 1.0, false,
    )
}

fn fig6c_composed_cap_residual_projection_with_weights(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
    tau_row_weight: f64,
    cap_row_weight: f64,
    include_validity_guards: bool,
) -> TestComposedCapProjectionReport {
    const FINITE_DIFFERENCE_STEP: f64 = 1.0e-6;
    const VALIDITY_GUARD_WEIGHT: f64 = 100.0;
    let variables = all_macro_patch_axis_variables();
    let tau_rows = fig6c_full_paper_direction_rows(map);
    let cap_residuals = fig6c_composed_cap_residual_vector(patch, guide, valence, sector);
    const VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let square_guards = if include_validity_guards {
        square_domain_guard_constraints(patch, valence, sector, &VALIDITY_SAMPLES, 1.0e-3)
    } else {
        Vec::new()
    };
    let orientation_guards = if include_validity_guards {
        orientation_guard_constraints(patch, valence, sector, &VALIDITY_SAMPLES, 1.0e-3)
    } else {
        Vec::new()
    };
    let row_count =
        tau_rows.len() + cap_residuals.len() + square_guards.len() + orientation_guards.len();
    let mut design = vec![vec![0.0_f64; variables.len()]; row_count];
    let mut rhs = vec![0.0_f64; row_count];

    for (row_index, row) in tau_rows.iter().enumerate() {
        rhs[row_index] = -tau_row_weight * scalar_residual(patch, row);
        for &(scalar_coeff, weight) in &row.terms {
            let variable_index =
                macro_patch_axis_variable_index(scalar_coeff.coeff, scalar_coeff.axis);
            design[row_index][variable_index] += tau_row_weight * weight;
        }
    }
    for (cap_index, &residual) in cap_residuals.iter().enumerate() {
        rhs[tau_rows.len() + cap_index] = -cap_row_weight * residual;
    }
    let square_start = tau_rows.len() + cap_residuals.len();
    let local_axes = sector_local_axis_components(valence, sector);
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    for (guard_index, guard) in square_guards.iter().copied().enumerate() {
        let row_index = square_start + guard_index;
        let current_value = sample_local_patch_axis(
            &local_patch,
            guard.piece_id,
            guard.u,
            guard.v,
            guard.local_axis,
        );
        rhs[row_index] = VALIDITY_GUARD_WEIGHT * (guard.target - current_value);
        let bu = test_bernstein_bi6_all(guard.u);
        let bv = test_bernstein_bi6_all(guard.v);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let coeff = Kp20PlanarTauCoeff::new(guard.piece_id, i, j);
                let bernstein_weight = bu_i * bv_j;
                for global_axis in 0..=1 {
                    let variable_index = macro_patch_axis_variable_index(coeff, global_axis);
                    design[row_index][variable_index] += VALIDITY_GUARD_WEIGHT
                        * bernstein_weight
                        * local_axes[guard.local_axis][global_axis];
                }
            }
        }
    }
    let orientation_start = square_start + square_guards.len();
    for (guard_index, guard) in orientation_guards.iter().copied().enumerate() {
        let row_index = orientation_start + guard_index;
        let local_piece = local_patch.piece(guard.piece_id);
        let (du, dv) = sample_local_patch_derivatives(local_piece, guard.u, guard.v);
        let det = du.x * dv.y - du.y * dv.x;
        rhs[row_index] = VALIDITY_GUARD_WEIGHT * (guard.target_det - det);
        add_orientation_active_design_row(
            &mut design[row_index],
            guard,
            du,
            dv,
            &local_axes,
            VALIDITY_GUARD_WEIGHT,
        );
    }

    for (variable_index, &(coeff, axis)) in variables.iter().enumerate() {
        let mut perturbed = *patch;
        let cp = &mut perturbed.piece_mut(coeff.piece).cps[coeff.i][coeff.j];
        if axis == 0 {
            cp.x += FINITE_DIFFERENCE_STEP;
        } else {
            cp.y += FINITE_DIFFERENCE_STEP;
        }
        let perturbed_residuals =
            fig6c_composed_cap_residual_vector(&perturbed, guide, valence, sector);
        for (cap_index, (&perturbed_residual, &base_residual)) in perturbed_residuals
            .iter()
            .zip(cap_residuals.iter())
            .enumerate()
        {
            design[tau_rows.len() + cap_index][variable_index] +=
                cap_row_weight * (perturbed_residual - base_residual) / FINITE_DIFFERENCE_STEP;
        }
    }

    let solution = solve_test_regularized_least_squares(&design, &rhs, 1.0e-10);
    let linear_residual = max_dense_row_residual(&design, &solution, &rhs);
    let mut full_step = *patch;
    let mut max_delta = 0.0_f64;
    for (&(coeff, axis), delta) in variables.iter().zip(solution.iter().copied()) {
        let cp = &mut full_step.piece_mut(coeff.piece).cps[coeff.i][coeff.j];
        if axis == 0 {
            cp.x += delta;
        } else {
            cp.y += delta;
        }
    }
    for chunk in solution.chunks(2) {
        let dx = chunk[0];
        let dy = chunk.get(1).copied().unwrap_or(0.0);
        max_delta = max_delta.max((dx * dx + dy * dy).sqrt());
    }
    let mut best_alpha = 0.0_f64;
    let mut best_patch = *patch;
    let mut best_cap_residual = max_abs(&cap_residuals);
    let mut valid_alpha = 0.0_f64;
    let mut valid_cap_residual = best_cap_residual;
    let mut valid_tau_full = scalar_row_set_residual(patch, &tau_rows);
    for alpha in [
        1.0_f64,
        0.5,
        0.25,
        0.125,
        0.0625,
        0.03125,
        0.015625,
        0.0078125,
        0.00390625,
        0.001953125,
        0.0009765625,
        0.00048828125,
        0.000244140625,
        0.0001220703125,
        0.00006103515625,
        0.000030517578125,
        0.0000152587890625,
    ] {
        let candidate = blend_macro_patch_points(patch, &full_step, alpha);
        let cap_residual = max_abs(&fig6c_composed_cap_residual_vector(
            &candidate, guide, valence, sector,
        ));
        if cap_residual < best_cap_residual {
            best_alpha = alpha;
            best_patch = candidate;
            best_cap_residual = cap_residual;
        }
        let candidate_tau_full = scalar_row_set_residual(&candidate, &tau_rows);
        let candidate_local = global_planar_patch_to_local(valence, sector, candidate);
        let candidate_square =
            test_square_domain_report(&candidate_local, &VALIDITY_SAMPLES, 1.0e-10);
        let candidate_orientation =
            test_orientation_detail_report(&candidate_local, &VALIDITY_SAMPLES, 1.0e-10);
        if candidate_tau_full <= 1.0e-5
            && candidate_square.outside == 0
            && candidate_orientation.negative == 0
            && cap_residual < valid_cap_residual
        {
            valid_alpha = alpha;
            valid_cap_residual = cap_residual;
            valid_tau_full = candidate_tau_full;
        }
    }
    TestComposedCapProjectionReport {
        patch: best_patch,
        alpha: best_alpha,
        valid_alpha,
        valid_cap_residual,
        valid_tau_full,
        max_delta: best_alpha * max_delta,
        linear_residual,
    }
}

fn fig6c_composed_cap_residual_vector(
    patch: &Kp20MacroPatch,
    guide: &Kp20Tri6GuideSector,
    valence: usize,
    sector: usize,
) -> Vec<f64> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    let sampled = sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(guide, &local_patch);
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let near_o = kp20_b_coeffs_near_o(c);
    let near_m = kp20_b_coeffs_near_m(c);
    let mut residuals = Vec::with_capacity((7 + 8) * 2 * 3);
    push_vec3_residuals(
        &mut residuals,
        &kp20_eq8_g1_residual(
            &sampled.piece(Kp20PieceId::NearOriginLeft).cps,
            &sampled.piece(Kp20PieceId::NearOriginRight).cps,
            near_o,
        ),
    );
    push_vec3_residuals(
        &mut residuals,
        &kp20_eq8_g1_residual(
            &sampled.piece(Kp20PieceId::NearMidpointLeft).cps,
            &sampled.piece(Kp20PieceId::NearMidpointRight).cps,
            near_m,
        ),
    );
    push_vec3_residuals(
        &mut residuals,
        &kp20_eq9_g2_residual(
            &sampled.piece(Kp20PieceId::NearOriginLeft).cps,
            &sampled.piece(Kp20PieceId::NearOriginRight).cps,
            near_o,
        ),
    );
    push_vec3_residuals(
        &mut residuals,
        &kp20_eq9_g2_residual(
            &sampled.piece(Kp20PieceId::NearMidpointLeft).cps,
            &sampled.piece(Kp20PieceId::NearMidpointRight).cps,
            near_m,
        ),
    );
    residuals
}

fn push_vec3_residuals<const N: usize>(out: &mut Vec<f64>, residuals: &[crate::geometry::Vec3; N]) {
    for residual in residuals {
        out.push(residual.x);
        out.push(residual.y);
        out.push(residual.z);
    }
}

fn max_abs(values: &[f64]) -> f64 {
    values
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
}

fn max_slice(values: &[f64]) -> f64 {
    values.iter().copied().fold(0.0_f64, f64::max)
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
struct TestRowSpaceAudit {
    row_count: usize,
    max_design_coeff: f64,
    max_rhs: f64,
    max_normal_target: f64,
    max_solution: f64,
    max_linear_residual: f64,
}

fn test_scalar_row_space_audit(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
    row_weight: f64,
) -> TestRowSpaceAudit {
    let variables = all_macro_patch_axis_variables();
    let mut design = vec![vec![0.0_f64; variables.len()]; rows.len()];
    let mut rhs = vec![0.0_f64; rows.len()];
    let mut max_design_coeff = 0.0_f64;
    for (row_index, row) in rows.iter().enumerate() {
        rhs[row_index] = -row_weight * scalar_residual(patch, row);
        for &(scalar_coeff, weight) in &row.terms {
            let variable_index =
                macro_patch_axis_variable_index(scalar_coeff.coeff, scalar_coeff.axis);
            let weighted = row_weight * weight;
            design[row_index][variable_index] += weighted;
            max_design_coeff = max_design_coeff.max(weighted.abs());
        }
    }
    let mut normal_target = vec![0.0_f64; variables.len()];
    for (row, &rhs_value) in design.iter().zip(rhs.iter()) {
        for (index, &coefficient) in row.iter().enumerate() {
            normal_target[index] += coefficient * rhs_value;
        }
    }
    let solution = solve_test_regularized_least_squares(&design, &rhs, 1.0e-10);
    let normalized_solution = solve_test_row_normalized_least_squares(&design, &rhs, 1.0e-10);
    TestRowSpaceAudit {
        row_count: rows.len(),
        max_design_coeff,
        max_rhs: max_abs(&rhs),
        max_normal_target: max_abs(&normal_target),
        max_solution: max_abs(&solution).max(max_abs(&normalized_solution)),
        max_linear_residual: max_dense_row_residual(&design, &solution, &rhs)
            .min(max_dense_row_residual(&design, &normalized_solution, &rhs)),
    }
}

fn add_cap_pair_v0_c0_rows(rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>) {
    add_cap_pair_v0_c0_rows_with_chart(rows, TestCapPairBoundaryChart::SameDirection);
}

fn add_cap_pair_v0_c0_rows_with_chart(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    chart: TestCapPairBoundaryChart,
) {
    for (left, right) in [
        (Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight),
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        ),
    ] {
        for i in 0..=6 {
            for axis in 0..=1 {
                rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                    (
                        Kp20PlanarTauScalarCoeff::new(Kp20PlanarTauCoeff::new(left, i, 0), axis),
                        1.0,
                    ),
                    (
                        Kp20PlanarTauScalarCoeff::new(chart.current_side_coeff(right, i, 0), axis),
                        -1.0,
                    ),
                ]));
            }
        }
    }
}

fn add_cap_pair_eq8_eq9_rows(rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>, valence: usize) {
    add_cap_pair_residual_rows(rows, valence, TestCapResidualFamily::Eq8);
    add_cap_pair_residual_rows(rows, valence, TestCapResidualFamily::Eq9);
}

#[derive(Debug, Clone, Copy)]
enum TestCapResidualFamily {
    Eq8,
    Eq9,
}

#[derive(Debug, Clone, Copy)]
struct TestCapResidualFamilySet {
    eq8: bool,
    eq9: bool,
}

impl TestCapResidualFamilySet {
    const fn new(eq8: bool, eq9: bool) -> Self {
        Self { eq8, eq9 }
    }
}

#[derive(Debug, Clone, Copy)]
enum TestCapEq8Join {
    NearOrigin,
    NearMidpoint,
}

#[derive(Debug, Clone, Copy)]
enum TestCapPairBoundaryChart {
    SameDirection,
    ReverseCurrentU,
}

impl TestCapPairBoundaryChart {
    fn current_side_coeff(self, piece_id: Kp20PieceId, i: usize, j: usize) -> Kp20PlanarTauCoeff {
        match self {
            Self::SameDirection => Kp20PlanarTauCoeff::new(piece_id, i, j),
            Self::ReverseCurrentU => Kp20PlanarTauCoeff::new(piece_id, 6 - i, j),
        }
    }
}

fn add_cap_pair_residual_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    family: TestCapResidualFamily,
) {
    add_cap_pair_residual_rows_with_chart(
        rows,
        valence,
        family,
        TestCapPairBoundaryChart::SameDirection,
    );
}

fn add_cap_pair_residual_rows_with_chart(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    family: TestCapResidualFamily,
    chart: TestCapPairBoundaryChart,
) {
    let c = (core::f64::consts::TAU / valence as f64).cos();
    for (left, right, b_coeffs) in [
        (
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            kp20_b_coeffs_near_o(c),
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
        ),
    ] {
        match family {
            TestCapResidualFamily::Eq8 => {
                let eq8 = kp20_eq8_jacobian_dense(b_coeffs);
                for dense_row in eq8 {
                    add_cap_pair_dense_residual_rows_with_chart(
                        rows, left, right, &dense_row, chart,
                    );
                }
            }
            TestCapResidualFamily::Eq9 => {
                let eq9 = kp20_eq9_jacobian_dense(b_coeffs);
                for dense_row in eq9 {
                    add_cap_pair_dense_residual_rows_with_chart(
                        rows, left, right, &dense_row, chart,
                    );
                }
            }
        }
    }
}

fn add_cap_pair_eq8_index_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    join: TestCapEq8Join,
    row_index: usize,
) {
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let (left, right, b_coeffs) = match join {
        TestCapEq8Join::NearOrigin => (
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            kp20_b_coeffs_near_o(c),
        ),
        TestCapEq8Join::NearMidpoint => (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
        ),
    };
    let eq8 = kp20_eq8_jacobian_dense(b_coeffs);
    add_cap_pair_dense_residual_rows(rows, left, right, &eq8[row_index]);
}

fn add_cap_pair_eq8_join_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    join: TestCapEq8Join,
) {
    for row_index in 0..7 {
        add_cap_pair_eq8_index_rows(rows, valence, join, row_index);
    }
}

fn log_near_origin_frame_summary(
    label: &str,
    valence: usize,
    sector: usize,
    patch: &Kp20MacroPatch,
) {
    let local = global_planar_patch_to_local(valence, sector, *patch);
    for piece_id in [
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let frame = near_origin_first_order_frame(&local, piece_id);
        eprintln!(
            "Appendix-A Fig.6c {label} {piece_id:?} first-order origin=({ox:e},{oy:e}) du=({dux:e},{duy:e}) dv=({dvx:e},{dvy:e}) det={det:e}",
            ox = frame.origin.x,
            oy = frame.origin.y,
            dux = frame.du.x,
            duy = frame.du.y,
            dvx = frame.dv.x,
            dvy = frame.dv.y,
            det = frame.determinant,
        );
    }
}

fn add_near_origin_right_transverse_matches_left_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
) {
    let left_origin = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 0);
    let left_transverse = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 1);
    let right_origin = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 0);
    let right_transverse = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, 1);
    for axis in 0..=1 {
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauScalarCoeff::new(right_transverse, axis), 1.0),
            (Kp20PlanarTauScalarCoeff::new(right_origin, axis), -1.0),
            (Kp20PlanarTauScalarCoeff::new(left_transverse, axis), -1.0),
            (Kp20PlanarTauScalarCoeff::new(left_origin, axis), 1.0),
        ]));
    }
}

fn add_fixed_near_origin_first_order_frame_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    source: &Kp20MacroPatch,
) {
    let left = source.piece(Kp20PieceId::NearOriginLeft);
    let targets = [
        ((0_usize, 0_usize), left.cps[0][0]),
        ((1_usize, 0_usize), left.cps[1][0]),
        ((0_usize, 1_usize), left.cps[0][1]),
    ];
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        for ((i, j), target) in targets {
            let coeff = Kp20PlanarTauCoeff::new(piece_id, i, j);
            for (axis, value) in [(0_usize, target.x), (1, target.y)] {
                rows.push(Kp20PlanarTauScalarLinearConstraint {
                    terms: vec![(Kp20PlanarTauScalarCoeff::new(coeff, axis), 1.0)],
                    target: value,
                });
            }
        }
    }
}

fn add_paper_quadrant_internal_c0_rows(rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>) {
    for edge in TestPaperQuadrantInternalEdge::ALL {
        add_paper_quadrant_internal_c0_edge_rows(rows, edge);
    }
}

#[derive(Debug, Clone, Copy)]
enum TestPaperQuadrantInternalEdge {
    LowerLeftToLowerRight,
    LowerLeftToUpperLeft,
    LowerRightToUpperRight,
    UpperLeftToUpperRight,
}

impl TestPaperQuadrantInternalEdge {
    const ALL: [Self; 4] = [
        Self::LowerLeftToLowerRight,
        Self::LowerLeftToUpperLeft,
        Self::LowerRightToUpperRight,
        Self::UpperLeftToUpperRight,
    ];
}

fn add_paper_quadrant_internal_c0_edge_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    edge: TestPaperQuadrantInternalEdge,
) {
    for k in 0_usize..=6 {
        let (a, b) = match edge {
            TestPaperQuadrantInternalEdge::LowerLeftToLowerRight => (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 6, k),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, 0, k),
            ),
            TestPaperQuadrantInternalEdge::LowerLeftToUpperLeft => (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, k, 6),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, k, 0),
            ),
            TestPaperQuadrantInternalEdge::LowerRightToUpperRight => (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginRight, k, 6),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, k, 0),
            ),
            TestPaperQuadrantInternalEdge::UpperLeftToUpperRight => (
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointLeft, 6, k),
                Kp20PlanarTauCoeff::new(Kp20PieceId::NearMidpointRight, 0, k),
            ),
        };
        add_axis_point_equality_rows(rows, a, b);
    }
}

fn add_axis_point_equality_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    a: Kp20PlanarTauCoeff,
    b: Kp20PlanarTauCoeff,
) {
    for axis in 0..=1 {
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauScalarCoeff::new(a, axis), 1.0),
            (Kp20PlanarTauScalarCoeff::new(b, axis), -1.0),
        ]));
    }
}

fn test_paper_quadrant_affine_macro_patch() -> Kp20MacroPatch {
    let mut patch = Kp20MacroPatch::zero();
    for (piece_id, maps) in kp20_paper_quadrant_sector_piece_maps() {
        let map = maps[0].1;
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                piece.cps[i][j] = Point3::new(
                    map.origin[0] + u * map.du[0] + v * map.dv[0],
                    map.origin[1] + u * map.du[1] + v * map.dv[1],
                    0.0,
                );
            }
        }
    }
    patch
}

fn add_cap_pair_eq8_rows_with_orientation(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    origin_b_sign: f64,
    midpoint_b_sign: f64,
    origin_swap: bool,
    midpoint_swap: bool,
) {
    let c = (core::f64::consts::TAU / valence as f64).cos();
    for (left, right, b_coeffs, swap) in [
        (
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            kp20_b_coeffs_near_o(c).map(|value| origin_b_sign * value),
            origin_swap,
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c).map(|value| midpoint_b_sign * value),
            midpoint_swap,
        ),
    ] {
        let (previous, current) = if swap { (right, left) } else { (left, right) };
        let eq8 = kp20_eq8_jacobian_dense(b_coeffs);
        for dense_row in eq8 {
            add_cap_pair_dense_residual_rows(rows, previous, current, &dense_row);
        }
    }
}

fn add_cap_pair_eq9_rows_with_terms(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    valence: usize,
    origin_terms: Kp20Eq9DiagnosticTerms,
    midpoint_terms: Kp20Eq9DiagnosticTerms,
) {
    let c = (core::f64::consts::TAU / valence as f64).cos();
    for (left, right, b_coeffs, terms) in [
        (
            Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearOriginRight,
            kp20_b_coeffs_near_o(c),
            origin_terms,
        ),
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
            kp20_b_coeffs_near_m(c),
            midpoint_terms,
        ),
    ] {
        let eq9 = kp20_eq9_jacobian_dense_with_diagnostic_terms(b_coeffs, terms);
        for dense_row in eq9 {
            add_cap_pair_dense_residual_rows(rows, left, right, &dense_row);
        }
    }
}

fn cap_pair_eq9_terms_residual(
    patch: &Kp20MacroPatch,
    valence: usize,
    origin_terms: Kp20Eq9DiagnosticTerms,
    midpoint_terms: Option<Kp20Eq9DiagnosticTerms>,
) -> f64 {
    let mut rows = Vec::new();
    add_cap_pair_eq9_rows_with_terms(
        &mut rows,
        valence,
        origin_terms,
        midpoint_terms.unwrap_or_else(Kp20Eq9DiagnosticTerms::paper),
    );
    scalar_row_set_residual(patch, &rows)
}

fn test_eq9_term_source_variants() -> Vec<Kp20Eq9DiagnosticTerms> {
    let mut variants = Vec::with_capacity(64);
    for signs in test_eq9_term_sign_variants() {
        variants.extend(test_eq9_source_variants_for_signs(signs));
    }
    variants
}

fn test_eq9_term_sign_variants() -> [Kp20Eq9TermSigns; 8] {
    let paper = Kp20Eq9TermSigns::paper();
    let mut variants = [paper; 8];
    let signs = [-1.0, 1.0];
    let mut index = 0;
    for two_b_duv in signs {
        for b_sq_duu in signs {
            for b_bprime_du in signs {
                variants[index] = Kp20Eq9TermSigns {
                    two_b_duv,
                    b_sq_duu,
                    b_bprime_du,
                    ..paper
                };
                index += 1;
            }
        }
    }
    variants
}

fn test_eq9_source_variants_for_signs(signs: Kp20Eq9TermSigns) -> [Kp20Eq9DiagnosticTerms; 8] {
    let paper = Kp20Eq9DiagnosticTerms {
        signs,
        ..Kp20Eq9DiagnosticTerms::paper()
    };
    let mut variants = [paper; 8];
    let sources = [Kp20Eq9TermSource::Previous, Kp20Eq9TermSource::Current];
    let mut index = 0;
    for duv_source in sources {
        for duu_source in sources {
            for du_source in sources {
                variants[index] = Kp20Eq9DiagnosticTerms {
                    signs,
                    duv_source,
                    duu_source,
                    du_source,
                };
                index += 1;
            }
        }
    }
    variants
}

fn add_cap_pair_dense_residual_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    left: Kp20PieceId,
    right: Kp20PieceId,
    dense_row: &[f64; 98],
) {
    add_cap_pair_dense_residual_rows_with_chart(
        rows,
        left,
        right,
        dense_row,
        TestCapPairBoundaryChart::SameDirection,
    );
}

fn add_cap_pair_dense_residual_rows_with_chart(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    left: Kp20PieceId,
    right: Kp20PieceId,
    dense_row: &[f64; 98],
    chart: TestCapPairBoundaryChart,
) {
    for axis in 0..=1 {
        let mut terms = Vec::new();
        for side in 0..=1 {
            let piece_id = if side == 0 { left } else { right };
            for i in 0..=6 {
                for j in 0..=6 {
                    let weight = dense_row[side * 49 + i * 7 + j];
                    if weight != 0.0 {
                        let coeff = if side == 0 {
                            Kp20PlanarTauCoeff::new(piece_id, i, j)
                        } else {
                            chart.current_side_coeff(piece_id, i, j)
                        };
                        terms.push((Kp20PlanarTauScalarCoeff::new(coeff, axis), weight));
                    }
                }
            }
        }
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(terms));
    }
}

#[derive(Debug, Clone, Copy)]
struct TestFoldSummary {
    fold_count: usize,
    min_jacobian_normal_cos: f64,
    piece_fold_counts: [usize; 4],
}

#[derive(Debug, Clone, Copy)]
struct TestCapPairC0Summary {
    near_origin: f64,
    near_midpoint: f64,
}

#[derive(Debug, Clone, Copy)]
struct TestReversedCurrentCapResidualSummary {
    c0: f64,
    eq8: f64,
    eq9: f64,
}

fn test_reversed_current_cap_residual_summary(
    cap: &Kp20Bi6MacroCapOutput,
    valence: usize,
) -> TestReversedCurrentCapResidualSummary {
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let mut summary = TestReversedCurrentCapResidualSummary {
        c0: 0.0,
        eq8: 0.0,
        eq9: 0.0,
    };
    for macro_patch in &cap.macro_patches {
        for (left, right, b_coeffs) in [
            (
                Kp20PieceId::NearOriginLeft,
                Kp20PieceId::NearOriginRight,
                kp20_b_coeffs_near_o(c),
            ),
            (
                Kp20PieceId::NearMidpointLeft,
                Kp20PieceId::NearMidpointRight,
                kp20_b_coeffs_near_m(c),
            ),
        ] {
            let left_grid = &macro_patch.piece(left).cps;
            let right_grid = reverse_bi6_grid_u(&macro_patch.piece(right).cps);
            summary.c0 = summary.c0.max(test_reversed_current_v0_pair_c0_distance(
                left_grid,
                &right_grid,
            ));
            summary.eq8 = summary.eq8.max(test_max_vec3(&kp20_eq8_g1_residual(
                left_grid,
                &right_grid,
                b_coeffs,
            )));
            summary.eq9 = summary.eq9.max(test_max_vec3(&kp20_eq9_g2_residual(
                left_grid,
                &right_grid,
                b_coeffs,
            )));
        }
    }
    summary
}

fn reverse_bi6_grid_u(grid: &[[Point3; 7]; 7]) -> [[Point3; 7]; 7] {
    let mut reversed = *grid;
    for i in 0..=6 {
        reversed[i] = grid[6 - i];
    }
    reversed
}

fn test_reversed_current_v0_pair_c0_distance(
    left: &[[Point3; 7]; 7],
    right_reversed_u: &[[Point3; 7]; 7],
) -> f64 {
    (0..=6)
        .map(|i| distance(left[i][0], right_reversed_u[i][0]))
        .fold(0.0_f64, f64::max)
}

fn test_max_vec3<const N: usize>(values: &[crate::geometry::Vec3; N]) -> f64 {
    values
        .iter()
        .map(|value| value.x.abs().max(value.y.abs()).max(value.z.abs()))
        .fold(0.0_f64, f64::max)
}

fn test_cap_v0_pair_c0_summary(cap: &Kp20Bi6MacroCapOutput) -> TestCapPairC0Summary {
    let mut summary = TestCapPairC0Summary {
        near_origin: 0.0,
        near_midpoint: 0.0,
    };
    for macro_patch in &cap.macro_patches {
        summary.near_origin = summary.near_origin.max(test_v0_pair_c0_distance(
            macro_patch.piece(Kp20PieceId::NearOriginLeft),
            macro_patch.piece(Kp20PieceId::NearOriginRight),
        ));
        summary.near_midpoint = summary.near_midpoint.max(test_v0_pair_c0_distance(
            macro_patch.piece(Kp20PieceId::NearMidpointLeft),
            macro_patch.piece(Kp20PieceId::NearMidpointRight),
        ));
    }
    summary
}

fn test_v0_pair_c0_distance(left: &Kp20Bi6Patch, right: &Kp20Bi6Patch) -> f64 {
    (0..=6)
        .map(|i| distance(left.cps[i][0], right.cps[i][0]))
        .fold(0.0_f64, f64::max)
}

fn test_macro_cap_from_sampled(
    n: usize,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Kp20Bi6MacroCapOutput {
    Kp20Bi6MacroCapOutput {
        n,
        patches_bi6: macro_patches
            .iter()
            .flat_map(|macro_patch| {
                Kp20PieceId::all()
                    .into_iter()
                    .map(|piece_id| macro_patch.piece(piece_id).cps)
            })
            .collect(),
        macro_patches,
        strip_projections: Vec::new(),
    }
}

fn test_summarize_cap_folds(cap: &Kp20Bi6MacroCapOutput) -> TestFoldSummary {
    const SELF_INTERSECTION_TOLERANCE: f64 = 1.0e-6;
    let mut fold_count = 0;
    let mut min_jacobian_normal_cos = f64::INFINITY;
    let mut piece_fold_counts = [0_usize; 4];
    for (patch_index, grid) in cap.patches_bi6.iter().enumerate() {
        let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
            grid: elevate_bi6_grid_to_bi9(*grid),
        }));
        let fold_report = patch.self_intersection_report(SELF_INTERSECTION_TOLERANCE);
        if fold_report.intersects() {
            fold_count += 1;
            piece_fold_counts[patch_index % 4] += 1;
            min_jacobian_normal_cos =
                min_jacobian_normal_cos.min(fold_report.min_jacobian_normal_cos);
        }
    }
    TestFoldSummary {
        fold_count,
        min_jacobian_normal_cos,
        piece_fold_counts,
    }
}

fn best_fig6c_internal_join_first_graph(
    patch: &Kp20MacroPatch,
    exact_pairs: &[(TestSection42StripMap, TestSection42StripMap)],
    previous_sector_rotation_radians: f64,
) -> Option<(
    f64,
    TestSection42StripMap,
    TestSection42StripMap,
    TestSection42StripMap,
    TestSection42StripMap,
)> {
    let two_piece_pairs = exact_pairs
        .iter()
        .copied()
        .filter(|&(left, right)| {
            matches!(
                fig6c_exact_pair_kind(left, right),
                Some(TestFig6cExactPairKind::TwoPieceJoin)
            )
        })
        .collect::<Vec<_>>();
    let mut best = None;
    for &(p, underlined_p) in &two_piece_pairs {
        for &(arrow_p, arrow_underlined_p) in &two_piece_pairs {
            if !fig6c_semantic_symbols_are_distinct([p, arrow_p, underlined_p, arrow_underlined_p])
            {
                continue;
            }
            let upper_sector_residual =
                strip_boundary_pair_residual(patch, p, arrow_p, previous_sector_rotation_radians);
            let lower_sector_residual = strip_boundary_pair_residual(
                patch,
                underlined_p,
                arrow_underlined_p,
                previous_sector_rotation_radians,
            );
            let sector_residual = upper_sector_residual.max(lower_sector_residual);
            if best
                .map(|(best_residual, _, _, _, _)| sector_residual < best_residual)
                .unwrap_or(true)
            {
                best = Some((
                    sector_residual,
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                ));
            }
        }
    }
    best
}

#[derive(Debug, Clone, Copy)]
struct TestFig6cSourceGraphReport {
    map: TestSection42StripSymbolMap,
    candidate_count: usize,
    full_residual: f64,
    c0_residual: f64,
    eq11_residual: f64,
    eq12_residual: f64,
    eq13_residual: f64,
}

fn best_fig6c_eq12_normalized_source_graph(
    patch: &Kp20MacroPatch,
    strips: &[TestSection42StripMap],
    exact_pairs: &[(TestSection42StripMap, TestSection42StripMap)],
    previous_sector_rotation_radians: f64,
) -> Option<(
    TestFig6cSourceGraphReport,
    TestFig6cSourceGraphReport,
    TestFig6cSourceGraphReport,
    TestFig6cSourceGraphReport,
    usize,
)> {
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
    let mut best = None;
    let mut best_c0 = None;
    let mut best_eq13 = None;
    let mut best_c0_closed = None;
    let mut c0_closed_count = 0_usize;
    let mut candidate_count = 0_usize;
    for &(p, arrow_p) in &sector_pairs {
        for &underlined_p in strips {
            for &arrow_underlined_p in strips {
                if !fig6c_semantic_symbols_are_distinct([
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                ]) {
                    continue;
                }
                candidate_count += 1;
                let map = TestSection42StripSymbolMap {
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                    previous_sector_rotation_radians,
                };
                let normalized = fig6c_eq12_eq13_normalized_source_patch(patch, map);
                let report = measure_fig6c_paper_direction_rows(&normalized, map);
                if report.c0_residual < 1.0e-12 {
                    c0_closed_count += 1;
                    if best_c0_closed
                        .map(|best_report: TestFig6cSourceGraphReport| {
                            report.full_residual < best_report.full_residual
                        })
                        .unwrap_or(true)
                    {
                        best_c0_closed = Some(TestFig6cSourceGraphReport {
                            candidate_count,
                            ..report
                        });
                    }
                }
                if best
                    .map(|best_report: TestFig6cSourceGraphReport| {
                        report.full_residual < best_report.full_residual
                    })
                    .unwrap_or(true)
                {
                    best = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                }
                if best_c0
                    .map(|best_report: TestFig6cSourceGraphReport| {
                        report.c0_residual < best_report.c0_residual
                    })
                    .unwrap_or(true)
                {
                    best_c0 = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                }
                if best_eq13
                    .map(|best_report: TestFig6cSourceGraphReport| {
                        report.eq13_residual < best_report.eq13_residual
                    })
                    .unwrap_or(true)
                {
                    best_eq13 = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                }
            }
        }
    }
    match (best, best_c0, best_eq13, best_c0_closed) {
        (
            Some(mut report),
            Some(mut c0_report),
            Some(mut eq13_report),
            Some(mut c0_closed_report),
        ) => {
            report.candidate_count = candidate_count;
            c0_report.candidate_count = candidate_count;
            eq13_report.candidate_count = candidate_count;
            c0_closed_report.candidate_count = candidate_count;
            Some((
                report,
                c0_report,
                eq13_report,
                c0_closed_report,
                c0_closed_count,
            ))
        }
        _ => None,
    }
}

fn best_fig6c_shared_row_projected_source_graph(
    patch: &Kp20MacroPatch,
    strips: &[TestSection42StripMap],
    exact_pairs: &[(TestSection42StripMap, TestSection42StripMap)],
    previous_sector_rotation_radians: f64,
) -> Option<(
    TestFig6cSourceGraphReport,
    f64,
    f64,
    TestFig6cSourceGraphReport,
    f64,
    f64,
    usize,
    usize,
    TestFig6cSourceGraphReport,
    f64,
    usize,
    TestFig6cSourceGraphReport,
)> {
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
    let mut best = None;
    let mut best_source_delta = 0.0_f64;
    let mut best_overlap = f64::INFINITY;
    let mut best_overlap_report = None;
    let mut best_overlap_source_delta = 0.0_f64;
    let mut minimum_overlap = f64::INFINITY;
    let mut overlap_closed_count = 0_usize;
    let mut independent_count = 0_usize;
    let mut independent_best = None;
    let mut independent_best_overlap = f64::INFINITY;
    let mut c0_eq11_eq12_closed_count = 0_usize;
    let mut best_c0_eq11_eq12_closed = None;
    let mut candidate_count = 0_usize;
    for &(p, arrow_p) in &sector_pairs {
        for &underlined_p in strips {
            for &arrow_underlined_p in strips {
                if !fig6c_semantic_symbols_are_distinct([
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                ]) {
                    continue;
                }
                candidate_count += 1;
                let map = TestSection42StripSymbolMap {
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                    previous_sector_rotation_radians,
                };
                let storage_is_independent =
                    fig6c_shared_row_projection_storage_is_independent(map);
                let (projected, source_delta) = fig6c_shared_row_projected_source_patch(patch, map);
                let overlap = fig6c_eq12_eq13_overlap_residual(&projected, map);
                let report = measure_fig6c_paper_direction_rows(&projected, map);
                if overlap < 1.0e-12 {
                    overlap_closed_count += 1;
                }
                if overlap < minimum_overlap {
                    best_overlap_report = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                    best_overlap_source_delta = source_delta;
                    minimum_overlap = overlap;
                }
                if storage_is_independent {
                    independent_count += 1;
                    if independent_best
                        .map(|best_report: TestFig6cSourceGraphReport| {
                            report.full_residual < best_report.full_residual
                        })
                        .unwrap_or(true)
                    {
                        independent_best = Some(TestFig6cSourceGraphReport {
                            candidate_count,
                            ..report
                        });
                        independent_best_overlap = overlap;
                    }
                }
                if report.c0_residual < 1.0e-12
                    && report.eq11_residual < 1.0e-12
                    && report.eq12_residual < 1.0e-12
                {
                    c0_eq11_eq12_closed_count += 1;
                    if best_c0_eq11_eq12_closed
                        .map(|best_report: TestFig6cSourceGraphReport| {
                            report.eq13_residual < best_report.eq13_residual
                        })
                        .unwrap_or(true)
                    {
                        best_c0_eq11_eq12_closed = Some(TestFig6cSourceGraphReport {
                            candidate_count,
                            ..report
                        });
                    }
                }
                if best
                    .map(|best_report: TestFig6cSourceGraphReport| {
                        report.full_residual < best_report.full_residual
                    })
                    .unwrap_or(true)
                {
                    best = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                    best_source_delta = source_delta;
                    best_overlap = overlap;
                }
            }
        }
    }
    match (best, best_overlap_report) {
        (Some(mut report), Some(mut overlap_report)) => {
            report.candidate_count = candidate_count;
            overlap_report.candidate_count = candidate_count;
            let mut independent_report = independent_best.unwrap_or(report);
            independent_report.candidate_count = candidate_count;
            let mut c0_eq11_eq12_closed_report = best_c0_eq11_eq12_closed.unwrap_or(report);
            c0_eq11_eq12_closed_report.candidate_count = candidate_count;
            Some((
                report,
                best_source_delta,
                best_overlap,
                overlap_report,
                best_overlap_source_delta,
                minimum_overlap,
                overlap_closed_count,
                independent_count,
                independent_report,
                independent_best_overlap,
                c0_eq11_eq12_closed_count,
                c0_eq11_eq12_closed_report,
            ))
        }
        _ => None,
    }
}

fn fig6c_shared_row_projection_storage_is_independent(map: TestSection42StripSymbolMap) -> bool {
    let source_p21 = map.p.stored_tau_coeff(2, 1);
    let source_arrow_p21 = map.arrow_p.stored_tau_coeff(2, 1);
    if source_p21 == source_arrow_p21 {
        return false;
    }

    let mut later_writes = Vec::new();
    for strip in [map.underlined_p, map.arrow_underlined_p] {
        for j in 0..=2 {
            later_writes.push(strip.stored_tau_coeff(0, j));
            if j >= 1 {
                later_writes.push(strip.stored_tau_coeff(5, j));
                later_writes.push(strip.stored_tau_coeff(1, j));
            }
        }
    }
    for strip in [map.underlined_p, map.arrow_underlined_p] {
        for i in [5_usize, 6] {
            later_writes.push(strip.stored_tau_coeff(i, 1));
        }
    }
    later_writes.push(map.p.stored_tau_coeff(4, 1));
    later_writes.push(map.arrow_p.stored_tau_coeff(4, 1));
    !later_writes.contains(&source_p21) && !later_writes.contains(&source_arrow_p21)
}

fn best_fig6c_source_row_basis_projected_graph(
    patch: &Kp20MacroPatch,
    strips: &[TestSection42StripMap],
    exact_pairs: &[(TestSection42StripMap, TestSection42StripMap)],
    previous_sector_rotation_radians: f64,
) -> Option<(TestFig6cSourceGraphReport, f64, f64)> {
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
    let mut best = None;
    let mut best_source_delta = 0.0_f64;
    let mut best_lsq_residual = f64::INFINITY;
    let mut candidate_count = 0_usize;
    for &(p, arrow_p) in &sector_pairs {
        for &underlined_p in strips {
            for &arrow_underlined_p in strips {
                if !fig6c_semantic_symbols_are_distinct([
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                ]) {
                    continue;
                }
                candidate_count += 1;
                let map = TestSection42StripSymbolMap {
                    p,
                    arrow_p,
                    underlined_p,
                    arrow_underlined_p,
                    previous_sector_rotation_radians,
                };
                let (projected, max_source_delta, row_lsq_residual) =
                    fig6c_source_row_basis_projected_source_patch(patch, map);
                let report = measure_fig6c_paper_direction_rows(&projected, map);
                if best
                    .map(|best_report: TestFig6cSourceGraphReport| {
                        report.full_residual < best_report.full_residual
                    })
                    .unwrap_or(true)
                {
                    best = Some(TestFig6cSourceGraphReport {
                        candidate_count,
                        ..report
                    });
                    best_source_delta = max_source_delta;
                    best_lsq_residual = row_lsq_residual;
                }
            }
        }
    }
    best.map(|mut report| {
        report.candidate_count = candidate_count;
        (report, best_source_delta, best_lsq_residual)
    })
}

#[derive(Debug, Clone, Copy)]
struct TestFig6cSourceRowVariable {
    strip: TestSection42StripMap,
    paper_i: usize,
    paper_j: usize,
    axis: usize,
}

#[derive(Debug, Clone, Copy)]
struct TestSquareDomainActiveConstraint {
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    local_axis: usize,
    target: f64,
}

#[derive(Debug, Clone, Copy)]
struct TestOrientationActiveConstraint {
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    target_det: f64,
}

#[derive(Debug, Clone, Copy)]
struct TestSquareActiveSetProjectionReport {
    patch: Kp20MacroPatch,
    iterations: usize,
    active_constraints: usize,
    orientation_constraints: usize,
    max_delta: f64,
    linear_residual: f64,
}

fn fig6c_source_row_basis_projected_source_patch(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> (Kp20MacroPatch, f64, f64) {
    let variables = fig6c_source_row_basis_variables(map);
    fig6c_project_source_patch_with_variables(patch, map, &variables)
}

fn fig6c_project_source_patch_with_variables(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    variables: &[TestFig6cSourceRowVariable],
) -> (Kp20MacroPatch, f64, f64) {
    let angle = map.previous_sector_rotation_radians;
    let mut base = *patch;
    apply_paper_eq12_strip_construction(&mut base, map.p, map.underlined_p, angle);
    apply_paper_eq12_strip_construction(&mut base, map.arrow_p, map.arrow_underlined_p, angle);
    apply_paper_eq13_strip_construction(&mut base, map);
    apply_paper_eq12_strip_construction(&mut base, map.p, map.underlined_p, angle);
    apply_paper_eq12_strip_construction(&mut base, map.arrow_p, map.arrow_underlined_p, angle);
    apply_paper_eq13_strip_construction(&mut base, map);

    let rows = fig6c_full_paper_direction_rows(map);
    let base_residuals = rows
        .iter()
        .map(|row| scalar_residual(&base, row))
        .collect::<Vec<_>>();
    let mut design = vec![vec![0.0_f64; variables.len()]; rows.len()];
    for (variable_index, variable) in variables.iter().copied().enumerate() {
        let mut perturbed = base;
        apply_fig6c_source_row_variable_delta(&mut perturbed, variable, 1.0, angle);
        apply_paper_eq12_strip_construction(&mut perturbed, map.p, map.underlined_p, angle);
        apply_paper_eq12_strip_construction(
            &mut perturbed,
            map.arrow_p,
            map.arrow_underlined_p,
            angle,
        );
        apply_paper_eq13_strip_construction(&mut perturbed, map);
        for (row_index, row) in rows.iter().enumerate() {
            design[row_index][variable_index] =
                scalar_residual(&perturbed, row) - base_residuals[row_index];
        }
    }
    let rhs = base_residuals
        .iter()
        .map(|residual| -residual)
        .collect::<Vec<_>>();
    let solution = solve_test_regularized_least_squares(&design, &rhs, 1.0e-12);
    let row_lsq_residual = max_dense_row_residual(&design, &solution, &rhs);

    let mut projected = base;
    for (variable, delta) in variables.iter().copied().zip(solution.iter().copied()) {
        apply_fig6c_source_row_variable_delta(&mut projected, variable, delta, angle);
    }
    apply_paper_eq12_strip_construction(&mut projected, map.p, map.underlined_p, angle);
    apply_paper_eq12_strip_construction(&mut projected, map.arrow_p, map.arrow_underlined_p, angle);
    apply_paper_eq13_strip_construction(&mut projected, map);

    let max_source_delta = solution
        .chunks(2)
        .map(|chunk| {
            let x = chunk[0];
            let y = chunk.get(1).copied().unwrap_or(0.0);
            (x * x + y * y).sqrt()
        })
        .fold(0.0_f64, f64::max);
    (projected, max_source_delta, row_lsq_residual)
}

fn fig6c_source_row_basis_variables(
    map: TestSection42StripSymbolMap,
) -> Vec<TestFig6cSourceRowVariable> {
    let mut variables = Vec::new();
    for strip in [map.p, map.arrow_p] {
        for paper_i in [2_usize, 4, 6] {
            for axis in [0_usize, 1] {
                variables.push(TestFig6cSourceRowVariable {
                    strip,
                    paper_i,
                    paper_j: 1,
                    axis,
                });
            }
        }
    }
    variables
}

fn apply_fig6c_source_row_variable_delta(
    patch: &mut Kp20MacroPatch,
    variable: TestFig6cSourceRowVariable,
    delta: f64,
    previous_sector_rotation_radians: f64,
) {
    let mut point = transformed_strip_point(
        patch,
        variable.strip,
        variable.paper_i,
        variable.paper_j,
        previous_sector_rotation_radians,
    );
    if variable.axis == 0 {
        point.x += delta;
    } else {
        point.y += delta;
    }
    set_transformed_strip_point(
        patch,
        variable.strip,
        variable.paper_i,
        variable.paper_j,
        point,
        previous_sector_rotation_radians,
    );
}

fn fig6c_boundary_row_basis_projected_source_patch(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> (Kp20MacroPatch, f64, f64) {
    let variables = fig6c_boundary_row_basis_variables(map);
    fig6c_project_source_patch_with_variables(patch, map, &variables)
}

fn fig6c_square_domain_active_set_projection(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    valence: usize,
    sector: usize,
    samples: &[f64],
) -> TestSquareActiveSetProjectionReport {
    fig6c_square_domain_active_set_projection_with_rows(patch, map, valence, sector, samples, &[])
}

fn fig6c_square_domain_active_set_projection_with_rows(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    valence: usize,
    sector: usize,
    samples: &[f64],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> TestSquareActiveSetProjectionReport {
    fig6c_square_domain_active_set_projection_with_rows_and_step_limit(
        patch, map, valence, sector, samples, extra_rows, None,
    )
}

fn fig6c_square_domain_active_set_projection_with_rows_and_step_limit(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    valence: usize,
    sector: usize,
    samples: &[f64],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
    max_step_delta: Option<f64>,
) -> TestSquareActiveSetProjectionReport {
    fig6c_square_domain_active_set_projection_with_rows_solver(
        patch,
        map,
        valence,
        sector,
        samples,
        extra_rows,
        max_step_delta,
        false,
    )
}

fn fig6c_square_domain_active_set_projection_with_rows_solver(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    valence: usize,
    sector: usize,
    samples: &[f64],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
    max_step_delta: Option<f64>,
    normalize_rows: bool,
) -> TestSquareActiveSetProjectionReport {
    let base_rows = fig6c_full_paper_direction_rows(map);
    square_domain_active_set_projection_with_base_and_extra_rows_solver(
        patch,
        &base_rows,
        valence,
        sector,
        samples,
        extra_rows,
        max_step_delta,
        normalize_rows,
    )
}

fn square_domain_active_set_projection_with_base_and_extra_rows_solver(
    patch: &Kp20MacroPatch,
    base_rows: &[Kp20PlanarTauScalarLinearConstraint],
    valence: usize,
    sector: usize,
    samples: &[f64],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
    max_step_delta: Option<f64>,
    normalize_rows: bool,
) -> TestSquareActiveSetProjectionReport {
    let mut projected = *patch;
    let mut total_max_delta = 0.0_f64;
    let mut last_active_constraints = 0_usize;
    let mut last_orientation_constraints = 0_usize;
    let mut last_linear_residual = 0.0_f64;
    const MAX_ACTIVE_SET_ITERATIONS: usize = 256;
    for iteration in 0..MAX_ACTIVE_SET_ITERATIONS {
        let active =
            square_domain_active_constraints(&projected, valence, sector, samples, 1.0e-10);
        let orientation_active =
            orientation_active_constraints(&projected, valence, sector, samples, 1.0e-10);
        last_active_constraints = active.len();
        last_orientation_constraints = orientation_active.len();
        let base_residual = scalar_row_set_residual(&projected, base_rows);
        let extra_residual = scalar_row_set_residual(&projected, extra_rows);
        if active.is_empty()
            && orientation_active.is_empty()
            && base_residual.max(extra_residual) < 1.0e-12
        {
            return TestSquareActiveSetProjectionReport {
                patch: projected,
                iterations: iteration,
                active_constraints: 0,
                orientation_constraints: 0,
                max_delta: total_max_delta,
                linear_residual: last_linear_residual,
            };
        }
        let (next, max_delta, linear_residual) = solve_square_domain_active_step_with_base_rows(
            &projected,
            valence,
            sector,
            &active,
            &orientation_active,
            base_rows,
            extra_rows,
            max_step_delta,
            normalize_rows,
        );
        projected = next;
        total_max_delta = total_max_delta.max(max_delta);
        last_linear_residual = linear_residual;
    }
    TestSquareActiveSetProjectionReport {
        patch: projected,
        iterations: MAX_ACTIVE_SET_ITERATIONS,
        active_constraints: last_active_constraints,
        orientation_constraints: last_orientation_constraints,
        max_delta: total_max_delta,
        linear_residual: last_linear_residual,
    }
}

fn square_domain_active_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
    tolerance: f64,
) -> Vec<TestSquareDomainActiveConstraint> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let piece = local_patch.piece(piece_id);
        for &u in samples {
            let bu = test_bernstein_bi6_all(u);
            for &v in samples {
                let bv = test_bernstein_bi6_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for (i, &bu_i) in bu.iter().enumerate() {
                    for (j, &bv_j) in bv.iter().enumerate() {
                        let weight = bu_i * bv_j;
                        let cp = piece.cps[i][j];
                        x += weight * cp.x;
                        y += weight * cp.y;
                    }
                }
                for (local_axis, value) in [(0_usize, x), (1_usize, y)] {
                    if value < -tolerance {
                        constraints.push(TestSquareDomainActiveConstraint {
                            piece_id,
                            u,
                            v,
                            local_axis,
                            target: 0.0,
                        });
                    } else if value > 1.0 + tolerance {
                        constraints.push(TestSquareDomainActiveConstraint {
                            piece_id,
                            u,
                            v,
                            local_axis,
                            target: 1.0,
                        });
                    }
                }
            }
        }
    }
    constraints
}

fn orientation_active_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
    tolerance: f64,
) -> Vec<TestOrientationActiveConstraint> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let piece = local_patch.piece(piece_id);
        for &u in samples {
            for &v in samples {
                let (du, dv) = sample_local_patch_derivatives(piece, u, v);
                let det = du.x * dv.y - du.y * dv.x;
                if det < -tolerance {
                    constraints.push(TestOrientationActiveConstraint {
                        piece_id,
                        u,
                        v,
                        target_det: 1.0e-3,
                    });
                }
            }
        }
    }
    constraints
}

fn square_domain_guard_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
    guard_distance: f64,
) -> Vec<TestSquareDomainActiveConstraint> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let piece = local_patch.piece(piece_id);
        for &u in samples {
            let bu = test_bernstein_bi6_all(u);
            for &v in samples {
                let bv = test_bernstein_bi6_all(v);
                let mut values = [0.0_f64; 2];
                for (i, &bu_i) in bu.iter().enumerate() {
                    for (j, &bv_j) in bv.iter().enumerate() {
                        let weight = bu_i * bv_j;
                        let cp = piece.cps[i][j];
                        values[0] += weight * cp.x;
                        values[1] += weight * cp.y;
                    }
                }
                for (local_axis, value) in values.into_iter().enumerate() {
                    if value <= guard_distance || value >= 1.0 - guard_distance {
                        constraints.push(TestSquareDomainActiveConstraint {
                            piece_id,
                            u,
                            v,
                            local_axis,
                            target: value.clamp(0.0, 1.0),
                        });
                    }
                }
            }
        }
    }
    constraints
}

fn orientation_guard_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
    determinant_floor: f64,
) -> Vec<TestOrientationActiveConstraint> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let piece = local_patch.piece(piece_id);
        for &u in samples {
            for &v in samples {
                let (du, dv) = sample_local_patch_derivatives(piece, u, v);
                let det = du.x * dv.y - du.y * dv.x;
                if det < determinant_floor {
                    constraints.push(TestOrientationActiveConstraint {
                        piece_id,
                        u,
                        v,
                        target_det: det.max(0.0),
                    });
                }
            }
        }
    }
    constraints
}

fn sample_local_patch_derivatives(patch: &Kp20Bi6Patch, u: f64, v: f64) -> (Point3, Point3) {
    let bu5 = test_bernstein_bi5_all(u);
    let bu6 = test_bernstein_bi6_all(u);
    let bv5 = test_bernstein_bi5_all(v);
    let bv6 = test_bernstein_bi6_all(v);
    let mut du = Point3::new(0.0, 0.0, 0.0);
    let mut dv = Point3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        for j in 0..=6 {
            let delta = sub_point(patch.cps[i + 1][j], patch.cps[i][j]);
            du = add_point(du, scale_point(delta, 6.0 * bu5[i] * bv6[j]));
        }
    }
    for i in 0..=6 {
        for j in 0..6 {
            let delta = sub_point(patch.cps[i][j + 1], patch.cps[i][j]);
            dv = add_point(dv, scale_point(delta, 6.0 * bu6[i] * bv5[j]));
        }
    }
    (du, dv)
}

fn add_orientation_active_design_row(
    design_row: &mut [f64],
    active_constraint: TestOrientationActiveConstraint,
    du: Point3,
    dv: Point3,
    local_axes: &[[f64; 2]; 2],
    row_weight: f64,
) {
    let bu5 = test_bernstein_bi5_all(active_constraint.u);
    let bu6 = test_bernstein_bi6_all(active_constraint.u);
    let bv5 = test_bernstein_bi5_all(active_constraint.v);
    let bv6 = test_bernstein_bi6_all(active_constraint.v);
    for i in 0..=6 {
        for j in 0..=6 {
            let mut local_sensitivity = [0.0_f64; 2];
            if i > 0 {
                let du_weight = 6.0 * bu5[i - 1] * bv6[j];
                local_sensitivity[0] += du_weight * dv.y;
                local_sensitivity[1] -= du_weight * dv.x;
            }
            if i < 6 {
                let du_weight = -6.0 * bu5[i] * bv6[j];
                local_sensitivity[0] += du_weight * dv.y;
                local_sensitivity[1] -= du_weight * dv.x;
            }
            if j > 0 {
                let dv_weight = 6.0 * bu6[i] * bv5[j - 1];
                local_sensitivity[0] -= dv_weight * du.y;
                local_sensitivity[1] += dv_weight * du.x;
            }
            if j < 6 {
                let dv_weight = -6.0 * bu6[i] * bv5[j];
                local_sensitivity[0] -= dv_weight * du.y;
                local_sensitivity[1] += dv_weight * du.x;
            }
            let coeff = Kp20PlanarTauCoeff::new(active_constraint.piece_id, i, j);
            for global_axis in 0..=1 {
                let variable_index = macro_patch_axis_variable_index(coeff, global_axis);
                design_row[variable_index] += row_weight
                    * (local_sensitivity[0] * local_axes[0][global_axis]
                        + local_sensitivity[1] * local_axes[1][global_axis]);
            }
        }
    }
}

fn solve_square_domain_active_step(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
    valence: usize,
    sector: usize,
    active: &[TestSquareDomainActiveConstraint],
    orientation_active: &[TestOrientationActiveConstraint],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
    max_step_delta: Option<f64>,
    normalize_rows: bool,
) -> (Kp20MacroPatch, f64, f64) {
    let base_rows = fig6c_full_paper_direction_rows(map);
    solve_square_domain_active_step_with_base_rows(
        patch,
        valence,
        sector,
        active,
        orientation_active,
        &base_rows,
        extra_rows,
        max_step_delta,
        normalize_rows,
    )
}

fn solve_square_domain_active_step_with_base_rows(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    active: &[TestSquareDomainActiveConstraint],
    orientation_active: &[TestOrientationActiveConstraint],
    base_rows: &[Kp20PlanarTauScalarLinearConstraint],
    extra_rows: &[Kp20PlanarTauScalarLinearConstraint],
    max_step_delta: Option<f64>,
    normalize_rows: bool,
) -> (Kp20MacroPatch, f64, f64) {
    const ROW_WEIGHT: f64 = 100.0;
    const DOMAIN_WEIGHT: f64 = 1.0;
    const ORIENTATION_WEIGHT: f64 = 1.0;
    let variables = all_macro_patch_axis_variables();
    let mut rows = base_rows.to_vec();
    rows.extend(extra_rows.iter().cloned());
    let row_count = rows.len() + active.len() + orientation_active.len();
    let mut design = vec![vec![0.0_f64; variables.len()]; row_count];
    let mut rhs = vec![0.0_f64; row_count];
    for (row_index, row) in rows.iter().enumerate() {
        rhs[row_index] = -ROW_WEIGHT * scalar_residual(patch, row);
        for &(scalar_coeff, weight) in &row.terms {
            let variable_index =
                macro_patch_axis_variable_index(scalar_coeff.coeff, scalar_coeff.axis);
            design[row_index][variable_index] += ROW_WEIGHT * weight;
        }
    }
    let local_axes = sector_local_axis_components(valence, sector);
    let local_patch = global_planar_patch_to_local(valence, sector, *patch);
    for (active_index, active_constraint) in active.iter().copied().enumerate() {
        let row_index = rows.len() + active_index;
        let current_value = sample_local_patch_axis(
            &local_patch,
            active_constraint.piece_id,
            active_constraint.u,
            active_constraint.v,
            active_constraint.local_axis,
        );
        rhs[row_index] = DOMAIN_WEIGHT * (active_constraint.target - current_value);
        let bu = test_bernstein_bi6_all(active_constraint.u);
        let bv = test_bernstein_bi6_all(active_constraint.v);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let coeff = Kp20PlanarTauCoeff::new(active_constraint.piece_id, i, j);
                let bernstein_weight = bu_i * bv_j;
                for global_axis in 0..=1 {
                    let variable_index = macro_patch_axis_variable_index(coeff, global_axis);
                    design[row_index][variable_index] += DOMAIN_WEIGHT
                        * bernstein_weight
                        * local_axes[active_constraint.local_axis][global_axis];
                }
            }
        }
    }
    let orientation_start = rows.len() + active.len();
    for (active_index, active_constraint) in orientation_active.iter().copied().enumerate() {
        let row_index = orientation_start + active_index;
        let local_piece = local_patch.piece(active_constraint.piece_id);
        let (du, dv) =
            sample_local_patch_derivatives(local_piece, active_constraint.u, active_constraint.v);
        let det = du.x * dv.y - du.y * dv.x;
        rhs[row_index] = ORIENTATION_WEIGHT * (active_constraint.target_det - det);
        add_orientation_active_design_row(
            &mut design[row_index],
            active_constraint,
            du,
            dv,
            &local_axes,
            ORIENTATION_WEIGHT,
        );
    }
    let solution = if normalize_rows {
        solve_test_row_normalized_least_squares(&design, &rhs, 1.0e-10)
    } else {
        solve_test_regularized_least_squares(&design, &rhs, 1.0e-10)
    };
    let linear_residual = max_dense_row_residual(&design, &solution, &rhs);
    let solution_max_delta = solution
        .chunks(2)
        .map(|chunk| {
            let dx = chunk[0];
            let dy = chunk.get(1).copied().unwrap_or(0.0);
            (dx * dx + dy * dy).sqrt()
        })
        .fold(0.0_f64, f64::max);
    let step_scale = max_step_delta
        .filter(|max_delta| solution_max_delta > *max_delta)
        .map(|max_delta| max_delta / solution_max_delta)
        .unwrap_or(1.0);
    let mut projected = *patch;
    let mut max_delta = 0.0_f64;
    for (&(coeff, axis), delta) in variables.iter().zip(solution.iter().copied()) {
        let delta = step_scale * delta;
        let cp = &mut projected.piece_mut(coeff.piece).cps[coeff.i][coeff.j];
        if axis == 0 {
            cp.x += delta;
        } else {
            cp.y += delta;
        }
    }
    for chunk in solution.chunks(2) {
        let dx = step_scale * chunk[0];
        let dy = step_scale * chunk.get(1).copied().unwrap_or(0.0);
        max_delta = max_delta.max((dx * dx + dy * dy).sqrt());
    }
    (projected, max_delta, linear_residual)
}

fn all_macro_patch_axis_variables() -> Vec<(Kp20PlanarTauCoeff, usize)> {
    let mut variables = Vec::with_capacity(4 * 7 * 7 * 2);
    for piece_id in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let coeff = Kp20PlanarTauCoeff::new(piece_id, i, j);
                variables.push((coeff, 0));
                variables.push((coeff, 1));
            }
        }
    }
    variables
}

fn macro_patch_axis_variable_index(coeff: Kp20PlanarTauCoeff, axis: usize) -> usize {
    (((coeff.piece.index() * 7 + coeff.i) * 7 + coeff.j) * 2) + axis
}

fn sample_local_patch_axis(
    patch: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    local_axis: usize,
) -> f64 {
    let piece = patch.piece(piece_id);
    let bu = test_bernstein_bi6_all(u);
    let bv = test_bernstein_bi6_all(v);
    let mut value = 0.0_f64;
    for (i, &bu_i) in bu.iter().enumerate() {
        for (j, &bv_j) in bv.iter().enumerate() {
            let cp = piece.cps[i][j];
            value += bu_i * bv_j * if local_axis == 0 { cp.x } else { cp.y };
        }
    }
    value
}

fn fig6c_boundary_row_basis_variables(
    map: TestSection42StripSymbolMap,
) -> Vec<TestFig6cSourceRowVariable> {
    let mut variables = fig6c_source_row_basis_variables(map);
    for strip in [map.p, map.arrow_p, map.underlined_p, map.arrow_underlined_p] {
        for paper_i in 0..=6 {
            for axis in [0_usize, 1] {
                variables.push(TestFig6cSourceRowVariable {
                    strip,
                    paper_i,
                    paper_j: 0,
                    axis,
                });
            }
        }
    }
    variables
}

fn fig6c_full_paper_direction_rows(
    map: TestSection42StripSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let mut rows = test_section42_strip_c0_rows(map);
    test_add_section42_strip_eq11_degree5_rows(&mut rows, map);
    test_add_section42_strip_paper_direction_eq12_rows(&mut rows, map);
    test_add_section42_strip_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn solve_test_regularized_least_squares(
    design: &[Vec<f64>],
    rhs: &[f64],
    regularization: f64,
) -> Vec<f64> {
    let variable_count = design.first().map(|row| row.len()).unwrap_or(0);
    let mut normal = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut target = vec![0.0_f64; variable_count];
    for (row, &rhs_value) in design.iter().zip(rhs) {
        // Sparse-support accumulation; same value-identity argument
        // as the production `solve_regularized_least_squares`
        // (PK-289): skipped terms are exact ±0.0 products.
        let support: Vec<usize> = (0..variable_count).filter(|&i| row[i] != 0.0).collect();
        for &i in &support {
            target[i] += row[i] * rhs_value;
            for &j in &support {
                normal[i][j] += row[i] * row[j];
            }
        }
    }
    for (i, row) in normal.iter_mut().enumerate() {
        row[i] += regularization;
    }
    solve_test_dense_linear_system(normal, target).unwrap_or_else(|| vec![0.0; variable_count])
}

fn solve_test_row_normalized_least_squares(
    design: &[Vec<f64>],
    rhs: &[f64],
    regularization: f64,
) -> Vec<f64> {
    let normalized_design = design
        .iter()
        .map(|row| {
            let norm = row.iter().map(|value| value * value).sum::<f64>().sqrt();
            if norm > 0.0 {
                row.iter().map(|value| value / norm).collect::<Vec<_>>()
            } else {
                row.clone()
            }
        })
        .collect::<Vec<_>>();
    let normalized_rhs = design
        .iter()
        .zip(rhs.iter().copied())
        .map(|(row, rhs_value)| {
            let norm = row.iter().map(|value| value * value).sum::<f64>().sqrt();
            if norm > 0.0 {
                rhs_value / norm
            } else {
                rhs_value
            }
        })
        .collect::<Vec<_>>();
    solve_test_regularized_least_squares(&normalized_design, &normalized_rhs, regularization)
}

fn solve_test_dense_linear_system(
    mut matrix: Vec<Vec<f64>>,
    mut rhs: Vec<f64>,
) -> Option<Vec<f64>> {
    let n = rhs.len();
    for pivot in 0..n {
        let mut best = pivot;
        let mut best_abs = matrix[pivot][pivot].abs();
        for (candidate, row) in matrix.iter().enumerate().skip(pivot + 1) {
            let value = row[pivot].abs();
            if value > best_abs {
                best = candidate;
                best_abs = value;
            }
        }
        if best_abs < 1.0e-18 {
            return None;
        }
        if best != pivot {
            matrix.swap(best, pivot);
            rhs.swap(best, pivot);
        }
        let pivot_value = matrix[pivot][pivot];
        for value in matrix[pivot].iter_mut().skip(pivot) {
            *value /= pivot_value;
        }
        rhs[pivot] /= pivot_value;
        let pivot_row = matrix[pivot].clone();
        for row_index in 0..n {
            if row_index == pivot {
                continue;
            }
            let factor = matrix[row_index][pivot];
            if factor == 0.0 {
                continue;
            }
            for column in pivot..n {
                matrix[row_index][column] -= factor * pivot_row[column];
            }
            rhs[row_index] -= factor * rhs[pivot];
        }
    }
    Some(rhs)
}

fn max_dense_row_residual(design: &[Vec<f64>], solution: &[f64], rhs: &[f64]) -> f64 {
    design
        .iter()
        .zip(rhs)
        .map(|(row, &target)| {
            let value = row
                .iter()
                .zip(solution)
                .map(|(&coeff, &variable)| coeff * variable)
                .sum::<f64>();
            (value - target).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn blend_macro_patch_points(
    left: &Kp20MacroPatch,
    right: &Kp20MacroPatch,
    right_weight: f64,
) -> Kp20MacroPatch {
    let left_weight = 1.0 - right_weight;
    let mut blended = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let left_point = left.piece(piece_id).cps[i][j];
                let right_point = right.piece(piece_id).cps[i][j];
                blended.piece_mut(piece_id).cps[i][j] = Point3::new(
                    left_weight * left_point.x + right_weight * right_point.x,
                    left_weight * left_point.y + right_weight * right_point.y,
                    left_weight * left_point.z + right_weight * right_point.z,
                );
            }
        }
    }
    blended
}

fn fig6c_eq12_normalized_source_patch(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> Kp20MacroPatch {
    let mut normalized = *patch;
    apply_paper_eq12_strip_construction(
        &mut normalized,
        map.p,
        map.underlined_p,
        map.previous_sector_rotation_radians,
    );
    apply_paper_eq12_strip_construction(
        &mut normalized,
        map.arrow_p,
        map.arrow_underlined_p,
        map.previous_sector_rotation_radians,
    );
    normalized
}

fn fig6c_eq12_eq13_normalized_source_patch(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> Kp20MacroPatch {
    let mut normalized = fig6c_eq12_normalized_source_patch(patch, map);
    apply_paper_eq13_strip_construction(&mut normalized, map);
    normalized
}

fn fig6c_shared_row_projected_source_patch(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> (Kp20MacroPatch, f64) {
    let mut projected = *patch;
    let angle = map.previous_sector_rotation_radians;
    apply_paper_eq12_strip_construction(&mut projected, map.p, map.underlined_p, angle);
    apply_paper_eq12_strip_construction(&mut projected, map.arrow_p, map.arrow_underlined_p, angle);
    apply_paper_eq13_strip_construction(&mut projected, map);

    let p21 = transformed_strip_point(&projected, map.p, 2, 1, angle);
    let arrow_p21 = transformed_strip_point(&projected, map.arrow_p, 2, 1, angle);
    let p41 = transformed_strip_point(&projected, map.p, 4, 1, angle);
    let arrow_p41 = transformed_strip_point(&projected, map.arrow_p, 4, 1, angle);
    let p61 = transformed_strip_point(&projected, map.p, 6, 1, angle);
    let arrow_p61 = transformed_strip_point(&projected, map.arrow_p, 6, 1, angle);
    let underlined_51 = transformed_strip_point(&projected, map.underlined_p, 5, 1, angle);
    let arrow_underlined_51 =
        transformed_strip_point(&projected, map.arrow_underlined_p, 5, 1, angle);

    let projected_p21 = sub_point(p41, scale_point(sub_point(underlined_51, p61), 4.0));
    let projected_arrow_p21 = sub_point(
        arrow_p41,
        scale_point(sub_point(arrow_underlined_51, arrow_p61), 4.0),
    );
    let max_source_delta =
        distance(p21, projected_p21).max(distance(arrow_p21, projected_arrow_p21));
    set_transformed_strip_point(&mut projected, map.p, 2, 1, projected_p21, angle);
    set_transformed_strip_point(
        &mut projected,
        map.arrow_p,
        2,
        1,
        projected_arrow_p21,
        angle,
    );

    apply_paper_eq12_strip_construction(&mut projected, map.p, map.underlined_p, angle);
    apply_paper_eq12_strip_construction(&mut projected, map.arrow_p, map.arrow_underlined_p, angle);
    apply_paper_eq13_strip_construction(&mut projected, map);
    (projected, max_source_delta)
}

fn apply_paper_eq12_strip_construction(
    patch: &mut Kp20MacroPatch,
    source: TestSection42StripMap,
    receiving: TestSection42StripMap,
    previous_sector_rotation_radians: f64,
) {
    for j in 0..=2 {
        let row6 = transformed_strip_point(patch, source, 6, j, previous_sector_rotation_radians);
        set_transformed_strip_point(
            patch,
            receiving,
            0,
            j,
            row6,
            previous_sector_rotation_radians,
        );
        if j >= 1 {
            let row4 =
                transformed_strip_point(patch, source, 4, j, previous_sector_rotation_radians);
            let row2 =
                transformed_strip_point(patch, source, 2, j, previous_sector_rotation_radians);
            let quarter_delta = scale_point(sub_point(row4, row2), 0.25);
            set_transformed_strip_point(
                patch,
                receiving,
                5,
                j,
                add_point(row6, quarter_delta),
                previous_sector_rotation_radians,
            );
            set_transformed_strip_point(
                patch,
                receiving,
                1,
                j,
                sub_point(row6, quarter_delta),
                previous_sector_rotation_radians,
            );
        }
    }
}

fn apply_paper_eq13_strip_construction(
    patch: &mut Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = (0.5 * b0_near_origin) / 30.0;
    let angle = map.previous_sector_rotation_radians;

    for i in [5_usize, 6] {
        let underlined_i0 = transformed_strip_point(patch, map.underlined_p, i, 0, angle);
        let underlined_i2 = transformed_strip_point(patch, map.underlined_p, i, 2, angle);
        let arrow_underlined_i2 =
            transformed_strip_point(patch, map.arrow_underlined_p, i, 2, angle);
        let quarter_delta = scale_point(sub_point(underlined_i2, arrow_underlined_i2), 0.25);
        set_transformed_strip_point(
            patch,
            map.underlined_p,
            i,
            1,
            add_point(underlined_i0, quarter_delta),
            angle,
        );
        set_transformed_strip_point(
            patch,
            map.arrow_underlined_p,
            i,
            1,
            sub_point(underlined_i0, quarter_delta),
            angle,
        );
    }

    let p40 = transformed_strip_point(patch, map.p, 4, 0, angle);
    let p42 = transformed_strip_point(patch, map.p, 4, 2, angle);
    let arrow_p42 = transformed_strip_point(patch, map.arrow_p, 4, 2, angle);
    let p50 = transformed_strip_point(patch, map.p, 5, 0, angle);
    let p60 = transformed_strip_point(patch, map.p, 6, 0, angle);
    let p52 = transformed_strip_point(patch, map.p, 5, 2, angle);
    let arrow_p52 = transformed_strip_point(patch, map.arrow_p, 5, 2, angle);
    let p62 = transformed_strip_point(patch, map.p, 6, 2, angle);
    let arrow_p62 = transformed_strip_point(patch, map.arrow_p, 6, 2, angle);

    let quarter = scale_point(sub_point(p42, arrow_p42), 0.25);
    let w3_term = scale_point(sub_point(p50, p60), w3);
    let w2_term = scale_point(sub_point(p52, arrow_p52), w2);
    let w1_term = scale_point(sub_point(p62, arrow_p62), w1);
    let correction = add_point(sub_point(add_point(quarter, w2_term), w3_term), w1_term);

    set_transformed_strip_point(patch, map.p, 4, 1, add_point(p40, correction), angle);
    set_transformed_strip_point(patch, map.arrow_p, 4, 1, sub_point(p40, correction), angle);
}

fn fig6c_eq12_eq13_overlap_residual(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> f64 {
    let angle = map.previous_sector_rotation_radians;
    let plain_eq12 = add_point(
        transformed_strip_point(patch, map.p, 6, 1, angle),
        scale_point(
            sub_point(
                transformed_strip_point(patch, map.p, 4, 1, angle),
                transformed_strip_point(patch, map.p, 2, 1, angle),
            ),
            0.25,
        ),
    );
    let arrow_eq12 = add_point(
        transformed_strip_point(patch, map.arrow_p, 6, 1, angle),
        scale_point(
            sub_point(
                transformed_strip_point(patch, map.arrow_p, 4, 1, angle),
                transformed_strip_point(patch, map.arrow_p, 2, 1, angle),
            ),
            0.25,
        ),
    );
    let underlined_base = transformed_strip_point(patch, map.underlined_p, 5, 0, angle);
    let quarter_delta = scale_point(
        sub_point(
            transformed_strip_point(patch, map.underlined_p, 5, 2, angle),
            transformed_strip_point(patch, map.arrow_underlined_p, 5, 2, angle),
        ),
        0.25,
    );
    let plain_eq13 = add_point(underlined_base, quarter_delta);
    let arrow_eq13 = sub_point(underlined_base, quarter_delta);

    distance(plain_eq12, plain_eq13).max(distance(arrow_eq12, arrow_eq13))
}

fn measure_fig6c_paper_direction_rows(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> TestFig6cSourceGraphReport {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let c0_rows = test_section42_strip_c0_rows(map);
    let mut eq11_rows = Vec::new();
    test_add_section42_strip_eq11_degree5_rows(&mut eq11_rows, map);
    let mut eq12_rows = Vec::new();
    test_add_section42_strip_paper_direction_eq12_rows(&mut eq12_rows, map);
    let mut eq13_rows = Vec::new();
    test_add_section42_strip_eq13_rows(&mut eq13_rows, map, b0_near_origin);
    let mut full_rows = Vec::new();
    full_rows.extend(c0_rows.iter().cloned());
    full_rows.extend(eq11_rows.iter().cloned());
    full_rows.extend(eq12_rows.iter().cloned());
    full_rows.extend(eq13_rows.iter().cloned());
    TestFig6cSourceGraphReport {
        map,
        candidate_count: 0,
        full_residual: scalar_row_set_residual(patch, &full_rows),
        c0_residual: scalar_row_set_residual(patch, &c0_rows),
        eq11_residual: scalar_row_set_residual(patch, &eq11_rows),
        eq12_residual: scalar_row_set_residual(patch, &eq12_rows),
        eq13_residual: scalar_row_set_residual(patch, &eq13_rows),
    }
}

fn measure_fig6c_endpoint_aware_section42_rows(
    patch: &Kp20MacroPatch,
    map: TestSection42StripSymbolMap,
) -> TestFig6cSourceGraphReport {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let c0_rows = test_section42_strip_c0_rows(map);
    let mut eq11_rows = Vec::new();
    test_add_section42_strip_eq11_degree5_rows(&mut eq11_rows, map);
    let mut eq12_rows = Vec::new();
    test_add_section42_strip_endpoint_aware_eq12_rows(&mut eq12_rows, map);
    let mut eq13_rows = Vec::new();
    test_add_section42_strip_eq13_rows(&mut eq13_rows, map, b0_near_origin);
    let mut full_rows = Vec::new();
    full_rows.extend(c0_rows.iter().cloned());
    full_rows.extend(eq11_rows.iter().cloned());
    full_rows.extend(eq12_rows.iter().cloned());
    full_rows.extend(eq13_rows.iter().cloned());
    TestFig6cSourceGraphReport {
        map,
        candidate_count: 0,
        full_residual: scalar_row_set_residual(patch, &full_rows),
        c0_residual: scalar_row_set_residual(patch, &c0_rows),
        eq11_residual: scalar_row_set_residual(patch, &eq11_rows),
        eq12_residual: scalar_row_set_residual(patch, &eq12_rows),
        eq13_residual: scalar_row_set_residual(patch, &eq13_rows),
    }
}

fn appendix_a_candidate_section42_strips() -> Vec<TestSection42StripMap> {
    Kp20PieceId::all()
        .into_iter()
        .flat_map(|piece| {
            [-1, 0, 1]
                .into_iter()
                .map(move |sector_rotation_steps| (piece, sector_rotation_steps))
        })
        .flat_map(|(piece, sector_rotation_steps)| {
            TestMacroEdge::ALL.into_iter().flat_map(move |edge| {
                [false, true]
                    .into_iter()
                    .map(move |reversed| TestSection42StripMap {
                        mapped: if sector_rotation_steps == 1 {
                            Kp20Section42MappedPiece::previous(piece)
                        } else {
                            Kp20Section42MappedPiece::current(piece)
                        },
                        edge,
                        reversed,
                        sector_rotation_steps,
                    })
            })
        })
        .collect()
}

fn exact_strip_boundary_pairs(
    patch: &Kp20MacroPatch,
    strips: &[TestSection42StripMap],
    previous_sector_rotation_radians: f64,
) -> Vec<(TestSection42StripMap, TestSection42StripMap)> {
    let mut pairs = Vec::new();
    for &left in strips {
        for &right in strips {
            if left.mapped.piece == right.mapped.piece
                && left.sector_rotation_steps == right.sector_rotation_steps
            {
                continue;
            }
            if strip_boundary_pair_residual(patch, left, right, previous_sector_rotation_radians)
                < 1.0e-12
            {
                pairs.push((left, right));
            }
        }
    }
    pairs
}

fn fig6c_semantic_symbols_are_distinct(strips: [TestSection42StripMap; 4]) -> bool {
    for left in 0..strips.len() {
        for right in (left + 1)..strips.len() {
            if strips[left].mapped.piece == strips[right].mapped.piece
                && strips[left].sector_rotation_steps == strips[right].sector_rotation_steps
            {
                return false;
            }
        }
    }
    true
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TestFig6cExactPairKind {
    SectorJoin,
    TwoPieceJoin,
}

fn fig6c_exact_pair_kind(
    left: TestSection42StripMap,
    right: TestSection42StripMap,
) -> Option<TestFig6cExactPairKind> {
    let same_sector = left.sector_rotation_steps == right.sector_rotation_steps;
    if same_sector
        && matches!(
            (left.edge, right.edge),
            (TestMacroEdge::J0, TestMacroEdge::J0)
        )
        && left.reversed != right.reversed
    {
        return Some(TestFig6cExactPairKind::SectorJoin);
    }
    if (left.sector_rotation_steps - right.sector_rotation_steps).abs() == 1
        && matches!(
            (left.edge, right.edge),
            (TestMacroEdge::J6, TestMacroEdge::I0) | (TestMacroEdge::I0, TestMacroEdge::J6)
        )
        && left.reversed != right.reversed
    {
        return Some(TestFig6cExactPairKind::SectorJoin);
    }
    if (left.sector_rotation_steps - right.sector_rotation_steps).abs() == 1
        && matches!(
            (left.edge, right.edge),
            (TestMacroEdge::J0, TestMacroEdge::I0) | (TestMacroEdge::I0, TestMacroEdge::J0)
        )
        && left.reversed == right.reversed
    {
        return Some(TestFig6cExactPairKind::SectorJoin);
    }
    if same_sector
        && matches!(
            (left.edge, right.edge),
            (TestMacroEdge::I0, TestMacroEdge::I6) | (TestMacroEdge::I6, TestMacroEdge::I0)
        )
        && left.reversed == right.reversed
    {
        return Some(TestFig6cExactPairKind::TwoPieceJoin);
    }
    if same_sector
        && matches!(
            (left.edge, right.edge),
            (TestMacroEdge::J0, TestMacroEdge::J6) | (TestMacroEdge::J6, TestMacroEdge::J0)
        )
        && left.reversed == right.reversed
    {
        return Some(TestFig6cExactPairKind::TwoPieceJoin);
    }
    None
}

#[derive(Debug, Clone, Copy)]
struct TestEq12Join {
    source_endpoint: usize,
    receiving_endpoint: usize,
    best_endpoint_residual: f64,
}

fn best_strip_eq12_join_residual(
    patch: &Kp20MacroPatch,
    source_closer_to_m: TestSection42StripMap,
    receiving_closer_to_o: TestSection42StripMap,
    previous_sector_rotation_radians: f64,
) -> TestEq12Join {
    let mut best = TestEq12Join {
        source_endpoint: 6,
        receiving_endpoint: 0,
        best_endpoint_residual: f64::INFINITY,
    };
    for source_endpoint in [0_usize, 6] {
        for receiving_endpoint in [0_usize, 6] {
            let residual = (0..=2)
                .map(|j| {
                    let source_point = transformed_strip_point(
                        patch,
                        source_closer_to_m,
                        source_endpoint,
                        j,
                        previous_sector_rotation_radians,
                    );
                    let receiving_point = transformed_strip_point(
                        patch,
                        receiving_closer_to_o,
                        receiving_endpoint,
                        j,
                        previous_sector_rotation_radians,
                    );
                    distance(source_point, receiving_point)
                })
                .fold(0.0_f64, f64::max);
            if residual < best.best_endpoint_residual {
                best = TestEq12Join {
                    source_endpoint,
                    receiving_endpoint,
                    best_endpoint_residual: residual,
                };
            }
        }
    }
    best
}

#[derive(Debug, Clone, Copy)]
struct TestSection42StripMap {
    mapped: Kp20Section42MappedPiece,
    edge: TestMacroEdge,
    reversed: bool,
    sector_rotation_steps: i32,
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
    if strip.sector_rotation_steps != 0 {
        rotate_xy(
            cp,
            previous_sector_rotation_radians * strip.sector_rotation_steps as f64,
        )
    } else {
        cp
    }
}

fn set_transformed_strip_point(
    patch: &mut Kp20MacroPatch,
    strip: TestSection42StripMap,
    paper_i: usize,
    paper_j: usize,
    transformed_point: Point3,
    previous_sector_rotation_radians: f64,
) {
    let stored_point = if strip.sector_rotation_steps != 0 {
        rotate_xy(
            transformed_point,
            -previous_sector_rotation_radians * strip.sector_rotation_steps as f64,
        )
    } else {
        transformed_point
    };
    let (i, j) = strip.stored_index(paper_i, paper_j);
    patch.piece_mut(strip.mapped.piece).cps[i][j] = stored_point;
}

fn add_point(left: Point3, right: Point3) -> Point3 {
    Point3::new(left.x + right.x, left.y + right.y, left.z + right.z)
}

fn sub_point(left: Point3, right: Point3) -> Point3 {
    Point3::new(left.x - right.x, left.y - right.y, left.z - right.z)
}

fn scale_point(point: Point3, scale: f64) -> Point3 {
    Point3::new(point.x * scale, point.y * scale, point.z * scale)
}

#[derive(Debug, Clone)]
struct TestSquareDomainReport {
    samples: usize,
    outside: usize,
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
    max_excess: f64,
    worst_sample: String,
}

#[derive(Debug, Clone)]
struct TestOrientationDetailReport {
    samples: usize,
    negative: usize,
    near_zero: usize,
    min_det: f64,
    max_det: f64,
    worst_sample: String,
}

#[derive(Debug, Clone)]
struct TestPieceDomainOrientationReport {
    piece: Kp20PieceId,
    domain_samples: usize,
    outside: usize,
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
    max_excess: f64,
    worst_domain_sample: String,
    orientation_samples: usize,
    negative: usize,
    near_zero: usize,
    min_det: f64,
    max_det: f64,
    worst_orientation_sample: String,
}

fn test_square_domain_report(
    tau: &Kp20MacroPatch,
    samples: &[f64],
    tolerance: f64,
) -> TestSquareDomainReport {
    let mut report = TestSquareDomainReport {
        samples: 0,
        outside: 0,
        min_u: f64::INFINITY,
        max_u: f64::NEG_INFINITY,
        min_v: f64::INFINITY,
        max_v: f64::NEG_INFINITY,
        max_excess: 0.0,
        worst_sample: String::new(),
    };
    for piece_id in Kp20PieceId::all() {
        let patch = tau.piece(piece_id);
        for &u in samples {
            let bu = test_bernstein_bi6_all(u);
            for &v in samples {
                let bv = test_bernstein_bi6_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for (i, &bu_i) in bu.iter().enumerate() {
                    for (j, &bv_j) in bv.iter().enumerate() {
                        let cp = patch.cps[i][j];
                        let weight = bu_i * bv_j;
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
                if total_excess > report.max_excess {
                    report.max_excess = total_excess;
                    report.worst_sample = format!("{piece_id:?}@({u:.2},{v:.2})=({x:.6},{y:.6})");
                }
            }
        }
    }
    report
}

fn test_piece_domain_orientation_reports(
    tau: &Kp20MacroPatch,
    samples: &[f64],
    tolerance: f64,
) -> Vec<TestPieceDomainOrientationReport> {
    let mut reports = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let patch = tau.piece(piece_id);
        let mut report = TestPieceDomainOrientationReport {
            piece: piece_id,
            domain_samples: 0,
            outside: 0,
            min_u: f64::INFINITY,
            max_u: f64::NEG_INFINITY,
            min_v: f64::INFINITY,
            max_v: f64::NEG_INFINITY,
            max_excess: 0.0,
            worst_domain_sample: String::new(),
            orientation_samples: 0,
            negative: 0,
            near_zero: 0,
            min_det: f64::INFINITY,
            max_det: f64::NEG_INFINITY,
            worst_orientation_sample: String::new(),
        };
        for &u in samples {
            let bu6 = test_bernstein_bi6_all(u);
            let bu5 = test_bernstein_bi5_all(u);
            for &v in samples {
                let bv6 = test_bernstein_bi6_all(v);
                let bv5 = test_bernstein_bi5_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for (i, &bu_i) in bu6.iter().enumerate() {
                    for (j, &bv_j) in bv6.iter().enumerate() {
                        let cp = patch.cps[i][j];
                        let weight = bu_i * bv_j;
                        x += weight * cp.x;
                        y += weight * cp.y;
                    }
                }
                let lower_excess = (-x).max(0.0).max((-y).max(0.0));
                let upper_excess = (x - 1.0).max(0.0).max((y - 1.0).max(0.0));
                let total_excess = lower_excess.max(upper_excess);
                report.domain_samples += 1;
                report.outside += usize::from(total_excess > tolerance);
                report.min_u = report.min_u.min(x);
                report.max_u = report.max_u.max(x);
                report.min_v = report.min_v.min(y);
                report.max_v = report.max_v.max(y);
                if total_excess > report.max_excess {
                    report.max_excess = total_excess;
                    report.worst_domain_sample =
                        format!("{piece_id:?}@({u:.2},{v:.2})=({x:.6},{y:.6})");
                }

                let mut du = Point3::new(0.0, 0.0, 0.0);
                let mut dv = Point3::new(0.0, 0.0, 0.0);
                for i in 0..6 {
                    for (j, &bv_j) in bv6.iter().enumerate() {
                        let delta = sub_point(patch.cps[i + 1][j], patch.cps[i][j]);
                        du = add_point(du, scale_point(delta, 6.0 * bu5[i] * bv_j));
                    }
                }
                for (i, &bu_i) in bu6.iter().enumerate() {
                    for j in 0..6 {
                        let delta = sub_point(patch.cps[i][j + 1], patch.cps[i][j]);
                        dv = add_point(dv, scale_point(delta, 6.0 * bu_i * bv5[j]));
                    }
                }
                let det = du.x * dv.y - du.y * dv.x;
                report.orientation_samples += 1;
                report.negative += usize::from(det < -tolerance);
                report.near_zero += usize::from(det.abs() <= tolerance);
                report.max_det = report.max_det.max(det);
                if det < report.min_det {
                    report.min_det = det;
                    report.worst_orientation_sample = format!("{piece_id:?}@({u:.2},{v:.2})");
                }
            }
        }
        reports.push(report);
    }
    reports
}

fn test_bernstein_bi6_all(t: f64) -> [f64; 7] {
    let u = 1.0 - t;
    [
        u.powi(6),
        6.0 * t * u.powi(5),
        15.0 * t.powi(2) * u.powi(4),
        20.0 * t.powi(3) * u.powi(3),
        15.0 * t.powi(4) * u.powi(2),
        6.0 * t.powi(5) * u,
        t.powi(6),
    ]
}

fn test_bernstein_bi5_all(t: f64) -> [f64; 6] {
    let u = 1.0 - t;
    [
        u.powi(5),
        5.0 * t * u.powi(4),
        10.0 * t.powi(2) * u.powi(3),
        10.0 * t.powi(3) * u.powi(2),
        5.0 * t.powi(4) * u,
        t.powi(5),
    ]
}

fn test_orientation_detail_report(
    tau: &Kp20MacroPatch,
    samples: &[f64],
    tolerance: f64,
) -> TestOrientationDetailReport {
    let mut report = TestOrientationDetailReport {
        samples: 0,
        negative: 0,
        near_zero: 0,
        min_det: f64::INFINITY,
        max_det: f64::NEG_INFINITY,
        worst_sample: String::new(),
    };
    for piece_id in Kp20PieceId::all() {
        let patch = tau.piece(piece_id);
        for &u in samples {
            let bu5 = test_bernstein_bi5_all(u);
            let bu6 = test_bernstein_bi6_all(u);
            for &v in samples {
                let bv5 = test_bernstein_bi5_all(v);
                let bv6 = test_bernstein_bi6_all(v);
                let mut du = Point3::new(0.0, 0.0, 0.0);
                let mut dv = Point3::new(0.0, 0.0, 0.0);
                for i in 0..6 {
                    for j in 0..=6 {
                        let delta = sub_point(patch.cps[i + 1][j], patch.cps[i][j]);
                        du = add_point(du, scale_point(delta, 6.0 * bu5[i] * bv6[j]));
                    }
                }
                for i in 0..=6 {
                    for j in 0..6 {
                        let delta = sub_point(patch.cps[i][j + 1], patch.cps[i][j]);
                        dv = add_point(dv, scale_point(delta, 6.0 * bu6[i] * bv5[j]));
                    }
                }
                let det = du.x * dv.y - du.y * dv.x;
                report.samples += 1;
                report.negative += usize::from(det < -tolerance);
                report.near_zero += usize::from(det.abs() <= tolerance);
                report.max_det = report.max_det.max(det);
                if det < report.min_det {
                    report.min_det = det;
                    report.worst_sample = format!("{piece_id:?}@({u:.2},{v:.2})");
                }
            }
        }
    }
    report
}

fn set_fig8_affine_square_piece(
    patch: &mut Kp20MacroPatch,
    valence: usize,
    sector: usize,
    piece_id: Kp20PieceId,
    outer_boundary_only: bool,
) {
    let origin = match piece_id {
        Kp20PieceId::NearOriginLeft => [0.0, 0.0],
        Kp20PieceId::NearOriginRight => [0.5, 0.0],
        Kp20PieceId::NearMidpointRight => [0.5, 0.5],
        Kp20PieceId::NearMidpointLeft => [0.0, 0.5],
    };
    for i in 0..=6 {
        for j in 0..=6 {
            if outer_boundary_only && i != 6 && j != 6 {
                continue;
            }
            patch.piece_mut(piece_id).cps[i][j] = local_to_global_planar(
                valence,
                sector,
                [origin[0] + i as f64 / 12.0, origin[1] + j as f64 / 12.0],
            );
        }
    }
}

fn fig8_affine_square_patch(valence: usize, sector: usize) -> Kp20MacroPatch {
    let mut patch = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        set_fig8_affine_square_piece(&mut patch, valence, sector, piece_id, false);
    }
    patch
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

fn test_add_section42_strip_endpoint_aware_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
) {
    // KP20 eq. (12) says the non-underlined strip closer to m feeds the
    // underlined strip closer to o. In the corrected oriented Fig. 6c chart
    // that is paper endpoint 6 on p/arrow-p to paper endpoint 0 on the
    // underlined receiving strips.
    test_add_section42_strip_eq12_variant_rows(rows, map, false, 6, 0);
}

fn test_add_section42_strip_paper_direction_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42StripSymbolMap,
) {
    for (source, receiving) in [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ] {
        for j in 0..=2 {
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 0, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 5, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, -0.25),
                    section42_term(source, 2, j, 0.25),
                ],
            );
            test_add_section42_strip_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 1, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, 0.25),
                    section42_term(source, 2, j, -0.25),
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
    let w3 = (0.5 * b0_near_origin) / 30.0;

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
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let strip = test_section42_strip_symbol_map(map, term.symbol);
        let (i, j) = strip.stored_index(term.i, term.j);
        let coeff = Kp20PlanarTauCoeff::new(strip.mapped.piece, i, j);
        if strip.sector_rotation_steps != 0 {
            let angle = map.previous_sector_rotation_radians * strip.sector_rotation_steps as f64;
            let (sin_theta, cos_theta) = angle.sin_cos();
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

    fn stored_tau_coeff(self, paper_i: usize, paper_j: usize) -> Kp20PlanarTauCoeff {
        let (i, j) = self.stored_index(paper_i, paper_j);
        Kp20PlanarTauCoeff::new(self.mapped.piece, i, j)
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

fn transposed_tensor_border_rows(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let split = split_kp20_tensor_border(tensor_border).expect("tensor-border splits");
    let previous = (sector + split.len() - 1) % split.len();
    let mut rows = Vec::with_capacity(42);
    for (piece, half_rows) in [
        (Kp20PieceId::NearMidpointLeft, split[previous].near_c.rows),
        (Kp20PieceId::NearMidpointRight, split[sector].near_c.rows),
    ] {
        for (depth, row) in half_rows.iter().enumerate() {
            for (u_index, cp) in row.iter().enumerate() {
                rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece, u_index, depth),
                    [cp.x, cp.y],
                ));
            }
        }
    }
    rows
}

fn tensor_border_half_assignment_rows(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    previous_uses_near_m: bool,
    current_uses_near_m: bool,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let split = split_kp20_tensor_border(tensor_border).expect("tensor-border splits");
    let previous = (sector + split.len() - 1) % split.len();
    let previous_rows = if previous_uses_near_m {
        split[previous].near_m.rows
    } else {
        split[previous].near_c.rows
    };
    let current_rows = if current_uses_near_m {
        split[sector].near_m.rows
    } else {
        split[sector].near_c.rows
    };
    let mut rows = Vec::with_capacity(42);
    for (piece, half_rows) in [
        (Kp20PieceId::NearMidpointLeft, previous_rows),
        (Kp20PieceId::NearMidpointRight, current_rows),
    ] {
        for (depth, row) in half_rows.iter().enumerate() {
            for (u_index, cp) in row.iter().enumerate() {
                rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece, depth, u_index),
                    [cp.x, cp.y],
                ));
            }
        }
    }
    rows
}

#[derive(Debug, Clone, Copy)]
struct TestTensorBorderChartVariant {
    previous_uses_near_m: bool,
    current_uses_near_m: bool,
    previous_reverse_u: bool,
    current_reverse_u: bool,
    transpose_indices: bool,
    swap_pieces: bool,
}

impl TestTensorBorderChartVariant {
    fn label(self) -> String {
        format!(
            "prev-{}{}-curr-{}{}-{}{}",
            if self.previous_uses_near_m {
                "near-m"
            } else {
                "near-c"
            },
            if self.previous_reverse_u { "-rev" } else { "" },
            if self.current_uses_near_m {
                "near-m"
            } else {
                "near-c"
            },
            if self.current_reverse_u { "-rev" } else { "" },
            if self.transpose_indices {
                "u-depth"
            } else {
                "depth-u"
            },
            if self.swap_pieces { "-swapped" } else { "" },
        )
    }
}

fn tensor_border_chart_variant_rows(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    variant: TestTensorBorderChartVariant,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let split = split_kp20_tensor_border(tensor_border).expect("tensor-border splits");
    let previous = (sector + split.len() - 1) % split.len();
    let previous_rows = if variant.previous_uses_near_m {
        split[previous].near_m.rows
    } else {
        split[previous].near_c.rows
    };
    let current_rows = if variant.current_uses_near_m {
        split[sector].near_m.rows
    } else {
        split[sector].near_c.rows
    };
    let (previous_piece, current_piece) = if variant.swap_pieces {
        (
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        )
    } else {
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        )
    };

    let mut rows = Vec::with_capacity(42);
    for (piece, half_rows, reverse_u) in [
        (previous_piece, previous_rows, variant.previous_reverse_u),
        (current_piece, current_rows, variant.current_reverse_u),
    ] {
        for (depth, row) in half_rows.iter().enumerate() {
            for (u_index, cp) in row.iter().enumerate() {
                let mapped_u = if reverse_u { 6 - u_index } else { u_index };
                let coeff = if variant.transpose_indices {
                    Kp20PlanarTauCoeff::new(piece, mapped_u, depth)
                } else {
                    Kp20PlanarTauCoeff::new(piece, depth, mapped_u)
                };
                rows.push(Kp20PlanarTauLinearConstraint::pin(coeff, [cp.x, cp.y]));
            }
        }
    }
    rows
}

fn tensor_border_chart_variant_rows_from_split_rows(
    previous_rows: [[Point3; 7]; 3],
    current_rows: [[Point3; 7]; 3],
    variant: TestTensorBorderChartVariant,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let (previous_piece, current_piece) = if variant.swap_pieces {
        (
            Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearMidpointLeft,
        )
    } else {
        (
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        )
    };

    let mut rows = Vec::with_capacity(42);
    for (piece, half_rows, reverse_u) in [
        (previous_piece, previous_rows, variant.previous_reverse_u),
        (current_piece, current_rows, variant.current_reverse_u),
    ] {
        for (depth, row) in half_rows.iter().enumerate() {
            for (u_index, cp) in row.iter().enumerate() {
                let mapped_u = if reverse_u { 6 - u_index } else { u_index };
                let coeff = if variant.transpose_indices {
                    Kp20PlanarTauCoeff::new(piece, mapped_u, depth)
                } else {
                    Kp20PlanarTauCoeff::new(piece, depth, mapped_u)
                };
                rows.push(Kp20PlanarTauLinearConstraint::pin(coeff, [cp.x, cp.y]));
            }
        }
    }
    rows
}

fn tensor_chart_coeff_depth(
    coeff: Kp20PlanarTauCoeff,
    variant: TestTensorBorderChartVariant,
) -> usize {
    if variant.transpose_indices {
        coeff.j
    } else {
        coeff.i
    }
}

fn tensor_chart_coeff_u(coeff: Kp20PlanarTauCoeff, variant: TestTensorBorderChartVariant) -> usize {
    if variant.transpose_indices {
        coeff.i
    } else {
        coeff.j
    }
}

fn without_tensor_chart_depth(
    rows: &[Kp20PlanarTauLinearConstraint],
    variant: TestTensorBorderChartVariant,
    omitted_depth: usize,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    rows.iter()
        .filter(|row| {
            let coeff = row.terms[0].0;
            tensor_chart_coeff_depth(coeff, variant) != omitted_depth
        })
        .cloned()
        .collect()
}

fn keep_tensor_chart_depths(
    rows: &[Kp20PlanarTauLinearConstraint],
    variant: TestTensorBorderChartVariant,
    kept_depths: &[usize],
) -> Vec<Kp20PlanarTauLinearConstraint> {
    rows.iter()
        .filter(|row| {
            let coeff = row.terms[0].0;
            kept_depths.contains(&tensor_chart_coeff_depth(coeff, variant))
        })
        .cloned()
        .collect()
}

fn reflect_tensor_chart_depth_across_boundary_for_piece(
    rows: &[Kp20PlanarTauLinearConstraint],
    variant: TestTensorBorderChartVariant,
    piece: Kp20PieceId,
    reflected_depth: usize,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        let coeff = row.terms[0].0;
        if coeff.piece == piece && tensor_chart_coeff_depth(coeff, variant) == reflected_depth {
            let u = tensor_chart_coeff_u(coeff, variant);
            let boundary_target = rows
                .iter()
                .find_map(|candidate| {
                    let candidate_coeff = candidate.terms[0].0;
                    if candidate_coeff.piece == piece
                        && tensor_chart_coeff_depth(candidate_coeff, variant) == 0
                        && tensor_chart_coeff_u(candidate_coeff, variant) == u
                    {
                        Some(candidate.target)
                    } else {
                        None
                    }
                })
                .expect("matching tensor boundary row");
            out.push(Kp20PlanarTauLinearConstraint::pin(
                coeff,
                [
                    2.0 * boundary_target[0] - row.target[0],
                    2.0 * boundary_target[1] - row.target[1],
                ],
            ));
        } else {
            out.push(row.clone());
        }
    }
    out
}

fn scale_tensor_chart_depth_from_boundary_for_piece(
    rows: &[Kp20PlanarTauLinearConstraint],
    variant: TestTensorBorderChartVariant,
    piece: Kp20PieceId,
    scaled_depth: usize,
    scale: f64,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        let coeff = row.terms[0].0;
        if coeff.piece == piece && tensor_chart_coeff_depth(coeff, variant) == scaled_depth {
            let u = tensor_chart_coeff_u(coeff, variant);
            let boundary_target = rows
                .iter()
                .find_map(|candidate| {
                    let candidate_coeff = candidate.terms[0].0;
                    if candidate_coeff.piece == piece
                        && tensor_chart_coeff_depth(candidate_coeff, variant) == 0
                        && tensor_chart_coeff_u(candidate_coeff, variant) == u
                    {
                        Some(candidate.target)
                    } else {
                        None
                    }
                })
                .expect("matching tensor boundary row");
            out.push(Kp20PlanarTauLinearConstraint::pin(
                coeff,
                [
                    boundary_target[0] + scale * (row.target[0] - boundary_target[0]),
                    boundary_target[1] + scale * (row.target[1] - boundary_target[1]),
                ],
            ));
        } else {
            out.push(row.clone());
        }
    }
    out
}

fn rotate_tensor_chart_piece_targets(
    rows: &[Kp20PlanarTauLinearConstraint],
    piece: Kp20PieceId,
    angle: f64,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    rows.iter()
        .map(|row| {
            let coeff = row.terms[0].0;
            if coeff.piece == piece {
                let target = rotate_xy(Point3::new(row.target[0], row.target[1], 0.0), angle);
                Kp20PlanarTauLinearConstraint::pin(coeff, [target.x, target.y])
            } else {
                row.clone()
            }
        })
        .collect()
}

fn test_fig8_midpoint_average_row(
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

fn test_fig8_indicated_average_variant_rows(
    include_diagonal_mirror: bool,
    include_rotated_quadrants: bool,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    if include_diagonal_mirror {
        rows.push(test_fig8_midpoint_average_row(
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
            rows.push(test_fig8_midpoint_average_row(piece_id, 1, 0, 0, 0, 2, 0));
            if include_diagonal_mirror {
                rows.push(test_fig8_midpoint_average_row(piece_id, 0, 1, 0, 0, 0, 2));
            }
        }
    }
    rows
}

#[derive(Debug, Clone, Copy)]
struct TestSection42ChartStrip {
    piece: Kp20PieceId,
    reverse_i: bool,
    reverse_j: bool,
    transpose: bool,
}

impl TestSection42ChartStrip {
    fn label(self) -> String {
        format!(
            "{:?}{}{}{}",
            self.piece,
            if self.reverse_i { "-rev-i" } else { "" },
            if self.reverse_j { "-rev-j" } else { "" },
            if self.transpose { "-transpose" } else { "" },
        )
    }

    fn coeff(self, paper_i: usize, paper_j: usize) -> Kp20PlanarTauCoeff {
        let i = if self.reverse_i { 6 - paper_i } else { paper_i };
        let j = if self.reverse_j { 6 - paper_j } else { paper_j };
        if self.transpose {
            Kp20PlanarTauCoeff::new(self.piece, j, i)
        } else {
            Kp20PlanarTauCoeff::new(self.piece, i, j)
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct TestSection42ChartMap {
    p: TestSection42ChartStrip,
    arrow_p: TestSection42ChartStrip,
    underlined_p: TestSection42ChartStrip,
    arrow_underlined_p: TestSection42ChartStrip,
}

fn test_section42_endpoint_aware_chart_rows(
    map: TestSection42ChartMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    test_add_section42_chart_c0_rows(&mut rows, map);
    test_add_section42_chart_eq11_rows(&mut rows, map);
    test_add_section42_chart_endpoint_aware_eq12_rows(&mut rows, map);
    test_add_section42_chart_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn test_add_section42_chart_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42ChartMap,
) {
    for i in 0..=6 {
        test_add_section42_chart_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        test_add_section42_chart_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

fn test_add_section42_chart_eq11_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42ChartMap,
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
        test_add_section42_chart_vector_rows(rows, map, &terms);
    }
}

fn test_add_section42_chart_endpoint_aware_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42ChartMap,
) {
    for (source, receiving) in [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ] {
        for j in 0..=2 {
            test_add_section42_chart_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 0, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            test_add_section42_chart_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 1, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, -0.25),
                    section42_term(source, 2, j, 0.25),
                ],
            );
            test_add_section42_chart_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 5, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, 0.25),
                    section42_term(source, 2, j, -0.25),
                ],
            );
        }
    }
}

fn test_add_section42_chart_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42ChartMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = 0.5 * b0_near_origin / 30.0;

    for i in [5, 6] {
        test_add_section42_chart_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        test_add_section42_chart_vector_rows(
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

    test_add_section42_chart_vector_rows(
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
    test_add_section42_chart_vector_rows(
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

fn test_add_section42_chart_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: TestSection42ChartMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let coeff = test_section42_chart_strip(map, term.symbol).coeff(term.i, term.j);
        x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), term.weight));
        y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), term.weight));
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

fn test_section42_chart_strip(
    map: TestSection42ChartMap,
    symbol: Kp20Section42PaperSymbol,
) -> TestSection42ChartStrip {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

fn vector_constraints_as_scalar_rows(
    rows: &[Kp20PlanarTauLinearConstraint],
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut scalar_rows = Vec::with_capacity(2 * rows.len());
    for row in rows {
        for axis in [0_usize, 1] {
            scalar_rows.push(Kp20PlanarTauScalarLinearConstraint {
                terms: row
                    .terms
                    .iter()
                    .map(|&(coeff, weight)| (Kp20PlanarTauScalarCoeff::new(coeff, axis), weight))
                    .collect(),
                target: row.target[axis],
            });
        }
    }
    scalar_rows
}

#[derive(Debug, Clone, Copy)]
struct TestComplex {
    re: f64,
    im: f64,
}

impl TestComplex {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn from_angle(angle: f64) -> Self {
        Self {
            re: angle.cos(),
            im: angle.sin(),
        }
    }

    const fn conjugate(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

fn complex_add(a: TestComplex, b: TestComplex) -> TestComplex {
    TestComplex::new(a.re + b.re, a.im + b.im)
}

fn complex_mul(a: TestComplex, b: TestComplex) -> TestComplex {
    TestComplex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

fn complex_from_point(point: Point3) -> TestComplex {
    TestComplex::new(point.x, point.y)
}

fn maybe_conjugate_complex(value: TestComplex, conjugate: bool) -> TestComplex {
    if conjugate { value.conjugate() } else { value }
}

fn complex_similarity_scale_for_pairs(pairs: &[(TestComplex, TestComplex)]) -> Option<TestComplex> {
    let mut numerator = TestComplex::new(0.0, 0.0);
    let mut denominator = 0.0_f64;
    for &(source, target) in pairs {
        numerator = complex_add(numerator, complex_mul(target, source.conjugate()));
        denominator += source.norm2();
    }
    if denominator <= 1.0e-20 {
        None
    } else {
        Some(TestComplex::new(
            numerator.re / denominator,
            numerator.im / denominator,
        ))
    }
}

fn unit_sector_endpoint_normalized_tensor_rows_for_test(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
    swap_targets: bool,
    conjugate_source: bool,
) -> Option<Kp20RegularTensorBorderCNet> {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = TestComplex::from_angle(alpha * sector as f64);
    let spoke1 = TestComplex::from_angle(alpha * (sector + 1) as f64);
    let target0 = if swap_targets { spoke1 } else { spoke0 };
    let target1 = if swap_targets { spoke0 } else { spoke1 };
    let source0 = maybe_conjugate_complex(
        complex_from_point(rows.sectors[sector][0][0]),
        conjugate_source,
    );
    let source1 = maybe_conjugate_complex(
        complex_from_point(rows.sectors[sector][0][3]),
        conjugate_source,
    );
    let scale = complex_similarity_scale_for_pairs(&[(source0, target0), (source1, target1)])?;
    Some(complex_scaled_tensor_border_rows_for_test(
        rows,
        scale,
        conjugate_source,
    ))
}

fn complex_scaled_tensor_border_rows_for_test(
    rows: &Kp20RegularTensorBorderCNet,
    scale: TestComplex,
    conjugate_source: bool,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let value = complex_mul(
                        maybe_conjugate_complex(
                            complex_from_point(sector[row][col]),
                            conjugate_source,
                        ),
                        scale,
                    );
                    out[row][col] = Point3::new(value.re, value.im, sector[row][col].z);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn reorder_regular_tensor_border_depth_rows_for_test(
    rows: &Kp20RegularTensorBorderCNet,
    order: [usize; 3],
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector_rows| {
            [
                sector_rows[order[0]],
                sector_rows[order[1]],
                sector_rows[order[2]],
            ]
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn kp15_beta_split_then_apply_for_test(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP15 beta split-then-apply diagnostic needs at least 3 sectors, got {n}"
        ));
    }
    let params = Bi6Reparams::for_valence(n);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        let mut near_c = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let mut near_m = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP15 beta split-then-apply sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
            let (left, right) = split_cubic_half_for_test(sector[depth]);
            near_c[depth] = left;
            near_m[depth] = right;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: Kp20BoundaryHalf {
                rows: apply_beta_reparameterization(&near_c, &params),
            },
            near_m: Kp20BoundaryHalf {
                rows: apply_beta_reparameterization(&near_m, &params),
            },
        });
    }
    Ok(out)
}

fn kp15_beta_full_then_split_for_test(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP15 beta full-then-split diagnostic needs at least 3 sectors, got {n}"
        ));
    }
    let params = Bi6Reparams::for_valence(n);
    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        for depth in 0..3 {
            for point in sector[depth] {
                if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                    return Err(format!(
                        "KP15 beta full-then-split sector {sector_index} has non-finite tensor-border point"
                    ));
                }
            }
        }
        let full_rows = apply_beta_reparameterization(sector, &params);
        let mut near_c = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];
        let mut near_m = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];
        for depth in 0..3 {
            let (left, right) = split_degree6_half_for_test(full_rows[depth]);
            near_c[depth] = left;
            near_m[depth] = right;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: Kp20BoundaryHalf { rows: near_c },
            near_m: Kp20BoundaryHalf { rows: near_m },
        });
    }
    Ok(out)
}

fn split_degree6_half_for_test(row: [Point3; 7]) -> ([Point3; 7], [Point3; 7]) {
    let mut levels = Vec::with_capacity(7);
    levels.push(row.to_vec());
    for depth in 1..=6 {
        let prev = &levels[depth - 1];
        let mut next = Vec::with_capacity(7 - depth);
        for i in 0..(7 - depth) {
            next.push(midpoint_for_test(prev[i], prev[i + 1]));
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

fn tensor_border_split_chart_variant_rows(
    split: &[Kp20SplitTensorBorder],
    sector: usize,
    variant: TestTensorBorderChartVariant,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let previous = (sector + split.len() - 1) % split.len();
    let previous_rows = if variant.previous_uses_near_m {
        split[previous].near_m.rows
    } else {
        split[previous].near_c.rows
    };
    let current_rows = if variant.current_uses_near_m {
        split[sector].near_m.rows
    } else {
        split[sector].near_c.rows
    };
    tensor_border_chart_variant_rows_from_split_rows(previous_rows, current_rows, variant)
}

fn split_regular_sector_rows_with_second_order_factor(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    previous_uses_near_m: bool,
    current_uses_near_m: bool,
    second_order_factor: f64,
) -> ([[Point3; 7]; 3], [[Point3; 7]; 3]) {
    let previous = (sector + tensor_border.sectors.len() - 1) % tensor_border.sectors.len();
    let previous_rows = split_sector_half_rows_with_second_order_factor(
        tensor_border.sectors[previous],
        previous_uses_near_m,
        second_order_factor,
    );
    let current_rows = split_sector_half_rows_with_second_order_factor(
        tensor_border.sectors[sector],
        current_uses_near_m,
        second_order_factor,
    );
    (previous_rows, current_rows)
}

fn split_sector_half_rows_with_second_order_factor(
    sector_rows: [[Point3; 4]; 3],
    use_near_m: bool,
    second_order_factor: f64,
) -> [[Point3; 7]; 3] {
    let c = (core::f64::consts::TAU / 5.0).cos();
    let d1_3 = 2.0 * c / 3.0;
    let delta = d1_3 / 4.0 - c / 12.0;
    let (d0, d1) = (
        [0.0, 0.0, 0.0, delta],
        [delta, 2.0 * delta, 4.0 * delta, d1_3],
    );
    let mut half_rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for depth in 0..3 {
        let (near_c, near_m) = split_cubic_half_for_test(sector_rows[depth]);
        half_rows[depth] = if use_near_m { near_m } else { near_c };
    }
    reparameterise_half_rows_with_second_order_factor(
        half_rows,
        if use_near_m { d1 } else { d0 },
        second_order_factor,
    )
}

fn split_cubic_half_for_test(row: [Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let p10 = midpoint_for_test(row[0], row[1]);
    let p11 = midpoint_for_test(row[1], row[2]);
    let p12 = midpoint_for_test(row[2], row[3]);
    let p20 = midpoint_for_test(p10, p11);
    let p21 = midpoint_for_test(p11, p12);
    let p30 = midpoint_for_test(p20, p21);
    ([row[0], p10, p20, p30], [p30, p21, p12, row[3]])
}

fn midpoint_for_test(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

fn reparameterise_half_rows_with_second_order_factor(
    rows: [[Point3; 4]; 3],
    d: [f64; 4],
    second_order_factor: f64,
) -> [[Point3; 7]; 3] {
    let source_row0 = elevate_cubic_to_degree6_for_test(rows[0]);
    let source_row1 = elevate_cubic_to_degree6_for_test(rows[1]);
    let source_row2 = elevate_cubic_to_degree6_for_test(rows[2]);
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut row2 = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0..7 {
        row1[i] = midpoint_for_test(source_row0[i], source_row1[i]);
        row2[i] = Point3::new(
            0.2 * source_row0[i].x + 0.6 * source_row1[i].x + 0.2 * source_row2[i].x,
            0.2 * source_row0[i].y + 0.6 * source_row1[i].y + 0.2 * source_row2[i].y,
            0.0,
        );
    }
    let diff = [
        sub_point(rows[1][0], rows[0][0]),
        sub_point(rows[1][1], rows[0][1]),
        sub_point(rows[1][2], rows[0][2]),
        sub_point(rows[1][3], rows[0][3]),
    ];
    let d_times_diff = multiply_scalar_cubic_by_point_cubic_for_test(d, diff);
    for i in 0..7 {
        row2[i] = add_point(row2[i], scale_point(d_times_diff[i], second_order_factor));
    }
    [source_row0, row1, row2]
}

fn elevate_cubic_to_degree6_for_test(row: [Point3; 4]) -> [Point3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut out = [Point3::new(0.0, 0.0, 0.0); 7];
    for k in 0_usize..=6 {
        let lo = k.saturating_sub(3);
        let hi = k.min(3);
        for i in lo..=hi {
            let weight = C3[i] * binomial3(k - i) / C6[k];
            out[k] = add_point(out[k], scale_point(row[i], weight));
        }
    }
    out
}

fn multiply_scalar_cubic_by_point_cubic_for_test(
    scalar: [f64; 4],
    points: [Point3; 4],
) -> [Point3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut out = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0..=3 {
        for j in 0..=3 {
            let k = i + j;
            let weight = scalar[i] * C3[i] * C3[j] / C6[k];
            out[k] = add_point(out[k], scale_point(points[j], weight));
        }
    }
    out
}

const fn binomial3(i: usize) -> f64 {
    match i {
        0 | 3 => 1.0,
        1 | 2 => 3.0,
        _ => 0.0,
    }
}

#[path = "kp20_bi6_tau_solve_assemble_tests_part1.rs"]
mod part1;
#[path = "kp20_bi6_tau_solve_assemble_tests_part2.rs"]
mod part2;
#[path = "kp20_bi6_tau_solve_assemble_tests_part3.rs"]
mod part3;
#[path = "kp20_bi6_tau_solve_assemble_tests_part4.rs"]
mod part4;
