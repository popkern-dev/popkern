// Paper-notation test module: loop indices, helper arities, and tuple shapes
// mirror KP15/KP20/PR08 subscripts and equation forms; keep the index form.
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::type_complexity
)]

use super::super::kp20_bi6_assemble::{
    Kp20Bi6MacroCapInput, Kp20Bi6MacroCapOutput,
    assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source,
};
use super::super::kp20_bi6_cap::{
    Kp20BoundaryHalf, Kp20SplitTensorBorder, split_kp20_tensor_border,
    split_kp20_tensor_border_with_d1_3,
};
use super::super::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, kp20_cap_net_figure_tensor_border_rows,
    kp20_tensor_border_from_regular_cnet_rows,
};
use super::super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedF3SolveOutput,
    Kp20PlanarTauAugmentedScalarCoeff, Kp20PlanarTauAugmentedScalarLinearConstraint,
    kp20_planar_tau_constraint_space_report,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
};
use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::kp20_bi6_cc_ring::{
    Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
    kp20_pr08_characteristic_tensor_border_candidate_rows,
};
use crate::geometry::gregory::kp20_bi6_guide_build::{Kp20GuideBuildInput, build_kp20_guide_ring};
use crate::geometry::gregory::kp20_bi6_local_constraint::{
    Kp20Eq9DiagnosticTerms, Kp20Eq9TermSigns, Kp20Eq9TermSource, kp20_eq8_g1_residual,
    kp20_eq8_jacobian_dense, kp20_eq9_g2_residual, kp20_eq9_jacobian_dense,
    kp20_eq9_jacobian_dense_with_diagnostic_terms, measure_kp20_assembly_constraints,
};
use crate::geometry::gregory::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use crate::geometry::gregory::{
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

#[test]
fn rotational_endpoint_aware_tensor_depth_policy_after_near_origin_fix_is_ranked() {
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
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
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
    let cases: [(&str, &[usize]); 7] = [
        ("depths-0-1-2", &[0, 1, 2]),
        ("depths-0-2-eq13-owns-1", &[0, 2]),
        ("depths-0-1", &[0, 1]),
        ("depths-1-2", &[1, 2]),
        ("depth-0-only", &[0]),
        ("depth-1-only", &[1]),
        ("depth-2-only", &[2]),
    ];
    let mut best = ("none", usize::MAX, usize::MAX, f64::INFINITY, 0.0_f64);
    for (label, kept_depths) in cases {
        let tensor_rows = keep_tensor_chart_depths(&source_rows, source_variant, kept_depths);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("tensor depth policy after near-origin fix solve succeeds");
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
            "tensor depth policy after near-origin fix {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
        "best tensor depth policy after near-origin fix: label={} neg={} outside={} residual={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "tensor depth policy after near-origin fix must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_row_space_seed_regularization_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const WEIGHTS: [f64; 9] = [
        0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0,
    ];
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
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
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
    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for weight in WEIGHTS {
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: weight,
            })
            .expect("corrected row-space seed-regularization solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "corrected row-space seed regularization weight={weight:e}: residual={residual:e} f3={f3:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            f3 = kp20_planar_bi6_tau_f3_energy(&local_tau),
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        let row = (
            weight,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected row-space seed regularization: weight={:e} neg={} outside={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(
        best.3.is_finite(),
        "seed-regularization ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_scale_after_near_origin_fix_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 7] = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5];
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
    let mut best = (
        1.0_f64,
        1.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("source cross-depth scale solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            if validity.orientation.negative_determinants <= 18
                || validity.domain.outside_triangle_samples <= 5
            {
                eprintln!(
                    "source cross-depth scale left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best source cross-depth scale after near-origin fix: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "source cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_scale_window_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 8] = [0.125, 0.175, 0.225, 0.25, 0.275, 0.325, 0.375, 0.5];
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
    let mut best = (
        1.0_f64,
        1.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("source cross-depth scale window solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            if validity.orientation.negative_determinants <= 17
                || validity.domain.outside_triangle_samples <= 1
            {
                eprintln!(
                    "source cross-depth scale window left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best source cross-depth scale window: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "source cross-depth scale window ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_source_cross_depth_spectral_scale_is_measured() {
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
    let scale_cases = [
        ("lambda", ring.lambda),
        ("one-minus-lambda", 1.0 - ring.lambda),
        ("half-one-minus-lambda", 0.5 * (1.0 - ring.lambda)),
        ("lambda-squared", ring.lambda * ring.lambda),
        ("quarter", 0.25),
        ("empirical-window-best", 0.225),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, scale) in scale_cases {
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &sector_rotated_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            scale,
        );
        let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            scale,
        );
        let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral source cross-depth scale solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "source cross-depth spectral scale {label}: scale={scale:e} lambda={lambda:e} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            lambda = ring.lambda,
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
            validity.orientation.min_determinant,
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
        "best source cross-depth spectral scale: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        (0.5 * (1.0 - ring.lambda) - 0.225).abs() < 1.0e-5,
        "Catmull-Clark spectral complement should explain the previous empirical cross-depth scale numerically"
    );
    assert!(
        best.3.is_finite(),
        "spectral source cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_spectral_source_depth1_projection_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let depth1_cases = [
        ("omit-depth1", None),
        ("same-spectral-depth1", Some(spectral_cross_depth_scale)),
        (
            "half-spectral-depth1",
            Some(0.5 * spectral_cross_depth_scale),
        ),
        (
            "sqrt-spectral-depth1",
            Some(spectral_cross_depth_scale.sqrt()),
        ),
        ("one-minus-lambda-depth1", Some(1.0 - ring.lambda)),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, depth1_scale) in depth1_cases {
        let mut scaled_rows = sector_rotated_rows.clone();
        if let Some(scale) = depth1_scale {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &scaled_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                1,
                scale,
            );
            scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                1,
                scale,
            );
        }
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &scaled_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            spectral_cross_depth_scale,
        );
        scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            spectral_cross_depth_scale,
        );
        let projected_tensor_rows = if depth1_scale.is_some() {
            scaled_rows.clone()
        } else {
            without_tensor_chart_depth(&scaled_rows, source_variant, 1)
        };
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral depth1 source projection solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "spectral source depth1 projection {label}: depth1_scale={depth1:e} depth2_scale={depth2:e} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            depth1 = depth1_scale.unwrap_or(f64::NAN),
            depth2 = spectral_cross_depth_scale,
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
            validity.orientation.min_determinant,
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
        "best spectral source depth1 projection: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite(),
        "spectral source depth1 projection ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_spectral_source_diagonal_symmetry_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let diagonal_rows = kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
        .expect("sector-local diagonal rows");
    let cases = [
        ("fig8-only", true, false),
        ("diagonal-only", false, true),
        ("fig8-and-diagonal", true, true),
    ];
    for (label, include_fig8, include_diagonal) in cases {
        let mut vector_constraints = Vec::new();
        if include_fig8 {
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        }
        if include_diagonal {
            let mut scalar_constraints_with_diagonal = scalar_constraints.clone();
            scalar_constraints_with_diagonal.extend(diagonal_rows.clone());
            vector_constraints.extend(projected_tensor_rows.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints_with_diagonal,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("spectral source diagonal-symmetry solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "spectral source diagonal symmetry {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                full_tensor = full_tensor_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            assert!(
                solved.max_constraint_residual.is_finite(),
                "spectral source diagonal case {label} must be finite"
            );
            continue;
        }
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("spectral source diagonal-symmetry solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "spectral source diagonal symmetry {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "spectral source diagonal case {label} must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_beta_source_under_corrected_chart_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let previous = valence - 1;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let beta_sources = [
        (
            "kp15-beta-split-then-apply",
            kp15_beta_split_then_apply_for_test(&tensor_border).expect("split beta source"),
        ),
        (
            "kp15-beta-full-then-split",
            kp15_beta_full_then_split_for_test(&tensor_border).expect("full beta source"),
        ),
    ];
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let previous_piece = Kp20PieceId::NearMidpointRight;
    let current_piece = Kp20PieceId::NearMidpointLeft;
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
    let mut best = (
        "none",
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (source_label, split) in beta_sources {
        let previous_rows = split[previous].near_m.rows;
        let current_rows = split[sector].near_c.rows;
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            previous_rows,
            current_rows,
            source_variant,
        );
        let rotated_rows = rotate_tensor_chart_piece_targets(&base_rows, previous_piece, -angle);
        for (scale_label, source_rows) in [
            ("raw-depth", rotated_rows.clone()),
            ("spectral-depth2", {
                let previous_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                    &rotated_rows,
                    source_variant,
                    previous_piece,
                    2,
                    spectral_cross_depth_scale,
                );
                scale_tensor_chart_depth_from_boundary_for_piece(
                    &previous_scaled,
                    source_variant,
                    current_piece,
                    2,
                    spectral_cross_depth_scale,
                )
            }),
        ] {
            let projected_tensor_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("corrected-chart beta source solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = source_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "corrected-chart beta source {source_label} {scale_label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
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
                source_label,
                scale_label,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.3 < best.3)
                || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best corrected-chart beta source: source={} scale={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "corrected-chart beta source ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_row_space_is_measured() {
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows.clone());
    let report = kp20_planar_tau_constraint_space_report(&vector_constraints, &scalar_constraints)
        .expect("corrected spectral row-space report");
    eprintln!(
        "corrected spectral source row-space: vars={} dense_rows={} independent_rows={} nullity={}",
        report.scalar_variable_count,
        report.dense_row_count,
        report.independent_row_count,
        report.nullity,
    );

    let mut vector_with_internal_c2 = kp20_planar_tau_internal_c2_constraints();
    vector_with_internal_c2.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_with_internal_c2.extend(projected_tensor_rows);
    let internal_c2_report =
        kp20_planar_tau_constraint_space_report(&vector_with_internal_c2, &scalar_constraints)
            .expect("corrected spectral internal-C2 row-space report");
    eprintln!(
        "corrected spectral source plus internal-C2 row-space: vars={} dense_rows={} independent_rows={} nullity={}",
        internal_c2_report.scalar_variable_count,
        internal_c2_report.dense_row_count,
        internal_c2_report.independent_row_count,
        internal_c2_report.nullity,
    );
    assert_eq!(report.scalar_variable_count, 392);
    assert!(
        report.nullity > 15,
        "corrected spectral source rows should still expose the missing sigma finite-basis reduction"
    );
    assert!(
        internal_c2_report.nullity <= report.nullity,
        "adding internal C2 should not increase the corrected spectral row-space nullity"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_fig8_average_variants_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (label, include_diagonal_mirror, include_rotated_quadrants) in [
        ("baseline", false, false),
        ("diagonal-mirror", true, false),
        ("rotated-quadrants", false, true),
        ("diagonal-mirror+rotated-quadrants", true, true),
    ] {
        let fig8_rows = test_fig8_indicated_average_variant_rows(
            include_diagonal_mirror,
            include_rotated_quadrants,
        );
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(fig8_rows);
        vector_constraints.extend(projected_tensor_rows.clone());
        let report =
            kp20_planar_tau_constraint_space_report(&vector_constraints, &scalar_constraints)
                .expect("Fig. 8 average variant row-space report");
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("Fig. 8 average variant solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "corrected spectral Fig8 average variant {label}: rows={rows} independent_rows={independent} nullity={nullity} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            rows = report.dense_row_count,
            independent = report.independent_row_count,
            nullity = report.nullity,
            residual = solved.max_constraint_residual,
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
            report.nullity,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
            || (row.1 == best.1 && row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected spectral Fig8 average variant: label={} neg={} outside={} nullity={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.4.is_finite(),
        "Fig. 8 average variant ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_internal_c2_solve_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    for (label, include_internal_c2) in [("without-internal-c2", false), ("with-internal-c2", true)]
    {
        let mut vector_constraints = Vec::new();
        if include_internal_c2 {
            vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        }
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("corrected spectral internal-C2 solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "corrected spectral source {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "corrected spectral internal-C2 case {label} must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_internal_c2_source_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let (projected_tensor_rows, projection_report) = kp20_project_rotational_paper_tensor_band(
        projected_source_input.clone(),
        &fixed_vector_constraints,
        &scalar_constraints,
    )
    .expect("corrected spectral internal-C2 source projection succeeds");
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_tensor_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("projected corrected spectral internal-C2 solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "corrected spectral internal-C2 source projection: projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        source_delta = projection_report.max_source_delta,
        adjusted = projection_report.adjusted_scalar_targets,
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    let mut deltas = projected_source_input
        .iter()
        .zip(&projected_tensor_rows)
        .map(|(source, projected)| {
            let coeff = source.terms[0].0;
            let dx = projected.target[0] - source.target[0];
            let dy = projected.target[1] - source.target[1];
            (dx.hypot(dy), coeff, dx, dy, source.target, projected.target)
        })
        .collect::<Vec<_>>();
    deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (rank, (norm, coeff, dx, dy, source, projected)) in deltas.iter().take(8).enumerate() {
        eprintln!(
            "  corrected spectral projection delta rank={rank} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} source=({sx:e},{sy:e}) projected=({px:e},{py:e})",
            piece = coeff.piece,
            i = coeff.i,
            j = coeff.j,
            sx = source[0],
            sy = source[1],
            px = projected[0],
            py = projected[1],
        );
    }
    assert!(
        solved.max_constraint_residual.is_finite()
            && projection_report.max_row_residual.is_finite(),
        "corrected spectral internal-C2 source projection must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_depth2_only_projection_is_measured() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let free_depth2_coeffs = projected_source_input
        .iter()
        .map(|row| row.terms[0].0)
        .filter(|&coeff| tensor_chart_coeff_depth(coeff, source_variant) == 2)
        .collect::<Vec<_>>();
    let (projected_tensor_rows, projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            projected_source_input.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &free_depth2_coeffs,
        )
        .expect("depth2-only corrected spectral projection succeeds");
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_tensor_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("depth2-only corrected spectral solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "corrected spectral internal-C2 depth2-only projection: projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        source_delta = projection_report.max_source_delta,
        adjusted = projection_report.adjusted_scalar_targets,
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite()
            && projection_report.max_row_residual.is_finite(),
        "depth2-only corrected spectral projection must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_corrected_spectral_depth2_piece_projection_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_source_input = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let mut best = ("none", f64::INFINITY, usize::MAX, usize::MAX, f64::INFINITY);
    for (label, free_piece) in [
        ("underlined-p-depth2", Kp20PieceId::NearMidpointLeft),
        ("arrow-underlined-p-depth2", Kp20PieceId::NearMidpointRight),
    ] {
        let free_coeffs = projected_source_input
            .iter()
            .map(|row| row.terms[0].0)
            .filter(|&coeff| {
                coeff.piece == free_piece && tensor_chart_coeff_depth(coeff, source_variant) == 2
            })
            .collect::<Vec<_>>();
        let (projected_tensor_rows, projection_report) =
            kp20_project_section42_source_tensor_rows_with_free_coeffs(
                projected_source_input.clone(),
                &fixed_vector_constraints,
                &scalar_constraints,
                &free_coeffs,
            )
            .expect("piece depth2 corrected spectral projection succeeds");
        let mut vector_constraints = fixed_vector_constraints.clone();
        vector_constraints.extend(projected_tensor_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("piece depth2 corrected spectral solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        eprintln!(
            "corrected spectral internal-C2 piece depth2 projection {label}: free={} projection_residual={projection_residual:e} source_delta={source_delta:e} adjusted={adjusted} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            free_coeffs.len(),
            projection_residual = projection_report.max_row_residual,
            source_delta = projection_report.max_source_delta,
            adjusted = projection_report.adjusted_scalar_targets,
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
            solved.max_constraint_residual,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            projection_report.max_source_delta,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best corrected spectral internal-C2 piece depth2 projection: label={} residual={:e} neg={} outside={} source_delta={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(best.1.is_finite(), "piece depth2 ranking must be finite");
}

#[test]
fn rotational_endpoint_aware_pr08_sources_without_projection_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: true,
        transpose_indices: false,
        swap_pieces: true,
    };
    let mut scalar_constraints = test_section42_endpoint_aware_chart_rows(
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
    scalar_constraints.extend(
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
            .expect("sector-local diagonal rows"),
    );

    let mut best = (
        String::from("none"),
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    let mut valid_candidates = 0_usize;
    for candidate in Kp20Pr08TensorBorderCandidate::all() {
        let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 candidate rows");
        let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
            &raw_rows, valence, sector, true, true,
        )
        .expect("candidate endpoint normalization");
        let tensor_border =
            kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
        let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
        let sector_rotated_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
            &sector_rotated_rows,
            source_variant,
            Kp20PieceId::NearMidpointLeft,
            2,
            spectral_cross_depth_scale,
        );
        let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
            &left_scaled,
            source_variant,
            Kp20PieceId::NearMidpointRight,
            2,
            spectral_cross_depth_scale,
        );

        for (depth_label, tensor_rows) in [
            (
                "depth0+2",
                without_tensor_chart_depth(&scaled_rows, source_variant, 1),
            ),
            ("depth0+1+2", scaled_rows.clone()),
        ] {
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(tensor_rows.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("unprojected PR08 source solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let is_valid = validity.is_strict_sector_bounded_reparameterization();
            if is_valid {
                valid_candidates += 1;
            }
            eprintln!(
                "unprojected PR08 source {candidate:?} {depth_label}: valid={is_valid} residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
                format!("{candidate:?}/{depth_label}"),
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.1 < best.1
                || (row.1 == best.1 && row.2 < best.2)
                || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
                || (row.1 == best.1 && row.2 == best.2 && row.3 == best.3 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best unprojected PR08 source: label={} valid_candidates={} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, valid_candidates, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite(),
        "unprojected PR08 source ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_promising_pr08_charts_are_ranked_with_live_and_diagnostic_diagonal() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let base_scalar_constraints = test_section42_endpoint_aware_chart_rows(
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
    let cases = [
        (
            "best-validity-before-diagonal",
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: true,
                previous_reverse_u: true,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: true,
            },
        ),
        (
            "best-domain-before-diagonal",
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: true,
            },
        ),
        (
            "scaled-top-domain-before-diagonal",
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows,
            TestTensorBorderChartVariant {
                previous_uses_near_m: true,
                current_uses_near_m: false,
                previous_reverse_u: false,
                current_reverse_u: true,
                transpose_indices: false,
                swap_pieces: false,
            },
        ),
    ];

    for include_diagnostic_diagonal in [false, true] {
        let mut scalar_constraints = base_scalar_constraints.clone();
        if include_diagnostic_diagonal {
            scalar_constraints.extend(
                kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)
                    .expect("sector-local diagonal rows"),
            );
        }
        let mode_label = if include_diagnostic_diagonal {
            "diagnostic-diagonal"
        } else {
            "live"
        };
        let mut best = (
            "none",
            usize::MAX,
            usize::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for (label, candidate, source_variant) in cases {
            let raw_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
                .expect("PR08 candidate rows");
            let normalized_rows = unit_sector_endpoint_normalized_tensor_rows_for_test(
                &raw_rows, valence, sector, true, true,
            )
            .expect("candidate endpoint normalization");
            let tensor_border =
                kp20_tensor_border_from_regular_cnet_rows(&normalized_rows).expect("tensor rows");
            let base_rows =
                tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
            let sector_rotated_rows = rotate_tensor_chart_piece_targets(
                &base_rows,
                Kp20PieceId::NearMidpointRight,
                -angle,
            );
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                spectral_cross_depth_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                spectral_cross_depth_scale,
            );
            for (depth_label, tensor_rows) in [
                (
                    "depth0+2",
                    without_tensor_chart_depth(&scaled_rows, source_variant, 1),
                ),
                ("depth0+1+2", scaled_rows.clone()),
            ] {
                let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
                vector_constraints.extend(tensor_rows);
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: planar_tau_seed_macro_patch(valence, sector),
                        vector_constraints,
                        scalar_constraints: scalar_constraints.clone(),
                        regularization_weight: 1.0e-10,
                    },
                )
                .expect("promising PR08 chart solve succeeds");
                let local_tau =
                    global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                let validity =
                    kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
                eprintln!(
                    "promising PR08 chart mode={mode_label} {label} candidate={} chart={} {depth_label}: valid={} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    candidate.name(),
                    source_variant.label(),
                    validity.is_strict_sector_bounded_reparameterization(),
                    residual = solved.max_constraint_residual,
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
                    validity.orientation.min_determinant,
                );
                if row.1 < best.1
                    || (row.1 == best.1 && row.2 < best.2)
                    || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
                {
                    best = row;
                }
            }
        }
        eprintln!(
            "best promising PR08 chart mode={mode_label}: label={} neg={} outside={} residual={:e} min_det={:e}",
            best.0, best.1, best.2, best.3, best.4,
        );
        assert!(
            best.3.is_finite(),
            "promising PR08 chart ranking must be finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_production_pr08_inner_q1_source_matches_live_domain_baseline() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let ring = kp20_catmull_clark_characteristic_ring(valence).expect("PR08 ring");
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(without_tensor_chart_depth(&scaled_rows, source_variant, 1));
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(
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
            ),
            regularization_weight: 1.0e-10,
        })
        .expect("production PR08 inner-Q1 source solve succeeds");
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    eprintln!(
        "production PR08 inner-Q1 live chart: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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

    assert_eq!(validity.domain.outside_triangle_samples, 0);
    assert_eq!(validity.orientation.negative_determinants, 37);
    assert!(solved.max_constraint_residual < 7.0e-2);
    assert!(full_tensor_residual < 1.0e-1);
}

#[test]
fn rotational_c2_edge_chart_production_pr08_inner_q1_source_is_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.125,
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(source_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge production PR08 inner-Q1 source solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    let piece_reports = test_piece_domain_orientation_reports(&local_solved, &SAMPLES, 1.0e-10);
    eprintln!(
        "C2-edge production PR08 inner-Q1 source: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    for report in piece_reports {
        eprintln!(
            "C2-edge production piece {:?}: domain outside={}/{} orientation neg={}/{} zero={} det=[{:.6},{:.6}] worst_det={}",
            report.piece,
            report.outside,
            report.domain_samples,
            report.negative,
            report.orientation_samples,
            report.near_zero,
            report.min_det,
            report.max_det,
            report.worst_orientation_sample,
        );
    }

    assert_eq!(validity.domain.outside_triangle_samples, 0);
    assert_eq!(validity.orientation.negative_determinants, 20);
    assert!(solved.max_constraint_residual < 5.5e-2);
    assert!(full_tensor_residual < 9.0e-2);
}

#[test]
fn rotational_c2_edge_eq14_reparameterized_pr08_source_factor_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const FACTORS: [f64; 13] = [
        -1.0, -0.5, -0.25, -0.2, -0.125, 0.0, 0.125, 0.2, 0.25, 0.5, 0.75, 1.0, 1.25,
    ];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border =
        kp20_pr08_stam_inner_q1_characteristic_tensor_border(valence).expect("PR08/Stam source");
    let source_variant = TestTensorBorderChartVariant {
        previous_uses_near_m: true,
        current_uses_near_m: false,
        previous_reverse_u: false,
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);

    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for second_order_factor in FACTORS {
        let (previous_rows, current_rows) = split_regular_sector_rows_with_second_order_factor(
            &tensor_border,
            sector,
            source_variant.previous_uses_near_m,
            source_variant.current_uses_near_m,
            second_order_factor,
        );
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            previous_rows,
            current_rows,
            source_variant,
        );
        let source_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let projected_source_rows = without_tensor_chart_depth(&source_rows, source_variant, 1);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(projected_source_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge eq14 reparameterized PR08 source solve succeeds");
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "C2-edge eq14 reparameterized PR08 source factor={second_order_factor:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            second_order_factor,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
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
        "best C2-edge eq14 reparameterized PR08 source factor={:e}: neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge eq14 reparameterized source ranking must remain finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_cross_depth_scale_window_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const SCALES: [f64; 12] = [
        -0.25,
        -0.125,
        0.0,
        0.125,
        0.175,
        0.2250058227408514,
        0.275,
        0.325,
        0.375,
        0.5,
        0.75,
        1.0,
    ];
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

    let mut best = (
        0.0_f64,
        0.0_f64,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in SCALES {
        for right_scale in SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("live cross-depth scale solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.4 < best.4)
                || (row.2 == best.2 && row.4 == best.4 && row.5 < best.5)
            {
                best = row;
            }
            if validity.orientation.negative_determinants <= 30
                || validity.domain.outside_triangle_samples == 0
            {
                eprintln!(
                    "live cross-depth scale left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
    eprintln!(
        "best live cross-depth scale: left={:e} right={:e} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7,
    );
    assert!(
        best.5.is_finite(),
        "live cross-depth scale ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_cross_depth_signed_refinement_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const LEFT_SCALES: [f64; 11] = [
        -0.35, -0.3, -0.25, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, 0.0,
    ];
    const RIGHT_SCALES: [f64; 9] = [
        0.125,
        0.175,
        0.2250058227408514,
        0.275,
        0.325,
        0.375,
        0.45,
        0.6,
        0.75,
    ];
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
    let mut best = (
        0.0_f64,
        0.0_f64,
        usize::MAX,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for left_scale in LEFT_SCALES {
        for right_scale in RIGHT_SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &sector_rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(projected_tensor_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("live signed cross-depth refinement solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.orientation.near_zero_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.2 < best.2
                || (row.2 == best.2 && row.4 < best.4)
                || (row.2 == best.2 && row.4 == best.4 && row.5 < best.5)
            {
                best = row;
            }
            if validity.orientation.negative_determinants <= 26 {
                eprintln!(
                    "live signed cross-depth left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
    eprintln!(
        "best live signed cross-depth refinement: left={:e} right={:e} neg={} zero={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6, best.7,
    );
    assert!(
        best.5.is_finite(),
        "live signed cross-depth refinement must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_cross_depth_best_is_localized() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(
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
            ),
            regularization_weight: 1.0e-10,
        })
        .expect("live signed cross-depth best localization solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "live signed cross-depth best: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
        "live signed cross-depth best localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_source_with_first_order_frame_is_measured() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
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
    let axes = sector_local_axis_components(valence, sector);
    for (label, include_right) in [("left-origin-frame", false), ("both-origin-frames", true)] {
        let mut augmented_scalar_constraints =
            first_order_centered_source_rows(Kp20PieceId::NearOriginLeft, axes, 0);
        if include_right {
            augmented_scalar_constraints.extend(first_order_centered_source_rows(
                Kp20PieceId::NearOriginRight,
                axes,
                0,
            ));
        }
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints: vector_constraints.clone(),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("live signed source first-order frame solve succeeds");
        let full_tensor_residual = scaled_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        let right = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginRight);
        eprintln!(
            "live signed source first-order {label}: residual={residual:e} scale={scale:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] left_det={left_det:e} right_det={right_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            full_tensor = full_tensor_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
            right_det = right.determinant,
        );
        assert!(
            solved.max_constraint_residual.is_finite(),
            "first-order frame solve must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_live_signed_source_first_order_frames_are_ranked() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(projected_tensor_rows);
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
    let axes = sector_local_axis_components(valence, sector);
    let frame_cases = [
        ("identity", [1.0, 0.0], [0.0, 1.0]),
        ("both-negated", [-1.0, 0.0], [0.0, -1.0]),
        ("swap", [0.0, 1.0], [1.0, 0.0]),
        ("swap-negated", [0.0, -1.0], [-1.0, 0.0]),
        ("fig8-left", [0.5, 0.5], [0.0, 1.0]),
        ("fig8-right-positive", [1.0, 0.0], [0.5, 0.5]),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0_f64,
    );
    for (label, du, dv) in frame_cases {
        let augmented_scalar_constraints = first_order_source_rows_for_local_frame(
            Kp20PieceId::NearOriginLeft,
            axes,
            [0.0, 0.0],
            du,
            dv,
            0,
        );
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints: vector_constraints.clone(),
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("live signed source first-order frame variant solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        let left = near_origin_first_order_frame(&local_tau, Kp20PieceId::NearOriginLeft);
        eprintln!(
            "live signed source first-order frame {label}: residual={residual:e} scale={scale:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] left_det={left_det:e}",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
            left_det = left.determinant,
        );
        let row = (
            label,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            solved.extra_values[0],
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best live signed source first-order frame: label={} neg={} outside={} residual={:e} min_det={:e} scale={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite(),
        "first-order frame variant ranking must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_signed_source_centered_projection_is_measured() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
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
    let axes = sector_local_axis_components(valence, sector);
    for (label, du, dv) in [
        ("identity", [1.0, 0.0], [0.0, 1.0]),
        ("fig8-left", [0.5, 0.5], [0.0, 1.0]),
    ] {
        let centered_rows = first_order_source_rows_for_local_frame(
            Kp20PieceId::NearOriginLeft,
            axes,
            [0.0, 0.0],
            du,
            dv,
            0,
        );
        let projected_rows = project_tensor_band_with_centered_first_order(
            source_rows.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &centered_rows,
        );
        let max_source_delta = source_rows
            .iter()
            .zip(&projected_rows)
            .map(|(source, projected)| {
                (projected.target[0] - source.target[0])
                    .hypot(projected.target[1] - source.target[1])
            })
            .fold(0.0_f64, f64::max);
        let mut deltas = source_rows
            .iter()
            .zip(&projected_rows)
            .map(|(source, projected)| {
                let coeff = source.terms[0].0;
                let dx = projected.target[0] - source.target[0];
                let dy = projected.target[1] - source.target[1];
                (dx.hypot(dy), coeff, dx, dy, source.target, projected.target)
            })
            .collect::<Vec<_>>();
        deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let mut vector_constraints = fixed_vector_constraints.clone();
        vector_constraints.extend(projected_rows);
        let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
            &Kp20PlanarTauAugmentedF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                augmented_scalar_constraints: centered_rows,
                extra_variable_count: 1,
                regularization_weight: 1.0e-10,
                extra_regularization_weight: 1.0e-10,
            },
        )
        .expect("centered projected source solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "live signed centered source projection {label}: residual={residual:e} scale={scale:e} max_source_delta={max_delta:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            scale = solved.extra_values[0],
            max_delta = max_source_delta,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        for (rank, (norm, coeff, dx, dy, source, projected)) in deltas.iter().take(8).enumerate() {
            eprintln!(
                "  centered projection {label} delta rank={rank} coeff={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e} source=({sx:e},{sy:e}) projected=({px:e},{py:e})",
                piece = coeff.piece,
                i = coeff.i,
                j = coeff.j,
                sx = source[0],
                sy = source[1],
                px = projected[0],
                py = projected[1],
            );
        }
        assert!(
            solved.max_constraint_residual.is_finite() && max_source_delta.is_finite(),
            "centered source projection must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_live_signed_source_appendix_a_refined_pins_are_ranked() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut best = (
        "none",
        0_usize,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, tau_family) in [
        ("paper-eq22", Kp20TauFamily::Paper),
        ("appendix-a-refined", Kp20TauFamily::PaperRefinedApex),
    ] {
        let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
            &local_identity_tri6_guide(),
            valence,
            sector,
            tau_family,
            kp20_paper_quadrant_sector_piece_maps(),
        )
        .expect("Appendix A tau samples with identity guide");
        let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
        for (source_label, include_source_rows) in
            [("without-source", false), ("with-source", true)]
        {
            for max_total_degree in [1_usize, 2, 3] {
                let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
                vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
                if include_source_rows {
                    vector_constraints.extend(source_rows.clone());
                }
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
                let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                    &Kp20PlanarTauF3ScalarSolveInput {
                        initial: planar_tau_seed_macro_patch(valence, sector),
                        vector_constraints,
                        scalar_constraints: scalar_constraints.clone(),
                        regularization_weight: 1.0e-10,
                    },
                )
                .expect("live signed source Appendix A central pin solve succeeds");
                let full_tensor_residual = scaled_rows
                    .iter()
                    .map(|row| {
                        let residual = vector_residual(&solved.tau_macro_patch, row);
                        residual[0].abs().max(residual[1].abs())
                    })
                    .fold(0.0_f64, f64::max);
                let local_solved =
                    global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
                let validity =
                    kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
                eprintln!(
                    "live signed Appendix-A pins {label} {source_label} degree={degree}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                    degree = max_total_degree,
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
                    max_total_degree,
                    validity.orientation.negative_determinants,
                    validity.domain.outside_triangle_samples,
                    solved.max_constraint_residual,
                    validity.orientation.min_determinant,
                    full_tensor_residual,
                );
                if include_source_rows
                    && (row.2 < best.2
                        || (row.2 == best.2 && row.3 < best.3)
                        || (row.2 == best.2 && row.3 == best.3 && row.4 < best.4))
                {
                    best = row;
                }
            }
        }
    }
    eprintln!(
        "best live signed Appendix-A pins: label={} degree={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite() && best.6.is_finite(),
        "Appendix A central pin ranking must remain finite"
    );
}

#[test]
fn rotational_appendix_a_refined_central_pins_row_stack_components_are_localized() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A refined tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mut central_pin_rows = Vec::new();
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        for i in 0..=2 {
            for j in 0..=(2 - i) {
                let target = paper_tau.piece(piece_id).cps[i][j];
                central_pin_rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece_id, i, j),
                    [target.x, target.y],
                ));
            }
        }
    }
    let scalar_rows = test_section42_endpoint_aware_chart_rows(
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
    let cases = [
        ("pins-only", false, false, false),
        ("pins-internal-c2", true, false, false),
        ("pins-fig8", false, true, false),
        ("pins-section42-scalar", false, false, true),
        ("pins-internal-c2-fig8", true, true, false),
        ("pins-internal-c2-section42", true, false, true),
        ("pins-fig8-section42", false, true, true),
        ("pins-full-live-no-source", true, true, true),
    ];
    for (label, include_internal_c2, include_fig8, include_scalar) in cases {
        let mut vector_constraints = central_pin_rows.clone();
        if include_internal_c2 {
            vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        }
        if include_fig8 {
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        }
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: if include_scalar {
                    scalar_rows.clone()
                } else {
                    Vec::new()
                },
                regularization_weight: 1.0e-10,
            })
            .expect("Appendix A refined row component solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "Appendix-A refined central pins row component {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            "Appendix A refined row component residual must be finite"
        );
    }

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
    let mut grouped_rows = Vec::new();
    let mut c0_rows = Vec::new();
    test_add_section42_chart_c0_rows(&mut c0_rows, map);
    grouped_rows.push(("section42-c0", c0_rows));
    let mut eq11_rows = Vec::new();
    test_add_section42_chart_eq11_rows(&mut eq11_rows, map);
    grouped_rows.push(("eq11-sector-curve", eq11_rows));
    let mut eq12_rows = Vec::new();
    test_add_section42_chart_endpoint_aware_eq12_rows(&mut eq12_rows, map);
    grouped_rows.push(("eq12-strip-join", eq12_rows));
    let mut eq13_rows = Vec::new();
    test_add_section42_chart_eq13_rows(&mut eq13_rows, map, 2.0 * c);
    grouped_rows.push(("eq13-lower-strip", eq13_rows));
    for (label, rows) in grouped_rows {
        let mut vector_constraints = central_pin_rows.clone();
        vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: rows,
                regularization_weight: 1.0e-10,
            })
            .expect("Appendix A refined scalar group solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "Appendix-A refined central pins with internal C2 scalar group {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            "Appendix A refined scalar group residual must be finite"
        );
    }

    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let c2_edge_scalar_rows = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut vector_constraints = central_pin_rows;
    vector_constraints.extend(kp20_planar_tau_internal_c2_constraints());
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: c2_edge_scalar_rows,
            regularization_weight: 1.0e-10,
        })
        .expect("Appendix A refined C2-edge chart solve succeeds");
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "Appendix-A refined central pins full rows c2-edge chart: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        "Appendix A refined C2-edge chart residual must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_live_source_c2_edge_section42_chart_is_measured() {
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
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        -0.125,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        0.5 * (1.0 - ring.lambda),
    );
    let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(source_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("live source C2-edge chart solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "live source C2-edge section42 chart: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "live source C2-edge chart solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_section42_source_phase_variants_are_ranked() {
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
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut ranked = Vec::new();
    for previous_uses_near_m in [false, true] {
        for current_uses_near_m in [false, true] {
            for previous_reverse_u in [false, true] {
                for current_reverse_u in [false, true] {
                    for transpose_indices in [false, true] {
                        for swap_pieces in [false, true] {
                            let source_variant = TestTensorBorderChartVariant {
                                previous_uses_near_m,
                                current_uses_near_m,
                                previous_reverse_u,
                                current_reverse_u,
                                transpose_indices,
                                swap_pieces,
                            };
                            let base_rows = tensor_border_chart_variant_rows(
                                &tensor_border,
                                sector,
                                source_variant,
                            );
                            for rotate_current in [false, true] {
                                let rotated_rows = if rotate_current {
                                    rotate_tensor_chart_piece_targets(
                                        &base_rows,
                                        Kp20PieceId::NearMidpointRight,
                                        -angle,
                                    )
                                } else {
                                    base_rows.clone()
                                };
                                for left_scale in [-0.125, 0.5 * (1.0 - ring.lambda)] {
                                    let left_scaled =
                                        scale_tensor_chart_depth_from_boundary_for_piece(
                                            &rotated_rows,
                                            source_variant,
                                            Kp20PieceId::NearMidpointLeft,
                                            2,
                                            left_scale,
                                        );
                                    let scaled_rows =
                                        scale_tensor_chart_depth_from_boundary_for_piece(
                                            &left_scaled,
                                            source_variant,
                                            Kp20PieceId::NearMidpointRight,
                                            2,
                                            0.5 * (1.0 - ring.lambda),
                                        );
                                    let source_rows =
                                        without_tensor_chart_depth(&scaled_rows, source_variant, 1);
                                    let mut vector_constraints =
                                        kp20_planar_tau_internal_c2_constraints();
                                    vector_constraints.extend(
                                        kp20_planar_tau_fig8_indicated_averaging_constraints(),
                                    );
                                    vector_constraints.extend(source_rows);
                                    let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                                        &Kp20PlanarTauF3ScalarSolveInput {
                                            initial: planar_tau_seed_macro_patch(valence, sector),
                                            vector_constraints,
                                            scalar_constraints: scalar_constraints.clone(),
                                            regularization_weight: 1.0e-10,
                                        },
                                    )
                                    .expect("C2-edge source phase solve succeeds");
                                    let full_tensor_residual = scaled_rows
                                        .iter()
                                        .map(|row| {
                                            let residual =
                                                vector_residual(&solved.tau_macro_patch, row);
                                            residual[0].abs().max(residual[1].abs())
                                        })
                                        .fold(0.0_f64, f64::max);
                                    let local_solved = global_planar_patch_to_local(
                                        valence,
                                        sector,
                                        solved.tau_macro_patch,
                                    );
                                    let validity = kp20_planar_bi6_tau_validity_report(
                                        &local_solved,
                                        &SAMPLES,
                                        1.0e-10,
                                        1.0e-10,
                                    );
                                    ranked.push((
                                        validity.orientation.negative_determinants,
                                        validity.domain.outside_triangle_samples,
                                        solved.max_constraint_residual,
                                        full_tensor_residual,
                                        validity.orientation.min_determinant,
                                        source_variant,
                                        rotate_current,
                                        left_scale,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ranked.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(a.2.partial_cmp(&b.2).unwrap())
    });
    for (rank, row) in ranked.iter().take(12).enumerate() {
        eprintln!(
            "C2-edge source phase rank={rank} neg={} outside={} residual={:e} full_tensor={:e} min_det={:e} rotate_current={} left_scale={:e} variant={:?}",
            row.0, row.1, row.2, row.3, row.4, row.6, row.7, row.5,
        );
    }
    let best = ranked.first().expect("source phase cases exist");
    assert!(
        best.2.is_finite() && best.3.is_finite(),
        "C2-edge source phase ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_section42_best_source_phase_depth_scales_are_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    const LEFT_SCALES: [f64; 11] = [
        -0.35, -0.3, -0.25, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, 0.0,
    ];
    const RIGHT_SCALES: [f64; 9] = [
        0.125,
        0.175,
        0.2250058227408514,
        0.275,
        0.325,
        0.375,
        0.45,
        0.6,
        0.75,
    ];
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let rotated_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut best = (
        0.0_f64,
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for left_scale in LEFT_SCALES {
        for right_scale in RIGHT_SCALES {
            let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
                &rotated_rows,
                source_variant,
                Kp20PieceId::NearMidpointLeft,
                2,
                left_scale,
            );
            let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
                &left_scaled,
                source_variant,
                Kp20PieceId::NearMidpointRight,
                2,
                right_scale,
            );
            let source_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
            let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
            vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
            vector_constraints.extend(source_rows);
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("C2-edge source depth scale solve succeeds");
            let full_tensor_residual = scaled_rows
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            let local_solved =
                global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
            if validity.domain.outside_triangle_samples == 0
                && (validity.orientation.negative_determinants <= 21
                    || solved.max_constraint_residual <= 5.5e-2)
            {
                eprintln!(
                    "C2-edge best phase depth left={left_scale:e} right={right_scale:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            let row = (
                left_scale,
                right_scale,
                validity.orientation.negative_determinants,
                validity.domain.outside_triangle_samples,
                solved.max_constraint_residual,
                validity.orientation.min_determinant,
                full_tensor_residual,
            );
            if row.3 < best.3
                || (row.3 == best.3 && row.2 < best.2)
                || (row.3 == best.3 && row.2 == best.2 && row.4 < best.4)
            {
                best = row;
            }
        }
    }
    eprintln!(
        "best C2-edge best phase depth scales: left={:e} right={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5, best.6,
    );
    assert!(
        best.4.is_finite() && best.6.is_finite(),
        "C2-edge best phase depth scale ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_without_depth_scaling_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(source_rows.clone());
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge paper-oriented source solve succeeds");
    let full_tensor_residual = source_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented source no depth scaling: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}] worst={worst:?}@({u:.2},{v:.2})",
        residual = solved.max_constraint_residual,
        full_tensor = full_tensor_residual,
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
    assert!(
        solved.max_constraint_residual.is_finite() && full_tensor_residual.is_finite(),
        "C2-edge paper-oriented source solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_projection_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut fixed_vector_constraints = kp20_planar_tau_internal_c2_constraints();
    fixed_vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let (projected_rows, projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            source_rows.clone(),
            &fixed_vector_constraints,
            &scalar_constraints,
            &source_rows
                .iter()
                .map(|row| row.terms[0].0)
                .collect::<Vec<_>>(),
        )
        .expect("C2-edge paper-oriented source projection succeeds");
    let max_source_delta = source_rows
        .iter()
        .zip(&projected_rows)
        .map(|(source, projected)| {
            (projected.target[0] - source.target[0]).hypot(projected.target[1] - source.target[1])
        })
        .fold(0.0_f64, f64::max);
    let mut vector_constraints = fixed_vector_constraints;
    vector_constraints.extend(projected_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("C2-edge projected paper-oriented source solve succeeds");
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented source projected: projection_residual={projection_residual:e} projection_delta={projection_delta:e} max_source_delta={max_delta:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        projection_residual = projection_report.max_row_residual,
        projection_delta = projection_report.max_source_delta,
        max_delta = max_source_delta,
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
        solved.max_constraint_residual.is_finite() && max_source_delta.is_finite(),
        "C2-edge paper-oriented source projection solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_affine_normalization_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let mut augmented_scalar_constraints = Vec::new();
    for row in &source_rows {
        let coeff = row.terms[0].0;
        let [x, y] = row.target;
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -x),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(1), -y),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(2), -1.0),
            ],
            target: x,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(3), -x),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(4), -y),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(5), -1.0),
            ],
            target: y,
        });
    }
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            augmented_scalar_constraints,
            extra_variable_count: 6,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
        },
    )
    .expect("C2-edge paper-oriented affine source solve succeeds");
    let d = &solved.extra_values;
    let affine_det = (1.0 + d[0]) * (1.0 + d[4]) - d[1] * d[3];
    let max_source_delta = source_rows
        .iter()
        .map(|row| {
            let [x, y] = row.target;
            (d[0] * x + d[1] * y + d[2]).hypot(d[3] * x + d[4] * y + d[5])
        })
        .fold(0.0_f64, f64::max);
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented affine source: residual={residual:e} affine_det={affine_det:e} max_source_delta={max_delta:e} extras=[{d0:e},{d1:e},{d2:e},{d3:e},{d4:e},{d5:e}] outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
        residual = solved.max_constraint_residual,
        max_delta = max_source_delta,
        d0 = d[0],
        d1 = d[1],
        d2 = d[2],
        d3 = d[3],
        d4 = d[4],
        d5 = d[5],
        outside = validity.domain.outside_triangle_samples,
        domain_samples = validity.domain.sample_count,
        neg = validity.orientation.negative_determinants,
        orientation_samples = validity.orientation.sample_count,
        zero = validity.orientation.near_zero_determinants,
        min_det = validity.orientation.min_determinant,
        max_det = validity.orientation.max_determinant,
    );
    assert!(
        solved.max_constraint_residual.is_finite() && affine_det.is_finite(),
        "C2-edge affine source normalization solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_translation_normalization_is_measured() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let mut augmented_scalar_constraints = Vec::new();
    for row in &source_rows {
        let coeff = row.terms[0].0;
        let [x, y] = row.target;
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 0)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(0), -1.0),
            ],
            target: x,
        });
        augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
            terms: vec![
                (
                    Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(coeff, 1)),
                    1.0,
                ),
                (Kp20PlanarTauAugmentedScalarCoeff::Extra(1), -1.0),
            ],
            target: y,
        });
    }
    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints: test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c),
            augmented_scalar_constraints,
            extra_variable_count: 2,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
        },
    )
    .expect("C2-edge paper-oriented translation source solve succeeds");
    let tx = solved.extra_values[0];
    let ty = solved.extra_values[1];
    let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let validity = kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
    eprintln!(
        "C2-edge paper-oriented translation source: residual={residual:e} translation=({tx:e},{ty:e}) outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
        "C2-edge translation source normalization solve must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_eq14_d1_3_sweep_is_ranked() {
    const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
    let valence = 5;
    let sector = 0;
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let default_d1_3 = 2.0 * c / 3.0;
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut best = (
        0.0_f64,
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for d1_3 in [
        -default_d1_3,
        -0.5 * default_d1_3,
        0.0,
        0.125,
        0.5 * default_d1_3,
        default_d1_3,
        1.5 * default_d1_3,
        2.0 * default_d1_3,
    ] {
        let split = split_kp20_tensor_border_with_d1_3(&tensor_border, d1_3).expect("eq14 split");
        let previous = (sector + split.len() - 1) % split.len();
        let base_rows = tensor_border_chart_variant_rows_from_split_rows(
            split[previous].near_m.rows,
            split[sector].near_c.rows,
            source_variant,
        );
        let source_rows =
            rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(source_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge d1_3 source solve succeeds");
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "C2-edge paper-oriented eq14 d1_3={d1_3:e}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            d1_3,
            validity.orientation.negative_determinants,
            validity.domain.outside_triangle_samples,
            solved.max_constraint_residual,
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge paper-oriented eq14 d1_3: d1_3={:e} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge d1_3 source ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_paper_oriented_source_depth_subsets_are_ranked() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, depths) in [
        ("depth-0", &[0_usize][..]),
        ("depth-1", &[1_usize][..]),
        ("depth-2", &[2_usize][..]),
        ("depth-0-1", &[0_usize, 1][..]),
        ("depth-0-2", &[0_usize, 2][..]),
        ("depth-1-2", &[1_usize, 2][..]),
        ("depth-0-1-2", &[0_usize, 1, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge source depth subset solve succeeds");
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "C2-edge paper-oriented source depth subset {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge paper-oriented source depth subset: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge source depth subset ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_appendix_a_central_pins_with_source_depths_are_ranked() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let local_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, local_tau);
    let mut central_pin_rows = Vec::new();
    for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        for i in 0..=2 {
            for j in 0..=(2 - i) {
                let target = paper_tau.piece(piece_id).cps[i][j];
                central_pin_rows.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece_id, i, j),
                    [target.x, target.y],
                ));
            }
        }
    }
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    );
    for (label, depths) in [
        ("depth-0", &[0_usize][..]),
        ("depth-1", &[1_usize][..]),
        ("depth-2", &[2_usize][..]),
        ("depth-0-2", &[0_usize, 2][..]),
        ("depth-1-2", &[1_usize, 2][..]),
        ("depth-0-1-2", &[0_usize, 1, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(central_pin_rows.clone());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge central source depth solve succeeds");
        let full_tensor_residual = source_rows
            .iter()
            .map(|row| {
                let residual = vector_residual(&solved.tau_macro_patch, row);
                residual[0].abs().max(residual[1].abs())
            })
            .fold(0.0_f64, f64::max);
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "C2-edge Appendix-A central pins source depth {label}: residual={residual:e} full_tensor={full_tensor:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
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
            validity.orientation.min_determinant,
            full_tensor_residual,
        );
        if row.2 < best.2
            || (row.2 == best.2 && row.1 < best.1)
            || (row.2 == best.2 && row.1 == best.1 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best C2-edge Appendix-A central pins source depth: label={} neg={} outside={} residual={:e} min_det={:e} full_tensor={:e}",
        best.0, best.1, best.2, best.3, best.4, best.5,
    );
    assert!(
        best.3.is_finite() && best.5.is_finite(),
        "C2-edge central source depth ranking must remain finite"
    );
}

#[test]
fn rotational_c2_edge_source_depth_candidates_are_compared_to_appendix_a_tau() {
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
        current_reverse_u: false,
        transpose_indices: true,
        swap_pieces: true,
    };
    let base_rows = tensor_border_chart_variant_rows(&tensor_border, sector, source_variant);
    let source_rows =
        rotate_tensor_chart_piece_targets(&base_rows, Kp20PieceId::NearMidpointRight, -angle);
    let c2_edge_map = TestSection42ChartMap {
        p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: TestSection42ChartStrip {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    };
    let scalar_constraints = test_section42_endpoint_aware_chart_rows(c2_edge_map, 2.0 * c);
    let paper_tau = sample_kp20_macro_patch_from_tau_family_maps(
        &local_identity_tri6_guide(),
        valence,
        sector,
        Kp20TauFamily::PaperRefinedApex,
        kp20_paper_quadrant_sector_piece_maps(),
    )
    .expect("Appendix A tau samples with identity guide");
    let paper_tau = local_planar_patch_to_global(valence, sector, paper_tau);
    for (label, depths) in [
        ("depth-1", &[1_usize][..]),
        ("depth-0-2", &[0_usize, 2][..]),
    ] {
        let kept_rows = keep_tensor_chart_depths(&source_rows, source_variant, depths);
        let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
        vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
        vector_constraints.extend(kept_rows);
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: scalar_constraints.clone(),
                regularization_weight: 1.0e-10,
            })
            .expect("C2-edge source depth candidate solve succeeds");
        let local_solved = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity =
            kp20_planar_bi6_tau_validity_report(&local_solved, &SAMPLES, 1.0e-10, 1.0e-10);
        let mut deltas = Vec::new();
        for piece in Kp20PieceId::all() {
            let solved_piece = solved.tau_macro_patch.piece(piece);
            let paper_piece = paper_tau.piece(piece);
            for i in 0..=6 {
                for j in 0..=6 {
                    let dx = solved_piece.cps[i][j].x - paper_piece.cps[i][j].x;
                    let dy = solved_piece.cps[i][j].y - paper_piece.cps[i][j].y;
                    deltas.push((dx.hypot(dy), piece, i, j, dx, dy));
                }
            }
        }
        deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        eprintln!(
            "C2-edge source candidate {label}: residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
            outside = validity.domain.outside_triangle_samples,
            domain_samples = validity.domain.sample_count,
            neg = validity.orientation.negative_determinants,
            orientation_samples = validity.orientation.sample_count,
            zero = validity.orientation.near_zero_determinants,
            min_det = validity.orientation.min_determinant,
            max_det = validity.orientation.max_determinant,
        );
        for (rank, (norm, piece, i, j, dx, dy)) in deltas.iter().take(10).enumerate() {
            eprintln!(
                "  vs Appendix-A rank={rank} piece={piece:?}[{i},{j}] delta=({dx:e},{dy:e}) norm={norm:e}"
            );
        }
        assert!(
            solved.max_constraint_residual.is_finite(),
            "C2-edge source candidate comparison residual must remain finite"
        );
    }
}

#[test]
fn rotational_endpoint_aware_best_scaled_source_frame_is_localized() {
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
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
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
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    vector_constraints.extend(projected_tensor_rows);
    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial: planar_tau_seed_macro_patch(valence, sector),
            vector_constraints,
            scalar_constraints,
            regularization_weight: 1.0e-10,
        })
        .expect("best scaled source frame localization solve succeeds");
    let full_tensor_residual = scaled_rows
        .iter()
        .map(|row| {
            let residual = vector_residual(&solved.tau_macro_patch, row);
            residual[0].abs().max(residual[1].abs())
        })
        .fold(0.0_f64, f64::max);
    let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
    let piece_reports = test_piece_domain_orientation_reports(&local_tau, &SAMPLES, 1.0e-10);
    eprintln!(
        "best scaled source frame localized residual={:e} full_tensor={:e}",
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
        "best scaled source frame localization must be finite"
    );
}

#[test]
fn rotational_endpoint_aware_eq13_b0_after_scaled_source_is_measured() {
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
    let spectral_cross_depth_scale = 0.5 * (1.0 - ring.lambda);
    let left_scaled = scale_tensor_chart_depth_from_boundary_for_piece(
        &sector_rotated_rows,
        source_variant,
        Kp20PieceId::NearMidpointLeft,
        2,
        spectral_cross_depth_scale,
    );
    let scaled_rows = scale_tensor_chart_depth_from_boundary_for_piece(
        &left_scaled,
        source_variant,
        Kp20PieceId::NearMidpointRight,
        2,
        spectral_cross_depth_scale,
    );
    let projected_tensor_rows = without_tensor_chart_depth(&scaled_rows, source_variant, 1);
    let base_map = TestSection42ChartMap {
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
    let b0_cases = [
        ("zero", 0.0),
        ("near-m-c", c),
        ("paper-2c", 2.0 * c),
        ("minus-c", -c),
        ("minus-2c", -2.0 * c),
    ];
    let mut best = (
        "none",
        usize::MAX,
        usize::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (label, b0) in b0_cases {
        let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
        vector_constraints.extend(projected_tensor_rows.clone());
        let solved =
            solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
                initial: planar_tau_seed_macro_patch(valence, sector),
                vector_constraints,
                scalar_constraints: test_section42_endpoint_aware_chart_rows(base_map, b0),
                regularization_weight: 1.0e-10,
            })
            .expect("eq13 b0 scaled source solve succeeds");
        let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
        let validity = kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
        eprintln!(
            "scaled source eq13 b0 {label}: b0={b0:e} residual={residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} det=[{min_det:e},{max_det:e}]",
            residual = solved.max_constraint_residual,
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
            validity.orientation.min_determinant,
        );
        if row.1 < best.1
            || (row.1 == best.1 && row.2 < best.2)
            || (row.1 == best.1 && row.2 == best.2 && row.3 < best.3)
        {
            best = row;
        }
    }
    eprintln!(
        "best scaled source eq13 b0: label={} neg={} outside={} residual={:e} min_det={:e}",
        best.0, best.1, best.2, best.3, best.4,
    );
    assert!(best.3.is_finite(), "eq13 b0 ranking must be finite");
}

#[test]
fn rotational_endpoint_aware_tensor_border_half_assignment_is_ranked() {
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
    let cases = [
        ("previous-near-c/current-near-c", false, false),
        ("previous-near-c/current-near-m", false, true),
        ("previous-near-m/current-near-c", true, false),
        ("previous-near-m/current-near-m", true, true),
    ];

    let mut best_fixed = ("none", f64::INFINITY, usize::MAX, usize::MAX);
    for (label, previous_uses_near_m, current_uses_near_m) in cases {
        let tensor_rows = tensor_border_half_assignment_rows(
            &tensor_border,
            sector,
            previous_uses_near_m,
            current_uses_near_m,
        );
        for project_targets in [false, true] {
            let tensor_constraints = if project_targets {
                let free_coeffs = tensor_rows
                    .iter()
                    .map(|row| row.terms[0].0)
                    .collect::<Vec<_>>();
                kp20_project_section42_source_tensor_rows_with_free_coeffs(
                    tensor_rows.clone(),
                    &fig8_rows,
                    &scalar_constraints,
                    &free_coeffs,
                )
                .expect("tensor row projection succeeds")
                .0
            } else {
                tensor_rows.clone()
            };
            let mut vector_constraints = fig8_rows.clone();
            vector_constraints.extend(tensor_constraints.clone());
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch(valence, sector),
                    vector_constraints,
                    scalar_constraints: scalar_constraints.clone(),
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("half-assignment solve succeeds");
            let local_tau = global_planar_patch_to_local(valence, sector, solved.tau_macro_patch);
            let validity =
                kp20_planar_bi6_tau_validity_report(&local_tau, &SAMPLES, 1.0e-10, 1.0e-10);
            let tensor_residual = tensor_constraints
                .iter()
                .map(|row| {
                    let residual = vector_residual(&solved.tau_macro_patch, row);
                    residual[0].abs().max(residual[1].abs())
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "endpoint-aware tensor-border half assignment {label} projected={project_targets}: residual={residual:e} tensor_residual={tensor_residual:e} outside={outside}/{domain_samples} neg={neg}/{orientation_samples} zero={zero} det=[{min_det:e},{max_det:e}]",
                residual = solved.max_constraint_residual,
                outside = validity.domain.outside_triangle_samples,
                domain_samples = validity.domain.sample_count,
                neg = validity.orientation.negative_determinants,
                orientation_samples = validity.orientation.sample_count,
                zero = validity.orientation.near_zero_determinants,
                min_det = validity.orientation.min_determinant,
                max_det = validity.orientation.max_determinant,
            );
            if !project_targets
                && (solved.max_constraint_residual < best_fixed.1
                    || (solved.max_constraint_residual == best_fixed.1
                        && validity.orientation.negative_determinants < best_fixed.2))
            {
                best_fixed = (
                    label,
                    solved.max_constraint_residual,
                    validity.orientation.negative_determinants,
                    validity.domain.outside_triangle_samples,
                );
            }
        }
    }
    eprintln!(
        "best fixed tensor-border half assignment: label={} residual={:e} neg={} outside={}",
        best_fixed.0, best_fixed.1, best_fixed.2, best_fixed.3
    );
    assert!(
        best_fixed.1.is_finite(),
        "half-assignment ranking must produce finite solves"
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
