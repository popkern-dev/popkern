use super::*;
use crate::geometry::gregory::kp_polar_boundary_input::TensorBorderBi3;
use crate::geometry::gregory::kp20_bi6_local_constraint::{
    kp20_eq8_g1_residual, kp20_eq9_g2_residual,
};
use crate::geometry::gregory::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use crate::geometry::gregory::kp20_planar_bi6_tau_f3_energy;
use crate::geometry::gregory::split_kp20_tensor_border;

fn p(x: f64, y: f64) -> Point3 {
    Point3::new(x, y, 0.0)
}

fn affine_tau_macro_patch() -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece in Kp20PieceId::all() {
        let offset = piece.index() as f64;
        for i in 0..=6 {
            for j in 0..=6 {
                out.piece_mut(piece).cps[i][j] = p(offset + i as f64 / 6.0, j as f64 / 6.0);
            }
        }
    }
    out
}

fn continuous_quadrant_affine_tau_macro_patch() -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece in Kp20PieceId::all() {
        let (offset_x, offset_y) = match piece {
            Kp20PieceId::NearOriginLeft => (0.0, 0.0),
            Kp20PieceId::NearOriginRight => (0.5, 0.0),
            Kp20PieceId::NearMidpointRight => (0.5, 0.5),
            Kp20PieceId::NearMidpointLeft => (0.0, 0.5),
        };
        for i in 0..=6 {
            for j in 0..=6 {
                out.piece_mut(piece).cps[i][j] = p(
                    offset_x + 0.5 * i as f64 / 6.0,
                    offset_y + 0.5 * j as f64 / 6.0,
                );
            }
        }
    }
    out
}

fn all_coeff_pins(tau: &Kp20MacroPatch) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut constraints = Vec::new();
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let cp = tau.piece(piece).cps[i][j];
                constraints.push(Kp20PlanarTauLinearConstraint::pin(
                    Kp20PlanarTauCoeff::new(piece, i, j),
                    [cp.x, cp.y],
                ));
            }
        }
    }
    constraints
}

fn boundary_pins(tau: &Kp20MacroPatch) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut constraints = Vec::new();
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                if i == 0 || i == 6 || j == 0 || j == 6 {
                    let cp = tau.piece(piece).cps[i][j];
                    constraints.push(Kp20PlanarTauLinearConstraint::pin(
                        Kp20PlanarTauCoeff::new(piece, i, j),
                        [cp.x, cp.y],
                    ));
                }
            }
        }
    }
    constraints
}

fn outer_boundary_pins(tau: &Kp20MacroPatch) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut constraints = Vec::new();
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let is_outer = match piece {
                    Kp20PieceId::NearOriginLeft => i == 0 || j == 0,
                    Kp20PieceId::NearOriginRight => i == 6 || j == 0,
                    Kp20PieceId::NearMidpointRight => i == 6 || j == 6,
                    Kp20PieceId::NearMidpointLeft => i == 0 || j == 6,
                };
                if is_outer {
                    let cp = tau.piece(piece).cps[i][j];
                    constraints.push(Kp20PlanarTauLinearConstraint::pin(
                        Kp20PlanarTauCoeff::new(piece, i, j),
                        [cp.x, cp.y],
                    ));
                }
            }
        }
    }
    constraints
}

fn planar_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::new();
    for sector in 0..n {
        let sector_offset = sector as f64;
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for (depth, row) in rows.iter_mut().enumerate() {
            for (i, cp) in row.iter_mut().enumerate() {
                let u = i as f64 / 3.0;
                let v = depth as f64 / 3.0;
                *cp = p(sector_offset + u, 0.25 * sector_offset + v);
            }
        }
        sectors.push(rows);
    }
    TensorBorderBi3 { sectors }
}

#[test]
fn planar_tau_f3_solver_preserves_fully_pinned_affine_tau() {
    let initial = affine_tau_macro_patch();

    let solved = solve_kp20_planar_tau_f3_diagnosed(&Kp20PlanarTauF3SolveInput {
        initial,
        constraints: all_coeff_pins(&initial),
        regularization_weight: 1.0e-12,
    })
    .expect("fully pinned affine solve succeeds");

    assert!(solved.max_constraint_residual <= 1.0e-10);
    assert!(kp20_planar_bi6_tau_f3_energy(&solved.tau_macro_patch) <= 1.0e-16);
    assert!(solved.kkt_condition_number.is_some_and(f64::is_finite));
    assert_macro_patch_close(&solved.tau_macro_patch, &initial, 1.0e-10);
}

#[test]
fn planar_tau_f3_solver_reduces_interior_cubic_perturbation_with_boundary_pins() {
    let boundary = affine_tau_macro_patch();
    let mut initial = boundary;
    for piece in Kp20PieceId::all() {
        for i in 1..6 {
            for j in 1..6 {
                let u = i as f64 / 6.0;
                initial.piece_mut(piece).cps[i][j].y += 0.1 * u * u * u;
            }
        }
    }
    let initial_energy = kp20_planar_bi6_tau_f3_energy(&initial);

    let solved = solve_kp20_planar_tau_f3(&Kp20PlanarTauF3SolveInput {
        initial,
        constraints: boundary_pins(&boundary),
        regularization_weight: 1.0e-12,
    })
    .expect("boundary-pinned F3 solve succeeds");
    let solved_energy = kp20_planar_bi6_tau_f3_energy(&solved.tau_macro_patch);

    assert!(solved.max_constraint_residual <= 1.0e-9);
    assert!(solved_energy < initial_energy * 0.1);
    for constraint in boundary_pins(&boundary) {
        let coeff = constraint.terms[0].0;
        let cp = solved.tau_macro_patch.piece(coeff.piece).cps[coeff.i][coeff.j];
        assert!((cp.x - constraint.target[0]).abs() <= 1.0e-9);
        assert!((cp.y - constraint.target[1]).abs() <= 1.0e-9);
    }
}

#[test]
fn planar_tau_tensor_border_constraints_match_kp20_split_rows() {
    let tensor_border = planar_tensor_border(5);
    let constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, 2).expect("tensor-border rows");
    let initialized = split_kp20_tensor_border(&tensor_border).expect("split succeeds");
    let mut tau = Kp20MacroPatch::zero();
    tau.piece_mut(Kp20PieceId::NearMidpointLeft).cps[0..3]
        .copy_from_slice(&initialized[1].near_c.rows);
    tau.piece_mut(Kp20PieceId::NearMidpointRight).cps[0..3]
        .copy_from_slice(&initialized[2].near_c.rows);

    assert_eq!(constraints.len(), 42);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_constraint_residual(&tau, constraint);
        assert!(
            residual[0].abs() <= 1.0e-12 && residual[1].abs() <= 1.0e-12,
            "row {row_index} residual {residual:?}"
        );
    }
}

#[test]
fn planar_tau_f3_solver_preserves_tensor_border_rows() {
    let tensor_border = planar_tensor_border(5);
    let constraints =
        kp20_planar_tau_tensor_border_constraints(&tensor_border, 0).expect("tensor-border rows");
    let mut initial = affine_tau_macro_patch();
    for piece in Kp20PieceId::all() {
        for i in 1..6 {
            for j in 1..6 {
                initial.piece_mut(piece).cps[i][j].x += 0.02 * (i * j) as f64;
            }
        }
    }

    let solved = solve_kp20_planar_tau_f3(&Kp20PlanarTauF3SolveInput {
        initial,
        constraints: constraints.clone(),
        regularization_weight: 1.0e-12,
    })
    .expect("tensor-border constrained F3 solve succeeds");

    assert!(solved.max_constraint_residual <= 1.0e-9);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_constraint_residual(&solved.tau_macro_patch, constraint);
        assert!(
            residual[0].abs() <= 1.0e-9 && residual[1].abs() <= 1.0e-9,
            "row {row_index} residual {residual:?}"
        );
    }
}

#[test]
fn planar_tau_sector_g1_g2_constraints_match_residual_extractors() {
    let mut tau = Kp20MacroPatch::zero();
    for piece in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
        let piece_offset = piece.index() as f64;
        for i in 0..=6 {
            for j in 0..=6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                tau.piece_mut(piece).cps[i][j] = p(
                    piece_offset + 0.3 * u + 0.2 * v + 0.07 * u * v,
                    -0.4 * piece_offset + 0.11 * u * u - 0.05 * v,
                );
            }
        }
    }
    let b = kp20_b_coeffs_near_o((core::f64::consts::TAU / 5.0).cos());
    let constraints = kp20_planar_tau_sector_g1_g2_constraints(
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearOriginRight,
        b,
    );
    let eq8 = kp20_eq8_g1_residual(
        &tau.piece(Kp20PieceId::NearOriginLeft).cps,
        &tau.piece(Kp20PieceId::NearOriginRight).cps,
        b,
    );
    let eq9 = kp20_eq9_g2_residual(
        &tau.piece(Kp20PieceId::NearOriginLeft).cps,
        &tau.piece(Kp20PieceId::NearOriginRight).cps,
        b,
    );

    assert_eq!(constraints.len(), 15);
    for row in 0..7 {
        let residual = evaluate_constraint_residual(&tau, &constraints[row]);
        assert!((residual[0] - eq8[row].x).abs() <= 1.0e-10);
        assert!((residual[1] - eq8[row].y).abs() <= 1.0e-10);
    }
    for row in 0..8 {
        let residual = evaluate_constraint_residual(&tau, &constraints[7 + row]);
        assert!((residual[0] - eq9[row].x).abs() <= 1.0e-10);
        assert!((residual[1] - eq9[row].y).abs() <= 1.0e-10);
    }
}

#[test]
fn planar_tau_sector_g1_g2_constraints_accept_constant_piece_pair() {
    let mut tau = Kp20MacroPatch::zero();
    for piece in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        for i in 0..=6 {
            for j in 0..=6 {
                tau.piece_mut(piece).cps[i][j] = p(0.25, -0.5);
            }
        }
    }
    let b = kp20_b_coeffs_near_m((core::f64::consts::TAU / 6.0).cos());
    let constraints = kp20_planar_tau_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        b,
    );

    assert_eq!(constraints.len(), 15);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_constraint_residual(&tau, constraint);
        assert!(
            residual[0].abs() <= 1.0e-12 && residual[1].abs() <= 1.0e-12,
            "row {row_index} residual {residual:?}"
        );
    }
}

#[test]
fn planar_tau_internal_c2_constraints_accept_continuous_quadrant_affine_tau() {
    let tau = continuous_quadrant_affine_tau_macro_patch();
    let constraints = kp20_planar_tau_internal_c2_constraints();

    assert_eq!(constraints.len(), 84);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_constraint_residual(&tau, constraint);
        assert!(
            residual[0].abs() <= 1.0e-12 && residual[1].abs() <= 1.0e-12,
            "row {row_index} residual {residual:?}"
        );
    }
}

#[test]
fn planar_tau_f3_solver_enforces_internal_c2_and_reduces_curvature() {
    let boundary = continuous_quadrant_affine_tau_macro_patch();
    let mut initial = boundary;
    for piece in Kp20PieceId::all() {
        for i in 1..6 {
            for j in 1..6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                initial.piece_mut(piece).cps[i][j].x += 0.04 * u * u * v;
                initial.piece_mut(piece).cps[i][j].y -= 0.03 * u * v * v;
            }
        }
    }
    let initial_energy = kp20_planar_bi6_tau_f3_energy(&initial);

    let mut constraints = outer_boundary_pins(&boundary);
    constraints.extend(kp20_planar_tau_internal_c2_constraints());
    let solved = solve_kp20_planar_tau_f3(&Kp20PlanarTauF3SolveInput {
        initial,
        constraints: constraints.clone(),
        regularization_weight: 1.0e-12,
    })
    .expect("outer-pinned internal-C2 F3 solve succeeds");
    let solved_energy = kp20_planar_bi6_tau_f3_energy(&solved.tau_macro_patch);

    assert!(solved.max_constraint_residual <= 1.0e-8);
    assert!(solved_energy < initial_energy * 0.1);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_constraint_residual(&solved.tau_macro_patch, constraint);
        assert!(
            residual[0].abs() <= 1.0e-8 && residual[1].abs() <= 1.0e-8,
            "row {row_index} residual {residual:?}"
        );
    }
}

#[test]
fn planar_tau_diagonal_symmetry_constraints_accept_symmetric_quadrant_affine_tau() {
    let tau = continuous_quadrant_affine_tau_macro_patch();
    let constraints = kp20_planar_tau_diagonal_symmetry_constraints();

    assert_eq!(constraints.len(), 196);
    for (row_index, constraint) in constraints.iter().enumerate() {
        let residual = evaluate_scalar_constraint_residual(&tau, constraint);
        assert!(
            residual.abs() <= 1.0e-12,
            "row {row_index} residual {residual}"
        );
    }
}

#[test]
fn planar_tau_f3_solver_enforces_internal_c2_and_diagonal_symmetry() {
    let boundary = continuous_quadrant_affine_tau_macro_patch();
    let mut initial = boundary;
    for piece in Kp20PieceId::all() {
        for i in 1..6 {
            for j in 1..6 {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                initial.piece_mut(piece).cps[i][j].x += 0.03 * u * u * (1.0 + v);
                initial.piece_mut(piece).cps[i][j].y -= 0.02 * v * v * (1.0 + u);
            }
        }
    }
    let initial_energy = kp20_planar_bi6_tau_f3_energy(&initial);

    let solved =
        solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
            initial,
            vector_constraints: {
                let mut constraints = outer_boundary_pins(&boundary);
                constraints.extend(kp20_planar_tau_internal_c2_constraints());
                constraints
            },
            scalar_constraints: kp20_planar_tau_diagonal_symmetry_constraints(),
            regularization_weight: 1.0e-12,
        })
        .expect("outer-pinned internal-C2 diagonal-symmetric F3 solve succeeds");
    let solved_energy = kp20_planar_bi6_tau_f3_energy(&solved.tau_macro_patch);

    assert!(solved.max_constraint_residual <= 1.0e-8);
    assert!(solved_energy < initial_energy * 0.1);
    for (row_index, constraint) in kp20_planar_tau_diagonal_symmetry_constraints()
        .iter()
        .enumerate()
    {
        let residual = evaluate_scalar_constraint_residual(&solved.tau_macro_patch, constraint);
        assert!(
            residual.abs() <= 1.0e-8,
            "row {row_index} residual {residual}"
        );
    }
}

fn assert_macro_patch_close(a: &Kp20MacroPatch, b: &Kp20MacroPatch, tol: f64) {
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let pa = a.piece(piece).cps[i][j];
                let pb = b.piece(piece).cps[i][j];
                assert!(
                    (pa.x - pb.x).abs() <= tol && (pa.y - pb.y).abs() <= tol,
                    "{piece:?}[{i}][{j}] mismatch: {pa:?} vs {pb:?}"
                );
            }
        }
    }
}

fn evaluate_scalar_constraint_residual(
    tau: &Kp20MacroPatch,
    constraint: &Kp20PlanarTauScalarLinearConstraint,
) -> f64 {
    let mut value = 0.0_f64;
    for &(coeff, weight) in &constraint.terms {
        let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
        let axis_value = if coeff.axis == 0 { cp.x } else { cp.y };
        value += weight * axis_value;
    }
    value - constraint.target
}

fn evaluate_constraint_residual(
    tau: &Kp20MacroPatch,
    constraint: &Kp20PlanarTauLinearConstraint,
) -> [f64; 2] {
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
