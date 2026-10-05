//! KP20 §4.4 step (iv): the sector-separating strip projection.
//!
//! Split from `kp20_paper_cap` for the file-size budget; see that
//! module for the assembly context and chart conventions.

use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_local_constraint::{kp20_eq8_jacobian_dense, kp20_eq9_jacobian_dense};
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use crate::geometry::Point3;
use popkern_math::PortableMath;

const P1: Kp20PieceId = Kp20PieceId::NearOriginLeft;
const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

/// How step (iv) treats the eq (8)/(9) rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20StripEnforcement {
    /// Paper-faithful: eq (8)/(9) are hard constraints (requires C²
    /// input; infeasible systems tear the strips on G¹-class data).
    HardG2,
    /// Containment for out-of-class input: C⁰/(11)/(12) stay hard,
    /// eq (8)/(9) become least-squares objectives.
    SoftG,
}

// ---- Step (iv): strip projection -----------------------------------------

/// One strip variable: (sector, piece, i, j).
type StripVar = (usize, Kp20PieceId, usize, usize);

/// One strip constraint row: variable terms and per-axis target.
type StripRow = (Vec<(usize, f64)>, [f64; 3]);

/// Strip variable set S for the curve between sectors `prev` and
/// `curr` (paper §4.4, mapped into the quadrant chart): the near-o
/// pieces contribute along-indices 4..6 at layers 0..2 plus the
/// (3,1)/(3,2) extras; the near-m pieces contribute along-indices
/// 0..3 at layers 0..2 (along 0 = the seam row, eq (12)-owned).
fn strip_variables(prev: usize, curr: usize) -> Vec<StripVar> {
    let mut vars = Vec::with_capacity(48);
    // Previous sector side: curve at j = 0 (along = i, layer = j).
    for i in 4..=6_usize {
        for j in 0..=2_usize {
            vars.push((prev, P1, i, j));
        }
    }
    // The curve coefficient at along-index 3 is determined through
    // the degree-5 rows (11) rather than sampled (the sampled value
    // is a corner-block average of a non-bi-6 function and breaks the
    // o-end G rows at ~1e-3 if frozen — measured).
    vars.push((prev, P1, 3, 0));
    vars.push((prev, P1, 3, 1));
    vars.push((prev, P1, 3, 2));
    for i in 0..=3_usize {
        for j in 0..=2_usize {
            vars.push((prev, P2, i, j));
        }
    }
    // Current sector side: curve at i = 0 (along = j, layer = i).
    for j in 4..=6_usize {
        for i in 0..=2_usize {
            vars.push((curr, P1, i, j));
        }
    }
    vars.push((curr, P1, 0, 3));
    vars.push((curr, P1, 1, 3));
    vars.push((curr, P1, 2, 3));
    for j in 0..=3_usize {
        for i in 0..=2_usize {
            vars.push((curr, P4, i, j));
        }
    }
    vars
}

/// Project the strip set of one curve onto the §4.1/§4.2 equation
/// manifold, minimizing distance to the sampled values. Returns the
/// largest coefficient movement, or `None` if the KKT solve fails.
pub(crate) fn project_strip(
    sectors: &mut [Kp20MacroPatch],
    n: usize,
    curve: usize,
    enforcement: Kp20StripEnforcement,
) -> Option<f64> {
    use super::direct_linalg::solve_dense;

    let prev = (curve + n - 1) % n;
    let curr = curve;
    let vars = strip_variables(prev, curr);

    // Constraint rows over [variables | constant]: each row is
    // Σ w·x = target, with context coefficients folded into target.
    let mut rows: Vec<StripRow> = Vec::new();

    // C0 along the curve: prev (i, 0) rows equal curr (0, i) rows for
    // the near-o pieces (i = 0..6) and prev P2 (i, 0) = curr P4
    // (0, i) for the near-m pieces.
    for i in 0..=6_usize {
        push_strip_row(
            &mut rows,
            &vars,
            sectors,
            vec![((prev, P1, i, 0), 1.0), ((curr, P1, 0, i), -1.0)],
        );
        push_strip_row(
            &mut rows,
            &vars,
            sectors,
            vec![((prev, P2, i, 0), 1.0), ((curr, P4, 0, i), -1.0)],
        );
    }
    // Degree-5 rows (11) on both halves of the shared curve
    // (expressed on the previous side; C0 carries them to the other).
    let binom6 = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    for piece in [P1, P2] {
        let mut terms = Vec::new();
        for (i, weight) in binom6.iter().enumerate() {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            terms.push(((prev, piece, i, 0), sign * weight));
        }
        push_strip_row(&mut rows, &vars, sectors, terms);
    }
    // Seam C2 rows (12), layers 0..2: previous side joins P1 (i = 6
    // end) to P2 (i = 0 end) in u; current side joins P1 (j = 6 end)
    // to P4 (j = 0 end) in v.
    for layer in 0..=2_usize {
        for (a_piece, b_piece, sector, transposed) in [(P1, P2, prev, false), (P1, P4, curr, true)]
        {
            let at = |i: usize, l: usize| -> StripVar {
                if transposed {
                    (sector, if i <= 6 { a_piece } else { b_piece }, l, i % 7)
                } else {
                    (sector, if i <= 6 { a_piece } else { b_piece }, i % 7, l)
                }
            };
            // Indices 0..6 = piece A, 7..13 = piece B (i % 7).
            push_strip_row(
                &mut rows,
                &vars,
                sectors,
                vec![(at(6, layer), 1.0), (at(7, layer), -1.0)],
            );
            push_strip_row(
                &mut rows,
                &vars,
                sectors,
                vec![
                    (at(6, layer), 1.0),
                    (at(5, layer), -1.0),
                    (at(8, layer), -1.0),
                    (at(7, layer), 1.0),
                ],
            );
            push_strip_row(
                &mut rows,
                &vars,
                sectors,
                vec![
                    (at(6, layer), 1.0),
                    (at(5, layer), -2.0),
                    (at(4, layer), 1.0),
                    (at(9, layer), -1.0),
                    (at(8, layer), 2.0),
                    (at(7, layer), -1.0),
                ],
            );
        }
    }
    // Full eq (8)/(9) rows for both half pairs: previous side natural
    // (along = i, across = j), current side transposed. Under SoftG
    // these become objective rows instead of constraints.
    let hard_row_count_before_g = rows.len();
    let c = (core::f64::consts::TAU / n as f64).pcos();
    for (prev_piece, curr_piece, b) in [
        (P1, P1, kp20_b_coeffs_near_o(c)),
        (P2, P4, kp20_b_coeffs_near_m(c)),
    ] {
        let jacobians: Vec<[f64; 98]> = kp20_eq8_jacobian_dense(b)
            .into_iter()
            .chain(kp20_eq9_jacobian_dense(b))
            .collect();
        for jac in jacobians {
            let mut terms = Vec::new();
            for side in 0..2_usize {
                for i in 0..7_usize {
                    for j in 0..7_usize {
                        let w = jac[side * 49 + i * 7 + j];
                        if w == 0.0 {
                            continue;
                        }
                        if side == 0 {
                            terms.push(((prev, prev_piece, i, j), w));
                        } else {
                            terms.push(((curr, curr_piece, j, i), w));
                        }
                    }
                }
            }
            push_strip_row(&mut rows, &vars, sectors, terms);
        }
    }

    // Constrained least squares: minimize distance to the sampled
    // values subject to the hard rows, per axis (shared KKT matrix).
    // Under SoftG, the G rows join the objective instead.
    let (hard_rows, soft_rows): (Vec<usize>, Vec<usize>) = match enforcement {
        Kp20StripEnforcement::HardG2 => ((0..rows.len()).collect(), Vec::new()),
        Kp20StripEnforcement::SoftG => (
            (0..hard_row_count_before_g).collect(),
            (hard_row_count_before_g..rows.len()).collect(),
        ),
    };
    let var_count = vars.len();
    let hard_row_refs: Vec<StripRow> = hard_rows.iter().map(|&i| rows[i].clone()).collect();
    let constraints_local = independent_strip_rows(var_count, &hard_row_refs);
    let constraints: Vec<usize> = constraints_local
        .into_iter()
        .map(|local| hard_rows[local])
        .collect();
    let size = var_count + constraints.len();
    let mut matrix = vec![vec![0.0_f64; size]; size];
    for (slot, row) in matrix.iter_mut().enumerate().take(var_count) {
        row[slot] = 1.0;
    }
    // SoftG: accumulate the G rows into the normal equations.
    for &row_index in &soft_rows {
        for &(vi, wi) in &rows[row_index].0 {
            for &(vj, wj) in &rows[row_index].0 {
                matrix[vi][vj] += wi * wj;
            }
        }
    }
    for (cidx, &row_index) in constraints.iter().enumerate() {
        for &(var, w) in &rows[row_index].0 {
            matrix[var][var_count + cidx] = w;
            matrix[var_count + cidx][var] = w;
        }
    }
    let mut max_movement = 0.0_f64;
    let mut solution = vec![Point3::new(0.0, 0.0, 0.0); var_count];
    for axis in 0..3_usize {
        let mut rhs = vec![0.0_f64; size];
        for (slot, var) in vars.iter().enumerate() {
            let value = strip_value(sectors, var);
            rhs[slot] = match axis {
                0 => value.x,
                1 => value.y,
                _ => value.z,
            };
        }
        for &row_index in &soft_rows {
            let target = rows[row_index].1[axis];
            for &(var, w) in &rows[row_index].0 {
                rhs[var] += w * target;
            }
        }
        for (cidx, &row_index) in constraints.iter().enumerate() {
            rhs[var_count + cidx] = rows[row_index].1[axis];
        }
        let solved = solve_dense(matrix.clone(), rhs)?;
        for (slot, value) in solution.iter_mut().enumerate() {
            match axis {
                0 => value.x = solved[slot],
                1 => value.y = solved[slot],
                _ => value.z = solved[slot],
            }
        }
    }
    for (slot, var) in vars.iter().enumerate() {
        let old = strip_value(sectors, var);
        let new = solution[slot];
        max_movement = max_movement
            .max((new.x - old.x).abs())
            .max((new.y - old.y).abs())
            .max((new.z - old.z).abs());
        sectors[var.0].piece_mut(var.1).cps[var.2][var.3] = new;
    }
    Some(max_movement)
}

fn strip_value(sectors: &[Kp20MacroPatch], var: &StripVar) -> Point3 {
    sectors[var.0].piece(var.1).cps[var.2][var.3]
}

fn push_strip_row(
    rows: &mut Vec<StripRow>,
    vars: &[StripVar],
    sectors: &[Kp20MacroPatch],
    terms: Vec<(StripVar, f64)>,
) {
    let mut var_terms = Vec::new();
    let mut target = [0.0_f64; 3];
    for (var, weight) in terms {
        match vars.iter().position(|v| *v == var) {
            Some(index) => var_terms.push((index, weight)),
            None => {
                let value = strip_value(sectors, &var);
                target[0] -= weight * value.x;
                target[1] -= weight * value.y;
                target[2] -= weight * value.z;
            }
        }
    }
    rows.push((var_terms, target));
}

/// Indices of a maximal linearly independent subset of the rows
/// (restricted to the variable columns).
fn independent_strip_rows(var_count: usize, rows: &[StripRow]) -> Vec<usize> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut keep = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        let mut dense = vec![0.0_f64; var_count];
        for &(var, w) in &row.0 {
            dense[var] += w;
        }
        for existing in &basis {
            let dot: f64 = dense.iter().zip(existing).map(|(a, b)| a * b).sum();
            for (value, e) in dense.iter_mut().zip(existing) {
                *value -= dot * e;
            }
        }
        let norm: f64 = dense.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm > 1.0e-9 {
            for value in dense.iter_mut() {
                *value /= norm;
            }
            basis.push(dense);
            keep.push(index);
        }
    }
    keep
}
