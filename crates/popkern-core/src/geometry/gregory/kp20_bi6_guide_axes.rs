//! The guide least-squares KKT solve shared by its three coordinates
//! (split from `kp20_bi6_guide_lsq.rs` under the 500-line budget).

use super::direct_linalg::{condition_number_dense, solve_dense_sides};
use super::kp20_bi6_guide_lsq::{CENTER_FLAT, LinearRow, ResidualRows, var_index};

/// One coordinate's guide solution and the shared KKT condition number.
pub(crate) struct AxisSolution {
    pub(crate) values: Vec<f64>,
    pub(crate) kkt_condition_number: f64,
}

/// The three coordinate solves share one KKT matrix (rows and
/// constraints are axis-independent; only the targets and the centre
/// pin differ), so it is assembled, conditioned and eliminated once
/// (PK-600); each axis's arithmetic is its single solve's exactly.
pub(crate) fn solve_axes(
    vars: usize,
    residuals: &ResidualRows,
    rhs: [&[f64]; 3],
    constraints: &[LinearRow],
    center: [f64; 3],
) -> Option<[AxisSolution; 3]> {
    let size = vars + constraints.len() + 1;
    let mut a = vec![vec![0.0_f64; size]; size];
    for row in &residuals.rows {
        accumulate_normal_matrix(&mut a, &row.terms, 1.0);
    }
    for (i, row) in a.iter_mut().enumerate().take(vars) {
        row[i] += 1.0e-12;
    }
    add_lagrange_terms(vars, &mut a, &[(var_index(0, CENTER_FLAT), 1.0)]);
    for (cidx, row) in constraints.iter().enumerate() {
        add_lagrange_terms(vars + 1 + cidx, &mut a, &row.terms);
    }
    let sides = (0..3)
        .map(|k| {
            let mut b = vec![0.0_f64; size];
            for (row, &target) in residuals.rows.iter().zip(rhs[k]) {
                accumulate_normal_side(&mut b, &row.terms, target, 1.0);
            }
            b[vars] = center[k];
            b
        })
        .collect();
    let kkt_condition_number = condition_number_dense(&a)?.condition_number;
    let mut x = solve_dense_sides(a, sides)?.into_iter();
    Some(core::array::from_fn(|_| AxisSolution {
        values: x.next().map(|x| x[..vars].to_vec()).unwrap_or_default(),
        kkt_condition_number,
    }))
}

fn accumulate_normal_matrix(a: &mut [Vec<f64>], terms: &[(usize, f64)], weight: f64) {
    for &(i, wi) in terms {
        for &(j, wj) in terms {
            a[i][j] += weight * wi * wj;
        }
    }
}

fn accumulate_normal_side(b: &mut [f64], terms: &[(usize, f64)], rhs: f64, weight: f64) {
    for &(i, wi) in terms {
        b[i] += weight * wi * rhs;
    }
}

#[cfg(test)]
pub(super) fn accumulate_normal_row(
    a: &mut [Vec<f64>],
    b: &mut [f64],
    terms: &[(usize, f64)],
    rhs: f64,
    weight: f64,
) {
    accumulate_normal_side(b, terms, rhs, weight);
    accumulate_normal_matrix(a, terms, weight);
}

fn add_lagrange_terms(row_idx: usize, a: &mut [Vec<f64>], terms: &[(usize, f64)]) {
    for &(var, weight) in terms {
        a[var][row_idx] += weight;
        a[row_idx][var] += weight;
    }
}

#[cfg(test)]
pub(super) fn add_lagrange_row(
    row_idx: usize,
    a: &mut [Vec<f64>],
    b: &mut [f64],
    terms: &[(usize, f64)],
    rhs: f64,
) {
    b[row_idx] = rhs;
    add_lagrange_terms(row_idx, a, terms);
}
