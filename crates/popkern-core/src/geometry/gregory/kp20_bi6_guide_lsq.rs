//! Direct KP18 guide solve for the KP20 total-degree-6 guide.
//!
//! KP18 describes the runtime guide as a linear least-squares fit
//! over the remaining guide degrees of freedom after fixing the
//! Catmull-Clark limit point and the sector-join equations. The paper
//! then suggests precomputed generating functions as the fast form.
//! This module implements the same linear construction directly:
//! unknowns are the scalar triangular BB coefficients of all sectors,
//! equality rows encode KP18 (2), (3), and the central cubic extension
//! rows from (4), (5), and residual rows fit the input tensor-border
//! degree-raised from bi-3 to bi-5.

use super::bernstein_tri6::{TRI6_CP_COUNT, flat_index};
use super::direct_linalg::{condition_number_dense, solve_dense};
use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp13_bi5_bi3_bi5_ops::degree_raise_1d_3_to_5;
use super::kp20_bi6_guide_continuity::kp18_guide_join_weights;
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use crate::geometry::Point3;

const DEGREE: usize = 6;
const CENTER_FLAT: usize = 0;

/// Error reported by the KP18 guide least-squares construction.
#[derive(Debug, Clone, PartialEq)]
pub enum Kp20GuideLsqError {
    /// KP18 excludes invalid or unsupported valences.
    UnsupportedValence(usize),
    /// Tensor-border shape does not match the valence.
    InvalidTensorBorder(String),
    /// The constrained normal system is singular for one coordinate.
    SingularSystem(usize),
}

impl std::fmt::Display for Kp20GuideLsqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedValence(n) => write!(f, "KP18 guide LSQ unsupported valence {n}"),
            Self::InvalidTensorBorder(s) => write!(f, "KP18 guide LSQ invalid tensor-border: {s}"),
            Self::SingularSystem(axis) => write!(f, "KP18 guide LSQ singular on axis {axis}"),
        }
    }
}

impl std::error::Error for Kp20GuideLsqError {}

/// Measured postconditions from a successful KP18 guide solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20GuideLsqReport {
    /// Valence of the guide ring.
    pub n: usize,
    /// Number of independent equality rows after rank filtering.
    pub independent_constraints: usize,
    /// Remaining scalar guide degrees of freedom after pinning the center.
    pub free_parameters_after_center_pin: usize,
    /// Maximum absolute residual of the fitted bi-5 tensor-border rows.
    pub max_fit_residual: f64,
    /// Root-mean-square residual of the fitted bi-5 tensor-border rows.
    pub rms_fit_residual: f64,
    /// Maximum absolute residual of the enforced equality constraints.
    pub max_constraint_residual: f64,
    /// Condition number of the KKT matrix used for the constrained solve.
    pub kkt_condition_number: f64,
}

/// Build a KP18/KP20 total-degree-6 guide by direct constrained LSQ.
pub fn solve_kp20_guide_lsq(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
) -> Result<Vec<Kp20Tri6GuideSector>, Kp20GuideLsqError> {
    solve_kp20_guide_lsq_with_report(tensor_border, star_point).map(|solved| solved.sectors)
}

/// Build a KP18/KP20 guide and return measured postconditions.
pub fn solve_kp20_guide_lsq_with_report(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
) -> Result<Kp20GuideLsqOutput, Kp20GuideLsqError> {
    solve_kp20_guide_lsq_inner(tensor_border, star_point, None, 0.0)
}

/// Build a KP18/KP20 guide with bounded least-squares prior rows.
pub fn solve_kp20_guide_lsq_with_prior(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    prior: &[Kp20Tri6GuideSector],
    prior_weight: f64,
) -> Result<Vec<Kp20Tri6GuideSector>, Kp20GuideLsqError> {
    solve_kp20_guide_lsq_inner(tensor_border, star_point, Some(prior), prior_weight)
        .map(|solved| solved.sectors)
}

fn solve_kp20_guide_lsq_inner(
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    prior: Option<&[Kp20Tri6GuideSector]>,
    prior_weight: f64,
) -> Result<Kp20GuideLsqOutput, Kp20GuideLsqError> {
    let n = tensor_border.sectors.len();
    if n < 3 || n == 4 || n > 10 {
        return Err(Kp20GuideLsqError::UnsupportedValence(n));
    }
    if !star_point.x.is_finite() || !star_point.y.is_finite() || !star_point.z.is_finite() {
        return Err(Kp20GuideLsqError::InvalidTensorBorder(
            "star point is non-finite".to_string(),
        ));
    }
    validate_tensor_border(tensor_border)?;

    let mut rows = residual_rows(tensor_border);
    if let Some(prior) = prior {
        if prior.len() != n || !prior_weight.is_finite() || prior_weight <= 0.0 {
            return Err(Kp20GuideLsqError::InvalidTensorBorder(
                "guide prior is invalid".to_string(),
            ));
        }
        append_prior_rows(&mut rows, prior, prior_weight);
    }
    let all_constraints = constraint_rows(n)?;
    let vars = n * TRI6_CP_COUNT;
    let constraints = independent_rows(vars, &all_constraints);
    let mut solved = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];

    let mut kkt_condition_number: f64 = 0.0;
    for (axis, solved_axis) in solved.iter_mut().enumerate() {
        let rhs = residual_rhs(&rows, axis);
        let fixed = axis_value(star_point, axis);
        let axis_solution = solve_axis(vars, &rows, &rhs, &constraints, fixed)
            .ok_or(Kp20GuideLsqError::SingularSystem(axis))?;
        kkt_condition_number = kkt_condition_number.max(axis_solution.kkt_condition_number);
        *solved_axis = axis_solution.values;
    }

    let mut sectors = vec![Kp20Tri6GuideSector::zero(); n];
    for s in 0..n {
        for cp in 0..TRI6_CP_COUNT {
            sectors[s].cps[cp] = Point3::new(
                solved[0][var_index(s, cp)],
                solved[1][var_index(s, cp)],
                solved[2][var_index(s, cp)],
            );
        }
        sectors[s].cps[CENTER_FLAT] = star_point;
        sectors[s].cps[flat_index(DEGREE, DEGREE, 0)] = tensor_border.sectors[s][0][0];
        sectors[s].cps[flat_index(DEGREE, 0, DEGREE)] = tensor_border.sectors[s][0][3];
    }
    let report = report(
        n,
        vars,
        &rows,
        &all_constraints,
        &solved,
        kkt_condition_number,
    );
    Ok(Kp20GuideLsqOutput { sectors, report })
}

/// Guide sectors plus measured solve report.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20GuideLsqOutput {
    /// One total-degree-6 guide sector per valence sector.
    pub sectors: Vec<Kp20Tri6GuideSector>,
    /// Measured postconditions from the direct solve.
    pub report: Kp20GuideLsqReport,
}

#[derive(Clone)]
struct LinearRow {
    terms: Vec<(usize, f64)>,
}

struct ResidualRows {
    rows: Vec<LinearRow>,
    targets: [Vec<f64>; 3],
}

fn validate_tensor_border(t: &TensorBorderBi3) -> Result<(), Kp20GuideLsqError> {
    for (s, sector) in t.sectors.iter().enumerate() {
        for (r, row) in sector.iter().enumerate() {
            for (i, p) in row.iter().enumerate() {
                if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
                    return Err(Kp20GuideLsqError::InvalidTensorBorder(format!(
                        "sector {s} row {r} CP {i} is non-finite"
                    )));
                }
            }
        }
    }
    Ok(())
}

fn residual_rows(t: &TensorBorderBi3) -> ResidualRows {
    let mut rows = Vec::new();
    let mut targets = [Vec::new(), Vec::new(), Vec::new()];
    for (s, sector) in t.sectors.iter().enumerate() {
        let raised = [
            degree_raise_1d_3_to_5(&sector[0]),
            degree_raise_1d_3_to_5(&sector[1]),
            degree_raise_1d_3_to_5(&sector[2]),
        ];
        for (q, raised_row) in raised.iter().enumerate() {
            for (j, target) in raised_row.iter().enumerate() {
                let u = j as f64 / 5.0;
                let v = 1.0 - u;
                rows.push(outer_edge_derivative_row(s, q, u, v));
                targets[0].push(target.x);
                targets[1].push(target.y);
                targets[2].push(target.z);
            }
        }
    }
    ResidualRows { rows, targets }
}

fn append_prior_rows(rows: &mut ResidualRows, prior: &[Kp20Tri6GuideSector], weight: f64) {
    for (s, sector) in prior.iter().enumerate() {
        for (cp, p) in sector.cps.iter().enumerate() {
            rows.rows.push(LinearRow {
                terms: vec![(var_index(s, cp), weight)],
            });
            rows.targets[0].push(weight * p.x);
            rows.targets[1].push(weight * p.y);
            rows.targets[2].push(weight * p.z);
        }
    }
}

fn outer_edge_derivative_row(sector: usize, order: usize, u: f64, v: f64) -> LinearRow {
    let mut terms = Vec::new();
    let row_degree = DEGREE - order;
    for a in 0..=row_degree {
        let b = row_degree - a;
        terms.push((
            var_index(sector, cp(a, b)),
            bernstein_1d(row_degree, a, u, v),
        ));
    }
    LinearRow { terms }
}

fn cp(a: usize, b: usize) -> usize {
    flat_index(DEGREE, a, b)
}

fn constraint_rows(n: usize) -> Result<Vec<LinearRow>, Kp20GuideLsqError> {
    let (_, w1, w2) =
        kp18_guide_join_weights(n).map_err(|_| Kp20GuideLsqError::UnsupportedValence(n))?;
    let mut rows = Vec::new();
    for s in 0..n {
        let r = (s + 1) % n;
        for i in 0..=6 {
            rows.push(eq(&[(r, 6 - i, i, 0, 1.0), (s, 6 - i, 0, i, -1.0)]));
        }
        for i in 1..=6 {
            rows.push(eq(&[
                (r, 6 - i, i - 1, 1, 1.0),
                (s, 6 - i, 1, i - 1, 1.0),
                (s, 6 - i, 0, i, -w1),
                (s, 6 - i + 1, 0, i - 1, -w2),
            ]));
        }
        for i in 2..=3 {
            rows.push(c2_constraint(s, r, i, n)?);
        }
        rows.push(c3_constraint(s, r, 3, n)?);
    }
    for s in 1..n {
        rows.push(LinearRow {
            terms: vec![
                (var_index(0, CENTER_FLAT), 1.0),
                (var_index(s, CENTER_FLAT), -1.0),
            ],
        });
    }
    Ok(rows)
}

fn c2_constraint(s: usize, r: usize, i: usize, n: usize) -> Result<LinearRow, Kp20GuideLsqError> {
    let (w0, w1, w2) =
        kp18_guide_join_weights(n).map_err(|_| Kp20GuideLsqError::UnsupportedValence(n))?;
    Ok(eq(&[
        (r, 6 - i, i - 2, 2, 1.0),
        (s, 6 - i, 2, i - 2, -w0 * w0),
        (s, 6 - i, 1, i - 1, -2.0 * w0 * w1),
        (s, 6 - i, 0, i, -w1 * w1),
        (s, 6 - i + 1, 1, i - 2, -2.0 * w0 * w2),
        (s, 6 - i + 1, 0, i - 1, -2.0 * w1 * w2),
        (s, 6 - i + 2, 0, i - 2, -w2 * w2),
    ]))
}

fn c3_constraint(s: usize, r: usize, i: usize, n: usize) -> Result<LinearRow, Kp20GuideLsqError> {
    let (w0, w1, w2) =
        kp18_guide_join_weights(n).map_err(|_| Kp20GuideLsqError::UnsupportedValence(n))?;
    Ok(eq(&[
        (r, 6 - i, i - 3, 3, 1.0),
        (s, 6 - i, 3, i - 3, -w0.powi(3)),
        (s, 6 - i, 2, i - 2, -3.0 * w0 * w0 * w1),
        (s, 6 - i, 1, i - 1, -3.0 * w0 * w1 * w1),
        (s, 6 - i, 0, i, -w1.powi(3)),
        (s, 6 - i + 1, 2, i - 3, -3.0 * w0 * w0 * w2),
        (s, 6 - i + 1, 1, i - 2, -6.0 * w0 * w1 * w2),
        (s, 6 - i + 1, 0, i - 1, -3.0 * w1 * w1 * w2),
        (s, 6 - i + 2, 1, i - 3, -3.0 * w0 * w2 * w2),
        (s, 6 - i + 2, 0, i - 2, -3.0 * w1 * w2 * w2),
        (s, 6 - i + 3, 0, i - 3, -w2.powi(3)),
    ]))
}

fn eq(terms: &[(usize, usize, usize, usize, f64)]) -> LinearRow {
    LinearRow {
        terms: terms
            .iter()
            .map(|&(s, _i, j, k, w)| (var_index(s, flat_index(DEGREE, j, k)), w))
            .collect(),
    }
}

struct AxisSolution {
    values: Vec<f64>,
    kkt_condition_number: f64,
}

fn solve_axis(
    vars: usize,
    residuals: &ResidualRows,
    rhs: &[f64],
    constraints: &[LinearRow],
    center: f64,
) -> Option<AxisSolution> {
    let size = vars + constraints.len() + 1;
    let mut a = vec![vec![0.0_f64; size]; size];
    let mut b = vec![0.0_f64; size];
    for (row, &target) in residuals.rows.iter().zip(rhs) {
        accumulate_normal_row(&mut a, &mut b, &row.terms, target, 1.0);
    }
    for (i, row) in a.iter_mut().enumerate().take(vars) {
        row[i] += 1.0e-12;
    }
    add_lagrange_row(
        vars,
        &mut a,
        &mut b,
        &[(var_index(0, CENTER_FLAT), 1.0)],
        center,
    );
    for (cidx, row) in constraints.iter().enumerate() {
        add_lagrange_row(vars + 1 + cidx, &mut a, &mut b, &row.terms, 0.0);
    }
    let kkt_condition_number = condition_number_dense(&a)?.condition_number;
    solve_dense(a, b).map(|x| AxisSolution {
        values: x[..vars].to_vec(),
        kkt_condition_number,
    })
}

fn accumulate_normal_row(
    a: &mut [Vec<f64>],
    b: &mut [f64],
    terms: &[(usize, f64)],
    rhs: f64,
    weight: f64,
) {
    for &(i, wi) in terms {
        b[i] += weight * wi * rhs;
        for &(j, wj) in terms {
            a[i][j] += weight * wi * wj;
        }
    }
}

fn add_lagrange_row(
    row_idx: usize,
    a: &mut [Vec<f64>],
    b: &mut [f64],
    terms: &[(usize, f64)],
    rhs: f64,
) {
    b[row_idx] = rhs;
    for &(var, weight) in terms {
        a[var][row_idx] += weight;
        a[row_idx][var] += weight;
    }
}

fn independent_rows(vars: usize, rows: &[LinearRow]) -> Vec<LinearRow> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut keep = Vec::new();
    for row in rows {
        let mut dense = vec![0.0_f64; vars];
        for &(i, w) in &row.terms {
            dense[i] += w;
        }
        for b in &basis {
            let dot = dense.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            for (x, bi) in dense.iter_mut().zip(b) {
                *x -= dot * bi;
            }
        }
        let norm = dense.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1.0e-10 {
            for x in &mut dense {
                *x /= norm;
            }
            basis.push(dense);
            keep.push(row.clone());
        }
    }
    keep
}

fn residual_rhs(rows: &ResidualRows, axis: usize) -> Vec<f64> {
    rows.targets[axis].clone()
}

fn report(
    n: usize,
    vars: usize,
    rows: &ResidualRows,
    constraints: &[LinearRow],
    solved: &[Vec<f64>; 3],
    kkt_condition_number: f64,
) -> Kp20GuideLsqReport {
    let mut max_fit: f64 = 0.0;
    let mut sum_fit_sq = 0.0;
    let mut count = 0_usize;
    for (row_idx, row) in rows.rows.iter().enumerate() {
        for (axis, solved_axis) in solved.iter().enumerate() {
            let residual = row_value(row, solved_axis) - rows.targets[axis][row_idx];
            max_fit = max_fit.max(residual.abs());
            sum_fit_sq += residual * residual;
            count += 1;
        }
    }
    let mut max_constraint: f64 = 0.0;
    for row in constraints {
        for solved_axis in solved {
            max_constraint = max_constraint.max(row_value(row, solved_axis).abs());
        }
    }
    let independent = independent_rows(vars, constraints).len();
    Kp20GuideLsqReport {
        n,
        independent_constraints: independent,
        free_parameters_after_center_pin: vars - independent - 1,
        max_fit_residual: max_fit,
        rms_fit_residual: (sum_fit_sq / count as f64).sqrt(),
        max_constraint_residual: max_constraint,
        kkt_condition_number,
    }
}

fn row_value(row: &LinearRow, values: &[f64]) -> f64 {
    row.terms.iter().map(|&(i, w)| w * values[i]).sum()
}

fn var_index(sector: usize, cp: usize) -> usize {
    sector * TRI6_CP_COUNT + cp
}

fn axis_value(p: Point3, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        2 => p.z,
        _ => unreachable!(),
    }
}

fn bernstein_1d(n: usize, i: usize, u: f64, v: f64) -> f64 {
    binom(n, i) as f64 * u.powi(i as i32) * v.powi((n - i) as i32)
}

const fn binom(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else if k == 0 || k == n {
        1
    } else {
        binom(n - 1, k - 1) + binom(n - 1, k)
    }
}

#[cfg(test)]
#[path = "kp20_bi6_guide_lsq_tests.rs"]
mod tests;
