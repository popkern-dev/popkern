//! Constrained degree-9 curve solver for chain-curve authoring.
//!
//! Solves one canonical curve per sub-chain: least-squares
//! against sampled row geometry, subject to exact interpolation
//! of pinned node positions and, when requested, tangent-plane
//! constraints `C'(x) · N = 0` at nodes. Position constraints
//! are axis-separable but tangent-plane constraints couple the
//! axes, so the solve uses one KKT system over all 30 control
//! coordinates.

use popkern_core::geometry::{Point3, Vec3};

use super::row_fit::bernstein9;

/// Tikhonov weight keeping the fit well-posed when samples are
/// sparse relative to the curve's degrees of freedom.
const FIT_REGULARIZATION: f64 = 1.0e-9;

/// Pivot floor below which the KKT system is declared singular.
const SOLVE_PIVOT_TOL: f64 = 1.0e-13;

/// Hard constraints for one sub-chain curve.
#[derive(Debug, Clone, Default)]
pub(super) struct ChainCurveConstraints {
    /// `(local param, position)` — the curve interpolates the
    /// pinned node position exactly.
    pub(super) positions: Vec<(f64, Point3)>,
    /// `(local param, plane normal)` — the curve tangent at the
    /// node lies in the plane perpendicular to the normal
    /// (`C'(x) · N = 0`). The normal's sign is irrelevant.
    pub(super) tangent_planes: Vec<(f64, Vec3)>,
}

impl ChainCurveConstraints {
    pub(super) fn row_count(&self) -> usize {
        3 * self.positions.len() + self.tangent_planes.len()
    }
}

/// Fit a degree-9 Bezier curve to `samples` under `constraints`.
/// Unknown layout: `[x_0..x_9, y_0..y_9, z_0..z_9]`.
pub(super) fn fit_constrained_degree9(
    samples: &[(f64, Point3)],
    constraints: &ChainCurveConstraints,
) -> Result<[Point3; 10], &'static str> {
    if samples.len() < 10 {
        return Err("sub_chain_sample_count_insufficient");
    }
    let constraint_rows = constraints.row_count();
    let n = 30 + constraint_rows;
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    // Objective normal block: A^T A per axis (identical blocks),
    // plus Tikhonov regularization.
    for &(t, target) in samples {
        let basis = bernstein9(t);
        for i in 0..10 {
            for j in 0..10 {
                let product = basis[i] * basis[j];
                matrix[i][j] += product;
                matrix[10 + i][10 + j] += product;
                matrix[20 + i][20 + j] += product;
            }
            rhs[i] += basis[i] * target.x;
            rhs[10 + i] += basis[i] * target.y;
            rhs[20 + i] += basis[i] * target.z;
        }
    }
    for i in 0..30 {
        matrix[i][i] += FIT_REGULARIZATION;
    }
    // Constraint rows (and their transposes in the KKT system).
    let mut row = 30;
    for &(t, position) in &constraints.positions {
        let basis = bernstein9(t);
        for (axis, value) in [(0, position.x), (1, position.y), (2, position.z)] {
            for i in 0..10 {
                matrix[row][axis * 10 + i] = basis[i];
                matrix[axis * 10 + i][row] = basis[i];
            }
            rhs[row] = value;
            row += 1;
        }
    }
    for &(t, normal) in &constraints.tangent_planes {
        let dbasis = dbernstein9(t);
        for i in 0..10 {
            matrix[row][i] = dbasis[i] * normal.x;
            matrix[i][row] = dbasis[i] * normal.x;
            matrix[row][10 + i] = dbasis[i] * normal.y;
            matrix[10 + i][row] = dbasis[i] * normal.y;
            matrix[row][20 + i] = dbasis[i] * normal.z;
            matrix[20 + i][row] = dbasis[i] * normal.z;
        }
        rhs[row] = 0.0;
        row += 1;
    }
    let solution = solve_dense(matrix, rhs).ok_or("sub_chain_fit_singular")?;
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = Point3::new(solution[index], solution[10 + index], solution[20 + index]);
    }
    Ok(out)
}

/// Derivative of the degree-9 Bernstein basis:
/// `B'_{i,9}(t) = 9 (B_{i-1,8}(t) - B_{i,8}(t))`.
pub(super) fn dbernstein9(t: f64) -> [f64; 10] {
    let b8 = bernstein8(t);
    let mut out = [0.0; 10];
    for (i, value) in out.iter_mut().enumerate() {
        let lower = if i > 0 { b8[i - 1] } else { 0.0 };
        let upper = if i < 9 { b8[i] } else { 0.0 };
        *value = 9.0 * (lower - upper);
    }
    out
}

fn bernstein8(t: f64) -> [f64; 9] {
    let omt = 1.0 - t;
    let coefficients = [1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0];
    let mut out = [0.0; 9];
    for (i, value) in out.iter_mut().enumerate() {
        *value = coefficients[i] * t.powi(i as i32) * omt.powi((8 - i) as i32);
    }
    out
}

pub(super) fn solve_dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for pivot in 0..n {
        let mut best = pivot;
        for row in pivot + 1..n {
            if matrix[row][pivot].abs() > matrix[best][pivot].abs() {
                best = row;
            }
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let scale = matrix[pivot][pivot];
        if scale.abs() <= SOLVE_PIVOT_TOL {
            return None;
        }
        for column in pivot..n {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..n {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for column in pivot..n {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    Some(rhs)
}
