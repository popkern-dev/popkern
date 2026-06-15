//! Constrained degree-9 curve solver for chain-curve authoring.
//!
//! Solves one canonical curve per sub-chain: least-squares
//! against sampled row geometry, subject to exact interpolation
//! of pinned node positions (hard KKT rows) and tangent-plane
//! constraints `C'(x) · N = 0` at nodes as SOFT-HEAVY objective
//! terms. Tangent planes used to be hard rows truncated to a
//! 28-row budget; on 8-9-node chains that silently dropped most
//! of them, leaving chained transverse legs up to 30% out of
//! their node planes — the dominant corner-zone G1 noise. Soft
//! tangent rows always apply (every node, no budget) at ~1e-6
//! relative accuracy, far inside the 1.25e-3 G1 gate. Position
//! constraints are axis-separable but tangent-plane terms couple
//! the axes, so the solve uses one system over all 30 control
//! coordinates.

use popkern_core::geometry::{Point3, Vec3};

use super::row_fit::bernstein9;

/// Tikhonov weight keeping the fit well-posed when samples are
/// sparse relative to the curve's degrees of freedom.
pub(super) const FIT_REGULARIZATION: f64 = 1.0e-9;

/// Oscillation penalty weight on the control polygon's THIRD
/// differences. Degree-9 least-squares through rough geometry
/// with many hard interpolation rows oscillates, and the
/// oscillation inflates the derivative coefficients that drive
/// the twist demands of adjacent boundaries apart. Third
/// differences vanish exactly for (degree-elevated) constant,
/// linear, and quadratic geometry, so genuine arc curvature is
/// unbiased — only the high-frequency wiggle is damped. Node
/// interpolation and tangent-plane rows stay HARD, so exact-G0
/// row identity and node coplanarity are unaffected.
pub(super) const FIT_SMOOTHING_WEIGHT: f64 = 2.0;

/// Weight on soft tangent-plane rows. Heavy enough that the
/// residual `C'(x)·N` is ~1e-6 of the geometry scale, light
/// enough not to swamp the position KKT pivots.
pub(super) const SOFT_TANGENT_PLANE_WEIGHT: f64 = 1.0e6;

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
    /// (`C'(x) · N = 0`), enforced as a soft-heavy objective
    /// term. The normal's sign is irrelevant.
    pub(super) tangent_planes: Vec<(f64, Vec3)>,
    /// `(cp index, target position, weight)` — soft position
    /// targets for individual control points, used to decompress
    /// end control spacing measured too thin on a first solve
    /// (magnitude retarget along the solved direction; the curve
    /// itself barely moves, only its speed profile).
    pub(super) spacing_targets: Vec<(usize, Point3, f64)>,
    /// `(local param, target derivative, weight)` — soft rows on
    /// `C'(t)`, used to decompress the speed profile at INTERIOR
    /// node parameters (T-junction cuts), which end-CP spacing
    /// targets cannot reach. Restriction rows inherit their leg
    /// magnitudes from `C'` at the cut; a compressed cut thins
    /// the crossing rows' frozen legs until their derived bands
    /// turn tangential and trim normals become noise.
    pub(super) derivative_targets: Vec<(f64, Vec3, f64)>,
}

impl ChainCurveConstraints {
    pub(super) fn row_count(&self) -> usize {
        3 * self.positions.len()
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
    for (i, matrix_row) in matrix.iter_mut().enumerate().take(30) {
        matrix_row[i] += FIT_REGULARIZATION;
    }
    for &(index, target, weight) in &constraints.spacing_targets {
        matrix[index][index] += weight;
        matrix[10 + index][10 + index] += weight;
        matrix[20 + index][20 + index] += weight;
        rhs[index] += weight * target.x;
        rhs[10 + index] += weight * target.y;
        rhs[20 + index] += weight * target.z;
    }
    for &(t, target, weight) in &constraints.derivative_targets {
        let dbasis = dbernstein9(t);
        for (axis, value) in [(0usize, target.x), (1, target.y), (2, target.z)] {
            for i in 0..10 {
                if dbasis[i] == 0.0 {
                    continue;
                }
                for j in 0..10 {
                    matrix[axis * 10 + i][axis * 10 + j] += weight * dbasis[i] * dbasis[j];
                }
                rhs[axis * 10 + i] += weight * dbasis[i] * value;
            }
        }
    }
    add_third_difference_penalty(&mut matrix, 0, 10, FIT_SMOOTHING_WEIGHT);
    add_third_difference_penalty(&mut matrix, 10, 10, FIT_SMOOTHING_WEIGHT);
    add_third_difference_penalty(&mut matrix, 20, 10, FIT_SMOOTHING_WEIGHT);
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
    debug_assert_eq!(row, n);
    // Soft-heavy tangent-plane terms: `weight · g gᵀ` over the 30
    // coefficient coordinates, target zero (`C'(x)·N = 0`).
    for &(t, normal) in &constraints.tangent_planes {
        let dbasis = dbernstein9(t);
        let mut gradient = [0.0_f64; 30];
        for i in 0..10 {
            gradient[i] = dbasis[i] * normal.x;
            gradient[10 + i] = dbasis[i] * normal.y;
            gradient[20 + i] = dbasis[i] * normal.z;
        }
        for i in 0..30 {
            if gradient[i] == 0.0 {
                continue;
            }
            for j in 0..30 {
                matrix[i][j] += SOFT_TANGENT_PLANE_WEIGHT * gradient[i] * gradient[j];
            }
        }
    }
    let solution = solve_dense(matrix, rhs).ok_or("sub_chain_fit_singular")?;
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = Point3::new(solution[index], solution[10 + index], solution[20 + index]);
    }
    Ok(out)
}

/// Add `weight * D^T D` for the third-difference operator `D`
/// over one axis block of `count` coefficients starting at
/// `offset` (penalty target is zero third difference; the
/// right-hand side is untouched, so polynomials up to quadratic
/// are unbiased).
pub(super) fn add_third_difference_penalty(
    matrix: &mut [Vec<f64>],
    offset: usize,
    count: usize,
    weight: f64,
) {
    for i in 0..count.saturating_sub(3) {
        let indices = [offset + i, offset + i + 1, offset + i + 2, offset + i + 3];
        let coefficients = [1.0, -3.0, 3.0, -1.0];
        for (row_position, &row_index) in indices.iter().enumerate() {
            for (column_position, &column_index) in indices.iter().enumerate() {
                matrix[row_index][column_index] +=
                    weight * coefficients[row_position] * coefficients[column_position];
            }
        }
    }
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

#[allow(clippy::needless_range_loop)] // Gaussian elimination is explicit row/column indexing.
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
