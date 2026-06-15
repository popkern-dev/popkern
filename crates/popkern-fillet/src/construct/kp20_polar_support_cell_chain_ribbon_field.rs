//! Degree-8 transverse ribbon fields along chain curves.
//!
//! A ribbon field `R(s)` spans, together with the chain curve
//! tangent `C'(s)`, the shared tangent plane of every face along
//! the chain. Row-1 controls are later DERIVED from `(C, R)`, so
//! the field's only requirements are: exact in-plane direction
//! at the nodes (where transverse chains must agree) and a
//! smooth transverse direction elsewhere. Degree 8 keeps the
//! derived product `β(t)·R(t)` (β linear) at degree 9 — exactly
//! representable in a bi-9 row.

use popkern_core::geometry::{Point3, Vec3};

use super::chain_solver::{add_third_difference_penalty, solve_dense};

/// Tikhonov weight for the field fit.
const RIBBON_FIELD_REGULARIZATION: f64 = 1.0e-9;

/// Oscillation (third-difference) penalty weight for the field
/// fit (see the chain-curve solver's smoothing weight).
const RIBBON_FIELD_SMOOTHING_WEIGHT: f64 = 0.5;

/// Maximum interpolated nodes: degree 8 has 9 coefficients per
/// axis; pinning more leaves no least-squares freedom.
const MAX_FIELD_NODES: usize = 8;

/// Weight of soft-heavy direction rows (authoritative chained
/// legs): near-exact where consistent, without consuming the
/// hard-constraint budget — node-rich sub-chains were silently
/// dropping leg constraints past the budget, leaving those rows'
/// frozen endpoints out of span.
const SOFT_DIRECTION_WEIGHT: f64 = 1.0e6;

/// Fit a degree-8 vector field over local `[0, 1]`: interpolate
/// `(param, direction)` node constraints exactly, weight the
/// `soft_directions` heavily, and least-squares the sampled
/// directions elsewhere. Returns the 9 control vectors.
pub(super) fn fit_ribbon_field(
    samples: &[(f64, Vec3)],
    nodes: &[(f64, Vec3)],
    soft_directions: &[(f64, Vec3)],
) -> Result<[Vec3; 9], &'static str> {
    if samples.len() < 9 {
        return Err("ribbon_field_sample_count_insufficient");
    }
    if nodes.len() > MAX_FIELD_NODES {
        return Err("ribbon_field_node_count_exceeds_degree");
    }
    let k = nodes.len();
    let n = 9 + k;
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = [vec![0.0_f64; n], vec![0.0_f64; n], vec![0.0_f64; n]];
    for &(t, direction) in samples {
        let basis = bernstein8(t);
        for i in 0..9 {
            for j in 0..9 {
                matrix[i][j] += basis[i] * basis[j];
            }
            rhs[0][i] += basis[i] * direction.x;
            rhs[1][i] += basis[i] * direction.y;
            rhs[2][i] += basis[i] * direction.z;
        }
    }
    for &(t, direction) in soft_directions {
        let basis = bernstein8(t);
        for (i, matrix_row) in matrix.iter_mut().enumerate().take(9) {
            for (j, entry) in matrix_row.iter_mut().enumerate().take(9) {
                *entry += SOFT_DIRECTION_WEIGHT * basis[i] * basis[j];
            }
            rhs[0][i] += SOFT_DIRECTION_WEIGHT * basis[i] * direction.x;
            rhs[1][i] += SOFT_DIRECTION_WEIGHT * basis[i] * direction.y;
            rhs[2][i] += SOFT_DIRECTION_WEIGHT * basis[i] * direction.z;
        }
    }
    for (i, matrix_row) in matrix.iter_mut().enumerate().take(9) {
        matrix_row[i] += RIBBON_FIELD_REGULARIZATION;
    }
    add_third_difference_penalty(&mut matrix, 0, 9, RIBBON_FIELD_SMOOTHING_WEIGHT);
    for (index, &(t, direction)) in nodes.iter().enumerate() {
        let basis = bernstein8(t);
        let row = 9 + index;
        for (i, &weight) in basis.iter().enumerate() {
            matrix[row][i] = weight;
            matrix[i][row] = weight;
        }
        rhs[0][row] = direction.x;
        rhs[1][row] = direction.y;
        rhs[2][row] = direction.z;
    }
    let xs = solve_dense(matrix.clone(), rhs[0].clone()).ok_or("ribbon_field_fit_singular")?;
    let ys = solve_dense(matrix.clone(), rhs[1].clone()).ok_or("ribbon_field_fit_singular")?;
    let zs = solve_dense(matrix, rhs[2].clone()).ok_or("ribbon_field_fit_singular")?;
    let mut out = [Vec3::new(0.0, 0.0, 0.0); 9];
    for (index, vector) in out.iter_mut().enumerate() {
        *vector = Vec3::new(xs[index], ys[index], zs[index]);
    }
    Ok(out)
}

/// Exact restriction of a degree-8 field to `[start, end]`
/// (oriented start → end), via two de Casteljau splits.
pub(super) fn ribbon_field_segment(field: &[Vec3; 9], start: f64, end: f64) -> [Vec3; 9] {
    let lo = start.min(end);
    let hi = start.max(end);
    let (_, tail) = de_casteljau_split_field(field, lo);
    let local = if (1.0 - lo).abs() < 1.0e-12 {
        0.0
    } else {
        (hi - lo) / (1.0 - lo)
    };
    let (mut segment, _) = de_casteljau_split_field(&tail, local);
    if start > end {
        segment.reverse();
    }
    segment
}

fn de_casteljau_split_field(field: &[Vec3; 9], t: f64) -> ([Vec3; 9], [Vec3; 9]) {
    let mut left = [Vec3::new(0.0, 0.0, 0.0); 9];
    let mut right = [Vec3::new(0.0, 0.0, 0.0); 9];
    let mut working = *field;
    left[0] = working[0];
    right[8] = working[8];
    for level in 1..9 {
        for index in 0..(9 - level) {
            working[index] = working[index] + (working[index + 1] - working[index]) * t;
        }
        left[level] = working[0];
        right[8 - level] = working[8 - level];
    }
    (left, right)
}

/// Evaluate the degree-8 field at `t`.
pub(super) fn eval_ribbon_field(field: &[Vec3; 9], t: f64) -> Vec3 {
    let basis = bernstein8(t);
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for (index, weight) in basis.iter().enumerate() {
        out = out + field[index] * *weight;
    }
    out
}

pub(super) fn bernstein8(t: f64) -> [f64; 9] {
    let omt = 1.0 - t;
    let coefficients = [1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0];
    let mut out = [0.0; 9];
    for (i, value) in out.iter_mut().enumerate() {
        *value = coefficients[i] * t.powi(i as i32) * omt.powi((8 - i) as i32);
    }
    out
}

/// Degree-8 Bezier derivative coefficients of a degree-9 row:
/// `d_i = 9 (P_{i+1} - P_i)`.
pub(super) fn row_derivative_coefficients(row: &[Point3; 10]) -> [Vec3; 9] {
    let mut out = [Vec3::new(0.0, 0.0, 0.0); 9];
    for (index, vector) in out.iter_mut().enumerate() {
        *vector = (row[index + 1] - row[index]) * 9.0;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_fit_interpolates_nodes_and_tracks_samples() {
        let mut samples = Vec::new();
        for index in 0..30 {
            let t = index as f64 / 29.0;
            samples.push((t, Vec3::new(1.0, t, t * t)));
        }
        let nodes = vec![
            (0.0, Vec3::new(1.0, 0.0, 0.0)),
            (0.5, Vec3::new(1.0, 0.5, 0.25)),
            (1.0, Vec3::new(1.0, 1.0, 1.0)),
        ];
        let field = fit_ribbon_field(&samples, &nodes, &[]).expect("fit must solve");
        for &(t, direction) in &nodes {
            let actual = eval_ribbon_field(&field, t);
            assert!(
                (actual - direction).dot(&(actual - direction)).sqrt() < 1.0e-9,
                "node not interpolated at t={t}"
            );
        }
        let mid = eval_ribbon_field(&field, 0.25);
        let expected = Vec3::new(1.0, 0.25, 0.0625);
        assert!(
            (mid - expected).dot(&(mid - expected)).sqrt() < 1.0e-6,
            "field should track the quadratic samples, got {mid:?}"
        );
    }

    #[test]
    fn field_segment_matches_source() {
        let mut field = [Vec3::new(0.0, 0.0, 0.0); 9];
        for (index, vector) in field.iter_mut().enumerate() {
            let t = index as f64;
            *vector = Vec3::new(t, t * t * 0.1, -t * 0.2);
        }
        let segment = ribbon_field_segment(&field, 0.25, 0.75);
        for sample in 0..=8 {
            let local = sample as f64 / 8.0;
            let global = 0.25 + 0.5 * local;
            let expected = eval_ribbon_field(&field, global);
            let actual = eval_ribbon_field(&segment, local);
            assert!(
                (expected - actual).dot(&(expected - actual)).sqrt() < 1.0e-12,
                "segment diverged at local {local}"
            );
        }
    }
}
