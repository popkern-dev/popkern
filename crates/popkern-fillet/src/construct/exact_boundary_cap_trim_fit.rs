//! Fixed-knot interpolation helpers for exact-boundary cap trim restoration.

use popkern_core::geometry::Point3;

pub(super) fn greville_parameters(
    expected: usize,
    degree: usize,
    knots: &[f64],
) -> Option<Vec<f64>> {
    if expected < 2 || knots.len() != expected + degree + 1 {
        return None;
    }
    let first = knots[degree];
    let last = knots[expected];
    let mut params = Vec::with_capacity(expected);
    for index in 0..expected {
        if index == 0 {
            params.push(first);
        } else if index + 1 == expected {
            params.push(last);
        } else {
            let sum: f64 = knots[(index + 1)..=(index + degree)].iter().sum();
            params.push(sum / degree as f64);
        }
    }
    Some(params)
}

#[allow(clippy::needless_range_loop)]
pub(super) fn interpolate_points_on_fixed_knots(
    points: &[Point3],
    params: &[f64],
    degree: usize,
    knots: &[f64],
) -> Option<Vec<Point3>> {
    let n = points.len();
    if params.len() != n || knots.len() != n + degree + 1 {
        return None;
    }
    let mut matrix = vec![vec![0.0; n]; n];
    for (row, &parameter) in params.iter().enumerate() {
        for col in 0..n {
            matrix[row][col] = basis_function(col, degree, parameter, knots);
        }
    }
    let x_values = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.x).collect::<Vec<_>>(),
    )?;
    let y_values = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.y).collect::<Vec<_>>(),
    )?;
    let z_values = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.z).collect::<Vec<_>>(),
    )?;
    Some(
        (0..n)
            .map(|index| Point3::new(x_values[index], y_values[index], z_values[index]))
            .collect(),
    )
}

fn basis_function(index: usize, degree: usize, parameter: f64, knots: &[f64]) -> f64 {
    if degree == 0 {
        if knots[index] <= parameter && parameter < knots[index + 1] {
            return 1.0;
        }
        if (parameter - knots[index + 1]).abs() < 1.0e-14 && knots[index] < knots[index + 1] {
            let last_span = ((index + 1)..knots.len().saturating_sub(1))
                .all(|knot_index| knots[knot_index] >= knots[knot_index + 1]);
            if last_span {
                return 1.0;
            }
        }
        return 0.0;
    }
    let mut result = 0.0;
    let left_denominator = knots[index + degree] - knots[index];
    if left_denominator.abs() > 1.0e-14 {
        result += ((parameter - knots[index]) / left_denominator)
            * basis_function(index, degree - 1, parameter, knots);
    }
    let right_denominator = knots[index + degree + 1] - knots[index + 1];
    if right_denominator.abs() > 1.0e-14 {
        result += ((knots[index + degree + 1] - parameter) / right_denominator)
            * basis_function(index + 1, degree - 1, parameter, knots);
    }
    result
}

#[allow(clippy::needless_range_loop)]
fn solve_linear_system(matrix: &[Vec<f64>], values: &[f64]) -> Option<Vec<f64>> {
    let n = matrix.len();
    if n == 0 || values.len() != n {
        return None;
    }
    let mut augmented = Vec::with_capacity(n);
    for row in 0..n {
        if matrix[row].len() != n {
            return None;
        }
        let mut augmented_row = matrix[row].clone();
        augmented_row.push(values[row]);
        augmented.push(augmented_row);
    }
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_value = augmented[col][col].abs();
        for row in (col + 1)..n {
            let candidate = augmented[row][col].abs();
            if candidate > pivot_value {
                pivot_value = candidate;
                pivot_row = row;
            }
        }
        if pivot_value < 1.0e-14 {
            return None;
        }
        if pivot_row != col {
            augmented.swap(col, pivot_row);
        }
        let pivot = augmented[col][col];
        for row in (col + 1)..n {
            let factor = augmented[row][col] / pivot;
            for entry in col..=n {
                augmented[row][entry] -= factor * augmented[col][entry];
            }
        }
    }
    let mut solution = vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = augmented[row][n];
        for col in (row + 1)..n {
            sum -= augmented[row][col] * solution[col];
        }
        solution[row] = sum / augmented[row][row];
    }
    Some(solution)
}

pub(super) fn knots_equivalent(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(left, right)| (*left - *right).abs() <= 1.0e-9)
}
