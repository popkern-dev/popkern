//! Derivative-constrained support-image fitting for boundary reauthoring.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use crate::apex_construction::exact_boundary_cap::pcurves::PreparedExactBoundaryPcurvePair;

const DERIVATIVE_ROW_WEIGHT: f64 = 32.0;
const ENDPOINT_ROW_WEIGHT: f64 = 64.0;
const TANGENT_PLANE_ROW_WEIGHT: f64 = 128.0;

pub(super) fn tangent_plane_constrained_sampled_curve(
    support_pair: &PreparedExactBoundaryPcurvePair,
    support_forward: bool,
    sample_count: usize,
) -> Option<NurbsCurve> {
    if sample_count < 2 {
        return None;
    }
    let degree = 5usize.min(sample_count.saturating_sub(1));
    let control_count = sample_count;
    let knots = open_uniform_knots(degree, control_count);
    let unknown_count = control_count * 3;
    let mut matrix = vec![vec![0.0; unknown_count]; unknown_count];
    let mut values = vec![0.0; unknown_count];

    for index in 0..sample_count {
        let parameter = index as f64 / (sample_count - 1) as f64;
        let (point, _) = support_image_sample(support_pair, parameter, support_forward)?;
        let endpoint_weight = if index == 0 || index + 1 == sample_count {
            ENDPOINT_ROW_WEIGHT
        } else {
            1.0
        };
        let basis = basis_row(degree, &knots, control_count, parameter);
        add_coordinate_row(
            &mut matrix,
            &mut values,
            &basis,
            0,
            endpoint_weight,
            point.x,
        );
        add_coordinate_row(
            &mut matrix,
            &mut values,
            &basis,
            1,
            endpoint_weight,
            point.y,
        );
        add_coordinate_row(
            &mut matrix,
            &mut values,
            &basis,
            2,
            endpoint_weight,
            point.z,
        );
        let normal = support_image_normal(support_pair, parameter, support_forward)?;
        let derivative_basis = derivative_basis_row(degree, &knots, control_count, parameter);
        add_tangent_plane_row(
            &mut matrix,
            &mut values,
            &derivative_basis,
            normal,
            TANGENT_PLANE_ROW_WEIGHT,
        );
    }

    let solution = solve_linear_system(&matrix, &values)?;
    let mut control_points = (0..control_count)
        .map(|index| {
            Point3::new(
                solution[index * 3],
                solution[index * 3 + 1],
                solution[index * 3 + 2],
            )
        })
        .collect::<Vec<_>>();
    control_points[0] = support_image_sample(support_pair, 0.0, support_forward)?.0;
    control_points[control_count - 1] = support_image_sample(support_pair, 1.0, support_forward)?.0;
    Some(NurbsCurve {
        degree,
        knots,
        control_points,
        weights: vec![1.0; control_count],
    })
}

pub(super) fn derivative_constrained_sampled_curve(
    support_pair: &PreparedExactBoundaryPcurvePair,
    support_forward: bool,
    sample_count: usize,
) -> Option<NurbsCurve> {
    if sample_count < 2 {
        return None;
    }
    let degree = 5usize.min(sample_count.saturating_sub(1));
    let control_count = sample_count.clamp(degree + 1, 256);
    let knots = open_uniform_knots(degree, control_count);
    let mut matrix = vec![vec![0.0; control_count]; control_count];
    let mut x_values = vec![0.0; control_count];
    let mut y_values = vec![0.0; control_count];
    let mut z_values = vec![0.0; control_count];

    for index in 0..sample_count {
        let parameter = index as f64 / (sample_count - 1) as f64;
        let (point, tangent) = support_image_sample(support_pair, parameter, support_forward)?;
        let endpoint_weight = if index == 0 || index + 1 == sample_count {
            ENDPOINT_ROW_WEIGHT
        } else {
            1.0
        };
        let basis = basis_row(degree, &knots, control_count, parameter);
        add_row(
            &mut matrix,
            &mut x_values,
            &mut y_values,
            &mut z_values,
            &basis,
            endpoint_weight,
            Vec3::new(point.x, point.y, point.z),
        );
        let derivative_basis = derivative_basis_row(degree, &knots, control_count, parameter);
        add_row(
            &mut matrix,
            &mut x_values,
            &mut y_values,
            &mut z_values,
            &derivative_basis,
            DERIVATIVE_ROW_WEIGHT,
            tangent,
        );
    }

    let x = solve_linear_system(&matrix, &x_values)?;
    let y = solve_linear_system(&matrix, &y_values)?;
    let z = solve_linear_system(&matrix, &z_values)?;
    let mut control_points = (0..control_count)
        .map(|index| Point3::new(x[index], y[index], z[index]))
        .collect::<Vec<_>>();
    control_points[0] = support_image_sample(support_pair, 0.0, support_forward)?.0;
    control_points[control_count - 1] = support_image_sample(support_pair, 1.0, support_forward)?.0;
    Some(NurbsCurve {
        degree,
        knots,
        control_points,
        weights: vec![1.0; control_count],
    })
}

pub(super) fn interpolated_sampled_curve(points: &[Point3], params: &[f64]) -> Option<NurbsCurve> {
    let degree = 5.min(points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(points, params, degree).ok()
}

fn support_image_sample(
    support_pair: &PreparedExactBoundaryPcurvePair,
    parameter: f64,
    support_forward: bool,
) -> Option<(Point3, Vec3)> {
    let support_parameter = if support_forward {
        parameter
    } else {
        1.0 - parameter
    };
    let (start, end) = support_pair.support_pcurve.parameter_range();
    let pcurve_t = start + (end - start) * support_parameter;
    let uv = support_pair.support_pcurve.point_at(pcurve_t);
    let uv_tangent = support_pair.support_pcurve.tangent_at(pcurve_t);
    let direction_sign = if support_forward { 1.0 } else { -1.0 };
    let tangent = support_pair.support_surface.partial_u(uv.x, uv.y) * uv_tangent.x
        + support_pair.support_surface.partial_v(uv.x, uv.y) * uv_tangent.y;
    Some((
        support_pair.support_surface.point_at(uv.x, uv.y),
        tangent * direction_sign * (end - start).abs(),
    ))
}

fn support_image_normal(
    support_pair: &PreparedExactBoundaryPcurvePair,
    parameter: f64,
    support_forward: bool,
) -> Option<Vec3> {
    let support_parameter = if support_forward {
        parameter
    } else {
        1.0 - parameter
    };
    let (start, end) = support_pair.support_pcurve.parameter_range();
    let pcurve_t = start + (end - start) * support_parameter;
    let uv = support_pair.support_pcurve.point_at(pcurve_t);
    support_pair
        .support_surface
        .partial_u(uv.x, uv.y)
        .cross(&support_pair.support_surface.partial_v(uv.x, uv.y))
        .normalized()
}

fn add_coordinate_row(
    matrix: &mut [Vec<f64>],
    values: &mut [f64],
    basis: &[f64],
    coordinate: usize,
    weight: f64,
    target: f64,
) {
    let row = basis
        .iter()
        .enumerate()
        .filter_map(|(index, basis_value)| {
            (basis_value.abs() > f64::EPSILON).then_some((index * 3 + coordinate, *basis_value))
        })
        .collect::<Vec<_>>();
    add_coupled_row(matrix, values, &row, weight, target);
}

fn add_tangent_plane_row(
    matrix: &mut [Vec<f64>],
    values: &mut [f64],
    derivative_basis: &[f64],
    normal: Vec3,
    weight: f64,
) {
    let mut row = Vec::new();
    for (index, basis_value) in derivative_basis.iter().enumerate() {
        if basis_value.abs() <= f64::EPSILON {
            continue;
        }
        row.push((index * 3, *basis_value * normal.x));
        row.push((index * 3 + 1, *basis_value * normal.y));
        row.push((index * 3 + 2, *basis_value * normal.z));
    }
    add_coupled_row(matrix, values, &row, weight, 0.0);
}

fn add_coupled_row(
    matrix: &mut [Vec<f64>],
    values: &mut [f64],
    row: &[(usize, f64)],
    weight: f64,
    target: f64,
) {
    for &(row_index, row_value) in row {
        values[row_index] += weight * row_value * target;
        for &(column_index, column_value) in row {
            matrix[row_index][column_index] += weight * row_value * column_value;
        }
    }
}

fn add_row(
    matrix: &mut [Vec<f64>],
    x_values: &mut [f64],
    y_values: &mut [f64],
    z_values: &mut [f64],
    row: &[f64],
    weight: f64,
    target: Vec3,
) {
    for r in 0..row.len() {
        x_values[r] += weight * row[r] * target.x;
        y_values[r] += weight * row[r] * target.y;
        z_values[r] += weight * row[r] * target.z;
        for c in 0..row.len() {
            matrix[r][c] += weight * row[r] * row[c];
        }
    }
}

fn basis_row(degree: usize, knots: &[f64], count: usize, parameter: f64) -> Vec<f64> {
    if parameter <= 1.0e-14 {
        let mut row = vec![0.0; count];
        row[0] = 1.0;
        return row;
    }
    if (1.0 - parameter).abs() <= 1.0e-14 {
        let mut row = vec![0.0; count];
        row[count - 1] = 1.0;
        return row;
    }
    (0..count)
        .map(|index| basis_function(index, degree, parameter, knots))
        .collect()
}

fn derivative_basis_row(degree: usize, knots: &[f64], count: usize, parameter: f64) -> Vec<f64> {
    let h = 1.0e-6;
    let lo = (parameter - h).max(0.0);
    let hi = (parameter + h).min(1.0);
    let span = hi - lo;
    if span <= 1.0e-12 {
        return vec![0.0; count];
    }
    let lo_basis = basis_row(degree, knots, count, lo);
    let hi_basis = basis_row(degree, knots, count, hi);
    hi_basis
        .into_iter()
        .zip(lo_basis)
        .map(|(hi, lo)| (hi - lo) / span)
        .collect()
}

fn open_uniform_knots(degree: usize, control_count: usize) -> Vec<f64> {
    let mut knots = Vec::with_capacity(control_count + degree + 1);
    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    let interior = control_count.saturating_sub(degree + 1);
    for index in 1..=interior {
        knots.push(index as f64 / (interior + 1) as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}

fn basis_function(index: usize, degree: usize, parameter: f64, knots: &[f64]) -> f64 {
    if degree == 0 {
        let is_last_span = index + 2 == knots.len();
        if (knots[index] <= parameter && parameter < knots[index + 1])
            || (is_last_span && (parameter - knots[index + 1]).abs() <= 1.0e-14)
        {
            return 1.0;
        }
        return 0.0;
    }
    let mut value = 0.0;
    let left = knots[index + degree] - knots[index];
    if left.abs() > 1.0e-14 {
        value +=
            (parameter - knots[index]) / left * basis_function(index, degree - 1, parameter, knots);
    }
    let right = knots[index + degree + 1] - knots[index + 1];
    if right.abs() > 1.0e-14 {
        value += (knots[index + degree + 1] - parameter) / right
            * basis_function(index + 1, degree - 1, parameter, knots);
    }
    value
}

#[allow(clippy::needless_range_loop)]
fn solve_linear_system(matrix: &[Vec<f64>], values: &[f64]) -> Option<Vec<f64>> {
    let count = matrix.len();
    if values.len() != count || matrix.iter().any(|row| row.len() != count) {
        return None;
    }
    let mut rows = matrix
        .iter()
        .zip(values)
        .map(|(row, value)| {
            let mut augmented = row.clone();
            augmented.push(*value);
            augmented
        })
        .collect::<Vec<_>>();
    for column in 0..count {
        let pivot = (column..count).max_by(|a, b| {
            rows[*a][column]
                .abs()
                .partial_cmp(&rows[*b][column].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if rows[pivot][column].abs() < 1.0e-12 {
            return None;
        }
        rows.swap(column, pivot);
        for row in (column + 1)..count {
            let factor = rows[row][column] / rows[column][column];
            for entry in column..=count {
                rows[row][entry] -= factor * rows[column][entry];
            }
        }
    }
    let mut solution = vec![0.0; count];
    for row in (0..count).rev() {
        let mut value = rows[row][count];
        for column in (row + 1)..count {
            value -= rows[row][column] * solution[column];
        }
        solution[row] = value / rows[row][row];
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_uniform_basis_keeps_clamped_endpoints_exact() {
        let knots = open_uniform_knots(5, 17);

        let start = basis_row(5, &knots, 17, 0.0);
        assert_eq!(start[0], 1.0);
        assert!(start[1..].iter().all(|value| *value == 0.0));

        let end = basis_row(5, &knots, 17, 1.0);
        assert_eq!(end[16], 1.0);
        assert!(end[..16].iter().all(|value| *value == 0.0));
    }
}
