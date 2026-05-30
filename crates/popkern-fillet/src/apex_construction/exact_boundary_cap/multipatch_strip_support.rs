//! Support-normal helper geometry for multi-patch strip construction.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;

const MIN_STRIP_WIDTH: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct BoundaryTangentCompatibilityReport {
    pub(super) min_reachable_deficit: f64,
    pub(super) sample_parameter: f64,
}

pub(super) fn boundary_tangent_compatibility(
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
) -> BoundaryTangentCompatibilityReport {
    let mut report = BoundaryTangentCompatibilityReport {
        min_reachable_deficit: 0.0,
        sample_parameter: samples
            .first()
            .map(|sample| sample.parameter)
            .unwrap_or(0.0),
    };
    for sample in samples {
        let Some(tangent) = boundary_tangent(boundary, sample.parameter) else {
            continue;
        };
        let Some(normal) = sample.support_normal.normalized() else {
            continue;
        };
        let dot = tangent.dot(&normal).abs().min(1.0);
        let reachable_alignment = (1.0 - dot * dot).max(0.0).sqrt();
        let min_reachable_deficit = 1.0 - reachable_alignment;
        if min_reachable_deficit > report.min_reachable_deficit {
            report = BoundaryTangentCompatibilityReport {
                min_reachable_deficit,
                sample_parameter: sample.parameter,
            };
        }
    }
    report
}

pub(super) fn endpoint_fallback_width_direction(
    samples: &[ExactBoundaryRibbonSample],
    parameter: f64,
) -> Option<Vec3> {
    if !is_endpoint_sample(parameter) {
        return None;
    }
    samples
        .iter()
        .filter(|sample| !is_endpoint_sample(sample.parameter))
        .filter_map(|sample| sample.cross_derivative.normalized())
        .next()
}

pub(super) fn strip_width_from_samples(samples: &[ExactBoundaryRibbonSample]) -> Option<f64> {
    let mut width_sum = 0.0;
    let mut count = 0;
    for sample in samples {
        let width = sample.cross_derivative.length();
        if width.is_finite() && width > MIN_STRIP_WIDTH {
            width_sum += width;
            count += 1;
        }
    }
    (count > 0).then_some(width_sum / count as f64)
}

pub(super) fn support_refined_boundary(
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
) -> Option<NurbsCurve> {
    let (start, end) = boundary.parameter_range();
    let tolerance = 1.0e-12 * (end - start).abs().max(1.0);
    let mut knots = samples
        .iter()
        .map(|sample| sample.parameter)
        .filter(|parameter| parameter.is_finite())
        .filter(|parameter| *parameter > start + tolerance && *parameter < end - tolerance)
        .collect::<Vec<_>>();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    knots.dedup_by(|a, b| (*a - *b).abs() <= tolerance);
    boundary.insert_knots(&knots).ok()
}

pub(super) fn support_frame_row_controls(
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
    strip_width: f64,
    width_scale: f64,
) -> Option<Vec<Point3>> {
    let fit_parameters = support_frame_fit_parameters(boundary, samples);
    let fit_targets =
        support_frame_targets(boundary, samples, strip_width, width_scale, &fit_parameters)?;
    least_squares_controls(
        &fit_targets,
        &fit_parameters,
        boundary.degree,
        &boundary.knots,
    )
    .or_else(|| {
        let parameters = (0..boundary.control_points.len())
            .map(|index| control_row_parameter(boundary, index))
            .collect::<Vec<_>>();
        let targets =
            support_frame_targets(boundary, samples, strip_width, width_scale, &parameters)?;
        fixed_knot_interpolation_controls(&targets, &parameters, boundary.degree, &boundary.knots)
    })
}

fn support_frame_fit_parameters(
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
) -> Vec<f64> {
    let (start, end) = boundary.parameter_range();
    let mut parameters = samples
        .iter()
        .map(|sample| sample.parameter)
        .chain((0..100).map(|index| start + (end - start) * index as f64 / 99.0))
        .collect::<Vec<_>>();
    parameters.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    parameters.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-12);
    parameters
}

fn support_frame_targets(
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
    strip_width: f64,
    width_scale: f64,
    parameters: &[f64],
) -> Option<Vec<Point3>> {
    parameters
        .iter()
        .map(|parameter| {
            let sample = sample_at_parameter(samples, *parameter);
            let direction = width_direction(boundary, *parameter, sample.support_normal)
                .or_else(|| {
                    width_direction_from_tangent(sample.boundary_tangent, sample.support_normal)
                })
                .or_else(|| {
                    width_direction_from_tangent(sample.cross_derivative, sample.support_normal)
                })
                .or_else(|| sample.cross_derivative.normalized())
                .or_else(|| endpoint_fallback_width_direction(samples, sample.parameter))?;
            Some(boundary.point_at(*parameter) + direction * strip_width * width_scale)
        })
        .collect()
}

pub(super) fn control_row_parameter(boundary: &NurbsCurve, index: usize) -> f64 {
    let (t_start, t_end) = boundary.parameter_range();
    if index == 0 || boundary.degree == 0 {
        return t_start;
    }
    if index + 1 == boundary.control_points.len() {
        return t_end;
    }
    let start = index + 1;
    let end = index + boundary.degree + 1;
    boundary.knots[start..end].iter().sum::<f64>() / boundary.degree as f64
}

pub(super) fn sample_at_parameter(
    samples: &[ExactBoundaryRibbonSample],
    parameter: f64,
) -> ExactBoundaryRibbonSample {
    let mut best = samples[0];
    let mut best_distance = (best.parameter - parameter).abs();
    for sample in &samples[1..] {
        let distance = (sample.parameter - parameter).abs();
        if distance < best_distance {
            best = *sample;
            best_distance = distance;
        }
    }
    best
}

pub(super) fn width_direction(
    boundary: &NurbsCurve,
    parameter: f64,
    support_normal: Vec3,
) -> Option<Vec3> {
    let tangent = boundary_tangent(boundary, parameter)?;
    width_direction_from_tangent(tangent, support_normal)
}

pub(super) fn width_direction_from_tangent(tangent: Vec3, support_normal: Vec3) -> Option<Vec3> {
    let normal = support_normal.normalized()?;
    normal.cross(&tangent).normalized()
}

pub(super) fn is_endpoint_sample(parameter: f64) -> bool {
    parameter <= 1.0e-12 || (1.0 - parameter).abs() <= 1.0e-12
}

fn boundary_tangent(boundary: &NurbsCurve, parameter: f64) -> Option<Vec3> {
    if let Some(tangent) = boundary.derivative_at(parameter, 1).normalized() {
        return Some(tangent);
    }
    let (start, end) = boundary.parameter_range();
    let span = end - start;
    if span.abs() <= 1.0e-12 {
        return None;
    }
    let inset = span.abs() * 1.0e-6;
    for fallback in [
        (parameter + inset).clamp(start, end),
        (parameter - inset).clamp(start, end),
        start + 0.5 * span,
    ] {
        if (fallback - parameter).abs() <= 1.0e-15 {
            continue;
        }
        if let Some(tangent) = boundary.derivative_at(fallback, 1).normalized() {
            return Some(tangent);
        }
    }
    None
}

fn fixed_knot_interpolation_controls(
    points: &[Point3],
    parameters: &[f64],
    degree: usize,
    knots: &[f64],
) -> Option<Vec<Point3>> {
    let count = points.len();
    if count == 0 || parameters.len() != count || knots.len() != count + degree + 1 {
        return None;
    }
    let mut matrix = vec![vec![0.0; count]; count];
    for (row, parameter) in parameters.iter().enumerate() {
        if (*parameter - knots[degree]).abs() <= 1.0e-14 {
            matrix[row][0] = 1.0;
            continue;
        }
        if (*parameter - knots[count]).abs() <= 1.0e-14 {
            matrix[row][count - 1] = 1.0;
            continue;
        }
        for column in 0..count {
            matrix[row][column] = basis_function(column, degree, *parameter, knots);
        }
    }
    let x = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.x).collect::<Vec<_>>(),
    )?;
    let y = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.y).collect::<Vec<_>>(),
    )?;
    let z = solve_linear_system(
        &matrix,
        &points.iter().map(|point| point.z).collect::<Vec<_>>(),
    )?;
    Some(
        (0..count)
            .map(|index| Point3::new(x[index], y[index], z[index]))
            .collect(),
    )
}

fn least_squares_controls(
    points: &[Point3],
    parameters: &[f64],
    degree: usize,
    knots: &[f64],
) -> Option<Vec<Point3>> {
    let count = knots.len().checked_sub(degree + 1)?;
    if count == 0 || points.len() != parameters.len() || points.len() < count {
        return None;
    }
    let mut normal_matrix = vec![vec![0.0; count]; count];
    let mut x_values = vec![0.0; count];
    let mut y_values = vec![0.0; count];
    let mut z_values = vec![0.0; count];
    for (point, parameter) in points.iter().zip(parameters) {
        let basis = (0..count)
            .map(|index| basis_function(index, degree, *parameter, knots))
            .collect::<Vec<_>>();
        for row in 0..count {
            x_values[row] += basis[row] * point.x;
            y_values[row] += basis[row] * point.y;
            z_values[row] += basis[row] * point.z;
            for column in 0..count {
                normal_matrix[row][column] += basis[row] * basis[column];
            }
        }
    }
    let x = solve_linear_system(&normal_matrix, &x_values)?;
    let y = solve_linear_system(&normal_matrix, &y_values)?;
    let z = solve_linear_system(&normal_matrix, &z_values)?;
    Some(
        (0..count)
            .map(|index| Point3::new(x[index], y[index], z[index]))
            .collect(),
    )
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
