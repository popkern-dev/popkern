//! Support-normal helper geometry for multi-patch strip construction.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Vec3};

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
