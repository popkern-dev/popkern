//! Support-boundary trace scoring for exact-boundary cap relock.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

pub(super) struct BoundaryTrace {
    pub(super) points: Vec<Point3>,
    pub(super) normals: Vec<Vec3>,
}

#[derive(Debug, Clone, Copy)]
struct TraceBoundaryCurveScore {
    point_gap: f64,
    tangent_deficit: f64,
}

pub(super) fn select_trace_boundary_curve(
    trace: &BoundaryTrace,
    interpolated: NurbsCurve,
    sampled: NurbsCurve,
) -> NurbsCurve {
    let interpolated_score = trace_boundary_curve_score(trace, &interpolated);
    let sampled_score = trace_boundary_curve_score(trace, &sampled);
    if sampled_score.point_gap > interpolated_score.point_gap * 2.0
        || sampled_score.point_gap > 1.0e-4
    {
        return interpolated;
    }
    if compare_trace_boundary_scores(sampled_score, interpolated_score).is_lt() {
        sampled
    } else {
        interpolated
    }
}

fn trace_boundary_curve_score(
    trace: &BoundaryTrace,
    curve: &NurbsCurve,
) -> TraceBoundaryCurveScore {
    let (start, end) = curve.parameter_range();
    let mut point_gap = 0.0_f64;
    let mut tangent_deficit = 0.0_f64;
    for index in 0..65 {
        let fraction = index as f64 / 64.0;
        let parameter = start + (end - start) * fraction;
        point_gap = point_gap.max(
            curve
                .point_at(parameter)
                .distance_to(&sample_trace(trace, fraction)),
        );
        if let (Some(curve_tangent), Some(trace_tangent)) = (
            curve_tangent(curve, parameter),
            trace_tangent(trace, fraction),
        ) {
            tangent_deficit = tangent_deficit.max(1.0 - curve_tangent.dot(&trace_tangent).abs());
        }
    }
    TraceBoundaryCurveScore {
        point_gap,
        tangent_deficit,
    }
}

fn compare_trace_boundary_scores(
    a: TraceBoundaryCurveScore,
    b: TraceBoundaryCurveScore,
) -> std::cmp::Ordering {
    let a_position_valid = a.point_gap <= 5.0e-3;
    let b_position_valid = b.point_gap <= 5.0e-3;
    match (a_position_valid, b_position_valid) {
        (true, false) => std::cmp::Ordering::Less,
        (false, true) => std::cmp::Ordering::Greater,
        (true, true) => a
            .tangent_deficit
            .total_cmp(&b.tangent_deficit)
            .then_with(|| a.point_gap.total_cmp(&b.point_gap)),
        (false, false) => a
            .point_gap
            .total_cmp(&b.point_gap)
            .then_with(|| a.tangent_deficit.total_cmp(&b.tangent_deficit)),
    }
}

pub(super) fn sample_trace(trace: &BoundaryTrace, fraction: f64) -> Point3 {
    let scaled = fraction.clamp(0.0, 1.0) * (trace.points.len() - 1) as f64;
    let lower = scaled.floor() as usize;
    let upper = scaled.ceil() as usize;
    if lower == upper {
        return trace.points[lower];
    }
    let local = scaled - lower as f64;
    let a = trace.points[lower];
    let b = trace.points[upper];
    Point3::new(
        a.x + (b.x - a.x) * local,
        a.y + (b.y - a.y) * local,
        a.z + (b.z - a.z) * local,
    )
}

fn trace_tangent(trace: &BoundaryTrace, fraction: f64) -> Option<Vec3> {
    if trace.points.len() < 2 {
        return None;
    }
    let h = (1.0 / (trace.points.len() - 1) as f64).max(1.0e-5);
    let before = (fraction - h).clamp(0.0, 1.0);
    let after = (fraction + h).clamp(0.0, 1.0);
    if (after - before).abs() <= f64::EPSILON {
        return None;
    }
    (sample_trace(trace, after) - sample_trace(trace, before)).normalized()
}

pub(super) fn sample_trace_normal(trace: &BoundaryTrace, fraction: f64) -> Option<Vec3> {
    if trace.normals.len() < 2 {
        return None;
    }
    let scaled = fraction.clamp(0.0, 1.0) * (trace.normals.len() - 1) as f64;
    let lower = scaled.floor() as usize;
    let upper = scaled.ceil() as usize;
    if lower == upper {
        return trace.normals[lower].normalized();
    }
    let local = scaled - lower as f64;
    (trace.normals[lower] * (1.0 - local) + trace.normals[upper] * local).normalized()
}

pub(super) fn reverse_nurbs_curve(nurbs: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = nurbs.parameter_range();
    let mut knots = nurbs
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect::<Vec<_>>();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = nurbs.control_points.clone();
    control_points.reverse();
    let mut weights = nurbs.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: nurbs.degree,
        knots,
        control_points,
        weights,
    }
}

fn curve_tangent(curve: &NurbsCurve, parameter: f64) -> Option<Vec3> {
    let (start, end) = curve.parameter_range();
    let span = (end - start).abs();
    if span <= f64::EPSILON {
        return None;
    }
    let h = (span * 1.0e-5).max(1.0e-6);
    let before = (parameter - h).clamp(start, end);
    let after = (parameter + h).clamp(start, end);
    if (after - before).abs() <= f64::EPSILON {
        return None;
    }
    let tangent = curve.point_at(after) - curve.point_at(before);
    tangent.normalized()
}
