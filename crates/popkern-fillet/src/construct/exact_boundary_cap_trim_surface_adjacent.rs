//! Adjacent normal-row transaction for exact-boundary cap relock.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use super::super::exact_boundary_cap_trim_fit::{
    greville_parameters, interpolate_points_on_fixed_knots,
};
use super::super::exact_boundary_cap_trim_restore::ExactBoundaryCapSurfaceSide;
use super::trace::{BoundaryTrace, sample_trace, sample_trace_normal};

pub(super) fn relock_nurbs_adjacent_row(
    surface: &mut NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    boundary: &NurbsCurve,
    support_trace: &BoundaryTrace,
    allow_wide_search: bool,
) -> Option<()> {
    let expected = side_control_count(surface, side);
    if expected < 2 || boundary.control_points.len() != expected {
        return None;
    }
    let degree = side_degree(surface, side).min(expected.saturating_sub(1));
    let knots = side_knots(surface, side);
    let params = greville_parameters(expected, degree, knots)?;
    let scale = cross_derivative_control_scale(surface, side)?;
    let mut adjacent_directions = Vec::with_capacity(expected);
    for parameter in &params {
        let fraction = parameter_fraction(*parameter, knots)?;
        let boundary_point = boundary.point_at(*parameter);
        let tangent = boundary_tangent(boundary, *parameter)?;
        let support_normal = sample_trace_normal(support_trace, fraction)?;
        let cross_direction = desired_cross_direction(side, tangent, support_normal)?;
        adjacent_directions.push((
            boundary_point,
            adjacent_offset_direction_for_side(side, cross_direction, scale),
        ));
    }
    let baseline = support_boundary_score_for_surface(surface, side, support_trace, &params, knots)
        .unwrap_or(AdjacentRowScore {
            point_gap: f64::INFINITY,
            normal_deficit: f64::INFINITY,
        });
    let mut best_surface = None;
    let mut best_score = baseline;
    let scale_factors: &[f64] = if allow_wide_search {
        &[0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
    } else {
        &[1.0]
    };
    for scale_factor in scale_factors {
        for polarity in [1.0, -1.0] {
            let adjacent_points = adjacent_directions
                .iter()
                .map(|(boundary_point, offset)| {
                    *boundary_point + *offset * *scale_factor * polarity
                })
                .collect::<Vec<_>>();
            let adjacent_curve =
                interpolate_points_on_fixed_knots(&adjacent_points, &params, degree, knots)?;
            let mut candidate = surface.clone();
            apply_adjacent_row_points(&mut candidate, side, &adjacent_curve)?;
            let Some(score) =
                support_boundary_score_for_surface(&candidate, side, support_trace, &params, knots)
            else {
                continue;
            };
            if compare_adjacent_row_scores(score, best_score).is_lt() {
                best_score = score;
                best_surface = Some(candidate);
            }
        }
    }
    if let Some(candidate) = best_surface {
        *surface = candidate;
    }
    Some(())
}

fn adjacent_offset_direction_for_side(
    side: ExactBoundaryCapSurfaceSide,
    cross_direction: Vec3,
    scale: f64,
) -> Vec3 {
    match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::UMin => {
            cross_direction * scale
        }
        ExactBoundaryCapSurfaceSide::VMax | ExactBoundaryCapSurfaceSide::UMax => {
            cross_direction * -scale
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct AdjacentRowScore {
    point_gap: f64,
    normal_deficit: f64,
}

fn compare_adjacent_row_scores(a: AdjacentRowScore, b: AdjacentRowScore) -> std::cmp::Ordering {
    let a_position_valid = a.point_gap <= 5.0e-3;
    let b_position_valid = b.point_gap <= 5.0e-3;
    match (a_position_valid, b_position_valid) {
        (true, false) => std::cmp::Ordering::Less,
        (false, true) => std::cmp::Ordering::Greater,
        (true, true) => a
            .normal_deficit
            .total_cmp(&b.normal_deficit)
            .then_with(|| a.point_gap.total_cmp(&b.point_gap)),
        (false, false) => a
            .point_gap
            .total_cmp(&b.point_gap)
            .then_with(|| a.normal_deficit.total_cmp(&b.normal_deficit)),
    }
}

fn support_boundary_score_for_surface(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    support_trace: &BoundaryTrace,
    params: &[f64],
    knots: &[f64],
) -> Option<AdjacentRowScore> {
    let mut point_gap = 0.0_f64;
    let mut normal_deficit = 0.0_f64;
    for parameter in params {
        let fraction = parameter_fraction(*parameter, knots)?;
        let support_normal = sample_trace_normal(support_trace, fraction)?;
        let (u, v) = side_normal_sample_uv(surface, side, *parameter)?;
        let support_point = sample_trace(support_trace, fraction);
        let cap_point = surface.point_at(u, v);
        point_gap = point_gap.max(cap_point.distance_to(&support_point));
        let cap_normal = surface.normal_at(u, v).normalized()?;
        let deficit = 1.0 - cap_normal.dot(&support_normal).abs();
        if !point_gap.is_finite() || !deficit.is_finite() {
            return None;
        }
        normal_deficit = normal_deficit.max(deficit);
    }
    Some(AdjacentRowScore {
        point_gap,
        normal_deficit,
    })
}

fn side_normal_sample_uv(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    parameter: f64,
) -> Option<(f64, f64)> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    Some(match side {
        ExactBoundaryCapSurfaceSide::VMin => (parameter, interior_parameter(v_min, v_min, v_max)),
        ExactBoundaryCapSurfaceSide::VMax => (parameter, interior_parameter(v_max, v_min, v_max)),
        ExactBoundaryCapSurfaceSide::UMin => (interior_parameter(u_min, u_min, u_max), parameter),
        ExactBoundaryCapSurfaceSide::UMax => (interior_parameter(u_max, u_min, u_max), parameter),
    })
}

fn apply_adjacent_row_points(
    surface: &mut NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    adjacent_points: &[Point3],
) -> Option<()> {
    // The adjacent row's END control points lie ON the transverse
    // boundary rows (a UMin-adjacent column starts in the VMin row and
    // ends in the VMax row). Writing them displaced committed neighbour
    // boundaries — the cap spokes and FilletEnd rows — by up to 5.4e-2;
    // that displacement was the historical source of the pyramid4 zoo
    // tessellation seam family (PK-124). The G1 correction pins them
    // and moves only the interior.
    let last = adjacent_points.len().checked_sub(1)?;
    for (index, point) in adjacent_points.iter().enumerate() {
        if index == 0 || index == last {
            continue;
        }
        let surface_index = adjacent_row_index(surface, side, index)?;
        surface.control_points[surface_index] = *point;
        surface.weights[surface_index] = 1.0;
    }
    Some(())
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}

fn cross_derivative_control_scale(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
) -> Option<f64> {
    let (degree, knots, at_min) = match side {
        ExactBoundaryCapSurfaceSide::VMin => (surface.degree_v, &surface.knots_v, true),
        ExactBoundaryCapSurfaceSide::VMax => (surface.degree_v, &surface.knots_v, false),
        ExactBoundaryCapSurfaceSide::UMin => (surface.degree_u, &surface.knots_u, true),
        ExactBoundaryCapSurfaceSide::UMax => (surface.degree_u, &surface.knots_u, false),
    };
    if degree == 0 {
        return None;
    }
    let boundary = if at_min {
        *knots.get(degree)?
    } else {
        *knots.get(knots.len().checked_sub(degree + 1)?)?
    };
    let neighbor = if at_min {
        knots
            .iter()
            .copied()
            .find(|knot| *knot > boundary + 1.0e-14)?
    } else {
        knots
            .iter()
            .rev()
            .copied()
            .find(|knot| *knot < boundary - 1.0e-14)?
    };
    Some((neighbor - boundary).abs() / degree as f64)
}

fn boundary_tangent(curve: &NurbsCurve, parameter: f64) -> Option<Vec3> {
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

fn desired_cross_direction(
    side: ExactBoundaryCapSurfaceSide,
    tangent: Vec3,
    support_normal: Vec3,
) -> Option<Vec3> {
    let normal = support_normal.normalized()?;
    let direction = match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::VMax => {
            normal.cross(&tangent)
        }
        ExactBoundaryCapSurfaceSide::UMin | ExactBoundaryCapSurfaceSide::UMax => {
            tangent.cross(&normal)
        }
    };
    direction.normalized()
}

fn adjacent_row_index(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    index: usize,
) -> Option<usize> {
    match side {
        ExactBoundaryCapSurfaceSide::VMin => Some(index * surface.num_v + 1),
        ExactBoundaryCapSurfaceSide::VMax => {
            Some(index * surface.num_v + surface.num_v.checked_sub(2)?)
        }
        ExactBoundaryCapSurfaceSide::UMin => Some(surface.num_v + index),
        ExactBoundaryCapSurfaceSide::UMax => {
            Some(surface.num_u.checked_sub(2)? * surface.num_v + index)
        }
    }
}

fn parameter_fraction(parameter: f64, knots: &[f64]) -> Option<f64> {
    let first = *knots.first()?;
    let last = *knots.last()?;
    if (last - first).abs() <= f64::EPSILON {
        return None;
    }
    Some(((parameter - first) / (last - first)).clamp(0.0, 1.0))
}

fn side_control_count(surface: &NurbsSurface, side: ExactBoundaryCapSurfaceSide) -> usize {
    match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::VMax => surface.num_u,
        ExactBoundaryCapSurfaceSide::UMin | ExactBoundaryCapSurfaceSide::UMax => surface.num_v,
    }
}

fn side_degree(surface: &NurbsSurface, side: ExactBoundaryCapSurfaceSide) -> usize {
    match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::VMax => surface.degree_u,
        ExactBoundaryCapSurfaceSide::UMin | ExactBoundaryCapSurfaceSide::UMax => surface.degree_v,
    }
}

fn side_knots(surface: &NurbsSurface, side: ExactBoundaryCapSurfaceSide) -> &[f64] {
    match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::VMax => &surface.knots_u,
        ExactBoundaryCapSurfaceSide::UMin | ExactBoundaryCapSurfaceSide::UMax => &surface.knots_v,
    }
}
