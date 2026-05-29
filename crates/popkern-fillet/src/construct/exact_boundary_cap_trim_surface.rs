//! Exact-boundary cap surface-row relock helpers.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::exact_boundary_cap_trim_fit::{
    greville_parameters, interpolate_points_on_fixed_knots, knots_equivalent,
};
use super::exact_boundary_cap_trim_restore::ExactBoundaryCapSurfaceSide;

pub(super) fn relock_exact_boundary_cap_side(
    brep: &mut BRep,
    half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    side: ExactBoundaryCapSurfaceSide,
    side_forward: bool,
) -> Option<()> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    let edge = brep.half_edges.get(half_edge)?.edge;
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let edge_forward = brep.edges.get(edge)?.half_edge == half_edge;
    let edge_direction = if side_forward == edge_forward {
        BoundarySampleDirection::Forward
    } else {
        BoundarySampleDirection::Reverse
    };
    let support_direction = if side_forward {
        BoundarySampleDirection::Forward
    } else {
        BoundarySampleDirection::Reverse
    };
    let support_trace = support_boundary_trace(brep, support_half_edge, support_direction);
    let surface = brep.faces.get_mut(face)?.surface.as_mut()?;
    let SurfaceGeometry::Nurbs(surface_nurbs) = &mut surface.geometry else {
        return None;
    };
    let boundary = support_trace
        .as_ref()
        .and_then(|trace| compatible_sampled_boundary_curve(surface_nurbs, side, trace))
        .or_else(|| compatible_edge_boundary_curve(surface_nurbs, side, curve, edge_direction));
    let boundary = boundary?;
    relock_nurbs_boundary_row(surface_nurbs, side, &boundary)?;
    if let Some(trace) = support_trace.as_ref() {
        relock_nurbs_adjacent_row(surface_nurbs, side, &boundary, trace);
    }
    Some(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BoundarySampleDirection {
    Forward,
    Reverse,
}

struct BoundaryTrace {
    points: Vec<Point3>,
    normals: Vec<Vec3>,
}

fn support_boundary_trace(
    brep: &BRep,
    support_half_edge: HalfEdgeHandle,
    direction: BoundarySampleDirection,
) -> Option<BoundaryTrace> {
    let half_edge = brep.half_edges.get(support_half_edge)?;
    let pcurve = half_edge.pcurve.as_ref()?;
    let face = brep.loops.get(half_edge.loop_)?.face;
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let (t_start, t_end) = pcurve.parameter_range();
    let mut points = Vec::with_capacity(65);
    let mut normals = Vec::with_capacity(65);
    for index in 0..65 {
        let fraction = index as f64 / 64.0;
        let support_fraction = match direction {
            BoundarySampleDirection::Forward => fraction,
            BoundarySampleDirection::Reverse => 1.0 - fraction,
        };
        let t = t_start + (t_end - t_start) * support_fraction;
        let uv = pcurve.point_at(t);
        points.push(surface.geometry.point_at(uv.x, uv.y));
        normals.push(surface.geometry.normal_at(uv.x, uv.y));
    }
    Some(BoundaryTrace { points, normals })
}

fn compatible_sampled_boundary_curve(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    trace: &BoundaryTrace,
) -> Option<NurbsCurve> {
    let expected = side_control_count(surface, side);
    let degree = side_degree(surface, side).min(expected.saturating_sub(1));
    if expected < 2 || degree < 1 {
        return None;
    }
    let side_knots = side_knots(surface, side);
    interpolate_trace_on_fixed_knots(trace, expected, degree, side_knots)
}

fn compatible_edge_boundary_curve(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    curve: &CurveOnEdge,
    direction: BoundarySampleDirection,
) -> Option<NurbsCurve> {
    let expected = side_control_count(surface, side);
    let degree = side_degree(surface, side).min(expected.saturating_sub(1));
    if expected < 2 || degree < 1 {
        return None;
    }
    let side_knots = side_knots(surface, side);
    let mut nurbs = curve
        .geometry
        .to_nurbs_in_range((curve.t_start, curve.t_end));
    if direction == BoundarySampleDirection::Reverse {
        nurbs = reverse_nurbs_curve(&nurbs);
    }
    if nurbs.degree == degree
        && nurbs.control_points.len() == expected
        && nurbs.weights.len() == expected
        && knots_equivalent(&nurbs.knots, side_knots)
    {
        return Some(nurbs);
    }
    interpolate_edge_curve_for_surface_side(curve, direction, expected, degree, side_knots)
}

fn interpolate_trace_on_fixed_knots(
    trace: &BoundaryTrace,
    expected: usize,
    degree: usize,
    knots: &[f64],
) -> Option<NurbsCurve> {
    let first = *knots.first()?;
    let last = *knots.last()?;
    if trace.points.len() < 2
        || knots.len() != expected + degree + 1
        || (last - first).abs() <= f64::EPSILON
    {
        return None;
    }
    let params = greville_parameters(expected, degree, knots)?;
    let points = params
        .iter()
        .map(|parameter| {
            let fraction = ((*parameter - first) / (last - first)).clamp(0.0, 1.0);
            sample_trace(trace, fraction)
        })
        .collect::<Vec<_>>();
    let control_points = interpolate_points_on_fixed_knots(&points, &params, degree, knots)?;
    Some(NurbsCurve {
        degree,
        knots: knots.to_vec(),
        weights: vec![1.0; control_points.len()],
        control_points,
    })
}

fn sample_trace(trace: &BoundaryTrace, fraction: f64) -> Point3 {
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

fn sample_trace_normal(trace: &BoundaryTrace, fraction: f64) -> Option<Vec3> {
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

fn parameter_fraction(parameter: f64, knots: &[f64]) -> Option<f64> {
    let first = *knots.first()?;
    let last = *knots.last()?;
    if (last - first).abs() <= f64::EPSILON {
        return None;
    }
    Some(((parameter - first) / (last - first)).clamp(0.0, 1.0))
}

fn interpolate_edge_curve_for_surface_side(
    curve: &CurveOnEdge,
    direction: BoundarySampleDirection,
    expected: usize,
    degree: usize,
    side_knots: &[f64],
) -> Option<NurbsCurve> {
    let first = *side_knots.first()?;
    let last = *side_knots.last()?;
    if side_knots.len() != expected + degree + 1 || (last - first).abs() <= f64::EPSILON {
        return None;
    }
    let mut params = Vec::with_capacity(expected);
    let mut points = Vec::with_capacity(expected);
    for index in 0..expected {
        let fraction = index as f64 / (expected - 1) as f64;
        let curve_fraction = match direction {
            BoundarySampleDirection::Forward => fraction,
            BoundarySampleDirection::Reverse => 1.0 - fraction,
        };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let side_t = first + (last - first) * fraction;
        params.push(side_t);
        points.push(curve.geometry.point_at(curve_t));
    }
    let mut interpolated =
        popkern_core::geometry::fitting::interpolate_curve(&points, &params, degree).ok()?;
    if !knots_equivalent(&interpolated.knots, side_knots) {
        interpolated = interpolate_edge_curve_on_fixed_knots(
            curve, direction, expected, degree, side_knots, first, last,
        )?;
    }
    Some(interpolated)
}

fn interpolate_edge_curve_on_fixed_knots(
    curve: &CurveOnEdge,
    direction: BoundarySampleDirection,
    expected: usize,
    degree: usize,
    knots: &[f64],
    first: f64,
    last: f64,
) -> Option<NurbsCurve> {
    let params = greville_parameters(expected, degree, knots)?;
    let points = params
        .iter()
        .map(|parameter| {
            let fraction = ((*parameter - first) / (last - first)).clamp(0.0, 1.0);
            let curve_fraction = match direction {
                BoundarySampleDirection::Forward => fraction,
                BoundarySampleDirection::Reverse => 1.0 - fraction,
            };
            let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
            curve.geometry.point_at(curve_t)
        })
        .collect::<Vec<_>>();
    let control_points = interpolate_points_on_fixed_knots(&points, &params, degree, knots)?;
    Some(NurbsCurve {
        degree,
        knots: knots.to_vec(),
        weights: vec![1.0; control_points.len()],
        control_points,
    })
}

fn relock_nurbs_boundary_row(
    surface: &mut NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    curve: &NurbsCurve,
) -> Option<()> {
    let expected = side_control_count(surface, side);
    if curve.control_points.len() != expected || curve.weights.len() != expected {
        return None;
    }
    for index in 0..expected {
        let surface_index = match side {
            ExactBoundaryCapSurfaceSide::VMin => index * surface.num_v,
            ExactBoundaryCapSurfaceSide::VMax => index * surface.num_v + (surface.num_v - 1),
            ExactBoundaryCapSurfaceSide::UMin => index,
            ExactBoundaryCapSurfaceSide::UMax => (surface.num_u - 1) * surface.num_v + index,
        };
        surface.control_points[surface_index] = curve.control_points[index];
        surface.weights[surface_index] = curve.weights[index];
    }
    Some(())
}

fn relock_nurbs_adjacent_row(
    surface: &mut NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
    boundary: &NurbsCurve,
    support_trace: &BoundaryTrace,
) -> Option<()> {
    let expected = side_control_count(surface, side);
    if expected < 2 || boundary.control_points.len() != expected {
        return None;
    }
    let degree = side_degree(surface, side).min(expected.saturating_sub(1));
    let knots = side_knots(surface, side);
    let params = greville_parameters(expected, degree, knots)?;
    let scale = cross_derivative_control_scale(surface, side)?;
    let mut adjacent_points = Vec::with_capacity(expected);
    for parameter in &params {
        let fraction = parameter_fraction(*parameter, knots)?;
        let boundary_point = boundary.point_at(*parameter);
        let tangent = boundary_tangent(boundary, *parameter)?;
        let support_normal = sample_trace_normal(support_trace, fraction)?;
        let cross_direction = desired_cross_direction(side, tangent, support_normal)?;
        adjacent_points.push(adjacent_point_for_side(
            side,
            boundary_point,
            cross_direction,
            scale,
        ));
    }
    let adjacent_curve =
        interpolate_points_on_fixed_knots(&adjacent_points, &params, degree, knots)?;
    for (index, point) in adjacent_curve.iter().enumerate() {
        let surface_index = adjacent_row_index(surface, side, index)?;
        surface.control_points[surface_index] = *point;
        surface.weights[surface_index] = 1.0;
    }
    Some(())
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

fn adjacent_point_for_side(
    side: ExactBoundaryCapSurfaceSide,
    boundary: Point3,
    cross_direction: Vec3,
    scale: f64,
) -> Point3 {
    match side {
        ExactBoundaryCapSurfaceSide::VMin | ExactBoundaryCapSurfaceSide::UMin => {
            boundary + cross_direction * scale
        }
        ExactBoundaryCapSurfaceSide::VMax | ExactBoundaryCapSurfaceSide::UMax => {
            boundary - cross_direction * scale
        }
    }
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

fn reverse_nurbs_curve(nurbs: &NurbsCurve) -> NurbsCurve {
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
