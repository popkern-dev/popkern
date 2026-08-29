//! Exact-boundary cap surface-row relock helpers.

#[path = "exact_boundary_cap_trim_surface_adjacent.rs"]
mod adjacent;
#[path = "exact_boundary_cap_trim_surface_trace.rs"]
mod trace;

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, NurbsSurface};
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::exact_boundary_cap_trim_fit::{
    greville_parameters, interpolate_points_on_fixed_knots, knots_equivalent,
};
use super::exact_boundary_cap_trim_restore::ExactBoundaryCapSurfaceSide;
use trace::{BoundaryTrace, reverse_nurbs_curve, sample_trace, select_trace_boundary_curve};

pub(super) fn relock_exact_boundary_cap_side(
    brep: &mut BRep,
    half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    side: ExactBoundaryCapSurfaceSide,
    side_forward: bool,
    allow_wide_adjacent_row_search: bool,
) -> Option<()> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    let edge = brep.half_edges.get(half_edge)?.edge;
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let edge_forward = brep.edges.get(edge)?.half_edge == half_edge;
    // Row writes replace control points in the side's INCREASING-parameter
    // storage order. `side_forward` converts half-edge order to storage
    // order for the edge-curve fallback; the support trace is oriented by
    // measurement against the current side endpoints below (the flag-based
    // orientation flipped rows whenever the partition reused a cap
    // half-edge running against the solved surface's parameter direction,
    // which scrambled committed-correct caps — PK-124).
    let edge_direction = if side_forward == edge_forward {
        BoundarySampleDirection::Forward
    } else {
        BoundarySampleDirection::Reverse
    };
    let side_endpoints = {
        let surface = brep.faces.get(face)?.surface.as_ref()?;
        let SurfaceGeometry::Nurbs(surface_nurbs) = &surface.geometry else {
            return None;
        };
        side_storage_order_endpoints(surface_nurbs, side)
    };
    let support_trace =
        support_boundary_trace(brep, support_half_edge, BoundarySampleDirection::Forward)
            .map(|trace| orient_trace_to_side_endpoints(trace, side_endpoints));
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
        adjacent::relock_nurbs_adjacent_row(
            surface_nurbs,
            side,
            &boundary,
            trace,
            allow_wide_adjacent_row_search,
        );
    }
    Some(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BoundarySampleDirection {
    Forward,
    Reverse,
}

/// Side endpoint images in increasing-parameter (storage) order.
fn side_storage_order_endpoints(
    surface: &NurbsSurface,
    side: ExactBoundaryCapSurfaceSide,
) -> (
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let (start_uv, end_uv) = match side {
        ExactBoundaryCapSurfaceSide::VMin => ((u_min, v_min), (u_max, v_min)),
        ExactBoundaryCapSurfaceSide::VMax => ((u_min, v_max), (u_max, v_max)),
        ExactBoundaryCapSurfaceSide::UMin => ((u_min, v_min), (u_min, v_max)),
        ExactBoundaryCapSurfaceSide::UMax => ((u_max, v_min), (u_max, v_max)),
    };
    (
        surface.point_at(start_uv.0, start_uv.1),
        surface.point_at(end_uv.0, end_uv.1),
    )
}

/// Orient a forward-sampled support trace into storage order by measuring
/// it against the current side endpoints. The committed cap side is within
/// commit tolerance of the support boundary, so endpoint distances decide
/// the direction robustly; a tie (degenerate side) keeps forward order.
fn orient_trace_to_side_endpoints(
    trace: BoundaryTrace,
    (side_start, side_end): (
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    ),
) -> BoundaryTrace {
    let (Some(first), Some(last)) = (trace.points.first(), trace.points.last()) else {
        return trace;
    };
    let forward_gap = first.distance_to(&side_start) + last.distance_to(&side_end);
    let reverse_gap = first.distance_to(&side_end) + last.distance_to(&side_start);
    if reverse_gap < forward_gap {
        BoundaryTrace {
            points: trace.points.into_iter().rev().collect(),
            normals: trace.normals.into_iter().rev().collect(),
        }
    } else {
        trace
    }
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
    let interpolated = NurbsCurve {
        degree,
        knots: knots.to_vec(),
        weights: vec![1.0; control_points.len()],
        control_points,
    };
    let sampled = NurbsCurve {
        degree,
        knots: knots.to_vec(),
        weights: vec![1.0; points.len()],
        control_points: points,
    };
    Some(select_trace_boundary_curve(trace, interpolated, sampled))
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
