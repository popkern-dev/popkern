//! Sampled pcurve/surface measurement helpers for the exact-boundary cap
//! trim restore pass: candidate support tangency scoring, pcurve
//! endpoint/trace gaps, and interior-clamped surface normal comparison.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::{half_edge_destination_vertex, surface_for_half_edge};

pub(super) fn support_normal_deficit_for_candidate(
    brep: &BRep,
    support_half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    cap_pcurve: &Pcurve,
) -> Option<f64> {
    let support_pcurve = brep.half_edges.get(support_half_edge)?.pcurve.as_ref()?;
    let support_surface = surface_for_half_edge(brep, support_half_edge)?;
    let (cap_start_t, cap_end_t) = cap_pcurve.parameter_range();
    let (support_start_t, support_end_t) = support_pcurve.parameter_range();
    let cap_start = pcurve_surface_point(cap_surface, cap_pcurve, cap_start_t);
    let cap_end = pcurve_surface_point(cap_surface, cap_pcurve, cap_end_t);
    let support_start = pcurve_surface_point(&support_surface, support_pcurve, support_start_t);
    let support_end = pcurve_surface_point(&support_surface, support_pcurve, support_end_t);
    let support_forward = cap_start.distance_to(&support_start) + cap_end.distance_to(&support_end)
        <= cap_start.distance_to(&support_end) + cap_end.distance_to(&support_start);
    let mut worst = 0.0_f64;
    for index in 0..17 {
        let parameter = 0.05 + 0.90 * index as f64 / 16.0;
        let support_parameter = if support_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let cap_uv = cap_pcurve.point_at(cap_start_t + (cap_end_t - cap_start_t) * parameter);
        let support_uv = support_pcurve
            .point_at(support_start_t + (support_end_t - support_start_t) * support_parameter);
        let deficit = interior_normal_deficit_between_surfaces(
            cap_surface,
            cap_uv.x,
            cap_uv.y,
            &support_surface,
            support_uv.x,
            support_uv.y,
        )?;
        worst = worst.max(deficit);
    }
    Some(worst)
}

pub(super) fn pcurve_surface_point(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

pub(super) fn interior_normal_deficit_between_surfaces(
    cap_surface: &SurfaceGeometry,
    cap_u: f64,
    cap_v: f64,
    support_surface: &SurfaceGeometry,
    support_u: f64,
    support_v: f64,
) -> Option<f64> {
    let (cap_u, cap_v) = crate::chart_sampling::periodic_interior_uv(cap_surface, cap_u, cap_v);
    let (support_u, support_v) =
        crate::chart_sampling::periodic_interior_uv(support_surface, support_u, support_v);
    let cap_normal = cap_surface
        .partial_u(cap_u, cap_v)
        .cross(&cap_surface.partial_v(cap_u, cap_v));
    let support_normal = support_surface
        .partial_u(support_u, support_v)
        .cross(&support_surface.partial_v(support_u, support_v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length <= 1.0e-300 || support_length <= 1.0e-300 {
        return None;
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

pub(super) fn pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let destination = brep
        .vertices
        .get(half_edge_destination_vertex(brep, half_edge)?)?
        .position;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    Some(
        start
            .distance_to(&origin)
            .max(end.distance_to(&destination)),
    )
}

pub(super) fn pcurve_curve_trace_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    curve: &CurveOnEdge,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    let forward = origin.distance_to(&curve_start) <= origin.distance_to(&curve_end);
    let (pc_start, pc_end) = pcurve.parameter_range();
    let mut worst = 0.0_f64;
    for index in 0..17 {
        let fraction = index as f64 / 16.0;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pc_t = pc_start + (pc_end - pc_start) * fraction;
        let uv = pcurve.point_at(pc_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        let gap = surface_point.distance_to(&curve.geometry.point_at(curve_t));
        if !gap.is_finite() {
            return None;
        }
        worst = worst.max(gap);
    }
    Some(worst)
}

pub(super) fn reverse_linear_pcurve(pcurve: &Pcurve) -> Pcurve {
    let (t_start, t_end) = pcurve.parameter_range();
    let start = pcurve.point_at(t_start);
    let end = pcurve.point_at(t_end);
    Pcurve::linear((end.x, end.y), (start.x, start.y), t_start, t_end)
}
