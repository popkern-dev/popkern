use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::SAMPLE_COUNT;

pub(super) fn samples_from_curve_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Vec<Point3> {
    let cap_start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let forward = p_start.distance_to(&cap_start) <= p_end.distance_to(&cap_start);
    (0..SAMPLE_COUNT)
        .map(|index| {
            let s = index as f64 / (SAMPLE_COUNT - 1) as f64;
            let local = if forward { s } else { 1.0 - s };
            let t = curve.t_start + local * (curve.t_end - curve.t_start);
            curve.geometry.point_at(t)
        })
        .collect()
}

pub(super) fn samples_from_support_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
) -> Vec<Point3> {
    let cap_start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let cap_end = brep.vertices[brep.half_edges[twin].origin].position;
    let (t0, t1) = support_pcurve.parameter_range();
    let p0 = support_surface_point(support_surface, support_pcurve, t0);
    let p1 = support_surface_point(support_surface, support_pcurve, t1);
    let forward = p0.distance_to(&cap_start) + p1.distance_to(&cap_end)
        <= p1.distance_to(&cap_start) + p0.distance_to(&cap_end);
    (0..SAMPLE_COUNT)
        .map(|index| {
            let s = index as f64 / (SAMPLE_COUNT - 1) as f64;
            let local = if forward { s } else { 1.0 - s };
            support_surface_point(support_surface, support_pcurve, t0 + local * (t1 - t0))
        })
        .collect()
}

pub(super) fn support_pcurve_matches_edge_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
    edge_curve: &CurveOnEdge,
    tolerance: f64,
) -> bool {
    let support_samples =
        samples_from_support_pcurve(brep, half_edge, twin, support_surface, support_pcurve);
    let edge_samples = samples_from_curve_for_half_edge(brep, half_edge, edge_curve);
    support_samples
        .iter()
        .zip(edge_samples.iter())
        .all(|(support, edge)| support.distance_to(edge) <= tolerance)
}

pub(super) fn support_pcurve_endpoints_match(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
    tolerance: f64,
) -> bool {
    support_pcurve_endpoint_gap(brep, half_edge, twin, support_surface, support_pcurve) <= tolerance
}

fn support_pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
) -> f64 {
    let cap_start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let cap_end = brep.vertices[brep.half_edges[twin].origin].position;
    let (t0, t1) = support_pcurve.parameter_range();
    let p0 = support_surface_point(support_surface, support_pcurve, t0);
    let p1 = support_surface_point(support_surface, support_pcurve, t1);
    let forward_gap = p0.distance_to(&cap_start).max(p1.distance_to(&cap_end));
    let reverse_gap = p1.distance_to(&cap_start).max(p0.distance_to(&cap_end));
    forward_gap.min(reverse_gap)
}

fn support_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}
