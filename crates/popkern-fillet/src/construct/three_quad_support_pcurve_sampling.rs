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

/// Last-resort sample source for a seam edge carrying no curve and no
/// pcurve on either side (measured: cube_extreme raw-corner cap seam).
/// The seam lies on the analytic support surface, so the endpoint
/// chord projected onto it is the best geometric authority available;
/// the downstream seam g0 acceptance gate remains the arbiter of the
/// authored geometry. Declines non-analytic supports (projection is
/// not exact there) and periodic branch jumps (an adjacent sample
/// step larger than the whole chord).
pub(super) fn chord_samples_projected_to_support(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
) -> Option<Vec<Point3>> {
    if !matches!(
        support_surface,
        SurfaceGeometry::Plane { .. } | SurfaceGeometry::Cylinder { .. }
    ) {
        return None;
    }
    let start = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let end = brep
        .vertices
        .get(brep.half_edges.get(twin)?.origin)?
        .position;
    let chord = start.distance_to(&end);
    if chord <= 0.0 {
        return None;
    }
    let samples: Vec<Point3> = (0..SAMPLE_COUNT)
        .map(|index| {
            let s = index as f64 / (SAMPLE_COUNT - 1) as f64;
            let p = Point3::new(
                start.x + (end.x - start.x) * s,
                start.y + (end.y - start.y) * s,
                start.z + (end.z - start.z) * s,
            );
            let (u, v) = super::helpers::project_to_uv(support_surface, &p);
            support_surface.point_at(u, v)
        })
        .collect();
    let max_step = samples
        .windows(2)
        .map(|pair| pair[0].distance_to(&pair[1]))
        .fold(0.0, f64::max);
    (max_step <= chord).then_some(samples)
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
