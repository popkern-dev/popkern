//! FaceEdge boundary curve recovery for DirectApex Gregory cap assembly.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::super::constructive_geometry::cylinder_arc_through_edge_endpoints;
use super::super::gregory_cap_boundary::half_edge_for_face_on_edge;

pub(super) fn oriented_support_face_edge_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    support_face: FaceHandle,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Option<NurbsCurve> {
    oriented_unique_planned_fillet_curve(brep, cap_face, edge, fillet_surfaces)
        .or_else(|| oriented_unique_cylinder_face_edge_curve(brep, cap_face, edge))
        .or_else(|| oriented_support_surface_curve(brep, cap_face, edge, support_face))
        .or_else(|| oriented_support_trim_curve(brep, cap_face, edge))
}

fn oriented_unique_planned_fillet_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Option<NurbsCurve> {
    let (start, end) = cap_edge_endpoints(brep, cap_face, edge)?;
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let mut candidates = fillet_surfaces
        .iter()
        .filter_map(|fillet| {
            let (geometry, t_start, t_end) =
                cylinder_arc_through_edge_endpoints(&fillet.surface, start, end)?;
            let curve = geometry.to_nurbs_in_range((t_start, t_end));
            curve_matches_endpoints(&curve, start, end, tolerance).then_some(curve)
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        candidates.pop()
    } else {
        None
    }
}

fn oriented_unique_cylinder_face_edge_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
) -> Option<NurbsCurve> {
    let (start, end) = cap_edge_endpoints(brep, cap_face, edge)?;
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let mut candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cap_face {
                return None;
            }
            let surface = &face_data.surface.as_ref()?.geometry;
            let (geometry, t_start, t_end) =
                cylinder_arc_through_edge_endpoints(surface, start, end)?;
            let curve = geometry.to_nurbs_in_range((t_start, t_end));
            curve_matches_endpoints(&curve, start, end, tolerance).then_some(curve)
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        candidates.pop()
    } else {
        None
    }
}

fn oriented_support_surface_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    support_face: FaceHandle,
) -> Option<NurbsCurve> {
    let (start, end) = cap_edge_endpoints(brep, cap_face, edge)?;
    let surface = &brep.faces.get(support_face)?.surface.as_ref()?.geometry;
    let (geometry, t_start, t_end) = cylinder_arc_through_edge_endpoints(surface, start, end)?;
    let curve = geometry.to_nurbs_in_range((t_start, t_end));
    curve_matches_endpoints(&curve, start, end, brep.tolerance.linear().max(1.0e-6))
        .then_some(curve)
}

fn cap_edge_endpoints(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
) -> Option<(Point3, Point3)> {
    let cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face)?;
    let support_half_edge = brep.half_edges.get(cap_half_edge)?.twin;
    let start = brep
        .vertices
        .get(brep.half_edges[cap_half_edge].origin)?
        .position;
    let end = brep
        .vertices
        .get(brep.half_edges.get(support_half_edge)?.origin)?
        .position;
    Some((start, end))
}

fn oriented_support_trim_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
) -> Option<NurbsCurve> {
    let cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face)?;
    let support_half_edge = brep.half_edges.get(cap_half_edge)?.twin;
    let support_half_edge_data = brep.half_edges.get(support_half_edge)?;
    let support_face = brep.loops.get(support_half_edge_data.loop_)?.face;
    let support_surface = brep.faces.get(support_face)?.surface.as_ref()?;
    let support_pcurve = support_half_edge_data.pcurve.as_ref()?;
    let start = brep
        .vertices
        .get(brep.half_edges[cap_half_edge].origin)?
        .position;
    let end = brep.vertices.get(support_half_edge_data.origin)?.position;
    fit_oriented_support_trim(
        support_surface,
        support_pcurve,
        start,
        end,
        brep.tolerance.linear().max(1.0e-6),
    )
}

fn fit_oriented_support_trim(
    support_surface: &popkern_core::geometry::surface::SurfaceOnFace,
    support_pcurve: &popkern_core::geometry::pcurve::Pcurve,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<NurbsCurve> {
    let (t0, t1) = support_pcurve.parameter_range();
    let endpoint_at = |t: f64| {
        let uv = support_pcurve.point_at(t);
        support_surface.geometry.point_at(uv.x, uv.y)
    };
    let p0 = endpoint_at(t0);
    let p1 = endpoint_at(t1);
    let forward = p0.distance_to(&start) + p1.distance_to(&end)
        <= p1.distance_to(&start) + p0.distance_to(&end);
    let n_samples = 40usize;
    let mut samples = (0..n_samples)
        .map(|index| {
            let s = index as f64 / (n_samples - 1) as f64;
            let local = if forward { s } else { 1.0 - s };
            endpoint_at(t0 + local * (t1 - t0))
        })
        .collect::<Vec<_>>();
    samples[0] = start;
    samples[n_samples - 1] = end;
    let params = (0..n_samples)
        .map(|index| index as f64 / (n_samples - 1) as f64)
        .collect::<Vec<_>>();
    let curve =
        popkern_core::geometry::fitting::fit_curve(&samples, &params, 3, samples.len()).ok()?;
    curve_matches_endpoints(&curve, start, end, tolerance).then_some(curve)
}

fn curve_matches_endpoints(curve: &NurbsCurve, start: Point3, end: Point3, tolerance: f64) -> bool {
    let (t0, t1) = curve.parameter_range();
    curve.point_at(t0).distance_to(&start) <= tolerance
        && curve.point_at(t1).distance_to(&end) <= tolerance
}
