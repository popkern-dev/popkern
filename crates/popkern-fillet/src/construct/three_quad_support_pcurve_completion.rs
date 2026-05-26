//! Cap/support seam completion for Gregory faces after surface assignment.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, nurbs_curve::NurbsCurve};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::helpers;

const SAMPLE_COUNT: usize = 12;

pub(crate) fn complete_missing_cap_support_seams(brep: &mut BRep, cap_face: FaceHandle) {
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return;
    };
    if !matches!(cap_surface, SurfaceGeometry::Gregory(_)) {
        return;
    }
    for half_edge in cap_face_half_edges(brep, cap_face) {
        complete_cap_support_half_edge(brep, cap_face, &cap_surface, half_edge);
    }
    complete_direct_cap_half_edges(brep, cap_face, &cap_surface);
}

fn complete_cap_support_half_edge(
    brep: &mut BRep,
    cap_face: FaceHandle,
    cap_surface: &SurfaceGeometry,
    half_edge: HalfEdgeHandle,
) {
    let Some(half_edge_data) = brep.half_edges.get(half_edge).cloned() else {
        return;
    };
    let twin = half_edge_data.twin;
    let Some(twin_data) = brep.half_edges.get(twin).cloned() else {
        return;
    };
    let Some(support_face) = brep
        .loops
        .get(twin_data.loop_)
        .map(|loop_data| loop_data.face)
    else {
        return;
    };
    if support_face == cap_face {
        return;
    }
    let Some(support_surface) = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return;
    };
    let edge = half_edge_data.edge;
    let needs_cap_pcurve = brep.half_edges[half_edge].pcurve.is_none();
    let needs_edge_curve = brep
        .edges
        .get(edge)
        .and_then(|edge| edge.curve.as_ref())
        .is_none();
    if !needs_cap_pcurve && !needs_edge_curve {
        return;
    }
    let Some(samples) = cap_oriented_samples(brep, half_edge, twin, &support_surface) else {
        return;
    };
    if needs_edge_curve {
        let mut edge_samples = samples.clone();
        if brep
            .edges
            .get(edge)
            .is_some_and(|edge_data| edge_data.half_edge != half_edge)
        {
            edge_samples.reverse();
        }
        let geometry = fit_nurbs_curve(&edge_samples)
            .map(CurveGeometry::Nurbs)
            .unwrap_or_else(|| line_geometry_from_endpoints(&edge_samples));
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry,
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    if needs_cap_pcurve {
        let uv_samples = samples
            .iter()
            .map(|point| cap_uv_point(cap_surface, point))
            .collect::<Vec<_>>();
        let pcurve = fit_nurbs_curve(&uv_samples)
            .map(Pcurve::from)
            .unwrap_or_else(|| linear_pcurve_from_endpoints(&uv_samples));
        brep.half_edges[half_edge].pcurve = Some(pcurve);
    }
}

fn cap_oriented_samples(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
) -> Option<Vec<Point3>> {
    if let Some(curve) = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.as_ref())
    {
        return Some(samples_from_curve_for_half_edge(brep, half_edge, curve));
    }
    let support_pcurve = brep.half_edges.get(twin)?.pcurve.as_ref()?;
    Some(samples_from_support_pcurve(
        brep,
        half_edge,
        twin,
        support_surface,
        support_pcurve,
    ))
}

fn samples_from_curve_for_half_edge(
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

fn samples_from_support_pcurve(
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

fn support_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}

fn cap_uv_point(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (mut u, mut v) = helpers::project_to_uv(surface, point);
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
    if u_lo.is_finite() && u_hi.is_finite() {
        u = u.clamp(u_lo, u_hi);
    }
    if v_lo.is_finite() && v_hi.is_finite() {
        v = v.clamp(v_lo, v_hi);
    }
    Point3::new(u, v, 0.0)
}

fn fit_nurbs_curve(points: &[Point3]) -> Option<NurbsCurve> {
    let params = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    fit_curve(points, &params, 3, points.len()).ok()
}

fn linear_pcurve_from_endpoints(points: &[Point3]) -> Pcurve {
    let start = points[0];
    let end = points[points.len() - 1];
    Pcurve::linear((start.x, start.y), (end.x, end.y), 0.0, 1.0)
}

fn line_geometry_from_endpoints(points: &[Point3]) -> CurveGeometry {
    let start = points[0];
    let end = points[points.len() - 1];
    CurveGeometry::Line {
        origin: start,
        direction: end - start,
    }
}

fn complete_direct_cap_half_edges(
    brep: &mut BRep,
    cap_face: FaceHandle,
    cap_surface: &SurfaceGeometry,
) {
    let half_edges = brep
        .half_edges
        .iter()
        .filter_map(|(half_edge, half_edge_data)| {
            let twin = half_edge_data.twin;
            let face = brep
                .loops
                .get(half_edge_data.loop_)
                .map(|loop_data| loop_data.face)?;
            let twin_face = brep
                .half_edges
                .get(twin)
                .and_then(|twin_data| brep.loops.get(twin_data.loop_))
                .map(|loop_data| loop_data.face)?;
            (face == cap_face && twin_face != cap_face).then_some(half_edge)
        })
        .collect::<Vec<_>>();
    for half_edge in half_edges {
        let half_edge_data = brep.half_edges[half_edge].clone();
        let twin = half_edge_data.twin;
        let start = brep.vertices[half_edge_data.origin].position;
        let end = brep.vertices[brep.half_edges[twin].origin].position;
        if brep.half_edges[half_edge].pcurve.is_none() {
            let uv_start = cap_uv_point(cap_surface, &start);
            let uv_end = cap_uv_point(cap_surface, &end);
            brep.half_edges[half_edge].pcurve = Some(Pcurve::linear(
                (uv_start.x, uv_start.y),
                (uv_end.x, uv_end.y),
                0.0,
                1.0,
            ));
        }
        if brep
            .edges
            .get(half_edge_data.edge)
            .and_then(|edge| edge.curve.as_ref())
            .is_none()
        {
            brep.edges[half_edge_data.edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: start,
                    direction: end - start,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
    }
}

fn cap_face_half_edges(brep: &BRep, cap_face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    for (_, edge_data) in &brep.edges {
        let first = edge_data.half_edge;
        let Some(first_data) = brep.half_edges.get(first) else {
            continue;
        };
        let second = first_data.twin;
        let Some(second_data) = brep.half_edges.get(second) else {
            continue;
        };
        let first_face = brep
            .loops
            .get(first_data.loop_)
            .map(|loop_data| loop_data.face);
        let second_face = brep
            .loops
            .get(second_data.loop_)
            .map(|loop_data| loop_data.face);
        match (first_face, second_face) {
            (Some(face), Some(other)) if face == cap_face && other != cap_face => out.push(first),
            (Some(other), Some(face)) if face == cap_face && other != cap_face => out.push(second),
            _ => {}
        };
    }
    out
}
