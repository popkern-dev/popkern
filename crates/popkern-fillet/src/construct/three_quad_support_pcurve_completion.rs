//! Cap/support seam completion for Gregory faces after surface assignment.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, nurbs_curve::NurbsCurve};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::helpers;

const SAMPLE_COUNT: usize = 12;
const STALE_CAP_PCURVE_TOLERANCE: f64 = 1.0e-5;

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
    let edge_curve = brep
        .edges
        .get(edge)
        .and_then(|edge| edge.curve.as_ref())
        .cloned();
    let support_pcurve = brep.half_edges.get(twin).and_then(|he| he.pcurve.as_ref());
    let support_owns_seam =
        matches!(support_surface, SurfaceGeometry::Cylinder { .. }) && support_pcurve.is_some();
    let Some(samples) =
        cap_oriented_samples(brep, half_edge, twin, &support_surface, support_owns_seam)
    else {
        return;
    };
    let mut support_samples = samples.clone();
    support_samples.reverse();
    let needs_edge_curve = edge_curve.as_ref().is_none_or(|curve| {
        support_owns_seam && edge_curve_sample_gap(brep, half_edge, curve, &samples) > 1.0e-5
    });
    let cap_boundary_owned = brep.surface_boundary_hint(half_edge).is_some();
    let needs_cap_pcurve = !cap_boundary_owned
        && brep.half_edges[half_edge]
            .pcurve
            .as_ref()
            .is_none_or(|pcurve| {
                cap_pcurve_sample_gap(brep, half_edge, cap_surface, pcurve, &samples)
                    > STALE_CAP_PCURVE_TOLERANCE
            });
    let needs_support_pcurve = !needs_edge_curve
        && match brep.half_edges[twin].pcurve.as_ref() {
            None => true,
            Some(_) if matches!(support_surface, SurfaceGeometry::Nurbs(_)) => false,
            Some(pcurve) => {
                cap_pcurve_sample_gap(brep, twin, &support_surface, pcurve, &support_samples)
                    > STALE_CAP_PCURVE_TOLERANCE
            }
        };
    if !needs_cap_pcurve && !needs_support_pcurve && !needs_edge_curve {
        return;
    }
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
        brep.half_edges[half_edge].pcurve = Some(best_pcurve_for_surface_samples(
            cap_surface,
            &uv_samples,
            &samples,
        ));
    }
    if needs_support_pcurve {
        let uv_samples = support_samples
            .iter()
            .map(|point| unclamped_uv_point(&support_surface, point))
            .collect::<Vec<_>>();
        brep.half_edges[twin].pcurve = Some(best_pcurve_for_surface_samples(
            &support_surface,
            &uv_samples,
            &support_samples,
        ));
    }
}

fn cap_oriented_samples(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    prefer_support_pcurve: bool,
) -> Option<Vec<Point3>> {
    if prefer_support_pcurve
        && let Some(support_pcurve) = brep.half_edges.get(twin)?.pcurve.as_ref()
    {
        return Some(samples_from_support_pcurve(
            brep,
            half_edge,
            twin,
            support_surface,
            support_pcurve,
        ));
    }
    if let Some(curve) = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.as_ref())
    {
        return Some(samples_from_curve_for_half_edge(brep, half_edge, curve));
    }
    if let Some(support_pcurve) = brep.half_edges.get(twin)?.pcurve.as_ref() {
        return Some(samples_from_support_pcurve(
            brep,
            half_edge,
            twin,
            support_surface,
            support_pcurve,
        ));
    }
    None
}

fn edge_curve_sample_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    target_samples: &[Point3],
) -> f64 {
    samples_from_curve_for_half_edge(brep, half_edge, curve)
        .iter()
        .zip(target_samples)
        .map(|(sample, target)| sample.distance_to(target))
        .fold(0.0, f64::max)
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

fn unclamped_uv_point(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (u, v) = helpers::project_to_uv(surface, point);
    Point3::new(u, v, 0.0)
}

fn fit_nurbs_curve(points: &[Point3]) -> Option<NurbsCurve> {
    let params = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    fit_curve(points, &params, 3, points.len()).ok()
}

fn best_pcurve_for_surface_samples(
    surface: &SurfaceGeometry,
    uv_samples: &[Point3],
    target_samples: &[Point3],
) -> Pcurve {
    let linear = linear_pcurve_from_endpoints(uv_samples);
    if uv_samples_are_collinear(uv_samples) {
        return linear;
    }
    let linear_gap = pcurve_surface_sample_gap(surface, &linear, target_samples);
    let Some(fitted) = fit_nurbs_curve(uv_samples).map(Pcurve::from) else {
        return linear;
    };
    let fitted_gap = pcurve_surface_sample_gap(surface, &fitted, target_samples);
    if fitted_gap <= linear_gap {
        fitted
    } else {
        linear
    }
}

fn uv_samples_are_collinear(samples: &[Point3]) -> bool {
    if samples.len() < 3 {
        return true;
    }
    let start = samples[0];
    let end = samples[samples.len() - 1];
    let dx = end.x - start.x;
    let dy = end.y - start.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len <= 1.0e-12 {
        return samples.iter().all(|sample| {
            let sx = sample.x - start.x;
            let sy = sample.y - start.y;
            (sx * sx + sy * sy).sqrt() <= 1.0e-10
        });
    }
    samples.iter().all(|sample| {
        let sx = sample.x - start.x;
        let sy = sample.y - start.y;
        (sx * dy - sy * dx).abs() / len <= 1.0e-8
    })
}

fn pcurve_surface_sample_gap(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    samples: &[Point3],
) -> f64 {
    if samples.len() < 2 {
        return f64::INFINITY;
    }
    let (t0, t1) = pcurve.parameter_range();
    samples
        .iter()
        .enumerate()
        .map(|(index, sample)| {
            let s = index as f64 / (samples.len() - 1) as f64;
            let point = cap_surface_point(surface, pcurve, t0 + s * (t1 - t0));
            if point_is_finite(&point) {
                point.distance_to(sample)
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}

fn cap_pcurve_sample_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    cap_pcurve: &Pcurve,
    samples: &[Point3],
) -> f64 {
    if samples.len() < 2 {
        return f64::INFINITY;
    }
    let cap_start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let (t0, t1) = cap_pcurve.parameter_range();
    let p0 = cap_surface_point(cap_surface, cap_pcurve, t0);
    let p1 = cap_surface_point(cap_surface, cap_pcurve, t1);
    if !point_is_finite(&p0) || !point_is_finite(&p1) {
        return f64::INFINITY;
    }
    let forward = p0.distance_to(&cap_start) <= p1.distance_to(&cap_start);
    samples
        .iter()
        .enumerate()
        .map(|(index, sample)| {
            let s = index as f64 / (samples.len() - 1) as f64;
            let local = if forward { s } else { 1.0 - s };
            let point = cap_surface_point(cap_surface, cap_pcurve, t0 + local * (t1 - t0));
            if point_is_finite(&point) {
                point.distance_to(sample)
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}

fn cap_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}

fn point_is_finite(point: &Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
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
        complete_cap_support_half_edge(brep, cap_face, cap_surface, half_edge);
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
