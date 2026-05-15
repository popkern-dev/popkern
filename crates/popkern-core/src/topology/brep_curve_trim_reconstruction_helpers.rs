//! Helper functions for curve-backed trim reconstruction.

use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use crate::geometry::{CurveEval, NurbsCurve, Point3, SurfaceEval};
use crate::topology::{BRep, HalfEdgeHandle};

pub(super) fn project_edge_curve_to_surface_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    samples: usize,
) -> Option<crate::geometry::pcurve::Pcurve> {
    let he = brep.half_edges.get(half_edge)?;
    let surface = &brep
        .faces
        .get(brep.loops.get(he.loop_)?.face)?
        .surface
        .as_ref()?
        .geometry;
    let curve = brep.edges.get(he.edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge)?;
    let samples = samples.max(2);
    let (period_u, period_v) = surface.periodicity();
    let mut prev_uv: Option<(f64, f64)> = None;
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let s = index as f64 / (samples - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let t = curve.t_start + (curve.t_end - curve.t_start) * curve_s;
        let p = curve.geometry.point_at(t);
        let (mut u, mut v) = if let Some((pu, pv)) = prev_uv {
            project_point_to_uv_near(surface, &p, pu, pv)?
        } else {
            project_point_to_uv(surface, &p)?
        };
        unwrap_periodic_uv(&mut u, &mut v, prev_uv, period_u, period_v);
        prev_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3.min(uv_points.len().saturating_sub(1));
    let curve = if uv_points.len() == 2 {
        linear_nurbs(uv_points[0], uv_points[1])
    } else {
        crate::geometry::fitting::fit_curve(&uv_points, &params, degree, uv_points.len()).ok()?
    };
    Some(crate::geometry::pcurve::Pcurve::from_nurbs(curve))
}

pub(super) fn sampled_edge_surface_disagreement(
    curve: &dyn CurveEval,
    curve_range: (f64, f64),
    half_edge: HalfEdgeHandle,
    brep: &BRep,
    pcurve: &crate::geometry::pcurve::Pcurve,
    samples: usize,
) -> Option<f64> {
    let surface = &brep
        .faces
        .get(brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face)?
        .surface
        .as_ref()?
        .geometry;
    let forward = half_edge_curve_forward(brep, half_edge)?;
    let samples = samples.max(2);
    let mut max = 0.0_f64;
    for index in 0..samples {
        let s = index as f64 / (samples - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let curve_t = curve_range.0 + (curve_range.1 - curve_range.0) * curve_s;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * s;
        let uv = pcurve.point_at(pcurve_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        max = max.max(curve.point_at(curve_t).distance_to(&surface_point));
    }
    Some(max)
}

pub(super) fn half_edge_endpoints(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let he = brep.half_edges.get(half_edge)?;
    let next = brep.half_edges.get(he.next)?;
    Some((
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(next.origin)?.position,
    ))
}

pub(super) fn half_edge_curve_forward(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<bool> {
    let he = brep.half_edges.get(half_edge)?;
    let curve = brep.edges.get(he.edge)?.curve.as_ref()?;
    let origin = brep.vertices.get(he.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&p0) <= origin.distance_to(&p1))
}

fn unwrap_periodic_uv(
    u: &mut f64,
    v: &mut f64,
    prev_uv: Option<(f64, f64)>,
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some((prev_u, prev_v)) = prev_uv else {
        return;
    };
    if let Some(period) = period_u {
        while *u - prev_u > period * 0.5 {
            *u -= period;
        }
        while prev_u - *u > period * 0.5 {
            *u += period;
        }
    }
    if let Some(period) = period_v {
        while *v - prev_v > period * 0.5 {
            *v -= period;
        }
        while prev_v - *v > period * 0.5 {
            *v += period;
        }
    }
}

fn linear_nurbs(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}
