//! STEP export pcurve synthesis for face-uses missing stored trims.

use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use popkern_core::geometry::{
    CurveEval, CurveOnEdge, NurbsCurve, Point3, SurfaceEval, SurfaceGeometry,
};
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle};

const EXPORT_PCURVE_SYNTHESIS_SAMPLES: usize = 17;
const STEP_EXPORT_PCURVE_ASSOCIATION_TOLERANCE: f64 = 1e-6;

pub(super) fn synthesize_export_half_edge_pcurve_on_surface(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    let he = brep.half_edges.get(half_edge)?;
    if he.edge != edge {
        return None;
    }
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge)?;
    let samples = EXPORT_PCURVE_SYNTHESIS_SAMPLES;
    let (period_u, period_v) = surface.periodicity();
    let mut prev_uv: Option<(f64, f64)> = None;
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);

    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let point = curve.geometry.point_at(curve_t);
        let (mut u, mut v) = if let Some((prev_u, prev_v)) = prev_uv {
            project_point_to_uv_near(surface, &point, prev_u, prev_v)?
        } else {
            project_point_to_uv(surface, &point)?
        };
        unwrap_periodic_uv(&mut u, &mut v, prev_uv, period_u, period_v);
        prev_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(fraction);
    }

    let nurbs = if matches!(surface, SurfaceGeometry::Plane { .. })
        || uv_points_are_linear(&uv_points, step_export_pcurve_association_tolerance(brep))
    {
        linear_nurbs(*uv_points.first()?, *uv_points.last()?)
    } else {
        let degree = 3.min(uv_points.len().saturating_sub(1));
        fit_curve(&uv_points, &params, degree, uv_points.len()).ok()?
    };
    let pcurve = Pcurve::from_nurbs(nurbs);
    (sampled_edge_surface_disagreement(brep, edge, half_edge, surface, &pcurve, samples)?
        <= step_export_pcurve_association_tolerance(brep))
    .then_some(pcurve)
}

pub(super) fn synthesize_export_face_use_curve_and_pcurve(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<(CurveOnEdge, Pcurve)> {
    let he = brep.half_edges.get(half_edge)?;
    if he.edge != edge {
        return None;
    }
    let next = brep.half_edges.get(he.next)?;
    let start = brep.vertices.get(he.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    let tolerance = step_export_pcurve_association_tolerance(brep);
    let mut start_uv = project_point_to_uv(surface, &start)?;
    let mut end_uv = project_point_to_uv_near(surface, &end, start_uv.0, start_uv.1)
        .or_else(|| project_point_to_uv(surface, &end))?;
    let (period_u, period_v) = surface.periodicity();
    unwrap_periodic_uv(
        &mut end_uv.0,
        &mut end_uv.1,
        Some(start_uv),
        period_u,
        period_v,
    );
    unwrap_periodic_uv(&mut start_uv.0, &mut start_uv.1, None, period_u, period_v);
    let surface_start = surface.point_at(start_uv.0, start_uv.1);
    let surface_end = surface.point_at(end_uv.0, end_uv.1);

    let samples = EXPORT_PCURVE_SYNTHESIS_SAMPLES;
    let mut points = Vec::with_capacity(samples);
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let u = start_uv.0 + (end_uv.0 - start_uv.0) * fraction;
        let v = start_uv.1 + (end_uv.1 - start_uv.1) * fraction;
        let point = if index == 0 {
            surface_start
        } else if index + 1 == samples {
            surface_end
        } else {
            surface.point_at(u, v)
        };
        points.push(point);
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(fraction);
    }
    let curve = CurveOnEdge {
        geometry: face_use_curve_geometry(&points, &params, tolerance)?,
        t_start: 0.0,
        t_end: 1.0,
    };
    let pcurve = Pcurve::from_nurbs(if uv_points_are_linear(&uv_points, tolerance) {
        linear_nurbs(*uv_points.first()?, *uv_points.last()?)
    } else {
        let degree = 3.min(uv_points.len().saturating_sub(1));
        fit_curve(&uv_points, &params, degree, uv_points.len()).ok()?
    });
    (sampled_edge_surface_disagreement_for_curve(&curve, surface, &pcurve, samples) <= tolerance)
        .then_some((curve, pcurve))
}

pub(super) fn step_export_pcurve_association_tolerance(brep: &BRep) -> f64 {
    // The gate decides whether a face use shares the global edge or
    // gets a face-local duplicate (which downgrades the shell to an
    // OPEN_SHELL surface model that consumers will not auto-sew).
    // Coherence is therefore judged at the model's DECLARED accuracy
    // — its accumulated imprecision — not at the raw linear
    // tolerance: an edge within the model tolerance is stitchable by
    // any consumer, and exporting it single-sided guarantees the
    // opposite (PK-096: Fusion showed unstitched faces at the cap
    // because 8.7e-6 rim residuals tripped a 1e-6 gate).
    brep.tolerance
        .linear()
        .max(brep.max_imprecision)
        .max(STEP_EXPORT_PCURVE_ASSOCIATION_TOLERANCE)
}

pub(super) fn export_pcurve_trace_matches_edge(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    if half_edge_data.edge != edge {
        return false;
    }
    let Some(curve) = brep.edges.get(edge).and_then(|edge| edge.curve.as_ref()) else {
        return false;
    };
    edge_curve_to_pcurve_trace_gap(curve, pcurve, surface, 17) <= tolerance
}

fn edge_curve_to_pcurve_trace_gap(
    curve: &CurveOnEdge,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    samples: usize,
) -> f64 {
    let samples = samples.max(2);
    let mut max_gap = 0.0_f64;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
        let point = curve.geometry.point_at(t);
        max_gap = max_gap.max(pcurve_trace_distance(
            pcurve,
            surface,
            &point,
            samples * 2 + 1,
        ));
    }
    max_gap
}

fn pcurve_trace_distance(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    point: &Point3,
    samples: usize,
) -> f64 {
    let samples = samples.max(2);
    let (t0, t1) = pcurve.parameter_range();
    if !t0.is_finite() || !t1.is_finite() {
        return f64::INFINITY;
    }
    let mut best = f64::INFINITY;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = t0 + (t1 - t0) * fraction;
        let uv = pcurve.point_at(t);
        let surface_point = surface.point_at(uv.x, uv.y);
        best = best.min(surface_point.distance_to(point));
    }
    best
}

fn half_edge_curve_forward(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<bool> {
    let he = brep.half_edges.get(half_edge)?;
    let curve = brep.edges.get(he.edge)?.curve.as_ref()?;
    let origin = brep.vertices.get(he.origin)?.position;
    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&start) <= origin.distance_to(&end))
}

fn sampled_edge_surface_disagreement(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    samples: usize,
) -> Option<f64> {
    let he = brep.half_edges.get(half_edge)?;
    if he.edge != edge {
        return None;
    }
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let forward = half_edge_curve_forward(brep, half_edge)?;
    let mut max_gap = 0.0_f64;
    for index in 0..samples.max(2) {
        let fraction = index as f64 / (samples.max(2) - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
        let uv = pcurve.point_at(pcurve_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        max_gap = max_gap.max(curve.geometry.point_at(curve_t).distance_to(&surface_point));
    }
    Some(max_gap)
}

fn sampled_edge_surface_disagreement_for_curve(
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    samples: usize,
) -> f64 {
    let samples = samples.max(2);
    let mut max_gap = 0.0_f64;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
        let uv = pcurve.point_at(pcurve_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        max_gap = max_gap.max(curve.geometry.point_at(curve_t).distance_to(&surface_point));
    }
    max_gap
}

fn face_use_curve_geometry(
    points: &[Point3],
    params: &[f64],
    tolerance: f64,
) -> Option<popkern_core::geometry::curve::CurveGeometry> {
    if points_are_linear(points, tolerance) {
        let start = *points.first()?;
        let end = *points.last()?;
        return Some(popkern_core::geometry::curve::CurveGeometry::Line {
            origin: start,
            direction: end - start,
        });
    }
    let degree = 3.min(points.len().saturating_sub(1));
    fit_curve(points, params, degree, points.len())
        .ok()
        .map(popkern_core::geometry::curve::CurveGeometry::Nurbs)
}

fn points_are_linear(points: &[Point3], tolerance: f64) -> bool {
    let Some(start) = points.first().copied() else {
        return false;
    };
    let Some(end) = points.last().copied() else {
        return false;
    };
    let direction = end - start;
    let denom = direction.dot(&direction);
    if denom <= tolerance * tolerance {
        return false;
    }
    points.iter().enumerate().all(|(index, point)| {
        let fraction = index as f64 / points.len().saturating_sub(1).max(1) as f64;
        let expected = start + direction * fraction;
        point.distance_to(&expected) <= tolerance
    })
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

fn uv_points_are_linear(points: &[Point3], tolerance: f64) -> bool {
    let Some(start) = points.first().copied() else {
        return false;
    };
    let Some(end) = points.last().copied() else {
        return false;
    };
    let direction = end - start;
    let denom = direction.dot(&direction);
    if denom <= tolerance * tolerance {
        return false;
    }
    points.iter().enumerate().all(|(index, point)| {
        let fraction = index as f64 / points.len().saturating_sub(1).max(1) as f64;
        let expected = start + direction * fraction;
        point.distance_to(&expected) <= tolerance
    })
}

fn linear_nurbs(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}
