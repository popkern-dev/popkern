//! STEP surface-curve normalization from validated pcurve/surface evidence.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, edge_segment};

use crate::step::diagnostics::{DiagnosticLog, Phase};

const STEP_SURFACE_CURVE_RECONCILE_SAMPLES: usize = 17;

pub(super) fn reconcile_step_surface_curve_edges(
    brep: &mut BRep,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) {
    let limit = tolerance.linear().max(1e-9);
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    let mut repaired_edges = 0usize;
    let mut declined_edges = 0usize;
    let mut max_before = 0.0_f64;
    let mut max_after = 0.0_f64;

    for edge in edge_handles {
        let Some(before) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if before.max_deviation <= limit {
            continue;
        }
        max_before = max_before.max(before.max_deviation);
        let Some(candidate) = best_surface_trim_edge_curve_candidate(brep, edge, limit) else {
            declined_edges += 1;
            continue;
        };
        let mut trial = brep.clone();
        trial.edges[edge].curve = Some(candidate);
        let Some(after) = trial.measure_edge_pcurve_deviation(edge) else {
            declined_edges += 1;
            continue;
        };
        let endpoint_gap = trial
            .edges
            .get(edge)
            .and_then(|edge_data| edge_data.curve.as_ref())
            .and_then(|curve| edge_curve_endpoint_gap(&trial, edge, curve))
            .unwrap_or(f64::INFINITY);
        if endpoint_gap <= limit
            && after.max_deviation <= limit
            && after.max_deviation < before.max_deviation
        {
            *brep = trial;
            repaired_edges += 1;
            max_after = max_after.max(after.max_deviation);
        } else {
            declined_edges += 1;
        }
    }

    if repaired_edges > 0 {
        log.warn(
            Phase::MapTopology,
            None,
            format!(
                "normalized {repaired_edges} STEP SURFACE_CURVE edge curve(s) from validated pcurve/surface evidence; max deviation {:.3e} -> {:.3e}",
                max_before, max_after
            ),
        );
    }
    if declined_edges > 0 {
        log.warn(
            Phase::MapTopology,
            None,
            format!(
                "declined {declined_edges} STEP SURFACE_CURVE edge normalization candidate(s); pcurve/surface evidence did not satisfy {:.3e}",
                limit
            ),
        );
    }
}

fn best_surface_trim_edge_curve_candidate(
    brep: &BRep,
    edge: EdgeHandle,
    tolerance: f64,
) -> Option<CurveOnEdge> {
    surface_trim_sources(brep, edge)
        .into_iter()
        .filter_map(|half_edge| surface_trim_edge_curve_candidate(brep, edge, half_edge, tolerance))
        .filter_map(|candidate| {
            let endpoint_gap = edge_curve_endpoint_gap(brep, edge, &candidate)?;
            let sampled_gap = sampled_surface_trim_disagreement(
                brep,
                edge,
                &candidate,
                STEP_SURFACE_CURVE_RECONCILE_SAMPLES,
            )?;
            (endpoint_gap <= tolerance && sampled_gap <= tolerance)
                .then_some((candidate, endpoint_gap.max(sampled_gap)))
        })
        .min_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate)
}

fn surface_trim_sources(brep: &BRep, edge: EdgeHandle) -> Vec<HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| {
            if data.edge != edge || data.pcurve.is_none() {
                return None;
            }
            let face = brep.faces.get(brep.loops.get(data.loop_)?.face)?;
            face.surface.is_some().then_some(half_edge)
        })
        .collect()
}

fn surface_trim_edge_curve_candidate(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    tolerance: f64,
) -> Option<CurveOnEdge> {
    let points = sample_surface_trim_in_edge_direction(
        brep,
        edge,
        half_edge,
        STEP_SURFACE_CURVE_RECONCILE_SAMPLES,
        tolerance,
    )?;
    let (start, end) = edge_segment(brep, edge)?;
    curve_from_surface_points(&points, start, end, tolerance)
}

fn sampled_surface_trim_disagreement(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
    samples: usize,
) -> Option<f64> {
    let mut max_gap = 0.0_f64;
    for half_edge in surface_trim_sources(brep, edge) {
        let surface_points =
            sample_surface_trim_in_edge_direction(brep, edge, half_edge, samples, f64::INFINITY)?;
        for (index, surface_point) in surface_points.iter().enumerate() {
            let fraction = index as f64 / (surface_points.len().saturating_sub(1).max(1)) as f64;
            let curve_t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
            max_gap = max_gap.max(curve.geometry.point_at(curve_t).distance_to(surface_point));
        }
    }
    Some(max_gap)
}

fn sample_surface_trim_in_edge_direction(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    samples: usize,
    endpoint_tolerance: f64,
) -> Option<Vec<Point3>> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let pcurve = half_edge_data.pcurve.as_ref()?;
    let face = brep.faces.get(brep.loops.get(half_edge_data.loop_)?.face)?;
    let surface = &face.surface.as_ref()?.geometry;
    let (start, end) = edge_segment(brep, edge)?;
    let forward_gap =
        pcurve_surface_endpoint_gap(pcurve, surface, pcurve.t_start, pcurve.t_end, start, end);
    let reverse_gap =
        pcurve_surface_endpoint_gap(pcurve, surface, pcurve.t_end, pcurve.t_start, start, end);
    if forward_gap.min(reverse_gap) > endpoint_tolerance {
        return None;
    }
    let reverse = reverse_gap < forward_gap;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let directed_fraction = if reverse { 1.0 - fraction } else { fraction };
        let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * directed_fraction;
        let uv = pcurve.point_at(t);
        points.push(surface.point_at(uv.x, uv.y));
    }
    if endpoint_tolerance.is_finite() {
        points[0] = start;
        let last = points.len() - 1;
        points[last] = end;
    }
    Some(points)
}

fn pcurve_surface_endpoint_gap(
    pcurve: &Pcurve,
    surface: &dyn SurfaceEval,
    start_t: f64,
    end_t: f64,
    start: Point3,
    end: Point3,
) -> f64 {
    let start_uv = pcurve.point_at(start_t);
    let end_uv = pcurve.point_at(end_t);
    surface
        .point_at(start_uv.x, start_uv.y)
        .distance_to(&start)
        .max(surface.point_at(end_uv.x, end_uv.y).distance_to(&end))
}

fn curve_from_surface_points(
    points: &[Point3],
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<CurveOnEdge> {
    if max_line_deviation(points, start, end) <= tolerance
        && let Some(line) = CurveOnEdge::line_segment(start, end, tolerance)
    {
        return Some(line);
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(polyline_nurbs(points)?),
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn edge_curve_endpoint_gap(brep: &BRep, edge: EdgeHandle, curve: &CurveOnEdge) -> Option<f64> {
    let (start, end) = edge_segment(brep, edge)?;
    Some(
        curve
            .geometry
            .point_at(curve.t_start)
            .distance_to(&start)
            .max(curve.geometry.point_at(curve.t_end).distance_to(&end)),
    )
}

fn max_line_deviation(points: &[Point3], start: Point3, end: Point3) -> f64 {
    let direction = end - start;
    let denom = direction.dot(&direction);
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    points
        .iter()
        .enumerate()
        .map(|(index, point)| {
            let fraction = index as f64 / (points.len().saturating_sub(1).max(1)) as f64;
            let line_point = start + direction * fraction;
            point.distance_to(&line_point)
        })
        .fold(0.0_f64, f64::max)
}

fn polyline_nurbs(points: &[Point3]) -> Option<NurbsCurve> {
    if points.len() < 2 {
        return None;
    }
    if points
        .windows(2)
        .map(|pair| pair[0].distance_to(&pair[1]))
        .fold(0.0_f64, f64::max)
        <= 0.0
    {
        return None;
    }
    let count = points.len();
    let mut knots = Vec::with_capacity(count + 2);
    knots.push(0.0);
    knots.push(0.0);
    for index in 1..count.saturating_sub(1) {
        knots.push(index as f64 / (count - 1) as f64);
    }
    knots.push(1.0);
    knots.push(1.0);
    Some(NurbsCurve {
        degree: 1,
        knots,
        control_points: points.to_vec(),
        weights: vec![1.0; count],
    })
}
