//! Geometry helpers for missing edge-curve completion.

use crate::geometry::{CurveEval, CurveGeometry, CurveOnEdge, NurbsCurve, Point3, SurfaceEval};
use crate::healing::deviation::edge_surface_deviation;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, edge_segment};

use super::{BrepEdgeCurveCompletionDeclineReason, BrepEdgeCurveCompletionOptions};

pub(super) struct SurfaceTrimSource {
    pub(super) half_edge: HalfEdgeHandle,
}

pub(super) fn curve_from_vertices(
    brep: &BRep,
    edge: EdgeHandle,
    options: BrepEdgeCurveCompletionOptions,
) -> Result<CurveOnEdge, BrepEdgeCurveCompletionDeclineReason> {
    if surfaced_uses(brep, edge) > 0 {
        return Err(BrepEdgeCurveCompletionDeclineReason::MissingSurfaceTrimSource);
    }
    if !options.allow_vertex_segments_without_surface {
        return Err(BrepEdgeCurveCompletionDeclineReason::VertexSegmentFallbackDisabled);
    }
    let (start, end) =
        edge_segment(brep, edge).ok_or(BrepEdgeCurveCompletionDeclineReason::MissingTopology)?;
    CurveOnEdge::line_segment(start, end, options.tolerance)
        .ok_or(BrepEdgeCurveCompletionDeclineReason::DegenerateEdge)
}

pub(super) fn curve_from_surface_trim(
    brep: &BRep,
    edge: EdgeHandle,
    source: &SurfaceTrimSource,
    options: BrepEdgeCurveCompletionOptions,
) -> Result<CurveOnEdge, BrepEdgeCurveCompletionDeclineReason> {
    let points = sample_surface_trim_in_edge_direction(brep, edge, source, options.samples)
        .ok_or(BrepEdgeCurveCompletionDeclineReason::MissingSurfaceTrimSource)?;
    let (start, end) =
        edge_segment(brep, edge).ok_or(BrepEdgeCurveCompletionDeclineReason::MissingTopology)?;
    curve_from_points(&points, start, end, options.tolerance)
        .ok_or(BrepEdgeCurveCompletionDeclineReason::DegenerateEdge)
}

pub(super) fn validate_candidate(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
    sources: &[SurfaceTrimSource],
    options: BrepEdgeCurveCompletionOptions,
) -> Option<(f64, f64)> {
    let endpoint_gap = validate_candidate_endpoints(brep, edge, curve)?;
    let mut sampled_disagreement = 0.0_f64;
    for source in sources {
        sampled_disagreement = sampled_disagreement.max(sampled_source_disagreement(
            brep,
            edge,
            curve,
            source,
            options.samples,
        )?);
    }
    Some((endpoint_gap, sampled_disagreement))
}

pub(super) fn validated_anchored_surface_trim_curve(
    brep: &BRep,
    edge: EdgeHandle,
    sources: &[SurfaceTrimSource],
    options: BrepEdgeCurveCompletionOptions,
) -> Option<(CurveOnEdge, f64)> {
    let source = sources.first()?;
    let mut points = sample_surface_trim_in_edge_direction(brep, edge, source, options.samples)?;
    let (start, end) = edge_segment(brep, edge)?;
    if points.len() < 2 || start.distance_to(&end) <= options.tolerance.max(1e-12) {
        return None;
    }
    points[0] = start;
    let last = points.len() - 1;
    points[last] = end;
    let candidate = curve_from_points(&points, start, end, options.tolerance)?;
    let endpoint_gap = validate_candidate_endpoints(brep, edge, &candidate)?;
    if endpoint_gap > options.tolerance {
        return None;
    }
    let max_surface_deviation =
        max_direct_surface_deviation(brep, &candidate, sources, options.samples)?;
    (max_surface_deviation <= options.tolerance).then_some((candidate, max_surface_deviation))
}

pub(super) fn validated_vertex_curve_for_surfaced_edge(
    brep: &BRep,
    edge: EdgeHandle,
    sources: &[SurfaceTrimSource],
    options: BrepEdgeCurveCompletionOptions,
) -> Option<CurveOnEdge> {
    if sources.is_empty() {
        return None;
    }
    let (start, end) = edge_segment(brep, edge)?;
    let candidate = CurveOnEdge::line_segment(start, end, options.tolerance)?;
    let endpoint_gap = validate_candidate_endpoints(brep, edge, &candidate)?;
    if endpoint_gap > options.tolerance {
        return None;
    }
    let max_surface_deviation =
        max_direct_surface_deviation(brep, &candidate, sources, options.samples)?;
    (max_surface_deviation <= options.tolerance).then_some(candidate)
}

pub(super) fn surface_trim_sources(brep: &BRep, edge: EdgeHandle) -> Vec<SurfaceTrimSource> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| {
            if data.edge != edge || data.pcurve.is_none() {
                return None;
            }
            let face = brep.faces.get(brep.loops.get(data.loop_)?.face)?;
            face.surface
                .is_some()
                .then_some(SurfaceTrimSource { half_edge })
        })
        .collect()
}

fn curve_from_points(
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

fn max_direct_surface_deviation(
    brep: &BRep,
    curve: &CurveOnEdge,
    sources: &[SurfaceTrimSource],
    samples: usize,
) -> Option<f64> {
    let mut max = 0.0_f64;
    for source in sources {
        let half_edge = brep.half_edges.get(source.half_edge)?;
        let face = brep.loops.get(half_edge.loop_)?.face;
        let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        let deviation = edge_surface_deviation(
            &curve.geometry,
            (curve.t_start, curve.t_end),
            surface,
            samples,
        );
        max = max.max(deviation.max_deviation);
    }
    Some(max)
}

fn validate_candidate_endpoints(brep: &BRep, edge: EdgeHandle, curve: &CurveOnEdge) -> Option<f64> {
    let (start, end) = edge_segment(brep, edge)?;
    Some(
        curve
            .geometry
            .point_at(curve.t_start)
            .distance_to(&start)
            .max(curve.geometry.point_at(curve.t_end).distance_to(&end)),
    )
}

fn sampled_source_disagreement(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
    source: &SurfaceTrimSource,
    samples: usize,
) -> Option<f64> {
    let samples = samples.max(2);
    let source_points = sample_surface_trim_in_edge_direction(brep, edge, source, samples)?;
    let mut max = 0.0_f64;
    for (index, surface_point) in source_points.iter().enumerate() {
        let f = index as f64 / (samples - 1) as f64;
        let t = curve.t_start + (curve.t_end - curve.t_start) * f;
        max = max.max(curve.geometry.point_at(t).distance_to(surface_point));
    }
    Some(max)
}

fn sample_surface_trim_in_edge_direction(
    brep: &BRep,
    edge: EdgeHandle,
    source: &SurfaceTrimSource,
    samples: usize,
) -> Option<Vec<Point3>> {
    let edge_data = brep.edges.get(edge)?;
    let half_edge = brep.half_edges.get(source.half_edge)?;
    let pcurve = half_edge.pcurve.as_ref()?;
    let face = brep.faces.get(brep.loops.get(half_edge.loop_)?.face)?;
    let surface = &face.surface.as_ref()?.geometry;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for i in 0..samples {
        let f = i as f64 / (samples - 1) as f64;
        let directed_f = if source.half_edge == edge_data.half_edge {
            f
        } else {
            1.0 - f
        };
        let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * directed_f;
        let uv = pcurve.point_at(t);
        points.push(surface.point_at(uv.x, uv.y));
    }
    Some(points)
}

fn surfaced_uses(brep: &BRep, edge: EdgeHandle) -> usize {
    brep.half_edges
        .iter()
        .filter(|(_, data)| {
            data.edge == edge
                && brep
                    .loops
                    .get(data.loop_)
                    .and_then(|loop_| brep.faces.get(loop_.face))
                    .is_some_and(|face| face.surface.is_some())
        })
        .count()
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
            let f = index as f64 / (points.len().saturating_sub(1).max(1)) as f64;
            let line_point = start + direction * f;
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
