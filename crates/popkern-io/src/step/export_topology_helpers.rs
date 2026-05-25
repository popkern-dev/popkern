//! Helper predicates and recovery utilities for STEP topology export.

use std::collections::HashSet;

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{
    BRep, BrepDegenerateEdgeKind, BrepDegenerateEdgeOptions, EdgeHandle, HalfEdge, HalfEdgeHandle,
    LoopHandle, VertexHandle, analyze_brep_degenerate_edges, edge_segment,
};

use crate::error::IoError;

const SURFACE_TRIM_EDGE_CURVE_SYNTHESIS_SAMPLES: usize = 17;
const STEP_EXPORT_SURFACE_TRIM_TOLERANCE: f64 = 1e-6;

pub(super) fn pcurve_endpoints_match_edge_segment(
    brep: &BRep,
    edge: EdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    let Some((start, end)) = edge_segment(brep, edge) else {
        return false;
    };
    let uv_start = pcurve.point_at(pcurve.t_start);
    let uv_end = pcurve.point_at(pcurve.t_end);
    let point_start = surface.point_at(uv_start.x, uv_start.y);
    let point_end = surface.point_at(uv_end.x, uv_end.y);
    let limit = tolerance.max(brep.tolerance.linear()).max(1e-9);
    let forward = point_start
        .distance_to(&start)
        .max(point_end.distance_to(&end));
    let reverse = point_start
        .distance_to(&end)
        .max(point_end.distance_to(&start));
    forward.min(reverse) <= limit
}

#[cfg(test)]
pub(super) fn edge_half_edges(brep: &BRep, edge_handle: EdgeHandle) -> Vec<HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| (data.edge == edge_handle).then_some(half_edge))
        .collect()
}

pub(super) fn synthesize_missing_edge_curve_from_surface_trims(
    brep: &BRep,
    edge: EdgeHandle,
) -> Option<CurveOnEdge> {
    let tolerance = brep
        .tolerance
        .linear()
        .max(STEP_EXPORT_SURFACE_TRIM_TOLERANCE);
    surface_trim_sources(brep, edge)
        .into_iter()
        .filter_map(|half_edge| surface_trim_edge_curve_candidate(brep, edge, half_edge, tolerance))
        .filter_map(|candidate| {
            let endpoint_gap = edge_curve_endpoint_gap(brep, edge, &candidate)?;
            let sampled_gap = sampled_surface_trim_disagreement(
                brep,
                edge,
                &candidate,
                SURFACE_TRIM_EDGE_CURVE_SYNTHESIS_SAMPLES,
            )?;
            (endpoint_gap <= tolerance && sampled_gap <= tolerance)
                .then_some((candidate, endpoint_gap.max(sampled_gap)))
        })
        .min_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate)
}

pub(super) fn edge_curve_surface_trim_disagreement(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
) -> Option<f64> {
    sampled_surface_trim_disagreement(brep, edge, curve, SURFACE_TRIM_EDGE_CURVE_SYNTHESIS_SAMPLES)
}

pub(super) fn collapsed_trim_missing_curve_edges(brep: &BRep) -> HashSet<EdgeHandle> {
    analyze_brep_degenerate_edges(
        brep,
        BrepDegenerateEdgeOptions {
            tolerance: brep.tolerance.linear(),
            samples: 5,
        },
    )
    .findings
    .into_iter()
    .filter_map(|finding| {
        let missing_curve = brep
            .edges
            .get(finding.edge)
            .is_some_and(|edge| edge.curve.is_none());
        (finding.kind == BrepDegenerateEdgeKind::CollapsedSurfaceTrim && missing_curve)
            .then_some(finding.edge)
    })
    .collect()
}

pub(super) fn edge_curve_same_sense(curve: &CurveOnEdge) -> bool {
    if matches!(
        curve.geometry,
        CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. }
    ) {
        curve.t_end >= curve.t_start
    } else {
        true
    }
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
        SURFACE_TRIM_EDGE_CURVE_SYNTHESIS_SAMPLES,
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

pub(super) fn edge_end_vertex(
    brep: &BRep,
    edge: EdgeHandle,
    he: &HalfEdge,
) -> Result<VertexHandle, IoError> {
    if brep.boundary_edges.contains(&edge) {
        return brep
            .half_edges
            .get(he.next)
            .map(|next| next.origin)
            .ok_or_else(|| {
                IoError::Export("Boundary edge references missing next half-edge".into())
            });
    }
    brep.half_edges
        .get(he.twin)
        .map(|twin| twin.origin)
        .ok_or_else(|| IoError::Export("Half-edge references missing twin".into()))
}

/// Synthesize a plane surface from the first three vertices of a loop.
pub(super) fn synthesize_plane(brep: &BRep, lh: LoopHandle) -> Result<SurfaceOnFace, IoError> {
    let loop_ = brep
        .loops
        .get(lh)
        .ok_or_else(|| IoError::Export("Loop not found for plane synthesis".into()))?;

    let mut positions = Vec::new();
    let start_he = loop_.half_edge;
    let mut current = start_he;
    loop {
        let he = brep
            .half_edges
            .get(current)
            .ok_or_else(|| IoError::Export("Half-edge not found during plane synthesis".into()))?;
        if let Some(v) = brep.vertices.get(he.origin) {
            positions.push(v.position);
        }
        if positions.len() >= 3 {
            break;
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if positions.len() < 3 {
        return Err(IoError::Export(
            "Cannot synthesize plane: loop has fewer than 3 vertices".into(),
        ));
    }

    let origin = positions[0];
    let u_vec = positions[1] - positions[0];
    let v_vec = positions[2] - positions[0];
    let normal = u_vec.cross(&v_vec);
    let normal = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let u_axis = u_vec.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v_axis = normal.cross(&u_axis);

    Ok(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        },
        orientation: true,
    })
}
