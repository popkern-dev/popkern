//! Edge-curve repair from a single authoritative face-side pcurve.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

pub(crate) fn repair_measured_edges_from_single_surface_pcurves(
    brep: &mut BRep,
    tolerance: f64,
    max_anchor_motion: f64,
) -> usize {
    let edges = brep.edges.keys().collect::<Vec<_>>();
    let mut repaired = 0;
    for edge in edges {
        let Some(report) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if report.max_deviation <= tolerance {
            continue;
        }
        let candidates = report
            .half_edges
            .iter()
            .filter(|half_edge| half_edge.max_deviation > tolerance)
            .map(|half_edge| half_edge.half_edge)
            .collect::<Vec<_>>();
        if candidates.into_iter().any(|half_edge| {
            repair_edge_from_single_surface_pcurve(
                brep,
                edge,
                half_edge,
                tolerance,
                max_anchor_motion,
            )
        }) {
            repaired += 1;
        }
    }
    repaired
}

fn repair_edge_from_single_surface_pcurve(
    brep: &mut BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    tolerance: f64,
    max_anchor_motion: f64,
) -> bool {
    let Some(points) = single_surface_points_in_edge_direction(brep, edge, half_edge, 257) else {
        return false;
    };
    let Some(first) = points.first().copied() else {
        return false;
    };
    let Some(last) = points.last().copied() else {
        return false;
    };
    let Some((start_vertex, end_vertex)) = edge_vertices(brep, edge) else {
        return false;
    };
    let current_start = brep.vertices[start_vertex].position;
    let current_end = brep.vertices[end_vertex].position;
    let anchor_motion = current_start
        .distance_to(&first)
        .max(current_end.distance_to(&last));
    if anchor_motion > max_anchor_motion {
        return false;
    }
    let incident_gap = incident_endpoint_gap_after_vertex_move(brep, edge, start_vertex, first)
        .max(incident_endpoint_gap_after_vertex_move(
            brep, edge, end_vertex, last,
        ));
    if incident_gap > tolerance {
        return false;
    }
    let Some(curve) = curve_from_surface_points(&points, tolerance) else {
        return false;
    };
    let snapshot = brep.clone();
    brep.vertices[start_vertex].position = first;
    brep.vertices[end_vertex].position = last;
    brep.edges[edge].curve = Some(curve);
    let accepted = brep
        .measure_edge_pcurve_deviation(edge)
        .is_some_and(|updated| updated.max_deviation <= tolerance)
        && brep.validate().is_ok();
    if accepted {
        true
    } else {
        *brep = snapshot;
        false
    }
}

fn single_surface_points_in_edge_direction(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    samples: usize,
) -> Option<Vec<Point3>> {
    let edge_data = brep.edges.get(edge)?;
    let half_edge_data = brep.half_edges.get(half_edge)?;
    if half_edge_data.edge != edge {
        return None;
    }
    let pcurve = half_edge_data.pcurve.as_ref()?;
    let surface = surface_for_half_edge(brep, half_edge)?;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let directed_fraction = if half_edge == edge_data.half_edge {
            fraction
        } else {
            1.0 - fraction
        };
        let parameter = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * directed_fraction;
        let uv = pcurve.point_at(parameter);
        points.push(surface.point_at(uv.x, uv.y));
    }
    Some(points)
}

fn curve_from_surface_points(points: &[Point3], tolerance: f64) -> Option<CurveOnEdge> {
    let start = *points.first()?;
    let end = *points.last()?;
    if max_line_deviation(points, start, end) <= tolerance
        && let Some(line) = CurveOnEdge::line_segment(start, end, tolerance)
    {
        return Some(line);
    }
    let parameters = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    let degree = 5.min(points.len().saturating_sub(1));
    if let Ok(curve) =
        popkern_core::geometry::fitting::interpolate_curve(points, &parameters, degree)
    {
        return Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve),
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(polyline_curve(points)?),
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn polyline_curve(points: &[Point3]) -> Option<NurbsCurve> {
    if points.len() < 2 {
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
            point.distance_to(&(start + direction * fraction))
        })
        .fold(0.0, f64::max)
}

fn edge_vertices(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    Some((half_edge.origin, next.origin))
}

fn incident_endpoint_gap_after_vertex_move(
    brep: &BRep,
    moved_edge: EdgeHandle,
    vertex: VertexHandle,
    proposed_position: Point3,
) -> f64 {
    brep.edges
        .iter()
        .filter_map(|(edge, edge_data)| {
            if edge == moved_edge {
                return None;
            }
            let curve = edge_data.curve.as_ref()?;
            let half_edge = brep.half_edges.get(edge_data.half_edge)?;
            let end_vertex = brep.half_edges.get(half_edge.next)?.origin;
            if half_edge.origin == vertex || end_vertex == vertex {
                Some(
                    curve
                        .geometry
                        .point_at(curve.t_start)
                        .distance_to(&proposed_position)
                        .min(
                            curve
                                .geometry
                                .point_at(curve.t_end)
                                .distance_to(&proposed_position),
                        ),
                )
            } else {
                None
            }
        })
        .fold(0.0, f64::max)
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<SurfaceGeometry> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    let surface = brep.faces.get(face)?.surface.as_ref()?.geometry.clone();
    Some(surface)
}
