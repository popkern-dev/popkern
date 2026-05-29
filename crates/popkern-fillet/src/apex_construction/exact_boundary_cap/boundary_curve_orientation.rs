//! Boundary curve endpoint orientation helpers.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, HalfEdgeHandle};

pub(super) fn cap_half_edge_vertex_span(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let origin = brep.vertices.get(half_edge_data.origin)?.position;
    let destination_vertex = brep
        .half_edges
        .get(half_edge_data.twin)
        .map(|twin| twin.origin)
        .or_else(|| {
            brep.half_edges
                .get(half_edge_data.next)
                .map(|next| next.origin)
        })?;
    let destination = brep.vertices.get(destination_vertex)?.position;
    Some((origin, destination))
}

pub(super) fn endpoint_orientation_gap(
    curve: &NurbsCurve,
    origin: Point3,
    destination: Point3,
) -> f64 {
    let (start, end) = curve.parameter_range();
    curve
        .point_at(start)
        .distance_to(&origin)
        .max(curve.point_at(end).distance_to(&destination))
}

pub(super) fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    let first = *curve.knots.first().unwrap_or(&0.0);
    let last = *curve.knots.last().unwrap_or(&1.0);
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|k| first + last - *k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points: Vec<Point3> = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}
