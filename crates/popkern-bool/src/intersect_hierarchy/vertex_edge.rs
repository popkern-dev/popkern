//! Vertex-edge contact detection between two BRep bodies.
//!
//! A vertex lies on an edge when the minimum distance from the vertex
//! to the edge curve is within the snap tolerance. The engine uses
//! `CurveEval::closest_point` to find the nearest parameter, then
//! checks the distance.

use popkern_core::geometry::CurveEval;
use popkern_core::topology::BRep;

use super::result::{Operand, VertexEdgeHit};

/// Find all vertex-edge contacts between A's vertices and B's edges,
/// and B's vertices and A's edges.
pub(crate) fn detect_vertex_edge(a: &BRep, b: &BRep, snap_tol: f64) -> Vec<VertexEdgeHit> {
    let mut hits = Vec::new();

    // A's vertices vs B's edges
    detect_one_direction(a, b, Operand::A, snap_tol, &mut hits);
    // B's vertices vs A's edges
    detect_one_direction(b, a, Operand::B, snap_tol, &mut hits);

    hits
}

/// Check vertices of `vert_brep` against edges of `edge_brep`.
fn detect_one_direction(
    vert_brep: &BRep,
    edge_brep: &BRep,
    vert_body: Operand,
    snap_tol: f64,
    hits: &mut Vec<VertexEdgeHit>,
) {
    for (vh, vert) in &vert_brep.vertices {
        let vp = vert.position;

        for (eh, edge) in &edge_brep.edges {
            let curve_on_edge = match &edge.curve {
                Some(c) => c,
                None => continue,
            };

            let curve = &curve_on_edge.geometry;
            let (t_start, t_end) = (curve_on_edge.t_start, curve_on_edge.t_end);

            // Use closest_point to find the nearest parameter
            let t = match curve.closest_point(&vp) {
                Some(t) => t,
                None => continue,
            };

            // Clamp to the edge's parameter range
            let t_clamped = t.clamp(t_start.min(t_end), t_start.max(t_end));

            let closest = curve.point_at(t_clamped);
            let dist = closest.distance_to(&vp);

            if dist < snap_tol {
                // Skip if the closest point is at an endpoint (that's a VV hit)
                let edge_len = (t_start - t_end).abs();
                let at_start = (t_clamped - t_start).abs() < edge_len * 0.01;
                let at_end = (t_clamped - t_end).abs() < edge_len * 0.01;
                if at_start || at_end {
                    continue;
                }

                hits.push(VertexEdgeHit {
                    vertex: vh,
                    vertex_body: vert_body,
                    edge: eh,
                    t_edge: t_clamped,
                    distance: dist,
                });
            }
        }
    }
}
