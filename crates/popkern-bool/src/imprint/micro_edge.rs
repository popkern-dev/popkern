//! Micro-edge detection and collapse.
//!
//! A micro-edge is shorter than the geometric tolerance — too short to
//! carry useful information. Collapse by merging endpoints.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// A detected micro-edge.
#[derive(Debug, Clone)]
pub struct MicroEdge {
    /// The edge handle.
    pub edge: EdgeHandle,
    /// Length of the edge.
    pub length: f64,
    /// Start vertex.
    pub start: VertexHandle,
    /// End vertex.
    pub end: VertexHandle,
}

/// Find all micro-edges in a BRep.
pub fn find_micro_edges(brep: &BRep, threshold: f64) -> Vec<MicroEdge> {
    let mut result = Vec::new();

    for (eh, edge) in &brep.edges {
        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let p_start = brep.vertices[he.origin].position;
        let p_end = brep.vertices[twin.origin].position;
        let length = p_start.distance_to(&p_end);

        if length < threshold {
            result.push(MicroEdge {
                edge: eh,
                length,
                start: he.origin,
                end: twin.origin,
            });
        }
    }

    result
}

/// Collapse micro-edges by merging their endpoints.
///
/// For each micro-edge, moves the end vertex to the start vertex's position
/// (or midpoint) and updates topology references. Returns the number of
/// edges collapsed.
///
/// Note: actual topological collapse (removing the edge and merging
/// vertices in the half-edge structure) requires KEF/KEMR which may
/// corrupt complex topology. This function only adjusts vertex positions
/// as a geometric repair, not a topological collapse.
pub fn collapse_micro_edge_positions(brep: &mut BRep, micro_edges: &[MicroEdge]) -> usize {
    let mut collapsed = 0;

    for me in micro_edges {
        let Some(v_start) = brep.vertices.get(me.start) else {
            continue;
        };
        let Some(v_end) = brep.vertices.get(me.end) else {
            continue;
        };

        // Move end vertex to midpoint
        let midpoint = Point3::new(
            (v_start.position.x + v_end.position.x) * 0.5,
            (v_start.position.y + v_end.position.y) * 0.5,
            (v_start.position.z + v_end.position.z) * 0.5,
        );

        brep.vertices[me.start].position = midpoint;
        brep.vertices[me.end].position = midpoint;
        collapsed += 1;
    }

    collapsed
}
