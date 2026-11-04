//! Micro-edge detection.
//!
//! A micro-edge is shorter than the geometric tolerance — too short to
//! carry useful information. The quality report lists them; the vertex
//! table resolves candidate vertices together so a boolean does not
//! make them.

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
