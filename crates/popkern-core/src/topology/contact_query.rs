//! Where a body touches itself, whatever its form (PK-1110, S11).
//!
//! A general body holds its contacts in its topology: a ring of more
//! than two face uses is an edge where solids meet, a vertex with more
//! than one use a point where they do (PK-1115: nothing is recorded).

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::topology::{BRep, ContactEntity, TopologyContact, TopologyContactKind};

impl BRep {
    /// Every contact of the body, named by the shared edge or vertex: each edge
    /// with more than two face uses, and each vertex with more than one
    /// use or ending such an edge (pieces meeting along an edge meet at
    /// its ends too, as the split form records them).
    pub fn contacts(&self) -> Vec<TopologyContact> {
        let mut out = Vec::new();
        let mut ends = std::collections::BTreeSet::new();
        for edge in self.edges.keys() {
            if self.ring_summary(edge).count <= 2 {
                continue;
            }
            if let Some((a, b)) = self.edge_vertices(edge) {
                ends.insert(a);
                ends.insert(b);
            }
            let location = match self.edges[edge].curve.as_ref() {
                Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
                None => self
                    .edge_vertices(edge)
                    .map_or(Point3::new(0.0, 0.0, 0.0), |(a, b)| {
                        let (a, b) = (self.vertices[a].position, self.vertices[b].position);
                        a + (b - a) * 0.5
                    }),
            };
            out.push(TopologyContact {
                kind: TopologyContactKind::Edge,
                location,
                entities: vec![ContactEntity::Edge(edge)],
            });
        }
        for (vertex, uses) in self.vertex_uses_all() {
            if uses.len() > 1 || ends.contains(&vertex) {
                out.push(TopologyContact {
                    kind: TopologyContactKind::Vertex,
                    location: self.vertices[vertex].position,
                    entities: vec![ContactEntity::Vertex(vertex)],
                });
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "contact_query_tests.rs"]
mod tests;
