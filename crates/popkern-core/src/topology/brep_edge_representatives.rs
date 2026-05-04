//! Edge representative normalization for topology rewrites.

use super::BRep;

impl BRep {
    /// Restore the invariant that each edge's representative half-edge points
    /// back to that same edge, and remove edge records with no half-edges.
    pub fn normalize_edge_representatives(&mut self) {
        let mut representative_by_edge = std::collections::BTreeMap::new();
        for (half_edge_handle, half_edge) in &self.half_edges {
            representative_by_edge
                .entry(half_edge.edge)
                .or_insert(half_edge_handle);
        }

        let edge_handles = self.edges.keys().collect::<Vec<_>>();
        for edge_handle in edge_handles {
            let Some(representative) = representative_by_edge.get(&edge_handle).copied() else {
                self.edges.remove(edge_handle);
                self.boundary_edges.remove(&edge_handle);
                continue;
            };
            let primary_is_valid = self
                .edges
                .get(edge_handle)
                .and_then(|edge| self.half_edges.get(edge.half_edge))
                .is_some_and(|half_edge| half_edge.edge == edge_handle);
            if !primary_is_valid && let Some(edge) = self.edges.get_mut(edge_handle) {
                edge.half_edge = representative;
            }
        }
    }
}
