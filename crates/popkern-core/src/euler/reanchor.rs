//! Vertex `outgoing`-anchor repair for half-edge removal surgery.
//!
//! `Vertex::outgoing` is DERIVED navigational state: any live incident
//! half-edge is a valid anchor, and every surgical site that removes
//! half-edges must re-anchor the affected vertices or their fan walks
//! silently skip them forever (measured: `BRep::torus`'s seam stitch
//! left one vertex dangling and the shell vertex recompute served a
//! stale position; the armed `validate_vertex_consistency` check now
//! rejects the class). `kve`/`mef` carry their own local repairs; the
//! remaining removal sites route through this shared helper.

use crate::topology::{BRep, HalfEdgeHandle};

impl BRep {
    /// Re-anchor every vertex whose `outgoing` points at one of the
    /// half-edges in `doomed`, choosing a surviving half-edge with the
    /// same origin (deterministic: iteration order of the slot map is
    /// slot order, stable for a given body). Call BEFORE removing the
    /// half-edges — or use [`BRep::remove_half_edge_reanchored`], which
    /// composes the repair with the removal. A vertex with no surviving
    /// incident half-edge keeps its stale anchor — it is topologically
    /// isolated and exempt from the consistency check.
    pub fn reanchor_vertex_outgoing(&mut self, doomed: &[HalfEdgeHandle]) {
        for &h in doomed {
            let Some(he) = self.half_edges.get(h) else {
                continue;
            };
            let v = he.origin;
            let Some(vertex) = self.vertices.get(v) else {
                continue;
            };
            if !doomed.contains(&vertex.outgoing) {
                continue;
            }
            if let Some((replacement, _)) = self
                .half_edges
                .iter()
                .find(|(hh, cand)| cand.origin == v && !doomed.contains(hh))
            {
                self.vertices[v].outgoing = replacement;
            }
        }
    }

    /// Remove a half-edge, first re-anchoring its origin vertex if its
    /// `outgoing` pointed at it. This is the ONLY sanctioned way to
    /// remove half-edges outside the anchor-aware Euler operators —
    /// enforced by `structural_half_edge_removal.rs`. Returns the
    /// removed half-edge like `SlotMap::remove`.
    pub fn remove_half_edge_reanchored(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Option<crate::topology::HalfEdge> {
        self.reanchor_vertex_outgoing(&[half_edge]);
        self.half_edges.remove(half_edge)
    }
}
