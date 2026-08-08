//! Whole-body orientation reversal (material-side mirror).
//!
//! Reverses every loop's traversal direction AND flips every stored
//! face orientation flag together, so the interior-left winding
//! invariant (outer loop CCW in UV ⇔ flag true) is preserved by
//! construction while the body's material side mirrors. This is the
//! honest form of the "clone the skin and flip the flags" mirror
//! that shell/offset builders need: flipping flags alone leaves the
//! cloned winding contradicting the flag on every face (the PK-244
//! dishonest-mirror class).
//!
//! Invariants and contract:
//! - Entity counts, handles, loop membership, and twin pairing are
//!   all unchanged — the operation is a pure re-orientation, an
//!   involution (applying it twice restores the original body).
//! - Per half-edge: `next`/`prev` swap and the origin moves to the
//!   old destination (the twin's origin). Twins remain opposite.
//! - Pcurves keep their UV image (the trim is identical); their
//!   parameterization direction is NOT rewritten, matching the
//!   repo-wide "consumers verify direction by endpoints" contract on
//!   [`HalfEdge::pcurve`](crate::topology::half_edge::HalfEdge).
//! - Vertex `outgoing` pointers are re-anchored (the old outgoing
//!   half-edge's origin moves away from the vertex; its twin's
//!   moves onto it).
//! - Edge curves are edge-level (undirected) and unchanged.
//! - Half-edges whose loop has been deleted (e.g. a skin builder
//!   removed a face but kept the shared edge ring) are reversed the
//!   same way — the operation is purely local, so dangling halves
//!   stay mutually consistent with their live twins.

use crate::topology::BRep;

impl BRep {
    /// Mirror the body's material side: reverse every loop's
    /// traversal and flip every face's orientation flag together.
    ///
    /// Preserves winding-flag honesty by construction (both the
    /// winding and the flag negate). Involution: calling twice
    /// restores the original connectivity and flags exactly.
    pub fn reverse_body_orientation(&mut self) {
        // 1. Per half-edge: swap next/prev, move origin to the old
        //    destination. Collect new origins first so the update
        //    does not observe partially-mutated state.
        let new_origins: Vec<(
            crate::topology::HalfEdgeHandle,
            crate::topology::VertexHandle,
        )> = self
            .half_edges
            .iter()
            .filter_map(|(handle, he)| {
                self.half_edges
                    .get(he.twin)
                    .map(|twin| (handle, twin.origin))
            })
            .collect();
        for (handle, origin) in new_origins {
            let he = &mut self.half_edges[handle];
            he.origin = origin;
            std::mem::swap(&mut he.next, &mut he.prev);
        }

        // 2. Re-anchor vertex outgoing pointers: the twin of the old
        //    outgoing half-edge now originates at this vertex.
        let new_outgoing: Vec<(
            crate::topology::VertexHandle,
            crate::topology::HalfEdgeHandle,
        )> = self
            .vertices
            .iter()
            .filter_map(|(vh, vertex)| self.half_edges.get(vertex.outgoing).map(|he| (vh, he.twin)))
            .collect();
        for (vh, outgoing) in new_outgoing {
            self.vertices[vh].outgoing = outgoing;
        }

        // 3. Flip every stored orientation flag in step with the
        //    winding so honesty is preserved.
        for (_, face) in &mut self.faces {
            if let Some(sof) = &mut face.surface {
                sof.orientation = !sof.orientation;
            }
        }
    }
}

#[cfg(test)]
#[path = "brep_reverse_orientation_tests.rs"]
mod brep_reverse_orientation_tests;
