//! Same-loop bridge dissolution for retained vertices.

use crate::error::CoreError;
use crate::topology::{BRep, HalfEdgeHandle};

impl BRep {
    /// Dissolve a same-loop bridge edge without removing either endpoint vertex.
    ///
    /// This is separate from `kev`: both endpoint vertices remain live and only
    /// the bridge edge pair is removed from the owning loop.
    pub fn dissolve_same_loop_bridge(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<(), CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let twin_handle = he_data.twin;
        let edge_handle = he_data.edge;
        let twin_data = self
            .half_edges
            .get(twin_handle)
            .ok_or_else(|| CoreError::InvalidHandle("twin handle is invalid".into()))?;

        if he_data.loop_ != twin_data.loop_ {
            return Err(CoreError::PreconditionFailed(
                "edge's half-edges are on different loops; use kef instead".into(),
            ));
        }
        if he_data.next == twin_handle && twin_data.next == half_edge {
            return Err(CoreError::PreconditionFailed(
                "dissolve would leave empty loop (bridge is loop's only edge)".into(),
            ));
        }

        let he_origin = he_data.origin;
        let twin_origin = twin_data.origin;
        let he_prev = he_data.prev;
        let he_next = he_data.next;
        let twin_prev = twin_data.prev;
        let twin_next = twin_data.next;
        let the_loop = he_data.loop_;

        let replacement_loop_half_edge = if he_next == twin_handle {
            self.half_edges[he_prev].next = twin_next;
            self.half_edges[twin_next].prev = he_prev;
            twin_next
        } else if twin_next == half_edge {
            self.half_edges[twin_prev].next = he_next;
            self.half_edges[he_next].prev = twin_prev;
            he_next
        } else {
            self.half_edges[he_prev].next = he_next;
            self.half_edges[he_next].prev = he_prev;
            self.half_edges[twin_prev].next = twin_next;
            self.half_edges[twin_next].prev = twin_prev;
            he_next
        };

        let loop_he = self.loops[the_loop].half_edge;
        if loop_he == half_edge || loop_he == twin_handle {
            self.loops[the_loop].half_edge = replacement_loop_half_edge;
        }
        self.repair_outgoing_after_same_loop_bridge_dissolve(he_origin, half_edge, twin_handle);
        self.repair_outgoing_after_same_loop_bridge_dissolve(twin_origin, twin_handle, half_edge);

        self.half_edges.remove(half_edge);
        self.half_edges.remove(twin_handle);
        self.edges.remove(edge_handle);

        Ok(())
    }

    fn repair_outgoing_after_same_loop_bridge_dissolve(
        &mut self,
        vertex: crate::topology::VertexHandle,
        removed: HalfEdgeHandle,
        removed_twin: HalfEdgeHandle,
    ) {
        if self.vertices[vertex].outgoing != removed {
            return;
        }
        self.vertices[vertex].outgoing = self
            .half_edges
            .iter()
            .find(|(half_edge, data)| {
                *half_edge != removed && *half_edge != removed_twin && data.origin == vertex
            })
            .map(|(half_edge, _)| half_edge)
            .unwrap_or_default();
    }
}
