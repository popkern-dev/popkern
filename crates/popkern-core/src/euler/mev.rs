//! mev — Make Edge, Vertex
//! kev — Kill Edge, Vertex

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::topology::{
    BRep, Edge, EdgeHandle, HalfEdge, HalfEdgeHandle, LoopHandle, Vertex, VertexHandle,
};

/// Result of the mev / mev_at_vertex operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MevResult {
    /// The newly created vertex.
    pub vertex: VertexHandle,
    /// The newly created edge.
    pub edge: EdgeHandle,
    /// Half-edge from existing vertex to new vertex.
    pub half_edge: HalfEdgeHandle,
    /// Half-edge from new vertex to existing vertex (twin).
    pub half_edge_twin: HalfEdgeHandle,
}

impl BRep {
    /// Make Edge, Vertex at a bare vertex — creates the first edge in a degenerate loop.
    ///
    /// Used after `mvfs` when the loop has no half-edges yet.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if `vertex` or `loop_` handles are invalid.
    /// - [`CoreError::PreconditionFailed`] if the loop already has half-edges.
    pub fn mev_at_vertex(
        &mut self,
        vertex: VertexHandle,
        loop_: LoopHandle,
        position: Point3,
    ) -> Result<MevResult, CoreError> {
        // Validate handles.
        if !self.vertices.contains_key(vertex) {
            return Err(CoreError::InvalidHandle("vertex handle is invalid".into()));
        }
        let loop_data = self
            .loops
            .get(loop_)
            .ok_or_else(|| CoreError::InvalidHandle("loop handle is invalid".into()))?;

        // Precondition: loop has no half-edges.
        if self.half_edges.contains_key(loop_data.half_edge) {
            return Err(CoreError::PreconditionFailed(
                "loop already has half-edges".into(),
            ));
        }

        // Create the new vertex.
        let new_vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // Create the edge.
        let edge = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
        });

        // Create half-edges with placeholder handles; we'll patch after insertion.
        let he = self.half_edges.insert(HalfEdge {
            origin: vertex,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });

        let he_twin = self.half_edges.insert(HalfEdge {
            origin: new_vertex,
            twin: he,
            next: he,
            prev: he,
            edge,
            loop_,
            pcurve: None,
        });

        // Patch he to point to its twin.
        let he_data = &mut self.half_edges[he];
        he_data.twin = he_twin;
        he_data.next = he_twin;
        he_data.prev = he_twin;

        // Update loop to point to the new half-edge.
        self.loops[loop_].half_edge = he;

        // Update edge to point to he.
        self.edges[edge].half_edge = he;

        // Update vertex outgoing pointers.
        self.vertices[vertex].outgoing = he;
        self.vertices[new_vertex].outgoing = he_twin;

        Ok(MevResult {
            vertex: new_vertex,
            edge,
            half_edge: he,
            half_edge_twin: he_twin,
        })
    }

    /// Make Edge, Vertex — general case that splices a new spur edge into an existing loop.
    ///
    /// Creates a new vertex at `position` and an edge from `origin(half_edge)` to it.
    /// The new half-edge pair is spliced before `half_edge` in the loop traversal.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the half-edge handle is invalid.
    pub fn mev(
        &mut self,
        half_edge: HalfEdgeHandle,
        position: Point3,
    ) -> Result<MevResult, CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let existing_vertex = he_data.origin;
        let the_loop = he_data.loop_;
        let old_prev = he_data.prev;

        // Create the new vertex.
        let new_vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // Create the edge.
        let edge = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
        });

        // Create he_out: existing_vertex -> new_vertex
        let he_out = self.half_edges.insert(HalfEdge {
            origin: existing_vertex,
            twin: Default::default(),
            next: Default::default(),
            prev: old_prev,
            edge,
            loop_: the_loop,
            pcurve: None,
        });

        // Create he_in: new_vertex -> existing_vertex
        let he_in = self.half_edges.insert(HalfEdge {
            origin: new_vertex,
            twin: he_out,
            next: half_edge,
            prev: he_out,
            edge,
            loop_: the_loop,
            pcurve: None,
        });

        // Patch he_out.
        let he_out_data = &mut self.half_edges[he_out];
        he_out_data.twin = he_in;
        he_out_data.next = he_in;

        // Splice into loop: old_prev.next = he_out
        self.half_edges[old_prev].next = he_out;

        // half_edge.prev = he_in
        self.half_edges[half_edge].prev = he_in;

        // Update edge.
        self.edges[edge].half_edge = he_out;

        // Update new vertex outgoing.
        self.vertices[new_vertex].outgoing = he_in;

        Ok(MevResult {
            vertex: new_vertex,
            edge,
            half_edge: he_out,
            half_edge_twin: he_in,
        })
    }

    /// Kill Edge, Vertex — inverse of mev.
    ///
    /// Removes the edge identified by `half_edge` and the destination vertex
    /// (the vertex at `twin.origin`), which must have valence 1.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the half-edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the destination vertex has valence > 1.
    pub fn kev(&mut self, half_edge: HalfEdgeHandle) -> Result<(), CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let he_twin_handle = he_data.twin;
        let edge_handle = he_data.edge;
        let he_prev = he_data.prev;
        let the_loop = he_data.loop_;

        let twin_data = self
            .half_edges
            .get(he_twin_handle)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge handle is invalid".into()))?;

        let dest_vertex = twin_data.origin;
        let twin_next = twin_data.next;

        // Validate that dest_vertex has valence 1: the only outgoing half-edge
        // from dest_vertex should be he_twin, and it should loop back (next == half_edge
        // means there's only the twin going out and coming back, so valence 1).
        // More precisely: count half-edges originating from dest_vertex.
        let valence = self
            .half_edges
            .values()
            .filter(|he| he.origin == dest_vertex)
            .count();

        if valence != 1 {
            return Err(CoreError::PreconditionFailed(format!(
                "destination vertex has valence {valence}, expected 1"
            )));
        }

        // Check if we're removing the last two half-edges in the loop (back to degenerate).
        let is_last_edge = he_prev == he_twin_handle && twin_next == half_edge;

        if is_last_edge {
            // The loop becomes empty again.
            self.loops[the_loop].half_edge = Default::default();
            // Reset the surviving vertex's outgoing.
            let surviving_vertex = he_data.origin;
            self.vertices[surviving_vertex].outgoing = Default::default();
        } else {
            // Unsplice: he_prev.next = twin_next, twin_next.prev = he_prev
            self.half_edges[he_prev].next = twin_next;
            self.half_edges[twin_next].prev = he_prev;

            // If the loop's half_edge was one of the removed ones, update it.
            let loop_he = self.loops[the_loop].half_edge;
            if loop_he == half_edge || loop_he == he_twin_handle {
                self.loops[the_loop].half_edge = twin_next;
            }
        }

        // Remove entities.
        self.half_edges.remove(half_edge);
        self.half_edges.remove(he_twin_handle);
        self.edges.remove(edge_handle);
        self.vertices.remove(dest_vertex);

        Ok(())
    }

    /// Dissolve a same-loop bridge edge — removes an edge whose
    /// two half-edges lie on the SAME loop, without merging faces
    /// (kef requires cross-face) and without removing a vertex
    /// (kev requires destination valence-1).
    ///
    /// Both endpoint vertices remain in the BRep with their valence
    /// reduced by 1.  The face's loop traversal is patched to skip
    /// the dissolved edge's two half-edges; the loop remains a
    /// single closed cycle (one fewer pair of half-edges).
    ///
    /// Used by Phase 6c-A2 §28.H Session 3.5 to clean up V's
    /// same-loop bridges that arise after `kef_cross_face_v_edges`
    /// merges per-plane V-corner triangles into their parent V-side
    /// sub-faces.  The first spur kef per plane creates a same-
    /// loop bridge for the second spur; the cube-edge stub
    /// becomes same-loop after the V-incident plane sub-faces
    /// merge.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if `half_edge` is invalid.
    /// - [`CoreError::PreconditionFailed`] if the edge's two
    ///   half-edges are on DIFFERENT loops (use `kef` instead).
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

        // Precondition: same loop (otherwise this is a kef job).
        if he_data.loop_ != twin_data.loop_ {
            return Err(CoreError::PreconditionFailed(
                "edge's half-edges are on different loops; use kef instead".into(),
            ));
        }

        // Precondition: bridge isn't the only edge in the loop —
        // otherwise dissolving would leave an empty loop.  Detect:
        // the bridge is degenerate when its half-edge's prev or
        // next is the twin (and vice versa), i.e., the loop
        // traversal goes he → twin → he with nothing between.
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

        // Splice out both half-edges from the loop.  For
        // ADJACENT slits (where he.next == twin or vice versa),
        // the splice has to skip both at once — splicing them
        // independently would leave one of the prev/next
        // pointers referring to a half-edge that's about to be
        // removed.  Handle the two adjacency cases explicitly.
        if he_next == twin_handle {
            // Loop: ... → he_prev → he → twin → twin_next → ...
            // Skip both: he_prev → twin_next.
            self.half_edges[he_prev].next = twin_next;
            self.half_edges[twin_next].prev = he_prev;
        } else if twin_next == half_edge {
            // Loop: ... → twin_prev → twin → he → he_next → ...
            self.half_edges[twin_prev].next = he_next;
            self.half_edges[he_next].prev = twin_prev;
        } else {
            // Non-adjacent: splice independently.
            self.half_edges[he_prev].next = he_next;
            self.half_edges[he_next].prev = he_prev;
            self.half_edges[twin_prev].next = twin_next;
            self.half_edges[twin_next].prev = twin_prev;
        }

        // If the loop's anchor half-edge was one of the dissolved
        // ones, retarget to a still-valid half-edge in the loop.
        let loop_he = self.loops[the_loop].half_edge;
        if loop_he == half_edge || loop_he == twin_handle {
            self.loops[the_loop].half_edge = he_next;
        }

        // Repair endpoint vertices' outgoing pointers if they
        // happened to reference the dissolved half-edges.
        if self.vertices[he_origin].outgoing == half_edge {
            // Find another outgoing half-edge for he_origin.
            let mut new_outgoing: Option<HalfEdgeHandle> = None;
            for (h, h_data) in self.half_edges.iter() {
                if h == half_edge || h == twin_handle {
                    continue;
                }
                if h_data.origin == he_origin {
                    new_outgoing = Some(h);
                    break;
                }
            }
            if let Some(new_out) = new_outgoing {
                self.vertices[he_origin].outgoing = new_out;
            } else {
                self.vertices[he_origin].outgoing = Default::default();
            }
        }
        if self.vertices[twin_origin].outgoing == twin_handle {
            let mut new_outgoing: Option<HalfEdgeHandle> = None;
            for (h, h_data) in self.half_edges.iter() {
                if h == half_edge || h == twin_handle {
                    continue;
                }
                if h_data.origin == twin_origin {
                    new_outgoing = Some(h);
                    break;
                }
            }
            if let Some(new_out) = new_outgoing {
                self.vertices[twin_origin].outgoing = new_out;
            } else {
                self.vertices[twin_origin].outgoing = Default::default();
            }
        }

        // Remove the dissolved entities.  Vertices stay (their
        // valence is just reduced by 1).
        self.half_edges.remove(half_edge);
        self.half_edges.remove(twin_handle);
        self.edges.remove(edge_handle);

        Ok(())
    }
}

#[cfg(test)]
#[path = "mev_tests.rs"]
mod tests;
