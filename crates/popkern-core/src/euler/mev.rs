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
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn mev_at_vertex_creates_first_edge() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

        let mev_r = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        assert_eq!(brep.vertices.len(), 2);
        assert_eq!(brep.edges.len(), 1);
        assert_eq!(brep.half_edges.len(), 2);

        // Half-edges are twins
        let he = &brep.half_edges[mev_r.half_edge];
        let twin = &brep.half_edges[mev_r.half_edge_twin];
        assert_eq!(he.twin, mev_r.half_edge_twin);
        assert_eq!(twin.twin, mev_r.half_edge);

        // Origins are correct
        assert_eq!(he.origin, mvfs_r.vertex);
        assert_eq!(twin.origin, mev_r.vertex);

        // Next/prev form a 2-cycle
        assert_eq!(he.next, mev_r.half_edge_twin);
        assert_eq!(twin.next, mev_r.half_edge);
        assert_eq!(he.prev, mev_r.half_edge_twin);
        assert_eq!(twin.prev, mev_r.half_edge);

        // Loop updated
        assert_eq!(brep.loops[mvfs_r.loop_].half_edge, mev_r.half_edge);

        // Vertex outgoing updated
        assert_eq!(brep.vertices[mvfs_r.vertex].outgoing, mev_r.half_edge);
        assert_eq!(brep.vertices[mev_r.vertex].outgoing, mev_r.half_edge_twin);
    }

    #[test]
    fn mev_extends_chain() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        // Extend from the new vertex: pass the twin half-edge (which originates at v1)
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();

        assert_eq!(brep.vertices.len(), 3);
        assert_eq!(brep.edges.len(), 2);
        assert_eq!(brep.half_edges.len(), 4);

        // The chain should be: v0 -> v1 -> v2 -> v1 -> v0
        // Following next from mev1.half_edge: he1 -> he2_out -> he2_in -> he1_twin -> he1
        let he1 = &brep.half_edges[mev1.half_edge];
        assert_eq!(he1.next, mev2.half_edge); // v0->v1 then v1->v2
    }

    #[test]
    fn mev_three_vertices_in_chain() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let _mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();

        assert_eq!(brep.vertices.len(), 4);
        assert_eq!(brep.edges.len(), 3);
        assert_eq!(brep.half_edges.len(), 6);
    }

    #[test]
    fn kev_removes_last_vertex() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_r = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        // kev on the half-edge whose twin's origin has valence 1
        // mev_r.half_edge goes v0->v1, so twin goes v1->v0. v1 has valence 1.
        brep.kev(mev_r.half_edge).unwrap();

        assert_eq!(brep.vertices.len(), 1);
        assert_eq!(brep.edges.len(), 0);
        assert_eq!(brep.half_edges.len(), 0);
    }

    #[test]
    fn kev_from_chain_end() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();

        // Remove v2 (end of chain). mev2.half_edge goes v1->v2, v2 has valence 1.
        brep.kev(mev2.half_edge).unwrap();

        assert_eq!(brep.vertices.len(), 2);
        assert_eq!(brep.edges.len(), 1);
        assert_eq!(brep.half_edges.len(), 2);
    }

    #[test]
    fn mev_at_vertex_invalid_vertex() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result =
            brep.mev_at_vertex(Default::default(), mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0));
        assert!(result.is_err());
    }

    #[test]
    fn mev_at_vertex_invalid_loop() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.mev_at_vertex(
            mvfs_r.vertex,
            Default::default(),
            Point3::new(1.0, 0.0, 0.0),
        );
        assert!(result.is_err());
    }

    #[test]
    fn mev_at_vertex_loop_already_has_half_edges() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Try again — loop already has half-edges
        let result = brep.mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(2.0, 0.0, 0.0));
        assert!(result.is_err());
    }

    #[test]
    fn mev_invalid_half_edge() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.mev(Default::default(), Point3::new(1.0, 0.0, 0.0));
        assert!(result.is_err());
    }

    #[test]
    fn kev_invalid_half_edge() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.kev(Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn kev_non_valence_1_vertex() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let _mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        // v1 has valence 2 now (connected to v0 and v2). kev should fail.
        let result = brep.kev(mev1.half_edge);
        assert!(result.is_err());
    }
}
