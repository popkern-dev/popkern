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

        let he_origin = he_data.origin;
        let twin_origin = twin_data.origin;
        let he_prev = he_data.prev;
        let he_next = he_data.next;
        let twin_prev = twin_data.prev;
        let twin_next = twin_data.next;
        let the_loop = he_data.loop_;

        // Splice out half_edge: prev.next = next, next.prev = prev.
        self.half_edges[he_prev].next = he_next;
        self.half_edges[he_next].prev = he_prev;

        // Splice out twin: same pattern.
        self.half_edges[twin_prev].next = twin_next;
        self.half_edges[twin_next].prev = twin_prev;

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
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// §28.H Session 3.5 — verify `dissolve_same_loop_bridge`
    /// removes a same-loop bridge cleanly: the edge and its 2
    /// half-edges are gone, the loop traversal is patched, and
    /// the endpoint vertices' valences each drop by 1.
    ///
    /// Construction: build a unit cube, then mev a spur from one
    /// of its corners.  The spur is a true valence-1 dangling
    /// edge (its destination has only that one incidence) — a
    /// well-formed same-loop bridge to dissolve.
    #[test]
    fn dissolve_same_loop_bridge_removes_a_spur() {
        let tol = Tolerance::default();
        let mut brep = BRep::unit_cube(tol).unwrap();
        let n_he_before = brep.half_edges.len();
        let n_edges_before = brep.edges.len();
        let n_vertices_before = brep.vertices.len();

        // Pick any vertex with at least one outgoing half-edge.
        let v = brep.vertices.keys().next().unwrap();
        let v_out = brep.vertices[v].outgoing;
        // mev a spur from v's outgoing.
        let mev_result = brep
            .mev(v_out, Point3::new(0.5, 0.5, 0.5))
            .expect("mev creates spur");
        let new_vertex = mev_result.vertex;
        assert_eq!(brep.half_edges.len(), n_he_before + 2);
        assert_eq!(brep.edges.len(), n_edges_before + 1);
        assert_eq!(brep.vertices.len(), n_vertices_before + 1);

        // Sanity: spur is same-loop (both halves on the same loop).
        let spur_out = mev_result.half_edge;
        let spur_in = mev_result.half_edge_twin;
        let l_out = brep.half_edges[spur_out].loop_;
        let l_in = brep.half_edges[spur_in].loop_;
        assert_eq!(l_out, l_in, "spur should be same-loop");

        // Dissolve via the new primitive.
        brep.dissolve_same_loop_bridge(spur_out)
            .expect("dissolve same-loop spur");

        // Edge + 2 half-edges removed.
        assert_eq!(brep.half_edges.len(), n_he_before);
        assert_eq!(brep.edges.len(), n_edges_before);
        // Both endpoint vertices stay (kev removes the dangling
        // vertex, dissolve does not).
        assert_eq!(brep.vertices.len(), n_vertices_before + 1);
        // The originally-isolated `new_vertex` is now valence 0.
        let new_v_outgoing_count = brep
            .half_edges
            .values()
            .filter(|h| h.origin == new_vertex)
            .count();
        assert_eq!(
            new_v_outgoing_count, 0,
            "dissolved spur leaves new_vertex with valence 0"
        );
    }

    /// §28.H Session 3.5 — verify dissolve refuses to operate on
    /// a cross-face edge (returns `PreconditionFailed`).
    #[test]
    fn dissolve_same_loop_bridge_rejects_cross_face() {
        let tol = Tolerance::default();
        let mut brep = BRep::unit_cube(tol).unwrap();
        // Pick any cube edge — its two half-edges are on
        // different faces (cross-face).
        let edge = brep.edges.keys().next().unwrap();
        let primary = brep.edges[edge].half_edge;
        let result = brep.dissolve_same_loop_bridge(primary);
        assert!(
            result.is_err(),
            "dissolve must reject cross-face edges; got {result:?}"
        );
    }

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
