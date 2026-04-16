//! kemr — Kill Edge, Make Ring
//! mekr — Make Edge, Kill Ring

use crate::error::CoreError;
use crate::topology::{BRep, Edge, EdgeHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle};

/// Result of the kemr (Kill Edge, Make Ring) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KemrResult {
    /// The newly created inner loop (hole).
    pub inner_loop: LoopHandle,
}

/// Result of the mekr (Make Edge, Kill Ring) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MekrResult {
    /// The newly created edge.
    pub edge: EdgeHandle,
    /// Half-edge from origin(he_outer) to origin(he_inner).
    pub half_edge: HalfEdgeHandle,
    /// Half-edge from origin(he_inner) to origin(he_outer) (twin).
    pub half_edge_twin: HalfEdgeHandle,
}

impl BRep {
    /// Kill Edge, Make Ring — removes an edge and splits a loop into outer + inner loop.
    ///
    /// Creates a hole (inner loop) in the face by removing the edge identified
    /// by `half_edge`. The removal must result in two disconnected cycles within
    /// the same face.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the half-edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the half-edge and its twin are in
    ///   different loops (use kef instead).
    pub fn kemr(&mut self, half_edge: HalfEdgeHandle) -> Result<KemrResult, CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let he_twin = he_data.twin;
        let edge_handle = he_data.edge;
        let he_prev = he_data.prev;
        let he_next = he_data.next;
        let he_loop = he_data.loop_;

        let twin_data = self
            .half_edges
            .get(he_twin)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge handle is invalid".into()))?;

        let twin_prev = twin_data.prev;
        let twin_next = twin_data.next;
        let twin_loop = twin_data.loop_;

        // Precondition: both half-edges must be in the same loop.
        // (If they're in different loops, this is a kef situation, not kemr.)
        if he_loop != twin_loop {
            return Err(CoreError::PreconditionFailed(
                "half_edge and twin must be in the same loop for kemr".into(),
            ));
        }

        let old_loop = he_loop;
        let face = self.loops[old_loop].face;

        // Unsplice the two half-edges from the loop.
        self.half_edges[he_prev].next = twin_next;
        self.half_edges[twin_next].prev = he_prev;
        self.half_edges[twin_prev].next = he_next;
        self.half_edges[he_next].prev = twin_prev;

        // Now there are two separate cycles:
        //   Cycle A: starting from twin_next (was after the twin)
        //   Cycle B: starting from he_next (was after he)
        // One of these keeps the old loop, the other becomes the new inner loop.

        // Determine which cycle is which. The old loop keeps whichever cycle
        // we assign to it. We'll assign the cycle starting from twin_next to
        // the old loop, and the cycle starting from he_next to the new inner loop.
        // (This is a convention choice; the key is consistency.)

        // Create the new inner loop.
        let inner_loop = self.loops.insert(Loop {
            half_edge: he_next,
            face,
        });

        // Walk cycle B (starting from he_next) and assign to inner_loop.
        let mut cursor = he_next;
        loop {
            self.half_edges[cursor].loop_ = inner_loop;
            cursor = self.half_edges[cursor].next;
            if cursor == he_next {
                break;
            }
        }

        // Update old loop's half_edge in case it pointed to one of the removed half-edges.
        self.loops[old_loop].half_edge = twin_next;

        // Add the inner loop to the face's inner_loops.
        self.faces[face].inner_loops.push(inner_loop);

        // Remove the edge and both half-edges.
        self.half_edges.remove(half_edge);
        self.half_edges.remove(he_twin);
        self.edges.remove(edge_handle);

        Ok(KemrResult { inner_loop })
    }

    /// Make Edge, Kill Ring — connects a vertex in the outer loop to a vertex
    /// in an inner loop, merging them.
    ///
    /// Inverse of kemr. Creates a new edge connecting `origin(he_outer)` to
    /// `origin(he_inner)`, merging the inner loop into the outer loop.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if either half-edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the half-edges are in the same loop,
    ///   or if `he_inner` is not in an inner loop of the same face.
    pub fn mekr(
        &mut self,
        he_outer: HalfEdgeHandle,
        he_inner: HalfEdgeHandle,
    ) -> Result<MekrResult, CoreError> {
        let outer_data = self
            .half_edges
            .get(he_outer)
            .ok_or_else(|| CoreError::InvalidHandle("he_outer handle is invalid".into()))?;
        let inner_data = self
            .half_edges
            .get(he_inner)
            .ok_or_else(|| CoreError::InvalidHandle("he_inner handle is invalid".into()))?;

        let outer_loop = outer_data.loop_;
        let inner_loop = inner_data.loop_;

        // Precondition: must be in different loops.
        if outer_loop == inner_loop {
            return Err(CoreError::PreconditionFailed(
                "he_outer and he_inner must be in different loops".into(),
            ));
        }

        // Precondition: both loops must be on the same face.
        let outer_face = self.loops[outer_loop].face;
        let inner_face = self.loops[inner_loop].face;
        if outer_face != inner_face {
            return Err(CoreError::PreconditionFailed(
                "he_outer and he_inner must be on the same face".into(),
            ));
        }

        let v_outer = outer_data.origin;
        let v_inner = inner_data.origin;
        let he_outer_prev = outer_data.prev;
        let he_inner_prev = inner_data.prev;

        // Create the edge.
        let edge = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
        });

        // Create new_he: v_outer -> v_inner (spliced into the merged loop)
        let new_he = self.half_edges.insert(HalfEdge {
            origin: v_outer,
            twin: Default::default(),
            next: he_inner,
            prev: he_outer_prev,
            edge,
            loop_: outer_loop,
            pcurve: None,
        });

        // Create new_he_twin: v_inner -> v_outer
        let new_he_twin = self.half_edges.insert(HalfEdge {
            origin: v_inner,
            twin: new_he,
            next: he_outer,
            prev: he_inner_prev,
            edge,
            loop_: outer_loop,
            pcurve: None,
        });

        // Patch twin pointer.
        self.half_edges[new_he].twin = new_he_twin;

        // Splice into the loops.
        self.half_edges[he_outer_prev].next = new_he;
        self.half_edges[he_inner_prev].next = new_he_twin;
        self.half_edges[he_outer].prev = new_he_twin;
        self.half_edges[he_inner].prev = new_he;

        // Walk the merged cycle and set all half-edges to the outer loop.
        let mut cursor = new_he;
        loop {
            self.half_edges[cursor].loop_ = outer_loop;
            cursor = self.half_edges[cursor].next;
            if cursor == new_he {
                break;
            }
        }

        // Update the outer loop's half_edge.
        self.loops[outer_loop].half_edge = new_he;

        // Update edge's half_edge.
        self.edges[edge].half_edge = new_he;

        // Remove the inner loop from the face's inner_loops.
        let face = &mut self.faces[outer_face];
        face.inner_loops.retain(|&l| l != inner_loop);

        // Remove the inner loop entity.
        self.loops.remove(inner_loop);

        Ok(MekrResult {
            edge,
            half_edge: new_he,
            half_edge_twin: new_he_twin,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Build a topology with a bridge edge suitable for kemr:
    ///
    /// v3 -- v0 -- v1 -- v2
    ///
    /// Three spur edges from v0: v0-v1, v1-v2 (chain), and v0-v3.
    /// The single loop traverses: v0->v1->v2->v1->v0->v3->v0
    /// The edge v0-v1 is a bridge. kemr on it splits into:
    ///   outer: v0->v3->v0
    ///   inner: v1->v2->v1
    #[test]
    fn kemr_creates_inner_loop() {
        let mut brep = BRep::new(Tolerance::default());

        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // v0->v1
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // v1->v2 (spur from v1)
        let _mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(2.0, 0.0, 0.0))
            .unwrap();
        // v0->v3 (spur from v0, spliced before mev1.half_edge)
        let _mev3 = brep
            .mev(mev1.half_edge, Point3::new(0.0, 1.0, 0.0))
            .unwrap();

        assert_eq!(brep.edges.len(), 3);
        assert_eq!(brep.half_edges.len(), 6);
        assert_eq!(brep.loops.len(), 1);

        // kemr on the bridge edge v0-v1
        let _kemr_r = brep.kemr(mev1.half_edge).unwrap();

        assert_eq!(brep.edges.len(), 2); // removed 1
        assert_eq!(brep.half_edges.len(), 4); // removed 2
        assert_eq!(brep.loops.len(), 2); // split into outer + inner

        // The face should have an inner loop.
        assert!(
            brep.faces.values().any(|f| !f.inner_loops.is_empty()),
            "Some face should have an inner loop"
        );

        // Verify both loops have valid half-edge cycles of length 2.
        for (_, l) in &brep.loops {
            let start = l.half_edge;
            let mut cursor = start;
            let mut count = 0;
            loop {
                count += 1;
                cursor = brep.half_edges[cursor].next;
                if cursor == start {
                    break;
                }
                assert!(count < 100, "infinite loop detected");
            }
            assert_eq!(count, 2, "each loop should have 2 half-edges");
        }
    }

    #[test]
    fn mekr_closes_inner_loop() {
        let mut brep = BRep::new(Tolerance::default());

        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let _mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(2.0, 0.0, 0.0))
            .unwrap();
        let _mev3 = brep
            .mev(mev1.half_edge, Point3::new(0.0, 1.0, 0.0))
            .unwrap();

        // kemr splits the loop
        let kemr_r = brep.kemr(mev1.half_edge).unwrap();

        assert_eq!(brep.edges.len(), 2);
        assert_eq!(brep.loops.len(), 2);

        // Now mekr: connect outer to inner, merging loops
        let inner_he = brep.loops[kemr_r.inner_loop].half_edge;
        let face_with_hole = brep
            .faces
            .values()
            .find(|f| !f.inner_loops.is_empty())
            .unwrap();
        let outer_he = brep.loops[face_with_hole.outer_loop].half_edge;

        let _mekr_r = brep.mekr(outer_he, inner_he).unwrap();

        assert_eq!(brep.edges.len(), 3); // restored
        assert_eq!(brep.loops.len(), 1); // merged back
        assert!(
            brep.faces.values().all(|f| f.inner_loops.is_empty()),
            "No face should have inner loops after mekr"
        );

        // Verify the single loop has 6 half-edges.
        let the_loop = brep.loops.values().next().unwrap();
        let start = the_loop.half_edge;
        let mut cursor = start;
        let mut count = 0;
        loop {
            count += 1;
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
            assert!(count < 100, "infinite loop detected");
        }
        assert_eq!(count, 6, "merged loop should have 6 half-edges");
    }

    #[test]
    fn kemr_precondition_different_loops() {
        let mut brep = BRep::new(Tolerance::default());

        // Build a square and close with mef — the mef edge has half-edges in different loops
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        // kemr on the mef edge should fail (different loops)
        let result = brep.kemr(mef1.half_edge);
        assert!(result.is_err());
    }

    #[test]
    fn kemr_invalid_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.kemr(Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn mekr_invalid_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.mekr(Default::default(), Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn mekr_same_loop_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Both half-edges are in the same loop
        let result = brep.mekr(mev1.half_edge, mev1.half_edge_twin);
        assert!(result.is_err());
    }

    #[test]
    fn mekr_different_face_fails() {
        let mut brep = BRep::new(Tolerance::default());
        // Build two separate solids — half-edges in loops on different faces
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        // mef created two loops on different faces
        // mekr between these two loops should fail (different faces)
        let he_new = mef.half_edge;
        let he_twin = mef.half_edge_twin;
        let result = brep.mekr(he_new, he_twin);
        assert!(result.is_err());
    }

    #[test]
    fn kemr_mekr_roundtrip() {
        let mut brep = BRep::new(Tolerance::default());

        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let _mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(2.0, 0.0, 0.0))
            .unwrap();
        let _mev3 = brep
            .mev(mev1.half_edge, Point3::new(0.0, 1.0, 0.0))
            .unwrap();

        let edges_before = brep.edges.len();
        let loops_before = brep.loops.len();

        // kemr then mekr should restore counts
        let kemr_r = brep.kemr(mev1.half_edge).unwrap();

        let inner_he = brep.loops[kemr_r.inner_loop].half_edge;
        let face_with_hole = brep
            .faces
            .values()
            .find(|f| !f.inner_loops.is_empty())
            .unwrap();
        let outer_he = brep.loops[face_with_hole.outer_loop].half_edge;

        brep.mekr(outer_he, inner_he).unwrap();

        assert_eq!(brep.edges.len(), edges_before);
        assert_eq!(brep.loops.len(), loops_before);
    }
}
