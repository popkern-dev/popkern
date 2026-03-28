//! mef — Make Edge, Face
//! kef — Kill Edge, Face

use crate::error::CoreError;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
};

/// Result of the mef (Make Edge, Face) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MefResult {
    /// The newly created edge.
    pub edge: EdgeHandle,
    /// The newly created face.
    pub face: FaceHandle,
    /// The newly created loop (on the new face).
    pub loop_: LoopHandle,
    /// Half-edge from origin(he1) to origin(he2).
    pub half_edge: HalfEdgeHandle,
    /// Half-edge from origin(he2) to origin(he1) (twin).
    pub half_edge_twin: HalfEdgeHandle,
}

impl BRep {
    /// Make Edge, Face — splits a face by connecting two vertices in the same loop.
    ///
    /// Creates a new edge connecting `origin(he1)` to `origin(he2)`, splitting
    /// the face into two. The new face inherits the parent face's surface
    /// geometry. The old face keeps the half-edges from `new_he` onward
    /// (through `he2`), and the new face gets the half-edges from `new_he_twin`
    /// onward (through `he1`).
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if `he1` or `he2` handles are invalid.
    /// - [`CoreError::PreconditionFailed`] if `he1` and `he2` are the same,
    ///   or not in the same loop.
    pub fn mef(
        &mut self,
        he1: HalfEdgeHandle,
        he2: HalfEdgeHandle,
    ) -> Result<MefResult, CoreError> {
        // Validate handles.
        let he1_data = self
            .half_edges
            .get(he1)
            .ok_or_else(|| CoreError::InvalidHandle("he1 handle is invalid".into()))?;
        let he2_data = self
            .half_edges
            .get(he2)
            .ok_or_else(|| CoreError::InvalidHandle("he2 handle is invalid".into()))?;

        // Precondition: he1 and he2 must be distinct.
        if he1 == he2 {
            return Err(CoreError::PreconditionFailed(
                "he1 and he2 must be distinct".into(),
            ));
        }

        // Precondition: he1 and he2 must be in the same loop.
        if he1_data.loop_ != he2_data.loop_ {
            return Err(CoreError::PreconditionFailed(
                "he1 and he2 must be in the same loop".into(),
            ));
        }

        let old_loop = he1_data.loop_;
        let old_face = self.loops[old_loop].face;
        let shell = self.faces[old_face].shell;

        let v1 = he1_data.origin;
        let v2 = he2_data.origin;
        let he1_prev = he1_data.prev;
        let he2_prev = he2_data.prev;

        // Create the edge.
        let edge = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
        });

        // Create new face (inherits parent surface).
        let new_face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: self.faces[old_face].surface.clone(),
        });

        // Create new loop for the new face.
        let new_loop = self.loops.insert(Loop {
            half_edge: Default::default(),
            face: new_face,
        });

        // Create new_he: v1 -> v2 (stays with old loop)
        let new_he = self.half_edges.insert(HalfEdge {
            origin: v1,
            twin: Default::default(),
            next: he2,
            prev: he1_prev,
            edge,
            loop_: old_loop,
            pcurve: None,
        });

        // Create new_he_twin: v2 -> v1 (goes to new loop)
        let new_he_twin = self.half_edges.insert(HalfEdge {
            origin: v2,
            twin: new_he,
            next: he1,
            prev: he2_prev,
            edge,
            loop_: new_loop,
            pcurve: None,
        });

        // Patch new_he's twin.
        self.half_edges[new_he].twin = new_he_twin;

        // Wire prev/next pointers.
        self.half_edges[he1_prev].next = new_he;
        self.half_edges[he2_prev].next = new_he_twin;
        self.half_edges[he1].prev = new_he_twin;
        self.half_edges[he2].prev = new_he;

        // Walk the cycle starting from new_he_twin, set all half-edges' loop to new_loop.
        let mut cursor = new_he_twin;
        loop {
            self.half_edges[cursor].loop_ = new_loop;
            cursor = self.half_edges[cursor].next;
            if cursor == new_he_twin {
                break;
            }
        }

        // Update old loop's half_edge to new_he.
        self.loops[old_loop].half_edge = new_he;

        // Update new loop's half_edge to new_he_twin.
        self.loops[new_loop].half_edge = new_he_twin;

        // Update new face's outer_loop.
        self.faces[new_face].outer_loop = new_loop;

        // Update edge's half_edge.
        self.edges[edge].half_edge = new_he;

        Ok(MefResult {
            edge,
            face: new_face,
            loop_: new_loop,
            half_edge: new_he,
            half_edge_twin: new_he_twin,
        })
    }

    /// Kill Edge, Face — inverse of mef.
    ///
    /// Removes the edge identified by `half_edge` and merges the two faces
    /// separated by it into one.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the half-edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the half-edge and its twin are
    ///   in the same loop (not separating two faces).
    pub fn kef(&mut self, half_edge: HalfEdgeHandle) -> Result<(), CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let he_twin_handle = he_data.twin;
        let edge_handle = he_data.edge;
        let he_prev = he_data.prev;
        let he_loop = he_data.loop_;

        let twin_data = self
            .half_edges
            .get(he_twin_handle)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge handle is invalid".into()))?;

        let twin_prev = twin_data.prev;
        let twin_next = twin_data.next;
        let twin_loop = twin_data.loop_;
        let he_next_via_twin = twin_next;

        // Precondition: they must be in different loops (different faces).
        if he_loop == twin_loop {
            return Err(CoreError::PreconditionFailed(
                "half_edge and twin must be in different loops/faces".into(),
            ));
        }

        // Determine which loop/face survives and which is killed.
        // The twin's loop/face is killed; the half_edge's loop/face survives.
        let surviving_loop = he_loop;
        let killed_loop = twin_loop;
        let killed_face = self.loops[killed_loop].face;

        // Unsplice: reconnect around the removed half-edges.
        // he_prev -> twin_next (skip over half_edge)
        // twin_prev -> he's next (skip over twin) — but he's next is he2 from mef context
        let he_next = self.half_edges[half_edge].next;

        self.half_edges[he_prev].next = he_next_via_twin;
        self.half_edges[he_next_via_twin].prev = he_prev;
        self.half_edges[twin_prev].next = he_next;
        self.half_edges[he_next].prev = twin_prev;

        // Walk the surviving loop and set all half-edges to surviving_loop.
        let start = he_next_via_twin;
        let mut cursor = start;
        loop {
            self.half_edges[cursor].loop_ = surviving_loop;
            cursor = self.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }

        // Update surviving loop's half_edge (in case it pointed to a removed one).
        self.loops[surviving_loop].half_edge = he_next_via_twin;

        // Remove entities.
        self.half_edges.remove(half_edge);
        self.half_edges.remove(he_twin_handle);
        self.edges.remove(edge_handle);
        self.loops.remove(killed_loop);
        self.faces.remove(killed_face);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn mef_splits_face_into_two() {
        let mut brep = BRep::new(Tolerance::default());
        // Build 4 vertices: v0(0,0) -> v1(1,0) -> v2(1,1) -> v3(0,1)
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

        // Before mef: 1 face, 1 loop
        assert_eq!(brep.faces.len(), 1);
        assert_eq!(brep.loops.len(), 1);

        // mef connects v3 to v0:
        // mev3.half_edge_twin originates at v3
        // mev1.half_edge originates at v0
        let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        assert_eq!(brep.faces.len(), 2);
        assert_eq!(brep.loops.len(), 2);
        assert_eq!(brep.edges.len(), 4); // 3 from mev + 1 from mef
        assert_eq!(brep.half_edges.len(), 8); // 6 from mev + 2 from mef
    }

    #[test]
    fn kef_merges_two_faces() {
        let mut brep = BRep::new(Tolerance::default());
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
        let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        assert_eq!(brep.faces.len(), 2);

        // kef on the mef edge merges back
        brep.kef(mef_r.half_edge).unwrap();

        assert_eq!(brep.faces.len(), 1);
        assert_eq!(brep.loops.len(), 1);
        assert_eq!(brep.edges.len(), 3);
        assert_eq!(brep.half_edges.len(), 6);
    }

    #[test]
    fn mef_precondition_same_loop() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        // Create a second, separate solid
        let mvfs2 = brep.mvfs(Point3::new(5.0, 0.0, 0.0)).unwrap();
        let mev2 = brep
            .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(6.0, 0.0, 0.0))
            .unwrap();

        // Try to mef between half-edges in different loops — should fail
        let result = brep.mef(mev1.half_edge, mev2.half_edge);
        assert!(result.is_err());
    }

    #[test]
    fn mef_same_half_edge_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // he1 == he2 should fail
        let result = brep.mef(mev.half_edge, mev.half_edge);
        assert!(result.is_err());
    }

    #[test]
    fn mef_invalid_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let result = brep.mef(Default::default(), mev.half_edge);
        assert!(result.is_err());
        let result = brep.mef(mev.half_edge, Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn kef_same_loop_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // This edge has both half-edges in the same loop, so kef should fail
        let result = brep.kef(mev.half_edge);
        assert!(result.is_err());
    }

    #[test]
    fn kef_invalid_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.kef(Default::default());
        assert!(result.is_err());
    }
}
