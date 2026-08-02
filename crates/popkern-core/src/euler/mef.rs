//! mef — Make Edge, Face
//! kef — Kill Edge, Face

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts};
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KefResult {
    pub removed_edge: EdgeHandle,
    pub removed_face: FaceHandle,
    pub removed_loop: LoopHandle,
    pub removed_half_edge: HalfEdgeHandle,
    pub removed_half_edge_twin: HalfEdgeHandle,
    pub surviving_loop: LoopHandle,
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
    #[track_caller]
    pub fn mef(
        &mut self,
        he1: HalfEdgeHandle,
        he2: HalfEdgeHandle,
    ) -> Result<MefResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let old_loop = self
            .half_edges
            .get(he1)
            .ok_or_else(|| CoreError::InvalidHandle("he1 handle is invalid".into()))?
            .loop_;
        let mut candidate = self.clone();
        let result = candidate.apply_mef(he1, he2)?;
        candidate.verify_mef_postconditions(before, result, old_loop)?;
        *self = candidate;
        Ok(result)
    }

    #[track_caller]
    fn apply_mef(
        &mut self,
        he1: HalfEdgeHandle,
        he2: HalfEdgeHandle,
    ) -> Result<MefResult, CoreError> {
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
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        // Create new face (inherits parent surface).
        let new_face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: self.faces[old_face].surface.clone(),
            name: None,
            tolerance: None,
            role: Default::default(),
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
    pub fn kef(&mut self, half_edge: HalfEdgeHandle) -> Result<KefResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_kef(half_edge)?;
        candidate.verify_kef_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_kef(&mut self, half_edge: HalfEdgeHandle) -> Result<KefResult, CoreError> {
        let he_data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?;

        let he_twin_handle = he_data.twin;
        let edge_handle = he_data.edge;
        let he_prev = he_data.prev;
        let he_loop = he_data.loop_;
        let he_origin = he_data.origin;

        let twin_data = self
            .half_edges
            .get(he_twin_handle)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge handle is invalid".into()))?;

        let twin_prev = twin_data.prev;
        let twin_next = twin_data.next;
        let twin_loop = twin_data.loop_;
        let twin_origin = twin_data.origin;
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

        // Both endpoint vertices survive kef; re-point any `outgoing`
        // that referenced a removed half-edge. `twin_next` originates
        // at `half_edge`'s origin and `he_next` at the twin's origin,
        // and both survive the unsplice (the different-loops
        // precondition rules out either being one of the removed
        // pair). Leaving these stale silently disables every
        // downstream vertex-star walk that starts at `outgoing`
        // (PK-232: tier-2 anchor collection returned false before
        // inspecting a single face).
        if let Some(v) = self.vertices.get_mut(he_origin)
            && (v.outgoing == half_edge || v.outgoing == he_twin_handle)
        {
            v.outgoing = he_next_via_twin;
        }
        if let Some(v) = self.vertices.get_mut(twin_origin)
            && (v.outgoing == half_edge || v.outgoing == he_twin_handle)
        {
            v.outgoing = he_next;
        }

        // Remove entities.
        self.half_edges.remove(half_edge);
        self.half_edges.remove(he_twin_handle);
        self.edges.remove(edge_handle);
        self.loops.remove(killed_loop);
        self.faces.remove(killed_face);

        Ok(KefResult {
            removed_edge: edge_handle,
            removed_face: killed_face,
            removed_loop: killed_loop,
            removed_half_edge: half_edge,
            removed_half_edge_twin: he_twin_handle,
            surviving_loop,
        })
    }

    fn verify_mef_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: MefResult,
        old_loop: LoopHandle,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "mef",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mef(),
        )?;
        let edge = self
            .edges
            .get(result.edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mef edge is not live".into()))?;
        if edge.half_edge != result.half_edge {
            return Err(CoreError::PreconditionFailed(
                "mef edge does not reference the created half-edge".into(),
            ));
        }
        let half_edge = self
            .half_edges
            .get(result.half_edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mef half-edge is not live".into()))?;
        let twin = self.half_edges.get(result.half_edge_twin).ok_or_else(|| {
            CoreError::PreconditionFailed("mef twin half-edge is not live".into())
        })?;
        if half_edge.twin != result.half_edge_twin
            || twin.twin != result.half_edge
            || half_edge.edge != result.edge
            || twin.edge != result.edge
            || half_edge.loop_ != old_loop
            || twin.loop_ != result.loop_
        {
            return Err(CoreError::PreconditionFailed(
                "mef created half-edge pair is not consistently wired".into(),
            ));
        }
        let new_face = self
            .faces
            .get(result.face)
            .ok_or_else(|| CoreError::PreconditionFailed("mef face is not live".into()))?;
        if new_face.outer_loop != result.loop_ {
            return Err(CoreError::PreconditionFailed(
                "mef face does not own the created loop".into(),
            ));
        }
        if self
            .loops
            .get(result.loop_)
            .is_none_or(|loop_data| loop_data.face != result.face)
        {
            return Err(CoreError::PreconditionFailed(
                "mef loop does not point to the created face".into(),
            ));
        }
        self.verify_mef_loop_contains(old_loop, &[result.half_edge])?;
        self.verify_mef_loop_contains(result.loop_, &[result.half_edge_twin])?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }

    fn verify_kef_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: KefResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "kef",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::kef(),
        )?;
        if self.edges.contains_key(result.removed_edge)
            || self.faces.contains_key(result.removed_face)
            || self.loops.contains_key(result.removed_loop)
            || self.half_edges.contains_key(result.removed_half_edge)
            || self.half_edges.contains_key(result.removed_half_edge_twin)
        {
            return Err(CoreError::PreconditionFailed(
                "kef left removed topology live".into(),
            ));
        }
        for face in self.faces.values() {
            if face.outer_loop == result.removed_loop
                || face.inner_loops.contains(&result.removed_loop)
            {
                return Err(CoreError::PreconditionFailed(
                    "kef left a live face referencing removed loop".into(),
                ));
            }
        }
        for half_edge in self.half_edges.values() {
            if half_edge.edge == result.removed_edge
                || half_edge.loop_ == result.removed_loop
                || half_edge.next == result.removed_half_edge
                || half_edge.next == result.removed_half_edge_twin
                || half_edge.prev == result.removed_half_edge
                || half_edge.prev == result.removed_half_edge_twin
            {
                return Err(CoreError::PreconditionFailed(
                    "kef left live half-edge references to removed topology".into(),
                ));
            }
        }
        let loop_anchor = self
            .loops
            .get(result.surviving_loop)
            .ok_or_else(|| CoreError::PreconditionFailed("kef surviving loop is not live".into()))?
            .half_edge;
        self.verify_mef_loop_contains(result.surviving_loop, &[loop_anchor])?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }

    fn verify_mef_loop_contains(
        &self,
        loop_: LoopHandle,
        required: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        let start = self
            .loops
            .get(loop_)
            .ok_or_else(|| CoreError::PreconditionFailed("loop is not live".into()))?
            .half_edge;
        let mut cursor = start;
        let mut visited = std::collections::BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed("loop references a missing half-edge".into())
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "loop walk left the expected loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "loop repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                if required.iter().all(|half_edge| visited.contains(half_edge)) {
                    return Ok(());
                }
                return Err(CoreError::PreconditionFailed(
                    "loop is missing required half-edge".into(),
                ));
            }
        }
        Err(CoreError::PreconditionFailed(
            "loop did not close during verification".into(),
        ))
    }
}

#[cfg(test)]
#[path = "mef_tests.rs"]
mod tests;
