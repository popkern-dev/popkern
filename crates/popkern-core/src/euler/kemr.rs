//! kemr — Kill Edge, Make Ring
//! mekr — Make Edge, Kill Ring

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts};
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

#[derive(Debug, Clone, Copy)]
struct KemrApplied {
    result: KemrResult,
    removed_edge: EdgeHandle,
    removed_half_edge: HalfEdgeHandle,
    removed_half_edge_twin: HalfEdgeHandle,
    outer_loop: LoopHandle,
    face: crate::topology::FaceHandle,
}

#[derive(Debug, Clone, Copy)]
struct MekrApplied {
    result: MekrResult,
    removed_inner_loop: LoopHandle,
    merged_loop: LoopHandle,
    face: crate::topology::FaceHandle,
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
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let applied = candidate.apply_kemr(half_edge)?;
        candidate.verify_kemr_postconditions(before, applied)?;
        *self = candidate;
        Ok(applied.result)
    }

    fn apply_kemr(&mut self, half_edge: HalfEdgeHandle) -> Result<KemrApplied, CoreError> {
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

        Ok(KemrApplied {
            result: KemrResult { inner_loop },
            removed_edge: edge_handle,
            removed_half_edge: half_edge,
            removed_half_edge_twin: he_twin,
            outer_loop: old_loop,
            face,
        })
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
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let applied = candidate.apply_mekr(he_outer, he_inner)?;
        candidate.verify_mekr_postconditions(before, applied)?;
        *self = candidate;
        Ok(applied.result)
    }

    fn apply_mekr(
        &mut self,
        he_outer: HalfEdgeHandle,
        he_inner: HalfEdgeHandle,
    ) -> Result<MekrApplied, CoreError> {
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
        if !self.faces[outer_face].inner_loops.contains(&inner_loop) {
            return Err(CoreError::PreconditionFailed(
                "he_inner loop must be listed as an inner loop of the face".into(),
            ));
        }

        let v_outer = outer_data.origin;
        let v_inner = inner_data.origin;
        let he_outer_prev = outer_data.prev;
        let he_inner_prev = inner_data.prev;

        // Create the edge.
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
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

        Ok(MekrApplied {
            result: MekrResult {
                edge,
                half_edge: new_he,
                half_edge_twin: new_he_twin,
            },
            removed_inner_loop: inner_loop,
            merged_loop: outer_loop,
            face: outer_face,
        })
    }

    fn verify_kemr_postconditions(
        &self,
        before: EulerTopologyCounts,
        applied: KemrApplied,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "kemr",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::kemr(),
        )?;
        if self.edges.contains_key(applied.removed_edge)
            || self.half_edges.contains_key(applied.removed_half_edge)
            || self.half_edges.contains_key(applied.removed_half_edge_twin)
        {
            return Err(CoreError::PreconditionFailed(
                "kemr left removed topology live".into(),
            ));
        }
        let face = self
            .faces
            .get(applied.face)
            .ok_or_else(|| CoreError::PreconditionFailed("kemr face is not live".into()))?;
        if !face.inner_loops.contains(&applied.result.inner_loop) {
            return Err(CoreError::PreconditionFailed(
                "kemr face does not own created inner loop".into(),
            ));
        }
        if self
            .loops
            .get(applied.result.inner_loop)
            .is_none_or(|loop_data| loop_data.face != applied.face)
        {
            return Err(CoreError::PreconditionFailed(
                "kemr created inner loop does not point to owning face".into(),
            ));
        }
        self.verify_kemr_loop(applied.outer_loop)?;
        self.verify_kemr_loop(applied.result.inner_loop)?;
        self.verify_no_removed_kemr_refs(applied)?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }

    fn verify_mekr_postconditions(
        &self,
        before: EulerTopologyCounts,
        applied: MekrApplied,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "mekr",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mekr(),
        )?;
        if self.loops.contains_key(applied.removed_inner_loop) {
            return Err(CoreError::PreconditionFailed(
                "mekr left removed inner loop live".into(),
            ));
        }
        let edge = self
            .edges
            .get(applied.result.edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mekr edge is not live".into()))?;
        if edge.half_edge != applied.result.half_edge {
            return Err(CoreError::PreconditionFailed(
                "mekr edge does not reference created half-edge".into(),
            ));
        }
        let half_edge = self
            .half_edges
            .get(applied.result.half_edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mekr half-edge is not live".into()))?;
        let twin = self
            .half_edges
            .get(applied.result.half_edge_twin)
            .ok_or_else(|| CoreError::PreconditionFailed("mekr twin is not live".into()))?;
        if half_edge.twin != applied.result.half_edge_twin
            || twin.twin != applied.result.half_edge
            || half_edge.loop_ != applied.merged_loop
            || twin.loop_ != applied.merged_loop
        {
            return Err(CoreError::PreconditionFailed(
                "mekr created half-edge pair is not consistently wired".into(),
            ));
        }
        if self.faces[applied.face]
            .inner_loops
            .contains(&applied.removed_inner_loop)
        {
            return Err(CoreError::PreconditionFailed(
                "mekr face still lists removed inner loop".into(),
            ));
        }
        self.verify_kemr_loop(applied.merged_loop)?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }

    fn verify_no_removed_kemr_refs(&self, applied: KemrApplied) -> Result<(), CoreError> {
        for half_edge in self.half_edges.values() {
            if half_edge.edge == applied.removed_edge
                || half_edge.next == applied.removed_half_edge
                || half_edge.next == applied.removed_half_edge_twin
                || half_edge.prev == applied.removed_half_edge
                || half_edge.prev == applied.removed_half_edge_twin
            {
                return Err(CoreError::PreconditionFailed(
                    "kemr left live half-edge references to removed topology".into(),
                ));
            }
        }
        Ok(())
    }

    fn verify_kemr_loop(&self, loop_: LoopHandle) -> Result<(), CoreError> {
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
                return Ok(());
            }
        }
        Err(CoreError::PreconditionFailed(
            "loop did not close during verification".into(),
        ))
    }
}

#[cfg(test)]
#[path = "kemr_tests.rs"]
mod tests;
