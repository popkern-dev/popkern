//! mev — Make Edge, Vertex
//! kev — Kill Edge, Vertex

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts, verify_loop_contains_half_edges};
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KevResult {
    pub vertex: VertexHandle,
    pub edge: EdgeHandle,
    pub half_edge: HalfEdgeHandle,
    pub half_edge_twin: HalfEdgeHandle,
    pub loop_: LoopHandle,
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
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_mev_at_vertex(vertex, loop_, position)?;
        candidate.verify_mev_postconditions("mev_at_vertex", before, result, loop_)?;
        *self = candidate;
        Ok(result)
    }

    pub(crate) fn apply_mev_at_vertex(
        &mut self,
        vertex: VertexHandle,
        loop_: LoopHandle,
        position: Point3,
    ) -> Result<MevResult, CoreError> {
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
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
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
    #[track_caller]
    pub fn mev(
        &mut self,
        half_edge: HalfEdgeHandle,
        position: Point3,
    ) -> Result<MevResult, CoreError> {
        let target_loop = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge handle is invalid".into()))?
            .loop_;
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_mev(half_edge, position)?;
        candidate.verify_mev_postconditions("mev", before, result, target_loop)?;
        *self = candidate;
        Ok(result)
    }

    /// [`BRep::mev`] unverified, for compound operators (extrude).
    #[track_caller]
    pub(crate) fn apply_mev(
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
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
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
    pub fn kev(&mut self, half_edge: HalfEdgeHandle) -> Result<KevResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_kev(half_edge)?;
        candidate.verify_kev_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_kev(&mut self, half_edge: HalfEdgeHandle) -> Result<KevResult, CoreError> {
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
        let is_last_edge = he_prev == he_twin_handle && twin_next == half_edge;
        if !is_last_edge
            && (he_data.next != he_twin_handle
                || twin_data.prev != half_edge
                || he_prev == he_twin_handle
                || twin_next == half_edge)
        {
            return Err(CoreError::PreconditionFailed(
                "kev requires an MEV spur: half_edge must be immediately followed by its twin"
                    .into(),
            ));
        }

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

        // Remove entities. The surviving source vertex's outgoing may
        // have pointed at the removed spur half-edge — re-anchor it to
        // a live half-edge (the unspliced `twin_next` qualifies).
        self.remove_half_edge_reanchored(half_edge);
        self.remove_half_edge_reanchored(he_twin_handle);
        self.edges.remove(edge_handle);
        self.vertices.remove(dest_vertex);

        Ok(KevResult {
            vertex: dest_vertex,
            edge: edge_handle,
            half_edge,
            half_edge_twin: he_twin_handle,
            loop_: the_loop,
        })
    }

    fn verify_mev_postconditions(
        &self,
        op: &'static str,
        before: EulerTopologyCounts,
        result: MevResult,
        expected_loop: LoopHandle,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            op,
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mev(),
        )?;
        let edge = self
            .edges
            .get(result.edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mev edge is not live".into()))?;
        if edge.half_edge != result.half_edge {
            return Err(CoreError::PreconditionFailed(
                "mev edge does not reference the created half-edge".into(),
            ));
        }
        let half_edge = self
            .half_edges
            .get(result.half_edge)
            .ok_or_else(|| CoreError::PreconditionFailed("mev half-edge is not live".into()))?;
        let twin = self.half_edges.get(result.half_edge_twin).ok_or_else(|| {
            CoreError::PreconditionFailed("mev twin half-edge is not live".into())
        })?;
        if half_edge.twin != result.half_edge_twin
            || twin.twin != result.half_edge
            || half_edge.edge != result.edge
            || twin.edge != result.edge
            || half_edge.loop_ != expected_loop
            || twin.loop_ != expected_loop
        {
            return Err(CoreError::PreconditionFailed(
                "mev created half-edge pair is not consistently wired".into(),
            ));
        }
        let vertex = self.vertices.get(result.vertex).ok_or_else(|| {
            CoreError::PreconditionFailed("mev created vertex is not live".into())
        })?;
        if vertex.outgoing != result.half_edge_twin || twin.origin != result.vertex {
            return Err(CoreError::PreconditionFailed(
                "mev created vertex does not own the created twin half-edge".into(),
            ));
        }
        verify_loop_contains_half_edges(
            self,
            expected_loop,
            &[result.half_edge, result.half_edge_twin],
        )?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }

    fn verify_kev_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: KevResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "kev",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::kev(),
        )?;
        if self.vertices.contains_key(result.vertex)
            || self.edges.contains_key(result.edge)
            || self.half_edges.contains_key(result.half_edge)
            || self.half_edges.contains_key(result.half_edge_twin)
        {
            return Err(CoreError::PreconditionFailed(
                "kev left removed topology live".into(),
            ));
        }
        if self.boundary_edges.contains(&result.edge) {
            return Err(CoreError::PreconditionFailed(
                "kev left removed edge in boundary set".into(),
            ));
        }
        for edge in self.edges.values() {
            if edge.half_edge == result.half_edge || edge.half_edge == result.half_edge_twin {
                return Err(CoreError::PreconditionFailed(
                    "kev left a live edge referencing removed half-edges".into(),
                ));
            }
        }
        for half_edge in self.half_edges.values() {
            let stale_reference = if half_edge.edge == result.edge {
                Some("edge")
            } else if half_edge.twin == result.half_edge || half_edge.twin == result.half_edge_twin
            {
                Some("twin")
            } else if half_edge.next == result.half_edge || half_edge.next == result.half_edge_twin
            {
                Some("next")
            } else if half_edge.prev == result.half_edge || half_edge.prev == result.half_edge_twin
            {
                Some("prev")
            } else if half_edge.origin == result.vertex {
                Some("origin")
            } else {
                None
            };
            if let Some(field) = stale_reference {
                return Err(CoreError::PreconditionFailed(format!(
                    "kev left live half-edge {field} reference to removed topology"
                )));
            }
        }
        let loop_anchor = self
            .loops
            .get(result.loop_)
            .ok_or_else(|| CoreError::PreconditionFailed("kev loop is not live".into()))?
            .half_edge;
        if self.half_edges.contains_key(loop_anchor) {
            verify_loop_contains_half_edges(self, result.loop_, &[loop_anchor])?;
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }
}

#[cfg(test)]
#[path = "mev_tests.rs"]
mod tests;
