//! split_vertex — Mäntylä's general MEV: split a vertex fan.
//!
//! Moves a contiguous fan sector of edges at a vertex to a new vertex and
//! splices a connector edge into the two cut wedges. The spur-only `mev`
//! cannot express this: it always inserts a zero-area detour into one
//! loop, while `split_vertex` re-homes existing edges (ΔV+1, ΔE+1, no
//! loop or face change). The shell gap-join corner choreography (design
//! §5 Stage 3a) peels per-face vertex runs with it; the connector edges
//! become the corner sphere arcs.

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts, verify_loop_contains_half_edges};
use crate::geometry::Point3;
use crate::topology::{BRep, Edge, HalfEdge, HalfEdgeHandle, Vertex};

/// Result of the split_vertex operator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitVertexResult {
    /// The newly created vertex (took the moved fan sector).
    pub vertex: crate::topology::VertexHandle,
    /// The newly created connector edge between the old and new vertex.
    pub edge: crate::topology::EdgeHandle,
    /// Connector half-edge from the old vertex to the new vertex. Lives in
    /// the loop of the wedge cut before `he_first`.
    pub half_edge: HalfEdgeHandle,
    /// Connector half-edge from the new vertex to the old vertex (twin).
    /// Lives in the loop of the wedge cut after `he_last`.
    pub half_edge_twin: HalfEdgeHandle,
    /// The outgoing half-edges whose origin moved to the new vertex, in
    /// fan order from `he_first` to `he_last`.
    pub moved: Vec<HalfEdgeHandle>,
}

impl BRep {
    /// Split a vertex fan: move the contiguous sector of outgoing
    /// half-edges from `he_first` to `he_last` (inclusive, walking the fan
    /// successor `twin(o).next`) to a new vertex at `position`, and
    /// connect old and new vertex with a new edge spliced into the two
    /// cut wedges.
    ///
    /// Both half-edges must originate at the same vertex, the fan walk
    /// between them must not cross a boundary (null-twin) edge, and the
    /// sector must leave at least one edge at the original vertex.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if either half-edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the half-edges originate at
    ///   different vertices, the fan wraps before reaching `he_last`, the
    ///   sector would consume the entire fan, or the fan is not walkable.
    pub fn split_vertex(
        &mut self,
        he_first: HalfEdgeHandle,
        he_last: HalfEdgeHandle,
        position: Point3,
    ) -> Result<SplitVertexResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_split_vertex(he_first, he_last, position)?;
        candidate.verify_split_vertex_postconditions(before, &result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_split_vertex(
        &mut self,
        he_first: HalfEdgeHandle,
        he_last: HalfEdgeHandle,
        position: Point3,
    ) -> Result<SplitVertexResult, CoreError> {
        let first_data = self
            .half_edges
            .get(he_first)
            .ok_or_else(|| CoreError::InvalidHandle("he_first handle is invalid".into()))?;
        let last_data = self
            .half_edges
            .get(he_last)
            .ok_or_else(|| CoreError::InvalidHandle("he_last handle is invalid".into()))?;

        let old_vertex = first_data.origin;
        if last_data.origin != old_vertex {
            return Err(CoreError::PreconditionFailed(
                "he_first and he_last must originate at the same vertex".into(),
            ));
        }

        // Collect the sector by walking the fan successor twin(o).next.
        let mut moved = vec![he_first];
        let bound = self.half_edges.len();
        while *moved.last().unwrap_or(&he_first) != he_last {
            let current = moved[moved.len() - 1];
            let twin = self.half_edges[current].twin;
            let twin_data = self.half_edges.get(twin).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "split_vertex fan walk crossed a boundary (null-twin) edge".into(),
                )
            })?;
            let next_out = twin_data.next;
            if next_out == he_first {
                return Err(CoreError::PreconditionFailed(
                    "he_last is not fan-reachable from he_first without wrapping".into(),
                ));
            }
            if self
                .half_edges
                .get(next_out)
                .is_none_or(|he| he.origin != old_vertex)
            {
                return Err(CoreError::PreconditionFailed(
                    "split_vertex fan walk left the vertex (inconsistent fan)".into(),
                ));
            }
            if moved.len() > bound {
                return Err(CoreError::PreconditionFailed(
                    "split_vertex fan walk did not terminate".into(),
                ));
            }
            moved.push(next_out);
        }

        // Sector must leave at least one edge at the old vertex: the
        // outgoing half-edge after he_last must not be he_first.
        let twin_last = self.half_edges[he_last].twin;
        let twin_last_data = self.half_edges.get(twin_last).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "split_vertex fan walk crossed a boundary (null-twin) edge".into(),
            )
        })?;
        let after_sector = twin_last_data.next;
        if after_sector == he_first {
            return Err(CoreError::PreconditionFailed(
                "split_vertex sector must leave at least one edge at the vertex".into(),
            ));
        }

        // Cut wedge 1: (prev(he_first), he_first) in he_first's loop.
        // Cut wedge 2: (twin(he_last), after_sector) in twin(he_last)'s loop.
        let cut1_in = self.half_edges[he_first].prev;
        let cut1_loop = self.half_edges[he_first].loop_;
        let cut2_loop = twin_last_data.loop_;

        // Create the new vertex and connector edge.
        let new_vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        // Connector halves: c_vu (old -> new) into cut wedge 1,
        // c_uv (new -> old) into cut wedge 2.
        let c_vu = self.half_edges.insert(HalfEdge {
            origin: old_vertex,
            twin: Default::default(),
            next: he_first,
            prev: cut1_in,
            edge,
            loop_: cut1_loop,
            pcurve: None,
        });
        let c_uv = self.half_edges.insert(HalfEdge {
            origin: new_vertex,
            twin: c_vu,
            next: after_sector,
            prev: twin_last,
            edge,
            loop_: cut2_loop,
            pcurve: None,
        });
        self.half_edges[c_vu].twin = c_uv;
        self.edges[edge].half_edge = c_vu;

        // Re-home the sector.
        for &he in &moved {
            self.half_edges[he].origin = new_vertex;
        }

        // Splice cut wedge 1: cut1_in -> c_vu -> he_first.
        self.half_edges[cut1_in].next = c_vu;
        self.half_edges[he_first].prev = c_vu;

        // Splice cut wedge 2: twin(he_last) -> c_uv -> after_sector.
        self.half_edges[twin_last].next = c_uv;
        self.half_edges[after_sector].prev = c_uv;

        // Anchors: the new vertex owns c_uv; the old vertex may have
        // pointed at a moved half-edge.
        self.vertices[new_vertex].outgoing = c_uv;
        let old_anchor = self.vertices[old_vertex].outgoing;
        if moved.contains(&old_anchor) {
            self.vertices[old_vertex].outgoing = c_vu;
        }

        Ok(SplitVertexResult {
            vertex: new_vertex,
            edge,
            half_edge: c_vu,
            half_edge_twin: c_uv,
            moved,
        })
    }

    fn verify_split_vertex_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: &SplitVertexResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "split_vertex",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mev(),
        )?;
        let edge = self
            .edges
            .get(result.edge)
            .ok_or_else(|| CoreError::PreconditionFailed("split_vertex edge is not live".into()))?;
        if edge.half_edge != result.half_edge {
            return Err(CoreError::PreconditionFailed(
                "split_vertex edge does not reference the connector half-edge".into(),
            ));
        }
        let c_vu = self.half_edges.get(result.half_edge).ok_or_else(|| {
            CoreError::PreconditionFailed("split_vertex connector half-edge is not live".into())
        })?;
        let c_uv = self.half_edges.get(result.half_edge_twin).ok_or_else(|| {
            CoreError::PreconditionFailed("split_vertex connector twin is not live".into())
        })?;
        if c_vu.twin != result.half_edge_twin
            || c_uv.twin != result.half_edge
            || c_vu.edge != result.edge
            || c_uv.edge != result.edge
            || c_uv.origin != result.vertex
        {
            return Err(CoreError::PreconditionFailed(
                "split_vertex connector pair is not consistently wired".into(),
            ));
        }
        for &he in &result.moved {
            if self
                .half_edges
                .get(he)
                .is_none_or(|data| data.origin != result.vertex)
            {
                return Err(CoreError::PreconditionFailed(
                    "split_vertex moved half-edge does not originate at the new vertex".into(),
                ));
            }
        }
        // A cut wedge may be a sheet BOUNDARY (loopless half-edges in
        // the materialized-twin rim form): the connector half spliced
        // into it inherits the dead loop handle and joins the free rim
        // ring — there is no live loop to verify (PK-326 rim corners).
        if self.loops.contains_key(c_vu.loop_) {
            verify_loop_contains_half_edges(self, c_vu.loop_, &[result.half_edge])?;
        }
        if self.loops.contains_key(c_uv.loop_) {
            verify_loop_contains_half_edges(self, c_uv.loop_, &[result.half_edge_twin])?;
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)
    }
}

#[cfg(test)]
#[path = "split_vertex_tests.rs"]
mod tests;
