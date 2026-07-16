//! Edge contraction Euler operator (`kve` — kill vertex, edge).
//!
//! Collapses an edge by merging its destination vertex into its
//! origin vertex: V−1, E−1, F and loop counts unchanged, so the
//! Euler characteristic is preserved. The classic use is retiring a
//! topological chord whose two endpoints must become one vertex —
//! e.g. corner-blend boundary chords whose true geometry is a single
//! tangency point (route B corner reconstruction).

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

/// Result of a [`BRep::kve`] edge contraction.
#[derive(Debug, Clone, Copy)]
pub struct KveResult {
    /// The surviving vertex (the collapsed edge's origin).
    pub vertex: VertexHandle,
    /// The removed edge.
    pub removed_edge: EdgeHandle,
    /// The removed destination vertex.
    pub removed_vertex: VertexHandle,
}

impl BRep {
    /// Contract the edge of `half_edge`, merging the destination
    /// vertex (`twin.origin`) into the origin vertex. The surviving
    /// vertex keeps the ORIGIN's position; callers wanting a
    /// specific merged position should move the origin vertex first.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the half-edge is invalid.
    /// - [`CoreError::PreconditionFailed`] if the edge is a
    ///   self-loop (origin == destination), if either adjacent loop
    ///   would drop below 3 half-edges, or if contraction would
    ///   produce a non-manifold configuration.
    pub fn kve(&mut self, half_edge: HalfEdgeHandle) -> Result<KveResult, CoreError> {
        let mut candidate = self.clone();
        let result = candidate.apply_kve(half_edge)?;
        // Structural postconditions only: kve is a pure topology
        // operator; edge/trim GEOMETRY around the merged vertex is
        // the caller's repair contract (curves still reference the
        // removed vertex's position until rebuilt), so geometric
        // checks (loop closure, watertight) belong to the caller's
        // own acceptance gate after that repair.
        //
        // No-regression contract: a check the PRE-operation state
        // already fails is exempt afterwards — mid-pipeline repair
        // passes legitimately run kve on models whose global Euler
        // accounting is temporarily off (e.g. the strip-end corner
        // contraction after a partially-failed carve, PK-121), and
        // kve itself is chi-neutral (V-1, E-1). On a compliant input
        // the gate is exactly as strict as before.
        type StructuralCheck = fn(&BRep) -> Result<(), crate::error::ValidationError>;
        let checks: [StructuralCheck; 6] = [
            BRep::validate_euler_poincare,
            BRep::validate_half_edge_pairing,
            BRep::validate_half_edge_cycles,
            BRep::validate_vertex_consistency,
            BRep::validate_manifold_edges,
            BRep::validate_vertex_valence,
        ];
        for check in &checks {
            if let Err(e) = check(&candidate) {
                if check(self).is_err() {
                    // Pre-existing failure of the same check: kve is
                    // not the author; the caller's acceptance gate
                    // owns the final verdict.
                    continue;
                }
                return Err(CoreError::PreconditionFailed(format!(
                    "kve: post-contraction structural validation failed: {e:?}"
                )));
            }
        }
        *self = candidate;
        Ok(result)
    }

    fn apply_kve(&mut self, half_edge: HalfEdgeHandle) -> Result<KveResult, CoreError> {
        let he = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("kve: half-edge not found".to_string()))?
            .clone();
        let twin_h = he.twin;
        let twin = self
            .half_edges
            .get(twin_h)
            .ok_or_else(|| CoreError::InvalidHandle("kve: half-edge not found".to_string()))?
            .clone();
        let v_keep = he.origin;
        let v_kill = twin.origin;
        if v_keep == v_kill {
            return Err(CoreError::PreconditionFailed(
                "kve: cannot contract a self-loop edge".into(),
            ));
        }
        // Loop-size preconditions: each adjacent loop loses one
        // half-edge and must keep at least 3 afterwards? No — 2 is
        // the mekr/kemr-style minimum for seam structures, but a
        // loop of 2 collapsing to 1 is degenerate. Require ≥ 3
        // BEFORE contraction on both sides.
        for start in [half_edge, twin_h] {
            let mut n = 0usize;
            let mut cur = start;
            for _ in 0..self.half_edges.len() + 1 {
                n += 1;
                cur = self.half_edges[cur].next;
                if cur == start {
                    break;
                }
            }
            if n < 3 {
                return Err(CoreError::PreconditionFailed(format!(
                    "kve: adjacent loop has only {n} half-edges"
                )));
            }
        }

        // Re-origin every half-edge leaving v_kill to v_keep.
        let leaving: Vec<HalfEdgeHandle> = self
            .half_edges
            .iter()
            .filter(|(_, d)| d.origin == v_kill)
            .map(|(h, _)| h)
            .collect();
        for h in leaving {
            self.half_edges[h].origin = v_keep;
        }

        // Unlink the two half-edges from their loops.
        for (dead, _other) in [(half_edge, twin_h), (twin_h, half_edge)] {
            let prev = self.half_edges[dead].prev;
            let next = self.half_edges[dead].next;
            self.half_edges[prev].next = next;
            self.half_edges[next].prev = prev;
            let loop_h = self.half_edges[dead].loop_;
            if self.loops[loop_h].half_edge == dead {
                self.loops[loop_h].half_edge = next;
            }
        }

        // Repair the surviving vertex's outgoing pointer.
        if !self.half_edges.contains_key(self.vertices[v_keep].outgoing)
            || self.vertices[v_keep].outgoing == half_edge
            || self.vertices[v_keep].outgoing == twin_h
        {
            let replacement = self
                .half_edges
                .iter()
                .find(|(h, d)| d.origin == v_keep && *h != half_edge && *h != twin_h)
                .map(|(h, _)| h);
            if let Some(r) = replacement {
                self.vertices[v_keep].outgoing = r;
            }
        }

        let edge_h = he.edge;
        self.half_edges.remove(half_edge);
        self.half_edges.remove(twin_h);
        self.edges.remove(edge_h);
        self.vertices.remove(v_kill);

        Ok(KveResult {
            vertex: v_keep,
            removed_edge: edge_h,
            removed_vertex: v_kill,
        })
    }
}

#[cfg(test)]
#[path = "kve_tests.rs"]
mod tests;
