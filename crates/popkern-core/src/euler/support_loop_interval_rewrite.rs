//! Support-loop interval rewrite — the constructive-winding bridge
//! primitive (PK-167 step 1).
//!
//! The face-family bridge authors a NEW face as the twin cycle of a
//! source support loop: a mirror overlay by construction, which is
//! why the V-free rewrite ships a mirror-wound shell plus an
//! orphaned source flap that the consolidation must repair. This
//! operator is the constructive replacement: it rewrites the SOURCE
//! support face's loop in place — a consecutive interval of its
//! half-edges (the near sides facing the replaced material) is
//! spliced out and replaced by fresh half-edges twinning already-
//! authored boundary sides (strip rails / cap edges). The face and
//! its winding are preserved, so no mirror shell, no flap, and no
//! orientation repair exist to begin with: `dF = 0`, `dV = 0`.
//!
//! Contract:
//! - `removed` is a non-empty consecutive directed sub-path of
//!   `loop_`, and a proper subset of it (the survivor sides keep the
//!   loop's identity and winding).
//! - `replacement` sides are twinless boundary half-edges (their
//!   edges are in `boundary_edges`), chained so their TWIN direction
//!   runs from the removed interval's start vertex to its end
//!   vertex: `twin(r).origin = destination(r)`, so the chain
//!   requirement is `destination(replacement[0]) == origin(removed[0])`,
//!   consecutive `origin(replacement[i]) == destination(replacement[i+1])`,
//!   and `origin(replacement.last()) == destination(removed.last())`.
//!   Vertex identity is required, not geometric nearness — the
//!   caller wires shared anchor vertices before requesting the
//!   rewrite.
//! - Removed sides whose edges lose their last half-edge retire the
//!   edge record; removed sides with a surviving twin expose the
//!   edge as a boundary side.
//!
//! Transactional: applied to a candidate clone; the rewritten loop
//! must close through the new chain and the Euler delta must match
//! (`dV = dF = dL = 0`, `dE = replacement.len() - retired`).

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdge, HalfEdgeHandle, LoopHandle};

/// Plan for one in-place support-loop interval rewrite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupportLoopIntervalRewritePlan {
    /// The support face's loop being rewritten.
    pub loop_: LoopHandle,
    /// Consecutive directed interval of `loop_` to splice out.
    pub removed: Vec<HalfEdgeHandle>,
    /// Twinless boundary sides to twin into the gap, ordered so the
    /// created twins run start-to-end along the loop direction (see
    /// module contract for the chaining rule).
    pub replacement: Vec<HalfEdgeHandle>,
}

/// Result of a committed interval rewrite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupportLoopIntervalRewriteResult {
    /// Fresh half-edges added to the rewritten loop, in loop order.
    pub created: Vec<HalfEdgeHandle>,
    /// Edge records retired because the removed side was their last.
    pub retired_edges: Vec<EdgeHandle>,
    /// Edges of removed sides whose twin survives — now boundary.
    pub exposed_edges: Vec<EdgeHandle>,
}

impl BRep {
    /// Rewrite a consecutive interval of a support loop in place,
    /// twinning already-authored boundary sides into the gap. See
    /// the module docs for the full contract.
    pub fn support_loop_interval_rewrite(
        &mut self,
        plan: &SupportLoopIntervalRewritePlan,
    ) -> Result<SupportLoopIntervalRewriteResult, CoreError> {
        let counts_before = (
            self.vertices.len(),
            self.edges.len(),
            self.faces.len(),
            self.loops.len(),
        );
        let mut candidate = self.clone();
        let result = candidate.apply_support_loop_interval_rewrite(plan)?;
        candidate.verify_support_loop_interval_rewrite(plan, counts_before, &result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_support_loop_interval_rewrite(
        &mut self,
        plan: &SupportLoopIntervalRewritePlan,
    ) -> Result<SupportLoopIntervalRewriteResult, CoreError> {
        let loop_data = self.loops.get(plan.loop_).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "interval rewrite: loop {:?} does not exist",
                plan.loop_
            ))
        })?;
        let face = loop_data.face;
        if plan.removed.is_empty() || plan.replacement.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite needs non-empty removed and replacement intervals".into(),
            ));
        }
        // `removed` must be a consecutive directed path of `loop_`,
        // and a proper subset (the survivors keep the loop alive).
        for &h in &plan.removed {
            let data = self.half_edges.get(h).ok_or_else(|| {
                CoreError::InvalidHandle(format!("interval rewrite: removed side {h:?} is stale"))
            })?;
            if data.loop_ != plan.loop_ {
                return Err(CoreError::PreconditionFailed(format!(
                    "interval rewrite: removed side {h:?} is not on the requested loop"
                )));
            }
        }
        for pair in plan.removed.windows(2) {
            if self.half_edges[pair[0]].next != pair[1] {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite: removed interval is not consecutive in loop order".into(),
                ));
            }
        }
        let before = self.half_edges[plan.removed[0]].prev;
        let last_removed = *plan
            .removed
            .last()
            .ok_or_else(|| CoreError::PreconditionFailed("empty removed interval".into()))?;
        let after = self.half_edges[last_removed].next;
        if plan.removed.contains(&before) || plan.removed.contains(&after) {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite: removed interval covers the whole loop".into(),
            ));
        }
        let start_vertex = self.half_edges[plan.removed[0]].origin;
        let end_vertex = self.half_edges[after].origin;

        // Replacement sides: twinless boundary sides whose twins will
        // chain start_vertex -> end_vertex along the loop direction.
        for &r in &plan.replacement {
            let data = self.half_edges.get(r).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "interval rewrite: replacement side {r:?} is stale"
                ))
            })?;
            if self.half_edges.contains_key(data.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "interval rewrite: replacement side {r:?} already has a live twin"
                )));
            }
            if !self.boundary_edges.contains(&data.edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "interval rewrite: replacement side {r:?} is not a boundary side"
                )));
            }
        }
        fn origin_of(brep: &BRep, h: HalfEdgeHandle) -> crate::topology::VertexHandle {
            brep.half_edges[h].origin
        }
        fn destination_of(brep: &BRep, h: HalfEdgeHandle) -> crate::topology::VertexHandle {
            let next = brep.half_edges[h].next;
            brep.half_edges[next].origin
        }
        if destination_of(self, plan.replacement[0]) != start_vertex {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite: replacement chain does not start at the interval's start \
                 vertex"
                    .into(),
            ));
        }
        for pair in plan.replacement.windows(2) {
            if origin_of(self, pair[0]) != destination_of(self, pair[1]) {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite: replacement sides do not chain vertex-to-vertex".into(),
                ));
            }
        }
        let last_replacement = *plan
            .replacement
            .last()
            .ok_or_else(|| CoreError::PreconditionFailed("empty replacement interval".into()))?;
        if origin_of(self, last_replacement) != end_vertex {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite: replacement chain does not end at the interval's end vertex"
                    .into(),
            ));
        }
        let _ = face;

        // --- Apply ---------------------------------------------------------
        // 1. Create the twin chain in loop order.
        let mut created: Vec<HalfEdgeHandle> = Vec::with_capacity(plan.replacement.len());
        for &r in &plan.replacement {
            let r_data = self.half_edges[r].clone();
            let origin = destination_of(self, r);
            let new_h = self.half_edges.insert(HalfEdge {
                origin,
                twin: r,
                next: HalfEdgeHandle::default(),
                prev: HalfEdgeHandle::default(),
                edge: r_data.edge,
                loop_: plan.loop_,
                pcurve: None,
            });
            self.half_edges[r].twin = new_h;
            self.boundary_edges.remove(&r_data.edge);
            created.push(new_h);
        }
        // 2. Splice the chain between `before` and `after`.
        let mut cursor = before;
        for &h in &created {
            self.half_edges[cursor].next = h;
            self.half_edges[h].prev = cursor;
            cursor = h;
        }
        self.half_edges[cursor].next = after;
        self.half_edges[after].prev = cursor;
        if self.loops[plan.loop_].half_edge == plan.removed[0]
            || plan.removed.contains(&self.loops[plan.loop_].half_edge)
        {
            self.loops[plan.loop_].half_edge = created[0];
        }
        // 3. Retire the removed sides.
        let mut retired_edges = Vec::new();
        let mut exposed_edges = Vec::new();
        for &h in &plan.removed {
            let data = self.half_edges.remove(h).ok_or_else(|| {
                CoreError::InvalidHandle(format!("interval rewrite: removed side {h:?} vanished"))
            })?;
            if self.half_edges.contains_key(data.twin) {
                self.boundary_edges.insert(data.edge);
                if let Some(edge) = self.edges.get_mut(data.edge)
                    && edge.half_edge == h
                {
                    edge.half_edge = data.twin;
                }
                exposed_edges.push(data.edge);
            } else if self.edges.remove(data.edge).is_some() {
                self.boundary_edges.remove(&data.edge);
                retired_edges.push(data.edge);
            }
        }
        Ok(SupportLoopIntervalRewriteResult {
            created,
            retired_edges,
            exposed_edges,
        })
    }

    fn verify_support_loop_interval_rewrite(
        &self,
        plan: &SupportLoopIntervalRewritePlan,
        counts_before: (usize, usize, usize, usize),
        result: &SupportLoopIntervalRewriteResult,
    ) -> Result<(), CoreError> {
        // Loop closes through the new chain.
        let start = self.loops[plan.loop_].half_edge;
        let mut cursor = start;
        let mut seen_created = 0usize;
        for _ in 0..self.half_edges.len() + 1 {
            let Some(data) = self.half_edges.get(cursor) else {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite postcondition: rewritten loop references a removed half-edge"
                        .into(),
                ));
            };
            if data.loop_ != plan.loop_ {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite postcondition: rewritten loop walks into another loop".into(),
                ));
            }
            if result.created.contains(&cursor) {
                seen_created += 1;
            }
            let next = data.next;
            if self.half_edges.get(next).map(|n| n.prev) != Some(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite: next/prev disagree on the rewritten loop".into(),
                ));
            }
            cursor = next;
            if cursor == start {
                break;
            }
        }
        if cursor != start {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite postcondition: rewritten loop does not close".into(),
            ));
        }
        if seen_created != result.created.len() {
            return Err(CoreError::PreconditionFailed(
                "interval rewrite postcondition: rewritten loop does not contain the full new chain".into(),
            ));
        }
        // Twin symmetry for the new chain.
        for &h in &result.created {
            let twin = self.half_edges[h].twin;
            if self.half_edges.get(twin).map(|t| t.twin) != Some(h) {
                return Err(CoreError::PreconditionFailed(
                    "interval rewrite postcondition: created side's twin does not point back"
                        .into(),
                ));
            }
        }
        // Euler deltas: faces, loops, vertices unchanged; edges may
        // only shrink by the retired count.
        let (v0, e0, f0, l0) = counts_before;
        let (v1, e1, f1, l1) = (
            self.vertices.len(),
            self.edges.len(),
            self.faces.len(),
            self.loops.len(),
        );
        if v1 != v0 || f1 != f0 || l1 != l0 {
            return Err(CoreError::PreconditionFailed(format!(
                "interval rewrite postcondition: dV={} dF={} dL={} (all must be 0)",
                v1 as i64 - v0 as i64,
                f1 as i64 - f0 as i64,
                l1 as i64 - l0 as i64
            )));
        }
        if e1 + result.retired_edges.len() != e0 {
            return Err(CoreError::PreconditionFailed(format!(
                "interval rewrite postcondition: dE={} but retired {}",
                e1 as i64 - e0 as i64,
                result.retired_edges.len()
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "support_loop_interval_rewrite_tests.rs"]
mod tests;
