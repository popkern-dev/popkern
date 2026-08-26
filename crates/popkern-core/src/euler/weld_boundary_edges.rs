//! Boundary-seam weld: fuse two geometrically coincident boundary
//! edges into one interior edge, merging their endpoint vertices.
//!
//! This is the re-close move of global self-intersection trimming
//! (PK-327): splitting two colliding faces along their shared
//! intersection curve authors one edge per face; removing the
//! unrealized region leaves the two copies as opposing boundary rims,
//! and the weld twins their surviving half-edges into a single
//! manifold edge.

use crate::error::CoreError;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve::CurveOnEdge;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

impl BRep {
    /// Weld boundary edge `kill` onto boundary edge `keep`: the two
    /// surviving half-edges become twins, `kill`'s endpoint vertices
    /// merge into `keep`'s, and `kill`'s edge record is deleted.
    /// `keep` keeps its curve geometry and leaves the boundary set.
    ///
    /// The two edges must traverse the same curve in OPPOSITE
    /// directions (manifold orientation): `kill`'s surviving
    /// half-edge must start where `keep`'s ends and end where it
    /// starts, within `tolerance`.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if either edge is missing.
    /// - [`CoreError::PreconditionFailed`] if either edge is not a
    ///   boundary edge, the endpoints do not oppose within tolerance,
    ///   or the curve midpoints disagree.
    /// - [`CoreError::Validation`] if the welded body fails structural
    ///   validation (the model is left untouched).
    pub fn weld_boundary_edge_pair(
        &mut self,
        keep: EdgeHandle,
        kill: EdgeHandle,
        tolerance: &Tolerance,
    ) -> Result<(), CoreError> {
        if keep == kill {
            return Err(CoreError::PreconditionFailed(
                "weld_boundary_edge_pair: keep and kill are the same edge".into(),
            ));
        }
        for eh in [keep, kill] {
            if !self.edges.contains_key(eh) {
                return Err(CoreError::InvalidHandle(format!(
                    "weld_boundary_edge_pair: edge {eh:?} is not in the body"
                )));
            }
            if !self.boundary_edges.contains(&eh) {
                return Err(CoreError::PreconditionFailed(format!(
                    "weld_boundary_edge_pair: edge {eh:?} is not a boundary edge"
                )));
            }
        }

        let tol = tolerance.linear() * 10.0;
        let hk = self.edges[keep].half_edge;
        let hx = self.edges[kill].half_edge;
        for (label, heh) in [("keep", hk), ("kill", hx)] {
            let Some(he) = self.half_edges.get(heh) else {
                return Err(CoreError::PreconditionFailed(format!(
                    "weld_boundary_edge_pair: {label} edge has no surviving half-edge"
                )));
            };
            if self.half_edges.contains_key(he.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "weld_boundary_edge_pair: {label} edge's half-edge still has a twin"
                )));
            }
        }

        // Endpoints: destination = origin of next-in-loop.
        let dest_of = |brep: &BRep, heh| {
            let next = brep.half_edges[heh].next;
            brep.half_edges[next].origin
        };
        let (k0, k1) = (self.half_edges[hk].origin, dest_of(self, hk));
        let (x0, x1) = (self.half_edges[hx].origin, dest_of(self, hx));
        let pos = |brep: &BRep, v: VertexHandle| brep.vertices[v].position;
        let d_opposed = (pos(self, x0) - pos(self, k1))
            .length()
            .max((pos(self, x1) - pos(self, k0)).length());
        if d_opposed > tol {
            return Err(CoreError::PreconditionFailed(format!(
                "weld_boundary_edge_pair: endpoints do not oppose (gap {d_opposed:.3e} > {tol:.3e})"
            )));
        }

        // Curve agreement at the midpoint (both curves must describe
        // the same 3D locus — the weld keeps `keep`'s geometry).
        if let (Some(ck), Some(cx)) = (
            self.edges[keep].curve.as_ref(),
            self.edges[kill].curve.as_ref(),
        ) {
            let mid = |c: &CurveOnEdge| c.geometry.point_at((c.t_start + c.t_end) * 0.5);
            let (mk, mx) = (mid(ck), mid(cx));
            let gap = (mk - mx).length();
            if gap > tol {
                return Err(CoreError::PreconditionFailed(format!(
                    "weld_boundary_edge_pair: curve midpoints disagree by {gap:.3e} > {tol:.3e}"
                )));
            }
        }

        let mut candidate = self.clone();

        // Merge kill-side endpoint vertices into keep's (x0 → k1,
        // x1 → k0), unless the split already shared them.
        for (from, to) in [(x0, k1), (x1, k0)] {
            if from == to || !candidate.vertices.contains_key(from) {
                continue;
            }
            candidate.weld_vertex_pair(to, from)?;
        }

        // Twin the survivors and retire the kill edge.
        candidate.half_edges[hk].twin = hx;
        candidate.half_edges[hx].twin = hk;
        candidate.half_edges[hx].edge = keep;
        candidate.edges.remove(kill);
        candidate.boundary_edges.remove(&keep);
        candidate.boundary_edges.remove(&kill);

        // Local postconditions only: a weld BATCH passes through
        // intermediate states (a slit welded half-way) whose global
        // Euler bookkeeping is legitimately off — the batch owner
        // validates the whole body after the last weld. Here: mutual
        // twins and both incident loops still walk closed.
        if candidate.half_edges[hk].twin != hx || candidate.half_edges[hx].twin != hk {
            return Err(CoreError::PreconditionFailed(
                "weld_boundary_edge_pair: twin rewire failed".into(),
            ));
        }
        for heh in [hk, hx] {
            let loop_h = candidate.half_edges[heh].loop_;
            if candidate.loops.contains_key(loop_h) {
                let start = candidate.loops[loop_h].half_edge;
                let mut cursor = start;
                let bound = candidate.half_edges.len().max(1);
                let mut closed = false;
                for _ in 0..=bound {
                    let Some(he) = candidate.half_edges.get(cursor) else {
                        break;
                    };
                    cursor = he.next;
                    if cursor == start {
                        closed = true;
                        break;
                    }
                }
                if !closed {
                    return Err(CoreError::PreconditionFailed(
                        "weld_boundary_edge_pair: an incident loop no longer walks closed".into(),
                    ));
                }
            }
        }

        *self = candidate;
        Ok(())
    }
}
