//! Identity queries and recovery operations for BRep precision tracking.
//!
//! Provides methods to check entity coincidence, imprecision overlap,
//! and to recover (tighten) entity tolerances.

use super::brep::BRep;
use crate::error::CoreError;
use crate::precision::EntityHandle;
use crate::topology::{EdgeHandle, VertexHandle};

impl BRep {
    /// Are these two entities the same within tracked imprecision?
    ///
    /// For vertex-vertex comparisons, returns `true` when the distance
    /// between the two positions is less than the maximum of their
    /// individual tolerances.
    ///
    /// Edge-edge and face-face comparisons are not yet implemented and
    /// return `false`.
    pub fn entities_coincident(&self, a: EntityHandle, b: EntityHandle) -> bool {
        match (a, b) {
            (EntityHandle::Vertex(va), EntityHandle::Vertex(vb)) => {
                let pos_a = self.vertices[va].position;
                let pos_b = self.vertices[vb].position;
                let tol_a = self.vertex_tolerance(va);
                let tol_b = self.vertex_tolerance(vb);
                let threshold = tol_a.max(tol_b);
                crate::adaptive_precision::adaptive_dist_sq_lt(
                    (pos_a.x, pos_a.y, pos_a.z),
                    (pos_b.x, pos_b.y, pos_b.z),
                    threshold * threshold,
                )
            }
            // Edge-edge and face-face: future work.
            _ => false,
        }
    }

    /// Do two distinct entities have overlapping imprecision regions?
    ///
    /// For vertex-vertex comparisons, returns `true` when the distance
    /// between the two positions is less than the sum of their
    /// individual tolerances.
    ///
    /// Edge-edge and face-face comparisons are not yet implemented and
    /// return `false`.
    pub fn imprecision_overlap(&self, a: EntityHandle, b: EntityHandle) -> bool {
        match (a, b) {
            (EntityHandle::Vertex(va), EntityHandle::Vertex(vb)) => {
                let pos_a = self.vertices[va].position;
                let pos_b = self.vertices[vb].position;
                let tol_a = self.vertex_tolerance(va);
                let tol_b = self.vertex_tolerance(vb);
                let threshold = tol_a + tol_b;
                crate::adaptive_precision::adaptive_dist_sq_lt(
                    (pos_a.x, pos_a.y, pos_a.z),
                    (pos_b.x, pos_b.y, pos_b.z),
                    threshold * threshold,
                )
            }
            // Edge-edge and face-face: future work.
            _ => false,
        }
    }

    /// Re-project edge curve onto adjacent surface, regenerating pcurve.
    ///
    /// Returns the new (tighter) imprecision after recovery.
    ///
    /// Current implementation is a stub: it reads the current edge tolerance,
    /// sets it via [`recover_entity_tolerance`], and returns the value.
    /// Future versions will perform actual re-intersection/re-projection.
    pub fn recover_edge_consistency(&mut self, eh: EdgeHandle) -> Result<f64, CoreError> {
        if !self.edges.contains_key(eh) {
            return Err(CoreError::InvalidHandle(format!("EdgeHandle {eh:?}")));
        }
        let tol = self.edge_tolerance(eh);
        self.recover_entity_tolerance(EntityHandle::Edge(eh), tol);
        Ok(tol)
    }

    /// Recompute vertex position from incident edge curves.
    ///
    /// Returns the new (tighter) imprecision after recovery.
    ///
    /// Current implementation is a stub: it reads the current vertex tolerance,
    /// sets it via [`recover_entity_tolerance`], and returns the value.
    /// Future versions will average incident edge endpoints and tighten.
    pub fn recover_vertex_precision(&mut self, vh: VertexHandle) -> Result<f64, CoreError> {
        if !self.vertices.contains_key(vh) {
            return Err(CoreError::InvalidHandle(format!("VertexHandle {vh:?}")));
        }
        let tol = self.vertex_tolerance(vh);
        self.recover_entity_tolerance(EntityHandle::Vertex(vh), tol);
        Ok(tol)
    }
}

#[cfg(test)]
#[path = "brep_precision_queries_tests.rs"]
mod brep_precision_queries_tests;
