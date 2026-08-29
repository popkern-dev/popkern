//! Vertex-welding primitives for [`BRep`].
//!
//! Two operations:
//!
//! - [`BRep::weld_coincident_vertices`] — globally merge every
//!   pair of vertices within `tolerance.linear()` into a single
//!   canonical representative. Used by the boolean engine to
//!   collapse vertices that two operands had at the same 3D
//!   position after concatenation.
//!
//! - [`BRep::weld_vertex_pair`] — targeted weld of one specific
//!   vertex into another. Used by callers (such as unfillet) that
//!   want to consolidate vertices they have moved to the same
//!   position, without touching other coincident pairs.
//!
//! Both functions update every half-edge whose `origin` was the
//! welded-away vertex to use the canonical representative, then
//! delete the welded-away vertex from the slotmap. Edge and
//! face structures are left intact — this is purely a vertex
//! operation. Callers that need to also collapse zero-length
//! edges resulting from welding should follow up with their own
//! edge-cleanup pass.

use std::collections::HashMap;

use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, VertexHandle};

impl BRep {
    /// Weld every pair of vertices whose 3D positions are within
    /// `tolerance.linear()` of each other into a single canonical
    /// representative.
    ///
    /// Walks all vertex pairs (O(n²)), unions coincident ones via
    /// union-find, rewrites every half-edge's `origin` to point at
    /// the canonical vertex of its class, and removes the
    /// non-canonical vertices from the slotmap.
    ///
    /// Returns the number of vertices welded away (i.e. removed).
    /// Edge / face / loop structures are not modified — only
    /// vertex references are updated.
    ///
    /// # Performance
    ///
    /// Runtime is O(n²) in vertex count, which is acceptable for
    /// typical sheet-body sizes seen by the boolean engine. A
    /// spatial index would only help for very large imprints.
    pub fn weld_coincident_vertices(&mut self, tolerance: &Tolerance) -> usize {
        let eps = tolerance.linear();
        let eps_sq = eps * eps;

        // Snapshot the (handle, position) pairs before mutating.
        let all: Vec<(VertexHandle, crate::geometry::Point3)> =
            self.vertices.iter().map(|(h, v)| (h, v.position)).collect();
        if all.len() < 2 {
            return 0;
        }

        // Union-find: each vertex starts as its own class.
        let mut parent: HashMap<VertexHandle, VertexHandle> =
            all.iter().map(|(h, _)| (*h, *h)).collect();

        // For each coincident pair, union.
        for (i, &(vh_i, pos_i)) in all.iter().enumerate() {
            for &(vh_j, pos_j) in all.iter().skip(i + 1) {
                let dx = pos_i.x - pos_j.x;
                let dy = pos_i.y - pos_j.y;
                let dz = pos_i.z - pos_j.z;
                if dx * dx + dy * dy + dz * dz <= eps_sq {
                    let root_i = uf_find(&mut parent, vh_i);
                    let root_j = uf_find(&mut parent, vh_j);
                    if root_i != root_j {
                        parent.insert(root_j, root_i);
                    }
                }
            }
        }

        // Resolve canonical handles after all unions.
        let keys: Vec<VertexHandle> = parent.keys().copied().collect();
        let mut canonical: HashMap<VertexHandle, VertexHandle> = HashMap::new();
        for vh in keys {
            let root = uf_find(&mut parent, vh);
            canonical.insert(vh, root);
        }

        // Rewrite half-edge origins to point at canonical reps.
        for (_, he) in self.half_edges.iter_mut() {
            if let Some(&canon) = canonical.get(&he.origin)
                && canon != he.origin
            {
                he.origin = canon;
            }
        }

        // Delete non-canonical vertices.
        let to_remove: Vec<VertexHandle> = canonical
            .iter()
            .filter_map(|(&vh, &canon)| if vh != canon { Some(vh) } else { None })
            .collect();
        let welds = to_remove.len();
        for vh in to_remove {
            self.vertices.remove(vh);
        }
        welds
    }

    /// Weld a specific `secondary` vertex into a `primary` vertex.
    ///
    /// After this call, every half-edge whose `origin` was
    /// `secondary` will use `primary` instead, and `secondary` is
    /// removed from the BRep. The primary's position is left
    /// unchanged — the caller is responsible for moving both
    /// vertices to the same 3D position before welding (or for
    /// accepting that the merged vertex's position is the
    /// primary's).
    ///
    /// This is the surgical counterpart to
    /// [`Self::weld_coincident_vertices`]: callers that want to
    /// merge specific vertices without affecting others use this.
    /// Used by `popkern_blend::unfillet` to collapse fillet
    /// corner vertices into restored sharp-edge corners.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if either handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if `primary == secondary`.
    pub fn weld_vertex_pair(
        &mut self,
        primary: VertexHandle,
        secondary: VertexHandle,
    ) -> Result<(), CoreError> {
        if primary == secondary {
            return Err(CoreError::PreconditionFailed(
                "weld_vertex_pair requires distinct primary and secondary vertices".into(),
            ));
        }
        if !self.vertices.contains_key(primary) {
            return Err(CoreError::InvalidHandle(format!(
                "weld_vertex_pair: primary {primary:?} does not exist"
            )));
        }
        if !self.vertices.contains_key(secondary) {
            return Err(CoreError::InvalidHandle(format!(
                "weld_vertex_pair: secondary {secondary:?} does not exist"
            )));
        }

        // Rewrite half-edge origins.
        for (_, he) in self.half_edges.iter_mut() {
            if he.origin == secondary {
                he.origin = primary;
            }
        }

        // Remove the secondary vertex.
        self.vertices.remove(secondary);

        Ok(())
    }
}

/// Union-find helper with path compression. Free function so the
/// recursive call doesn't require capturing `&mut self`.
fn uf_find(parent: &mut HashMap<VertexHandle, VertexHandle>, vh: VertexHandle) -> VertexHandle {
    let p = parent[&vh];
    if p == vh {
        return vh;
    }
    let root = uf_find(parent, p);
    parent.insert(vh, root);
    root
}

#[cfg(test)]
#[path = "brep_weld_tests.rs"]
mod tests;
