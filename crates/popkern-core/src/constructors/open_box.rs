//! Open-box constructor — a cube with one face removed, producing
//! the first true open-shell BRep in popkern. This is the test
//! vehicle for §8.3 cap insertion on free-boundary fillet endpoints
//! and for task #62 (open-shell support).
//!
//! The implementation starts from a `unit_cube`, identifies one
//! face, removes it and its half-edges, then marks the 4 former
//! boundary edges as `boundary_edges` (twin set to null sentinel).
//! The result has 5 faces, 12 edges (4 boundary + 8 interior),
//! and passes `validate()` with the relaxed manifold checks.

use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, HalfEdgeHandle};

impl BRep {
    /// Build a unit open box — a `[0, 1]³` cube with one face
    /// removed. The removed face is the first face in slotmap
    /// order (deterministic).
    ///
    /// **Topology after removal:**
    /// - 8 vertices, 12 edges, 5 faces
    /// - 4 **boundary edges** (the edges that bordered the
    ///   removed face; their half-edge twin is the null sentinel
    ///   `HalfEdgeHandle::default()`)
    /// - 8 interior edges (each still has 2 half-edges in
    ///   different faces, manifold)
    ///
    /// Passes `validate()` with the relaxed manifold / Euler
    /// checks enabled by `boundary_edges`.
    ///
    /// # Errors
    ///
    /// - [`CoreError::Validation`] if downstream validation fails
    ///   (shouldn't happen for well-formed inputs).
    pub fn open_box(tolerance: Tolerance) -> Result<Self, CoreError> {
        let mut brep = Self::unit_cube(tolerance)?;

        // Pick the first face to remove.
        let face_to_remove = brep
            .faces
            .keys()
            .next()
            .ok_or_else(|| CoreError::PreconditionFailed("unit_cube has no faces".into()))?;

        // Collect the half-edges on the doomed face's outer loop.
        let loop_h = brep.faces[face_to_remove].outer_loop;
        let start_he = brep.loops[loop_h].half_edge;
        let mut doomed_hes = Vec::new();
        let mut cursor = start_he;
        loop {
            doomed_hes.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }

        // For each doomed half-edge, its twin survives but
        // becomes a boundary half-edge (twin = null). Mark the
        // parent edge as boundary.
        for &heh in &doomed_hes {
            let twin_h = brep.half_edges[heh].twin;
            let edge_h = brep.half_edges[heh].edge;
            // Point the edge's canonical half-edge at the surviving
            // twin (it may currently point at the doomed one).
            brep.edges[edge_h].half_edge = twin_h;
            // Null out the twin pointer on the surviving half-edge.
            brep.half_edges[twin_h].twin = HalfEdgeHandle::default();
            brep.boundary_edges.insert(edge_h);
        }

        // Fix vertex outgoing pointers that reference doomed
        // half-edges. For each doomed he, if its origin vertex's
        // outgoing pointer IS that he, switch to the twin (which
        // has a different origin but... actually we need ANY
        // surviving he with this origin).
        let doomed_set: std::collections::HashSet<_> = doomed_hes.iter().copied().collect();
        for &heh in &doomed_hes {
            let origin = brep.half_edges[heh].origin;
            if brep.vertices[origin].outgoing == heh {
                // Find any surviving half-edge with this origin.
                let replacement = brep
                    .half_edges
                    .iter()
                    .find(|(h, he)| he.origin == origin && !doomed_set.contains(h))
                    .map(|(h, _)| h);
                if let Some(rep) = replacement {
                    brep.vertices[origin].outgoing = rep;
                }
            }
        }

        // Remove the doomed half-edges.
        for &heh in &doomed_hes {
            brep.half_edges.remove(heh);
        }

        // Remove the loop and face.
        brep.loops.remove(loop_h);
        brep.faces.remove(face_to_remove);

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_box_topology_counts() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        assert_eq!(ob.vertices.len(), 8);
        assert_eq!(ob.edges.len(), 12);
        assert_eq!(ob.faces.len(), 5);
    }

    #[test]
    fn open_box_has_four_boundary_edges() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        assert_eq!(ob.boundary_edges.len(), 4);
    }

    #[test]
    fn open_box_passes_validation() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        ob.validate().unwrap();
    }

    #[test]
    fn open_box_boundary_edges_have_null_twin() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        for &edge_h in &ob.boundary_edges {
            let he = ob.edges[edge_h].half_edge;
            let twin = ob.half_edges[he].twin;
            assert_eq!(
                twin,
                HalfEdgeHandle::default(),
                "boundary edge's surviving half-edge should have null twin"
            );
        }
    }

    #[test]
    fn open_box_interior_edges_are_manifold() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        for (edge_h, edge) in &ob.edges {
            if ob.boundary_edges.contains(&edge_h) {
                continue;
            }
            let he = &ob.half_edges[edge.half_edge];
            let twin = &ob.half_edges[he.twin];
            assert_eq!(
                twin.edge, edge_h,
                "interior edge twin should reference back"
            );
        }
    }
}
