//! Half-edge structural checks: pairing, cycles, and vertex consistency.

use crate::error::ValidationError;
use crate::topology::{BRep, TopologyRole};

impl BRep {
    /// Check 2: Every half-edge's twin must exist and the twin's twin must point back.
    pub fn validate_half_edge_pairing(&self) -> Result<(), ValidationError> {
        for (he_handle, he_data) in &self.half_edges {
            // Boundary half-edges have twin = default (null) — skip.
            if self.boundary_edges.contains(&he_data.edge) {
                continue;
            }
            let twin = self
                .half_edges
                .get(he_data.twin)
                .ok_or(ValidationError::BrokenTwin {
                    half_edge: he_handle,
                })?;

            if twin.twin != he_handle {
                return Err(ValidationError::BrokenTwin {
                    half_edge: he_handle,
                });
            }
        }
        Ok(())
    }

    /// Check 3: Following `next` from any half-edge must return to the start
    /// within a bounded number of steps.
    pub fn validate_half_edge_cycles(&self) -> Result<(), ValidationError> {
        let bound = self.half_edges.len();

        for (he_handle, he_data) in &self.half_edges {
            if self.half_edge_is_construction_residual_for_cycle_validation(he_data) {
                continue;
            }
            let mut cursor = he_data.next;
            let mut steps = 0;

            while cursor != he_handle {
                steps += 1;
                if steps > bound {
                    return Err(ValidationError::BrokenCycle {
                        half_edge: he_handle,
                        loop_: he_data.loop_,
                    });
                }
                match self.half_edges.get(cursor) {
                    Some(next_he) => {
                        if self.half_edge_is_construction_residual_for_cycle_validation(next_he) {
                            return Err(ValidationError::BrokenCycle {
                                half_edge: he_handle,
                                loop_: he_data.loop_,
                            });
                        }
                        cursor = next_he.next;
                    }
                    None => {
                        return Err(ValidationError::BrokenCycle {
                            half_edge: he_handle,
                            loop_: he_data.loop_,
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn half_edge_is_construction_residual_for_cycle_validation(
        &self,
        half_edge: &crate::topology::HalfEdge,
    ) -> bool {
        self.edges
            .get(half_edge.edge)
            .is_some_and(|edge| edge.role == TopologyRole::ConstructionResidual)
            || self
                .loops
                .get(half_edge.loop_)
                .and_then(|loop_data| self.faces.get(loop_data.face))
                .is_some_and(|face| face.role == TopologyRole::ConstructionResidual)
    }

    /// Check 4: For every vertex with a valid outgoing half-edge, that half-edge's
    /// origin must be the vertex. A vertex whose stored outgoing handle is
    /// non-null but MISSING from the arena is also inconsistent: it
    /// points at a REMOVED half-edge, and the dangling anchor silently
    /// exempts the vertex from every outgoing-fan walk (measured on
    /// `BRep::torus`: the seam stitch left one vertex dangling and the
    /// shell vertex recompute skipped it, serving a stale position; all
    /// removal sites now route through `remove_half_edge_reanchored`,
    /// enforced structurally). A NULL (never-set) anchor stays exempt —
    /// numerous construction paths wire half-edge origins directly and
    /// anchor lazily; arming that class is its own audit (filed).
    pub fn validate_vertex_consistency(&self) -> Result<(), ValidationError> {
        for (v_handle, v_data) in &self.vertices {
            let outgoing = v_data.outgoing;
            if !self.half_edges.contains_key(outgoing) {
                if outgoing != crate::topology::HalfEdgeHandle::default() {
                    return Err(ValidationError::InconsistentVertex {
                        vertex: v_handle,
                        half_edge: outgoing,
                    });
                }
                continue;
            }
            let he_data = &self.half_edges[outgoing];
            if he_data.origin != v_handle {
                return Err(ValidationError::InconsistentVertex {
                    vertex: v_handle,
                    half_edge: outgoing,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    fn make_cube() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();
        brep
    }

    #[test]
    fn cube_passes_half_edge_pairing() {
        let brep = make_cube();
        assert!(brep.validate_half_edge_pairing().is_ok());
    }

    #[test]
    fn cube_passes_half_edge_cycles() {
        let brep = make_cube();
        assert!(brep.validate_half_edge_cycles().is_ok());
    }

    #[test]
    fn cube_passes_vertex_consistency() {
        let brep = make_cube();
        assert!(brep.validate_vertex_consistency().is_ok());
    }

    #[test]
    fn corrupted_twin_fails_pairing() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt the twin pointer.
        brep.half_edges[mev.half_edge].twin = Default::default();
        assert!(brep.validate_half_edge_pairing().is_err());
    }

    #[test]
    fn corrupted_next_fails_cycle() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Make next point to itself instead of twin — cycle never returns.
        brep.half_edges[mev.half_edge].next = mev.half_edge;
        // This creates a cycle of length 1 that doesn't include the twin.
        // The twin's cycle will never reach itself in bounded steps.
        assert!(brep.validate_half_edge_cycles().is_err());
    }

    #[test]
    fn corrupted_outgoing_fails_vertex_consistency() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Point the origin vertex's outgoing to the twin (which originates from new vertex).
        brep.vertices[mvfs.vertex].outgoing = mev.half_edge_twin;
        assert!(brep.validate_vertex_consistency().is_err());
    }

    #[test]
    fn broken_cycle_next_points_to_invalid() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Make next point to a non-existent handle
        brep.half_edges[mev.half_edge].next = Default::default();
        let result = brep.validate_half_edge_cycles();
        assert!(result.is_err());
    }

    #[test]
    fn twin_twin_mismatch_fails_pairing() {
        // Create a valid structure, then corrupt twin's twin to NOT point back
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        // Make mev1's twin's twin point to something else (mev2's half_edge)
        brep.half_edges[mev1.half_edge_twin].twin = mev2.half_edge;
        let result = brep.validate_half_edge_pairing();
        assert!(result.is_err());
    }
}
