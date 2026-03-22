//! Manifold edge and vertex valence checks.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Check 5: Every edge's half-edge must have a twin, and that twin must
    /// reference the same edge.
    pub fn validate_manifold_edges(&self) -> Result<(), ValidationError> {
        for (edge_handle, edge_data) in &self.edges {
            let he = self
                .half_edges
                .get(edge_data.half_edge)
                .ok_or(ValidationError::NonManifoldEdge { edge: edge_handle })?;

            let twin = self
                .half_edges
                .get(he.twin)
                .ok_or(ValidationError::NonManifoldEdge { edge: edge_handle })?;

            if twin.edge != edge_handle {
                return Err(ValidationError::NonManifoldEdge { edge: edge_handle });
            }
        }
        Ok(())
    }

    /// Check 6: Every vertex must have at least 2 edges incident to it.
    ///
    /// Vertices that were just created by `mvfs` and don't have edges yet
    /// (their outgoing handle is the default/null handle) are skipped.
    pub fn validate_vertex_valence(&self) -> Result<(), ValidationError> {
        for (v_handle, v_data) in &self.vertices {
            // Skip vertices without any outgoing half-edge (mvfs-only vertices).
            if !self.half_edges.contains_key(v_data.outgoing) {
                continue;
            }

            // Count distinct edges incident to this vertex by walking outgoing
            // half-edges. We walk the "fan" around the vertex:
            // start at outgoing, follow twin.next to get the next outgoing he.
            let start = v_data.outgoing;
            let mut count = 0usize;
            let mut cursor = start;
            let bound = self.half_edges.len();
            let mut steps = 0;

            loop {
                count += 1;
                steps += 1;
                if steps > bound {
                    // Safety: avoid infinite loop on corrupted data.
                    break;
                }
                let he = &self.half_edges[cursor];
                let twin = he.twin;
                match self.half_edges.get(twin) {
                    Some(twin_data) => {
                        cursor = twin_data.next;
                        if cursor == start {
                            break;
                        }
                    }
                    None => break,
                }
            }

            if count < 2 {
                return Err(ValidationError::LowValence {
                    vertex: v_handle,
                    valence: count,
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
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0)
            .unwrap();
        brep
    }

    #[test]
    fn cube_passes_manifold_edges() {
        let brep = make_cube();
        assert!(brep.validate_manifold_edges().is_ok());
    }

    #[test]
    fn cube_passes_vertex_valence() {
        let brep = make_cube();
        assert!(brep.validate_vertex_valence().is_ok());
    }

    #[test]
    fn mvfs_only_skips_vertex_valence() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // mvfs vertex has no outgoing half-edge, so it's skipped
        assert!(brep.validate_vertex_valence().is_ok());
    }

    #[test]
    fn single_edge_fails_vertex_valence() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Each vertex has valence 1 — should fail.
        assert!(brep.validate_vertex_valence().is_err());
    }

    #[test]
    fn corrupted_edge_he_fails_manifold() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt the edge's half_edge pointer.
        brep.edges[mev.edge].half_edge = Default::default();
        assert!(brep.validate_manifold_edges().is_err());
    }

    #[test]
    fn corrupted_twin_edge_fails_manifold() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt the twin's edge pointer to not match.
        brep.half_edges[mev.half_edge_twin].edge = Default::default();
        assert!(brep.validate_manifold_edges().is_err());
    }

    #[test]
    fn vertex_valence_with_corrupted_twin() {
        // Create a cube then corrupt a twin to trigger the safety break in validate_vertex_valence
        let mut brep = make_cube();
        // Find a vertex and corrupt its fan by removing a twin
        let he_handle = brep.half_edges.keys().next().unwrap();
        let twin = brep.half_edges[he_handle].twin;
        // Make twin point to invalid handle so the fan walk breaks early
        brep.half_edges[he_handle].twin = Default::default();
        // Should not panic; may return Ok or Err depending on resulting valence
        let _ = brep.validate_vertex_valence();
    }
}
