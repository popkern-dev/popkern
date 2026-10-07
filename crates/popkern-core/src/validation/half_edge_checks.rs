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
    ///
    /// Each cycle is walked once: a walk that returns to its start proves
    /// every half-edge on it, and those are not walked again (PK-709 — a
    /// walk from every half-edge cost O(Σ loop length²), and the Euler
    /// operators run this check after every step, which made extruding a
    /// 1000-edge profile cubic). `next` is a function, so a walk that
    /// meets a half-edge already proven to lie on another cycle can never
    /// return to its start: that start is reported, as the exhaustive
    /// walk would report it once past the bound. The first failing
    /// half-edge in arena order is the one reported, as before.
    pub fn validate_half_edge_cycles(&self) -> Result<(), ValidationError> {
        let bound = self.half_edges.len();
        let mut on_cycle: std::collections::HashSet<crate::topology::HalfEdgeHandle> =
            std::collections::HashSet::with_capacity(bound);
        let mut members = Vec::new();

        for (he_handle, he_data) in &self.half_edges {
            if self.half_edge_is_construction_residual_for_cycle_validation(he_data)
                || on_cycle.contains(&he_handle)
            {
                continue;
            }
            let broken = || ValidationError::BrokenCycle {
                half_edge: he_handle,
                loop_: he_data.loop_,
            };
            members.clear();
            members.push(he_handle);
            let mut cursor = he_data.next;
            let mut steps = 0;

            while cursor != he_handle {
                steps += 1;
                if steps > bound || on_cycle.contains(&cursor) {
                    return Err(broken());
                }
                match self.half_edges.get(cursor) {
                    Some(next_he) => {
                        if self.half_edge_is_construction_residual_for_cycle_validation(next_he) {
                            return Err(broken());
                        }
                        members.push(cursor);
                        cursor = next_he.next;
                    }
                    None => return Err(broken()),
                }
            }
            on_cycle.extend(members.iter().copied());
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
    /// MISSING from the arena is also inconsistent: a non-null handle
    /// points at a REMOVED half-edge (all removal sites route through
    /// `remove_half_edge_reanchored`, enforced structurally), and a
    /// NULL (never-set) anchor on a vertex that half-edges originate at
    /// means a construction path wired origins without anchoring
    /// (PK-337) — both silently exempt the vertex from every
    /// outgoing-fan walk (measured on `BRep::torus`: one dangling
    /// anchor and the shell vertex recompute served a stale position).
    /// Only genuinely isolated vertices (no incident half-edges at
    /// all) are exempt.
    pub fn validate_vertex_consistency(&self) -> Result<(), ValidationError> {
        let mut has_incident: std::collections::BTreeSet<crate::topology::VertexHandle> =
            std::collections::BTreeSet::new();
        for (_, he) in &self.half_edges {
            has_incident.insert(he.origin);
        }
        for (v_handle, v_data) in &self.vertices {
            let outgoing = v_data.outgoing;
            if !self.half_edges.contains_key(outgoing) {
                if outgoing != crate::topology::HalfEdgeHandle::default()
                    || has_incident.contains(&v_handle)
                {
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

    /// PK-324/PK-337 detection verification: both dangling-anchor
    /// classes fire — a stale non-null handle (points at a removed
    /// half-edge) and a never-set null handle on a vertex that live
    /// half-edges originate at.
    #[test]
    fn dangling_and_never_set_anchors_fail_vertex_consistency() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        assert!(brep.validate_vertex_consistency().is_ok());

        // Never-set: null anchor while a half-edge originates here.
        let saved = brep.vertices[mvfs.vertex].outgoing;
        brep.vertices[mvfs.vertex].outgoing = Default::default();
        assert!(brep.validate_vertex_consistency().is_err());
        brep.vertices[mvfs.vertex].outgoing = saved;

        // Stale: non-null anchor pointing at a removed half-edge.
        let clone = brep.half_edges[mev.half_edge].clone();
        let spare_he = brep.half_edges.insert(clone);
        brep.remove_half_edge_reanchored(spare_he);
        brep.vertices[mvfs.vertex].outgoing = spare_he;
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

    /// The exhaustive walk the one-walk-per-cycle check replaced: the
    /// first half-edge in arena order whose `next` walk does not return
    /// within the bound (residual half-edges aside — the cube has none).
    fn first_broken_exhaustively(brep: &BRep) -> Option<crate::topology::HalfEdgeHandle> {
        let bound = brep.half_edges.len();
        brep.half_edges.iter().find_map(|(h, data)| {
            let mut cursor = data.next;
            let mut steps = 0;
            while cursor != h {
                steps += 1;
                let next = brep.half_edges.get(cursor)?;
                if steps > bound {
                    return Some(h);
                }
                cursor = next.next;
            }
            None
        })
    }

    #[test]
    fn a_walk_leading_into_another_cycle_reports_the_same_half_edge_as_the_exhaustive_walk() {
        use crate::error::ValidationError;
        let cube = make_cube();
        assert!(cube.validate_half_edge_cycles().is_ok());
        let handles: Vec<_> = cube.half_edges.keys().collect();
        // Redirect each half-edge in turn into a cycle of another loop
        // (a rho shape: its walk never comes back) and compare.
        for &h in &handles {
            let mut brep = cube.clone();
            let own_loop = brep.half_edges[h].loop_;
            let Some(other) = handles
                .iter()
                .copied()
                .find(|&o| brep.half_edges[o].loop_ != own_loop)
            else {
                continue;
            };
            brep.half_edges[h].next = other;
            let expected = first_broken_exhaustively(&brep);
            match brep.validate_half_edge_cycles() {
                Err(ValidationError::BrokenCycle { half_edge, .. }) => {
                    assert_eq!(Some(half_edge), expected, "redirected {h:?}");
                }
                other => panic!("redirected {h:?}: expected a broken cycle, got {other:?}"),
            }
        }
    }

    #[test]
    fn a_dangling_next_is_a_broken_cycle() {
        use crate::error::ValidationError;
        let mut brep = make_cube();
        let h = brep.half_edges.keys().next().expect("a half-edge");
        brep.half_edges[h].next = crate::topology::HalfEdgeHandle::default();
        assert!(matches!(
            brep.validate_half_edge_cycles(),
            Err(ValidationError::BrokenCycle { .. })
        ));
    }
}
