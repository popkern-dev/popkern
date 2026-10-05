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
use popkern_math::PortableMath;

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
        brep.remove_face_to_boundary(face_to_remove)?;
        Ok(brep)
    }

    /// Build an open n-gonal prism — [`BRep::prism`] with one
    /// n-gonal cap removed, leaving an `n`-sided hole. The N ≥ 5
    /// fixture for general-N hole healing (PK-248): `open_box`
    /// covers the 4-sided class, this covers everything above it.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if no n-sided cap face
    ///   exists (degenerate constructor input).
    /// - [`CoreError::Validation`] if the opened shell fails
    ///   validation.
    pub fn open_prism(
        n: usize,
        radius: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        let profile: Vec<crate::geometry::Point3> = (0..n)
            .map(|k| {
                let angle = std::f64::consts::TAU * (k as f64) / (n as f64);
                crate::geometry::Point3::new(radius * angle.pcos(), radius * angle.psin(), 0.0)
            })
            .collect();
        let mut brep = Self::prism(
            &profile,
            crate::geometry::Vec3::new(0.0, 0.0, 1.0),
            height,
            tolerance,
        )?;
        // Remove one n-gonal cap: the face whose outer loop has n
        // half-edges and whose vertices share one z (side walls are
        // quads; both caps qualify — take the first).
        let cap = brep
            .faces
            .keys()
            .find(|&face| {
                let loop_h = brep.faces[face].outer_loop;
                let start = brep.loops[loop_h].half_edge;
                let mut count = 0usize;
                let mut cursor = start;
                loop {
                    count += 1;
                    cursor = brep.half_edges[cursor].next;
                    if cursor == start || count > n + 1 {
                        break;
                    }
                }
                count == n
            })
            .ok_or_else(|| {
                CoreError::PreconditionFailed(format!("prism({n}) has no {n}-sided cap face"))
            })?;
        brep.remove_face_to_boundary(cap)?;
        Ok(brep)
    }

    /// Remove `face_to_remove`, registering its former edges as
    /// boundary edges (null-sentinel twins) — the shared open-shell
    /// surgery behind [`BRep::open_box`] and [`BRep::open_prism`].
    fn remove_face_to_boundary(
        &mut self,
        face_to_remove: crate::topology::FaceHandle,
    ) -> Result<(), CoreError> {
        let brep = self;
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
            // twin (it may currently point at the doomed one) — and
            // keep the canonical-forward convention true: the twin
            // traverses the edge curve BACKWARD, so the stored curve
            // range reverses with the repoint (PK-274: without this,
            // the forward-oriented pcurve trace measure pairs the
            // ends crosswise and its direction-free fallback bottoms
            // at half its sampling spacing — the four rim edges
            // measured permanently "stale" at exactly 3.125e-3 and
            // every downstream fix_stale pass re-cleared and
            // reverted them forever).
            if brep.edges[edge_h].half_edge == heh
                && let Some(curve) = brep.edges[edge_h].curve.as_mut()
            {
                std::mem::swap(&mut curve.t_start, &mut curve.t_end);
            }
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
            brep.remove_half_edge_reanchored(heh);
        }

        // Remove the loop and face.
        brep.loops.remove(loop_h);
        brep.faces.remove(face_to_remove);
        if let Some((replacement_face, _)) = brep.faces.iter().next() {
            let shell = brep.faces[replacement_face].shell;
            brep.shells[shell].face = replacement_face;
        }

        brep.validate().map_err(CoreError::Validation)?;
        Ok(())
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

    /// PK-274 pin: the rim repoint keeps the canonical-forward curve
    /// convention (the stored curve reverses with the canonical), so
    /// the builder ships zero stale-measuring pcurves. Regression
    /// mode: the four rim edges measure exactly half the trace
    /// sampling spacing (3.125e-3 on the unit box) forever, and
    /// every downstream fix_stale pass re-clears and reverts them.
    #[test]
    fn open_box_ships_no_stale_pcurves() {
        let ob = BRep::open_box(Tolerance::default()).unwrap();
        assert_eq!(
            ob.stale_pcurve_handles(),
            vec![],
            "open_box rim repoint must preserve the canonical-forward curve convention"
        );
    }

    #[test]
    fn open_prism_ships_no_stale_pcurves() {
        let prism = BRep::open_prism(5, 1.0, 2.0, Tolerance::default()).unwrap();
        assert_eq!(prism.stale_pcurve_handles(), vec![]);
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

#[cfg(test)]
mod open_prism_tests {
    use super::*;

    #[test]
    fn open_prism_has_n_boundary_edges_and_validates() {
        for n in [5usize, 6] {
            let prism = BRep::open_prism(n, 1.0, 1.0, Tolerance::default()).unwrap();
            assert_eq!(prism.boundary_edges.len(), n, "n={n}");
            assert_eq!(prism.faces.len(), n + 1, "n={n}: walls + one cap");
            prism.validate().unwrap();
        }
    }
}
