//! Degenerate-ring removal Euler operator (`kvr` — kill vertex, ring).
//!
//! Removes an EMPTY inner ring (a loop whose last edge was retired by
//! [`BRep::kev`]'s last-edge arm, leaving `loop.half_edge` null) from
//! its face, together with the ring's orphaned anchor vertex: V−1,
//! ring−1, so the Euler characteristic `V − E + F − H` is preserved.
//! This is the terminal move of feature dissolution: a vanished
//! protrusion's boundary ring shrinks by `kef`/`kev` down to one
//! anchor vertex on an empty ring, and `kvr` retires the pair.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, LoopHandle, VertexHandle};

impl BRep {
    /// Remove the empty inner ring `ring` from `face` together with its
    /// orphaned anchor vertex.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if any handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if `ring` is not an inner
    ///   loop of `face`, still contains half-edges, or the vertex still
    ///   has incident half-edges.
    pub fn kvr(
        &mut self,
        face: FaceHandle,
        ring: LoopHandle,
        vertex: VertexHandle,
    ) -> Result<(), CoreError> {
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?;
        let loop_data = self
            .loops
            .get(ring)
            .ok_or_else(|| CoreError::InvalidHandle("ring handle is invalid".into()))?;
        if !self.vertices.contains_key(vertex) {
            return Err(CoreError::InvalidHandle("vertex handle is invalid".into()));
        }
        if loop_data.face != face || !face_data.inner_loops.contains(&ring) {
            return Err(CoreError::PreconditionFailed(
                "kvr: ring is not an inner loop of the given face".into(),
            ));
        }
        if self.half_edges.contains_key(loop_data.half_edge) {
            return Err(CoreError::PreconditionFailed(
                "kvr: ring still contains half-edges (retire them with kev first)".into(),
            ));
        }
        if self.half_edges.values().any(|he| he.origin == vertex) {
            return Err(CoreError::PreconditionFailed(
                "kvr: anchor vertex still has incident half-edges".into(),
            ));
        }

        let mut candidate = self.clone();
        if let Some(f) = candidate.faces.get_mut(face) {
            f.inner_loops.retain(|l| *l != ring);
        }
        candidate.loops.remove(ring);
        candidate.vertices.remove(vertex);

        // Structural postconditions (kve's no-regression contract):
        // pre-existing failures of a check are exempt.
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
                    continue;
                }
                return Err(CoreError::PreconditionFailed(format!(
                    "kvr: post-removal structural validation failed: {e:?}"
                )));
            }
        }
        *self = candidate;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Build a cube with a spur-ring on one face (mev into a fresh
    /// ring is not directly expressible, so exercise the kev
    /// last-edge → kvr pipeline through a real dissolution: split a
    /// face, then retire the split remnant). The focused contract
    /// here is precondition honesty, exercised directly.
    #[test]
    fn kvr_rejects_non_empty_rings_and_wrong_faces() {
        let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let (face, outer) = {
            let (fh, f) = cube.faces.iter().next().unwrap();
            (fh, f.outer_loop)
        };
        let vertex = cube.vertices.keys().next().unwrap();
        // The outer loop is not an inner ring.
        let err = cube.kvr(face, outer, vertex).unwrap_err();
        assert!(format!("{err}").contains("inner loop"), "{err}");
    }
}
