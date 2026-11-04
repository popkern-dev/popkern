//! Vertex `outgoing`-anchor repair for half-edge removal surgery.
//!
//! `Vertex::outgoing` is DERIVED navigational state: any live incident
//! half-edge is a valid anchor, and every surgical site that removes
//! half-edges must re-anchor the affected vertices or their fan walks
//! silently skip them forever (measured: `BRep::torus`'s seam stitch
//! left one vertex dangling and the shell vertex recompute served a
//! stale position; the armed `validate_vertex_consistency` check now
//! rejects the class). `kve`/`mef` carry their own local repairs; the
//! remaining removal sites route through this shared helper.

use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

impl BRep {
    /// Re-anchor every vertex whose `outgoing` points at one of the
    /// half-edges in `doomed`, choosing a surviving half-edge with the
    /// same origin (deterministic: iteration order of the slot map is
    /// slot order, stable for a given body). Call BEFORE removing the
    /// half-edges — or use [`BRep::remove_half_edge_reanchored`], which
    /// composes the repair with the removal. A vertex with no surviving
    /// incident half-edge keeps its stale anchor — it is topologically
    /// isolated and exempt from the consistency check.
    pub fn reanchor_vertex_outgoing(&mut self, doomed: &[HalfEdgeHandle]) {
        for &h in doomed {
            let Some(he) = self.half_edges.get(h) else {
                continue;
            };
            let v = he.origin;
            let Some(vertex) = self.vertices.get(v) else {
                continue;
            };
            if !doomed.contains(&vertex.outgoing) {
                continue;
            }
            if let Some((replacement, _)) = self
                .half_edges
                .iter()
                .find(|(hh, cand)| cand.origin == v && !doomed.contains(hh))
            {
                self.vertices[v].outgoing = replacement;
            }
        }
    }

    /// Remove a half-edge, first re-anchoring its origin vertex if its
    /// `outgoing` pointed at it. This is the ONLY sanctioned way to
    /// remove half-edges outside the anchor-aware Euler operators —
    /// enforced by `structural_half_edge_removal.rs`. Returns the
    /// removed half-edge like `SlotMap::remove`.
    pub fn remove_half_edge_reanchored(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Option<crate::topology::HalfEdge> {
        self.reanchor_vertex_outgoing(&[half_edge]);
        self.half_edges.remove(half_edge)
    }

    /// Remove a face whose loops and half-edges are already gone, and
    /// keep its shell's anchors true. `Shell::face` and `Solid::shell`
    /// are anchors like `Vertex::outgoing`: a shell anchored on the
    /// removed face moves to a surviving face of the same shell; a
    /// shell left with no face is removed (it bounds nothing — the
    /// `validate_all` reference check rejects it, PK-966); a solid
    /// anchored on that shell moves to another of its shells, or is
    /// removed when it has none. Returns whether the shell went too.
    pub fn remove_face_reanchored(&mut self, face: FaceHandle) -> bool {
        let Some(removed) = self.faces.remove(face) else {
            return false;
        };
        let shell = removed.shell;
        let Some(shell_data) = self.shells.get(shell) else {
            return false;
        };
        let solid = shell_data.solid;
        if let Some((survivor, _)) = self.faces.iter().find(|(_, f)| f.shell == shell) {
            if self.shells[shell].face == face {
                self.shells[shell].face = survivor;
            }
            return false;
        }
        self.shells.remove(shell);
        if self.solids.get(solid).is_some_and(|s| s.shell == shell) {
            match self.shells.iter().find(|(_, s)| s.solid == solid) {
                Some((other, _)) => self.solids[solid].shell = other,
                None => {
                    self.solids.remove(solid);
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn removing_the_anchor_face_moves_the_shell_anchor() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let shell = brep.shells.keys().next().expect("a shell");
        let anchor = brep.shells[shell].face;
        assert!(!brep.remove_face_reanchored(anchor));
        assert_ne!(brep.shells[shell].face, anchor);
        assert!(brep.faces.contains_key(brep.shells[shell].face));
    }

    #[test]
    fn removing_the_last_face_removes_the_shell_and_its_solid() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces: Vec<_> = brep.faces.keys().collect();
        let gone: Vec<bool> = faces
            .iter()
            .map(|&f| brep.remove_face_reanchored(f))
            .collect();
        assert_eq!(gone.iter().filter(|&&g| g).count(), 1);
        assert_eq!(gone.last(), Some(&true));
        assert!(brep.shells.is_empty());
        assert!(brep.solids.is_empty());
    }
}
