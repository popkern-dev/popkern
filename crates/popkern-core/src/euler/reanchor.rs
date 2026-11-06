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

use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, ShellHandle};

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
        self.unlink_fin(half_edge);
        self.half_edges.remove(half_edge)
    }

    /// Remove a face whose loops and half-edges are already gone, and
    /// keep its shell's anchors true. `Shell::face` and `Region::shell`
    /// are anchors like `Vertex::outgoing`: a shell anchored on the
    /// removed face moves to a surviving face of the same shell; a
    /// shell left with no face is removed (it bounds nothing — the
    /// `validate_all` reference check rejects it, PK-966); a region
    /// anchored on that shell moves to another of its shells, or is
    /// removed when it has none. Returns whether the shell went too.
    pub fn remove_face_reanchored(&mut self, face: FaceHandle) -> bool {
        let Some(removed) = self.faces.remove(face) else {
            return false;
        };
        let back_gone = self.release_side(removed.back, face);
        if removed.front != removed.back {
            self.release_side(removed.front, face);
        }
        back_gone
    }

    /// After `face` went, keep the shell one of its sides lay in true:
    /// re-anchor it on a face still in it, or remove it when none is,
    /// re-pointing or removing its region (never the outside) if the
    /// shell was the region's peripheral one. Returns whether it went.
    fn release_side(&mut self, shell: ShellHandle, face: FaceHandle) -> bool {
        let Some(shell_data) = self.shells.get(shell) else {
            return false;
        };
        let region = shell_data.region;
        if let Some(&survivor) = self.shell_faces(shell).first() {
            if self.shells[shell].face == face {
                self.set_shell_anchor(shell, survivor);
            }
            return false;
        }
        // A shell still holding wire edges or an acorn outlives its last
        // face side.
        if !self.shells[shell].wire_edges.is_empty() || self.shells[shell].acorn.is_some() {
            self.shells[shell].face = FaceHandle::default();
            return false;
        }
        self.shells.remove(shell);
        if region != self.outside && self.regions.get(region).is_some_and(|r| r.shell == shell) {
            match self.shells.iter().find(|(_, s)| s.region == region) {
                Some((other, _)) => self.regions[region].shell = other,
                None => {
                    self.regions.remove(region);
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
        assert!(brep.solid_regions().is_empty());
        assert_eq!(brep.regions.keys().collect::<Vec<_>>(), vec![brep.outside]);
    }
}
