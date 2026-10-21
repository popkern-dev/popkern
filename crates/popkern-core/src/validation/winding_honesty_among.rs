//! Region form of the winding-honesty audit and the orientation-flag
//! reconcile (PK-815, docs/architecture/local-boolean-scope.md "the
//! touched region"): a boolean that carried most of its body unchanged
//! re-derives the flags of the faces it touched without auditing the
//! rest.

use std::collections::VecDeque;

use slotmap::SecondaryMap;

use super::winding_honesty::FaceWindingHonesty;
use crate::topology::{BRep, FaceHandle};

impl BRep {
    /// [`BRep::audit_winding_honesty`]'s verdicts for `faces` only, the
    /// same ones the whole-body audit gives them. The second tier
    /// resolves an indeterminate face from its outer-loop neighbours and
    /// chains to a fixpoint, so the audit also measures, transitively,
    /// every neighbour of an indeterminate face it reaches: the
    /// determinate ones anchor the chain, the indeterminate ones extend
    /// it. Faces outside that closure cannot change a verdict here.
    /// Unknown, dead, duplicate or surfaceless handles are skipped.
    pub fn audit_winding_honesty_among(&self, faces: &[FaceHandle]) -> Vec<FaceWindingHonesty> {
        let mut wanted: SecondaryMap<FaceHandle, ()> = SecondaryMap::new();
        let mut seen: SecondaryMap<FaceHandle, ()> = SecondaryMap::new();
        let mut queue: VecDeque<FaceHandle> = VecDeque::new();
        for &face in faces {
            if self.faces.contains_key(face) && seen.insert(face, ()).is_none() {
                wanted.insert(face, ());
                queue.push_back(face);
            }
        }
        let mut verdicts = Vec::new();
        while let Some(face) = queue.pop_front() {
            let Some(verdict) = self.face_winding_tier1(face) else {
                continue;
            };
            if verdict.honest.is_none() {
                for neighbour in self.outer_loop_neighbours(face) {
                    if seen.insert(neighbour, ()).is_none() {
                        queue.push_back(neighbour);
                    }
                }
            }
            verdicts.push(verdict);
        }
        self.resolve_by_tangency_witness(&mut verdicts);
        verdicts.retain(|v| wanted.contains_key(v.face));
        verdicts
    }

    /// [`BRep::reconcile_orientation_flags`] over `faces` only: their
    /// flags are re-derived from the audit (rounds as there); every
    /// other face keeps its flag. Returns the number of flags flipped.
    pub fn reconcile_orientation_flags_among(&mut self, faces: &[FaceHandle]) -> usize {
        let mut flipped = 0;
        for _ in 0..4 {
            let dishonest: Vec<FaceHandle> = self
                .audit_winding_honesty_among(faces)
                .into_iter()
                .filter(|v| v.honest == Some(false))
                .map(|v| v.face)
                .collect();
            if dishonest.is_empty() {
                break;
            }
            for face in dishonest {
                if let Some(sof) = self.faces.get_mut(face).and_then(|f| f.surface.as_mut()) {
                    sof.orientation = !sof.orientation;
                    flipped += 1;
                }
            }
        }
        flipped
    }

    /// The faces across `face`'s outer-loop edges (the second tier's
    /// witness candidates), itself excluded.
    fn outer_loop_neighbours(&self, face: FaceHandle) -> Vec<FaceHandle> {
        let mut out = Vec::new();
        let Some(start) = self
            .faces
            .get(face)
            .and_then(|f| self.loops.get(f.outer_loop))
            .map(|l| l.half_edge)
        else {
            return out;
        };
        let mut cur = start;
        for _ in 0..self.half_edges.len() + 1 {
            let Some(he) = self.half_edges.get(cur) else {
                break;
            };
            if let Some(n) = self
                .half_edges
                .get(he.twin)
                .and_then(|t| self.loops.get(t.loop_))
                .map(|l| l.face)
                && n != face
            {
                out.push(n);
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "winding_honesty_among_tests.rs"]
mod tests;
