//! Free-boundary mirror loops of sheet faces (PK-978).
//!
//! A sheet's face holds both half-edges of each boundary edge: the
//! outer loop, and the same edges reversed as an inner loop. That inner
//! loop is no hole — it traces the boundary the outer loop already
//! bounds — and every containment and hole test must leave it out, or
//! even-odd containment calls the whole face outside. A sheet's real
//! hole carries a mirror too (its rim's reversed copy, also an inner
//! loop): exactly one loop of each mirrored pair counts.

use crate::topology::{BRep, FaceHandle, LoopHandle};

impl BRep {
    /// Whether inner loop `ring` of `face` is a free-boundary mirror:
    /// every half-edge's twin lies in ONE other loop of `face`, and that
    /// partner is the face's outer loop or an inner loop with a lower
    /// handle (so one loop of each mirrored pair of inner loops counts).
    pub fn is_mirror_loop(&self, face: FaceHandle, ring: LoopHandle) -> bool {
        (|| -> Option<bool> {
            let face_data = self.faces.get(face)?;
            if ring == face_data.outer_loop {
                return Some(false);
            }
            let start = self.loops.get(ring)?.half_edge;
            let mut partner: Option<LoopHandle> = None;
            let mut cursor = start;
            for _ in 0..=self.half_edges.len() {
                let he = self.half_edges.get(cursor)?;
                let twin_loop = he.mate().and_then(|m| self.half_edges.get(m))?.loop_;
                if twin_loop == ring || self.loops.get(twin_loop)?.face != face {
                    return Some(false);
                }
                match partner {
                    None => partner = Some(twin_loop),
                    Some(p) if p == twin_loop => {}
                    Some(_) => return Some(false),
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
            let partner = partner?;
            Some(partner == face_data.outer_loop || partner < ring)
        })()
        .unwrap_or(false)
    }
}

#[cfg(test)]
#[path = "brep_mirror_loops_tests.rs"]
mod tests;
