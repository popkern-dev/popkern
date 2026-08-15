//! Region-scoped pcurve recompute (PK-274).
//!
//! Nested inside `pcurve_recompute` (shares its module context, like
//! `pcurve_stale_repair`): the scoped variants reuse the parent's
//! private per-face recompute and periodic-winding close.
//!
//! Whole-model recompute at operator tails is priced by model size —
//! 23 Newton point-inversions per pcurve × every half-edge of every
//! face, including the ~90% of faces the operation never touched. An
//! operation only invalidates pcurves in the neighbourhood it
//! rewired, so tails pass their dirty faces dilated by one ring (see
//! `BRep::validation_region`); the whole-model forms remain for
//! serving boundaries and the debug cross-check
//! (`BRep::stale_pcurve_handles`).

use crate::topology::{BRep, FaceHandle};

impl BRep {
    /// As [`BRep::recompute_all_pcurves`], restricted to `faces`
    /// (deduplicated; unknown or dead handles ignored). Strict abort
    /// semantics: the first face whose projection fails abandons the
    /// remaining faces, exactly like the whole-model form.
    pub fn recompute_pcurves_for(
        &mut self,
        faces: &[FaceHandle],
    ) -> Result<(), crate::error::CoreError> {
        let mut seen = std::collections::BTreeSet::new();
        for &fh in faces {
            if !self.faces.contains_key(fh) || !seen.insert(fh) {
                continue;
            }
            self.recompute_face_pcurves(fh)?;
        }
        Ok(())
    }

    /// Recompute pcurves for every face with at least one half-edge
    /// MISSING its pcurve (strict abort semantics, like
    /// [`BRep::recompute_all_pcurves`]). The boolean tails use this:
    /// their whole-model recompute existed to fill the pcurves the
    /// merge seams and collinear-vertex dissolutions left empty, and
    /// a face with a complete pcurve inventory was not touched by
    /// either (PK-274).
    pub fn recompute_pcurves_where_missing(&mut self) -> Result<(), crate::error::CoreError> {
        let mut faces = std::collections::BTreeSet::new();
        for half_edge in self.half_edges.values() {
            if half_edge.pcurve.is_none()
                && let Some(loop_data) = self.loops.get(half_edge.loop_)
            {
                faces.insert(loop_data.face);
            }
        }
        let faces: Vec<FaceHandle> = faces.into_iter().collect();
        self.recompute_pcurves_for(&faces)
    }

    /// As [`BRep::recompute_all_pcurves_lenient`], restricted to
    /// `faces`: one face's projection failure does not abandon the
    /// faces after it.
    pub fn recompute_pcurves_for_lenient(&mut self, faces: &[FaceHandle]) {
        let mut seen = std::collections::BTreeSet::new();
        for &fh in faces {
            if !self.faces.contains_key(fh) || !seen.insert(fh) {
                continue;
            }
            let _ = self.recompute_face_pcurves(fh);
        }
    }
}

#[cfg(test)]
#[path = "pcurve_region_recompute_tests.rs"]
mod tests;
