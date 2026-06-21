//! Coherence validation (L1+L3 checks).
//!
//! Stage A2 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
//! coherence-invariants.md`. Lives in a sibling file to `brep.rs` so
//! adding L1/L3 (and L2 in Stage A3) does not push the parent over
//! the 500-line cap.
//!
//! L1 and L3 are presence-only checks — they walk the BRep's edge
//! and half-edge slot maps and collect entities that lack their
//! coherence-required representation (curve / pcurve) AND have no
//! [`PendingCoherenceMarker`] excusing the absence. L2 sampling
//! lands in Stage A3.

use crate::topology::coherence::{
    CoherenceReport, CoherenceValidationScope, PendingCoherenceMarker,
};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

impl BRep {
    /// Run the coherence invariants L1+L3 (and L2 once Stage A3
    /// lands) against the requested scope.
    ///
    /// Returns an empty [`CoherenceReport`] when every tracked
    /// invariant holds. Populated fields name the violating
    /// entities. The BRep is not mutated.
    ///
    /// Stage A2 covers L1 (missing edge curves) and L3 (missing
    /// pcurves on half-edges whose face has a surface). Stage A3
    /// extends to L2 (adjacent-face surface disagreement at the
    /// shared edge).
    pub fn validate_coherence(&self, scope: &CoherenceValidationScope) -> CoherenceReport {
        let mut report = CoherenceReport::default();
        match scope {
            CoherenceValidationScope::PendingOnly => self.check_pending_only(&mut report),
            CoherenceValidationScope::AllSteadyState => self.check_all_steady_state(&mut report),
            CoherenceValidationScope::Faces(face_set) => {
                self.check_face_scope(face_set, &mut report)
            }
        }
        report
    }

    fn check_pending_only(&self, report: &mut CoherenceReport) {
        for marker in &self.pending_coherence {
            match marker {
                PendingCoherenceMarker::EdgeCurvePending { edge, .. } => {
                    if self.edge_lacks_curve(*edge) {
                        report.l1_missing_edge_curves.push(*edge);
                    }
                }
                PendingCoherenceMarker::HalfEdgePcurvePending { half_edge, .. } => {
                    if self.half_edge_lacks_pcurve_on_surfaced_face(*half_edge) {
                        report.l3_missing_pcurves.push(*half_edge);
                    }
                }
            }
        }
    }

    fn check_all_steady_state(&self, report: &mut CoherenceReport) {
        for (edge_handle, edge) in &self.edges {
            if edge.curve.is_none() && !self.has_pending_curve_marker(edge_handle) {
                report.l1_missing_edge_curves.push(edge_handle);
            }
        }
        for (he_handle, _half_edge) in &self.half_edges {
            if self.half_edge_lacks_pcurve_on_surfaced_face(he_handle)
                && !self.has_pending_pcurve_marker(he_handle)
            {
                report.l3_missing_pcurves.push(he_handle);
            }
        }
    }

    fn check_face_scope(
        &self,
        face_set: &std::collections::BTreeSet<crate::topology::FaceHandle>,
        report: &mut CoherenceReport,
    ) {
        for (he_handle, he) in &self.half_edges {
            let face = self
                .loops
                .get(he.loop_)
                .map(|l| l.face)
                .filter(|f| face_set.contains(f));
            if face.is_none() {
                continue;
            }
            let edge_handle = he.edge;
            if let Some(edge) = self.edges.get(edge_handle)
                && edge.curve.is_none()
                && !self.has_pending_curve_marker(edge_handle)
                && !report.l1_missing_edge_curves.contains(&edge_handle)
            {
                report.l1_missing_edge_curves.push(edge_handle);
            }
            if self.half_edge_lacks_pcurve_on_surfaced_face(he_handle)
                && !self.has_pending_pcurve_marker(he_handle)
            {
                report.l3_missing_pcurves.push(he_handle);
            }
        }
    }

    fn edge_lacks_curve(&self, edge: EdgeHandle) -> bool {
        self.edges.get(edge).is_some_and(|e| e.curve.is_none())
    }

    fn half_edge_lacks_pcurve_on_surfaced_face(&self, half_edge: HalfEdgeHandle) -> bool {
        let Some(he) = self.half_edges.get(half_edge) else {
            return false;
        };
        if he.pcurve.is_some() {
            return false;
        }
        let face = match self.loops.get(he.loop_) {
            Some(l) => l.face,
            None => return false,
        };
        self.faces.get(face).is_some_and(|f| f.surface.is_some())
    }

    fn has_pending_curve_marker(&self, edge: EdgeHandle) -> bool {
        self.pending_coherence.iter().any(|m| {
            matches!(
                m,
                PendingCoherenceMarker::EdgeCurvePending { edge: e, .. } if *e == edge
            )
        })
    }

    fn has_pending_pcurve_marker(&self, half_edge: HalfEdgeHandle) -> bool {
        self.pending_coherence.iter().any(|m| {
            matches!(
                m,
                PendingCoherenceMarker::HalfEdgePcurvePending { half_edge: h, .. } if *h == half_edge
            )
        })
    }
}

#[cfg(test)]
#[path = "validate_coherence_tests.rs"]
mod tests;
