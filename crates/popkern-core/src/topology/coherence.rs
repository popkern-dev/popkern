//! Edge-curve / face-boundary / pcurve coherence invariants.
//!
//! Stage A1 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
//! coherence-invariants.md`. Ships the data types and the
//! `BRep::record_pending_authoring` API for tracking transient L1/L3
//! violations across operator chains. The actual `validate_coherence`
//! method lands in Stage A2.
//!
//! The three invariants tracked here:
//!
//! - **L1** — every edge in a steady-state BRep has a curve.
//! - **L2** — both incident face surfaces interpolate the edge
//!   curve at the half-edges' pcurve images within tolerance.
//! - **L3** — every half-edge has a pcurve mapping the edge curve
//!   into its face's parameter space within tolerance.
//!
//! Transient violations (mid-surgery state) are recorded via
//! [`PendingCoherenceMarker`] so the operator chain that introduced
//! the violation can resolve it before returning to a steady-state
//! boundary.

use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Transient marker carrying an obligation to author a curve or
/// pcurve before a BRep returns to a steady-state boundary. Recorded
/// by an operator that legitimately cannot author at the moment it
/// creates the entity (e.g. because the canonical curve must be
/// derived from a face surface that has not yet been attached).
///
/// NOT serialized: markers are operator-local state, MUST be empty
/// when a BRep is serialized at a steady-state boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PendingCoherenceMarker {
    /// Edge created without a curve. Owner is a `&'static str`
    /// identifying the operator (e.g. `"split_edge"`) so the
    /// downstream consumer can attribute the obligation.
    EdgeCurvePending {
        /// Edge missing its curve representation.
        edge: EdgeHandle,
        /// Operator label identifying who created the edge.
        owner_operator: &'static str,
    },
    /// Half-edge created without a pcurve on a face that has a
    /// surface (pcurve is the missing piece, not the surface).
    HalfEdgePcurvePending {
        /// Half-edge missing its pcurve.
        half_edge: HalfEdgeHandle,
        /// Operator label identifying who created the half-edge.
        owner_operator: &'static str,
    },
}

/// Scope of a `validate_coherence` call.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Stage A2 wires validate_coherence
pub enum CoherenceValidationScope {
    /// Check only entities flagged with a [`PendingCoherenceMarker`].
    /// Used as a focused gate at operator boundaries.
    PendingOnly,
    /// Check every edge / half-edge / adjacent-face pair in the
    /// BRep. Used by validation passes and end-of-pipeline gates.
    AllSteadyState,
    /// Restrict the check to a known face set (e.g. cap faces for
    /// the cap/support boundary solve).
    Faces(std::collections::BTreeSet<FaceHandle>),
}

/// Per-edge L2 disagreement record. Emitted by L2 sampling when
/// two adjacent face surfaces evaluate to different 3D points at
/// the same edge parameter beyond `coherence_tolerance`.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)] // Stage A3 populates these via L2 sampling
pub struct L2Disagreement {
    /// Edge shared by the two disagreeing faces.
    pub edge: EdgeHandle,
    /// First face (the one whose half-edge sits on `edge.half_edge`).
    pub face_a: FaceHandle,
    /// Adjacent face (the twin half-edge's face).
    pub face_b: FaceHandle,
    /// Worst sampled distance between the two faces' surface
    /// boundary evaluations along the edge.
    pub worst_gap: f64,
    /// Parameter `t ∈ [0, 1]` of the half-edge on `face_a` where
    /// the worst gap was sampled.
    pub worst_parameter: f64,
}

/// Report returned by `BRep::validate_coherence`. Empty fields mean
/// the corresponding invariant holds at the requested scope.
#[derive(Debug, Default, Clone, PartialEq)]
#[allow(dead_code)] // Stage A2/A3 populate these
pub struct CoherenceReport {
    /// L1 violations: edges with no curve and no pending marker.
    pub l1_missing_edge_curves: Vec<EdgeHandle>,
    /// L2 violations: face pairs whose surfaces disagree at their
    /// shared edge beyond `coherence_tolerance`.
    pub l2_face_boundary_disagreements: Vec<L2Disagreement>,
    /// L3 violations: half-edges with no pcurve on a face whose
    /// surface is present, and no pending marker.
    pub l3_missing_pcurves: Vec<HalfEdgeHandle>,
}

impl CoherenceReport {
    /// True if every tracked invariant held at the requested scope.
    #[allow(dead_code)] // first caller lands in Stage A2 tests
    pub fn is_empty(&self) -> bool {
        self.l1_missing_edge_curves.is_empty()
            && self.l2_face_boundary_disagreements.is_empty()
            && self.l3_missing_pcurves.is_empty()
    }
}

impl BRep {
    /// Record a transient L1/L3 obligation. Idempotent — recording
    /// the same marker twice has no effect. The owning operator (or
    /// a downstream consumer) MUST clear the marker before the
    /// BRep returns to a steady-state boundary.
    ///
    /// The marker is stored in [`BRep::pending_coherence`] for
    /// later consumption by `validate_coherence(PendingOnly)`
    /// (Stage A2).
    pub fn record_pending_authoring(&mut self, marker: PendingCoherenceMarker) {
        if !self.pending_coherence.contains(&marker) {
            self.pending_coherence.push(marker);
        }
    }

    /// Clear all pending markers for an edge or half-edge. Called
    /// when the obligation has been satisfied (curve / pcurve
    /// authored).
    #[allow(dead_code)] // first caller lands in Stage C
    pub fn clear_pending_authoring_for_edge(&mut self, edge: EdgeHandle) {
        self.pending_coherence.retain(|m| {
            !matches!(
                m,
                PendingCoherenceMarker::EdgeCurvePending { edge: e, .. } if *e == edge
            )
        });
    }

    /// Clear all pending markers for a half-edge.
    #[allow(dead_code)] // first caller lands in Stage C
    pub fn clear_pending_authoring_for_half_edge(&mut self, half_edge: HalfEdgeHandle) {
        self.pending_coherence.retain(|m| !matches!(
            m,
            PendingCoherenceMarker::HalfEdgePcurvePending { half_edge: h, .. } if *h == half_edge
        ));
    }
}

#[cfg(test)]
#[path = "coherence_tests.rs"]
mod tests;
