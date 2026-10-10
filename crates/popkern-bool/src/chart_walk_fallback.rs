//! Structured record of the chart walk → 3D-Newell walk fallback.
//!
//! Each face of a cellular boolean is walked on its parameter chart
//! (`merge::walk::chart_subdivision`); when the chart substrate cannot
//! be built (a seam balance it cannot close, an orbit collision) the
//! face walks on the 3D-Newell walker instead. The reasons used to
//! surface only when the boolean then FAILED; a served boolean dropped
//! them, so a fallback that happened to serve correctly was invisible
//! and its removal untestable (PK-750: a rod wall whose seam a cut line
//! took over walked on the fallback with every gate green).
//!
//! The registry is process-global and drained by
//! [`take_chart_walk_fallback_events`], like
//! [`crate::merge_fallback`].

use std::sync::Mutex;

/// One face that walked on the 3D-Newell fallback.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChartWalkFallbackEvent {
    /// Whether the face belongs to the first operand (else the second).
    pub operand_a: bool,
    /// The face, in its operand's BRep.
    pub face: popkern_core::topology::FaceHandle,
    /// Why the chart substrate declined.
    pub reason: String,
}

static EVENTS: Mutex<Vec<ChartWalkFallbackEvent>> = Mutex::new(Vec::new());

/// Record one face's chart-walk fallback.
pub(crate) fn record_chart_walk_fallback(face: &crate::merge::graph::FaceRef, reason: &str) {
    if let Ok(mut events) = EVENTS.lock() {
        events.push(ChartWalkFallbackEvent {
            operand_a: face.operand == crate::merge::graph::Operand::A,
            face: face.face,
            reason: reason.to_string(),
        });
    }
}

/// Drain and return every chart-walk fallback recorded by this process
/// since the last call.
pub fn take_chart_walk_fallback_events() -> Vec<ChartWalkFallbackEvent> {
    EVENTS
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}
