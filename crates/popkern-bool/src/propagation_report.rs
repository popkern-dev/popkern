//! Structured record of cellular verdict-propagation declines.
//!
//! The cellular classifier (PK-496) classifies a component of
//! operand-A fragments from one seed ray plus flip parity, checked by
//! audit rays. A component that declines falls back to one ray per
//! fragment. The result is the same body, but the cost is quadratic on
//! multi-wrap walls. A decline therefore changes nothing a validator
//! sees, and before this registry it was visible only under
//! `POPKERN_TRACE=cellular`: an inconclusive audit declined the whole
//! wall of an 80 mm thread for two days, taking its merge stage from
//! 1.5 s to 30 s, and only the release perf gate noticed (PK-786).
//!
//! The registry is process-global and drained by
//! [`take_propagation_declines`], like `merge_fallback`.

use std::sync::Mutex;

/// Why a propagation component declined to the per-fragment path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropagationDeclineReason {
    /// The seed's ray verdict was not plain Inside or Outside.
    SeedNotPlain,
    /// An audit's decisive verdict contradicted the seed under flip
    /// parity.
    AuditDisagrees,
    /// Audits and reserves ran out before enough of them were
    /// decisive: every candidate sample sat on operand B's boundary.
    AuditsInconclusive,
}

/// One component that declined verdict propagation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PropagationDecline {
    /// Fragments in the component, all now classified one by one.
    pub component_size: usize,
    /// Why it declined.
    pub reason: PropagationDeclineReason,
}

static DECLINES: Mutex<Vec<PropagationDecline>> = Mutex::new(Vec::new());

/// Record a decline: append to the process registry and say so once
/// on stderr. A decline is a cost cliff, so it is never silent.
pub(crate) fn record_propagation_decline(component_size: usize, reason: PropagationDeclineReason) {
    eprintln!(
        "[boolean] verdict propagation declined for a component of {component_size} fragments: {reason:?}"
    );
    if let Ok(mut declines) = DECLINES.lock() {
        declines.push(PropagationDecline {
            component_size,
            reason,
        });
    }
}

/// Drain and return every decline recorded by this process since the
/// last call. Empty when every planned component propagated.
pub fn take_propagation_declines() -> Vec<PropagationDecline> {
    DECLINES
        .lock()
        .map(|mut declines| std::mem::take(&mut *declines))
        .unwrap_or_default()
}
