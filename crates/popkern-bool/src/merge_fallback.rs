//! Structured record of the cellular→graph-merge fallback boundary.
//!
//! The boolean pipeline's primary selection is the order-independent
//! cellular decomposition; the legacy graph-merge selection is its
//! explicit fallback. Every time the fallback fires, the pipeline
//! records WHY (structured, process-queryable) and says so once on
//! stderr — a silent route change here previously cost a whole
//! campaign of misdirected pcurve archaeology, because the cellular
//! error was discarded with `.ok()` on all builds.
//!
//! The registry is process-global and drained by
//! [`take_merge_fallback_events`]; tests use it to pin the fallback
//! fire-rate over a corpus.

use std::sync::Mutex;

use crate::error::BoolOperation;

/// Why the boolean pipeline fell back from cellular selection to the
/// legacy graph-merge selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeFallbackReason {
    /// `POPKERN_CELLULAR=0` forced graph-merge (diagnostic
    /// kill-switch; never set in production).
    CellularDisabled,
    /// The cellular assembly returned an error before producing a
    /// body. Carries the error's display form.
    CellularError(String),
    /// The cellular assembly produced a body that failed validation.
    /// Carries the validation error's display form.
    CellularInvalid(String),
}

/// One firing of the cellular→graph-merge fallback.
#[derive(Debug, Clone)]
pub struct MergeFallbackEvent {
    /// The boolean operation being executed when the fallback fired.
    pub operation: BoolOperation,
    /// Why cellular selection did not serve.
    pub reason: MergeFallbackReason,
}

static EVENTS: Mutex<Vec<MergeFallbackEvent>> = Mutex::new(Vec::new());

/// Record a fallback firing: append to the process registry and emit
/// one unconditional stderr line. The stderr line is deliberate — a
/// fallback at this boundary changes which 10k-line assembly serves
/// the result, and that must never be invisible.
pub(crate) fn record_merge_fallback(operation: BoolOperation, reason: MergeFallbackReason) {
    eprintln!("[boolean] cellular->graph-merge fallback ({operation:?}): {reason:?}");
    if let Ok(mut events) = EVENTS.lock() {
        events.push(MergeFallbackEvent { operation, reason });
    }
}

/// Drain and return every fallback event recorded by this process
/// since the last call. Empty when the cellular path served every
/// boolean — the expected steady state.
pub fn take_merge_fallback_events() -> Vec<MergeFallbackEvent> {
    EVENTS
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}
