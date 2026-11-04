//! Structured record of the cellular stage's declines.
//!
//! The boolean pipeline's selection is the order-independent cellular
//! decomposition (the complex primary, the cellular classifier its
//! fallback). Graph-merge selection once caught its declines; plan
//! Stage 7d retired it, and a decline now goes to the complex rescue
//! or declines the boolean. Every time the cellular stage declines,
//! the pipeline records WHY (structured, process-queryable) and says
//! so once on stderr — a silent route change here previously cost a
//! whole campaign of misdirected pcurve archaeology, because the
//! cellular error was discarded with `.ok()` on all builds.
//!
//! The registry is process-global and drained by
//! [`take_merge_fallback_events`]; tests use it to pin the fallback
//! fire-rate over a corpus.

use std::sync::Mutex;

use crate::error::BoolOperation;

/// Why the cellular stage declined.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeFallbackReason {
    /// The cellular assembly returned an error before producing a
    /// body. Carries the error's display form.
    CellularError(String),
    /// The cellular assembly produced a body that failed validation.
    /// Carries the validation error's display form.
    CellularInvalid(String),
}

/// One cellular decline.
#[derive(Debug, Clone)]
pub struct MergeFallbackEvent {
    /// The boolean operation being executed when the fallback fired.
    pub operation: BoolOperation,
    /// Why cellular selection did not serve.
    pub reason: MergeFallbackReason,
}

static EVENTS: Mutex<Vec<MergeFallbackEvent>> = Mutex::new(Vec::new());

/// Record a decline: append to the process registry and emit one
/// unconditional stderr line. The stderr line is deliberate — a
/// decline here changes what serves the result (the complex rescue,
/// or nothing), and that must never be invisible.
pub(crate) fn record_merge_fallback(operation: BoolOperation, reason: MergeFallbackReason) {
    eprintln!("[boolean] cellular declined ({operation:?}): {reason:?}");
    crate::complex_shadow::append_log_line(&format!("MergeFallback {operation:?} {reason:?}"));
    if let Ok(mut events) = EVENTS.lock() {
        events.push(MergeFallbackEvent { operation, reason });
    }
}

/// Drain and return every decline recorded by this process since the
/// last call. Empty when the cellular path served every boolean — the
/// expected steady state.
pub fn take_merge_fallback_events() -> Vec<MergeFallbackEvent> {
    EVENTS
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}
