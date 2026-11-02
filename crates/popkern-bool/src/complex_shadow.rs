//! Shadow comparison of the cell complex against the cellular
//! classifier (PK-976, plan Stage 7a).
//!
//! With `POPKERN_COMPLEX_SHADOW=1` every cellular boolean also
//! classifies its fragments through the cell complex and records how
//! the two verdicts compare. Nothing served changes. The registry is
//! process-global and drained by [`take_complex_shadow_events`]; the
//! suite-wide measurement (walk completeness, agreement, ray counts)
//! reads it.

use std::sync::Mutex;

use crate::error::BoolOperation;

/// How the complex verdict compared with the cellular one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComplexShadowOutcome {
    /// Identical keep/flip for every fragment.
    Agree,
    /// The verdicts differ on these fragments.
    Disagree {
        /// Fragment indices whose keep or flip differ, each with how.
        fragments: Vec<(usize, ShadowDifference)>,
    },
    /// The complex declined; the decline's debug form.
    Declined(String),
}

/// How the complex's verdict on one fragment differs from the cellular
/// classifier's. Only [`ShadowDifference::Membership`] says the two
/// would serve different point sets; the others are the complex
/// deciding by construction what the cellular path leaves to the
/// assembler and post-processing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadowDifference {
    /// Both keep the fragment and differ on its flip. The complex
    /// orients every kept face by its component's folded-normal
    /// convention; the cellular assembler re-wires a ring it was handed
    /// the wrong way reversed (measured: every such difference in the
    /// workspace suite was a cellular ring re-wired, the complex's
    /// none).
    Orientation,
    /// A coincident patch both keep, by a different copy: the complex
    /// keeps the copy whose loops bound the whole patch, the cellular
    /// classifier the workpiece's, which its own seam may split.
    CoincidentCopy,
    /// The fragment fills a hole of its own face. The complex never
    /// lifts it; the cellular classifier kept it, and post-processing
    /// removes it with its self-twinned rim.
    HoleFiller,
    /// The verdicts disagree on whether the fragment bounds the result.
    Membership,
}

/// One shadow comparison.
#[derive(Debug, Clone)]
pub struct ComplexShadowEvent {
    /// The operation.
    pub operation: BoolOperation,
    /// Agreement.
    pub outcome: ComplexShadowOutcome,
    /// Fragment count.
    pub fragments: usize,
    /// Point classifications the complex needed (the cellular path
    /// casts at least one per sampled fragment).
    pub classifications: usize,
    /// Non-manifold edges and vertices the complex found in the
    /// result of its own selection.
    pub nonmanifold: (usize, usize),
}

static EVENTS: Mutex<Vec<ComplexShadowEvent>> = Mutex::new(Vec::new());

/// What the complex did where the cellular classifier could not serve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComplexRescueOutcome {
    /// It served the boolean.
    Served,
    /// Its selected boundary is non-manifold: the boolean declines
    /// unless graph-merge produces a body `validate_all` accepts.
    NonManifold,
    /// It declined (the decline's display form); graph-merge follows.
    Declined(String),
    /// Its body failed the cellular acceptance; graph-merge follows.
    Invalid(String),
    /// Its split of a non-manifold boundary served, but the pieces
    /// failed `validate_all` with their contacts recorded; graph-merge's
    /// coincident copies follow.
    SplitRejected(String),
}

/// One complex rescue attempt.
#[derive(Debug, Clone)]
pub struct ComplexRescueEvent {
    /// The operation.
    pub operation: BoolOperation,
    /// What happened.
    pub outcome: ComplexRescueOutcome,
}

static RESCUES: Mutex<Vec<ComplexRescueEvent>> = Mutex::new(Vec::new());

/// Record a rescue attempt (and say it once on stderr: a rescue changes
/// which classifier served the result).
pub(crate) fn record_complex_rescue(operation: BoolOperation, outcome: ComplexRescueOutcome) {
    eprintln!("[boolean] complex rescue ({operation:?}): {outcome:?}");
    if let Ok(mut events) = RESCUES.lock() {
        events.push(ComplexRescueEvent { operation, outcome });
    }
}

/// Record that the complex, serving as the primary classifier, declined
/// and the cellular classifier decided instead (plan Stage 7c; its
/// counter must read 0 before the cellular classifier retires).
pub(crate) fn record_complex_primary_fallback(
    operation: BoolOperation,
    diag: &crate::error::BoolDiagnostics,
) {
    let reason = diag
        .error
        .as_ref()
        .map_or_else(|| "no error".to_string(), |e| e.to_string());
    if popkern_core::trace::trace_enabled("complex") {
        eprintln!("[complex] primary declined ({operation:?}), cellular decides: {reason}");
    }
    append_log_line(&format!("PrimaryFallback {operation:?} {reason}"));
}

/// Drain every rescue attempt recorded since the last call.
pub fn take_complex_rescue_events() -> Vec<ComplexRescueEvent> {
    RESCUES
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}

/// Whether shadow comparison is on for this process.
pub(crate) fn shadow_enabled() -> bool {
    std::env::var_os("POPKERN_COMPLEX_SHADOW").is_some_and(|v| v == "1")
}

/// The label a log line carries: the test thread's name, or, on a
/// thread the test spawned (a bigger stack, a worker pool), the test
/// the process was started for: nextest runs each test in its own
/// process, naming it as the first non-flag argument.
fn log_label() -> String {
    if let Some(name) = std::thread::current().name()
        && name != "main"
    {
        return name.to_string();
    }
    std::env::args()
        .skip(1)
        .find(|a| !a.starts_with('-'))
        .unwrap_or_else(|| "?".to_string())
}

/// Append `line` to `POPKERN_COMPLEX_SHADOW_LOG`, when set, labelled
/// with the test it ran under, so a whole test-suite run (one process
/// per test) can be measured at once.
pub(crate) fn append_log_line(line: &str) {
    let Some(path) = std::env::var_os("POPKERN_COMPLEX_SHADOW_LOG") else {
        return;
    };
    use std::io::Write;
    if let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    {
        let _ = writeln!(f, "{}\t{line}", log_label());
    }
}

/// Record one comparison. With `POPKERN_COMPLEX_SHADOW_LOG=<path>` the
/// event is also appended to that file as one line, so a whole
/// test-suite run (one process per test) can be measured at once.
pub(crate) fn record_complex_shadow(event: ComplexShadowEvent) {
    if popkern_core::trace::trace_enabled("complex") {
        eprintln!("[complex] shadow {event:?}");
    }
    append_log_line(&format!("{event:?}"));
    if let Ok(mut events) = EVENTS.lock() {
        events.push(event);
    }
}

/// Drain every shadow comparison recorded since the last call.
pub fn take_complex_shadow_events() -> Vec<ComplexShadowEvent> {
    EVENTS
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}
