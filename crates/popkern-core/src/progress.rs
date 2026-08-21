//! Op-agnostic progress reporting for long-running kernel operations.
//!
//! Design: `docs/architecture/interactive-operation-control.md`
//! (PK-276). One concrete event shape serves shell, sweep,
//! tessellation, and smooth-mesh progress, so a UI subscribes to a
//! single type across operations, and round-based pipelines can report
//! honest per-run stage totals (which a static phase enum cannot).
//!
//! Boundary note: `popkern-fillet` and `popkern-bool` keep their
//! shipped, requirement-numbered progress APIs (fillet's carries §18
//! per-edge telemetry; bool's is typed on its error `Stage`), and
//! `popkern_core::healing::ProgressCallback` keeps its polling `&self`
//! shape. What is unified across every operation is the
//! [`crate::cancellation::CancellationToken`].

use std::time::Duration;

use crate::cancellation::{CancellationToken, EventSink};

/// One progress event emitted by an operation pipeline, after the
/// named stage completes.
#[derive(Debug, Clone, Copy)]
pub struct ProgressEvent {
    /// Stable semantic label of the stage that just completed, e.g.
    /// `"shell.offset_surfaces"` or `"tessellate.encroachment_round"`.
    /// Labels are diagnostic names — stable across releases, unique
    /// within an operation.
    pub stage: &'static str,
    /// 0-based position of this stage in the run.
    pub stage_index: usize,
    /// Honest per-run stage count (dynamic rounds included). An
    /// early-exiting pipeline may emit fewer events than announced,
    /// never more.
    pub total_stages: usize,
    /// Wall-clock duration of the stage that just completed.
    pub duration: Duration,
    /// Stage-specific entity count (faces offset, rounds inserted,
    /// surfaces checked, …).
    pub entity_count: usize,
}

impl ProgressEvent {
    /// Percentage complete `0..=100` (integer, rounded down), computed
    /// from the stage position. Returns 100 for a degenerate zero
    /// total.
    pub fn percent(&self) -> u8 {
        if self.total_stages == 0 {
            return 100;
        }
        (((self.stage_index + 1) * 100 / self.total_stages).min(100)) as u8
    }
}

/// A progress callback an operation pipeline invokes after each stage.
///
/// Implementors typically drive a UI progress bar or write metrics to
/// a log. Keep implementations cheap — the pipeline calls the reporter
/// inline (or drains it on the calling thread of a detached run), so a
/// slow callback adds to perceived wall time.
pub trait ProgressReporter: Send {
    /// Called immediately after a pipeline stage completes.
    fn on_stage(&mut self, event: ProgressEvent);
}

/// A no-op reporter that ignores all events.
#[derive(Debug, Default, Clone, Copy)]
pub struct NullReporter;

impl ProgressReporter for NullReporter {
    fn on_stage(&mut self, _event: ProgressEvent) {}
}

/// A reporter that records every event, for tests and telemetry.
#[derive(Debug, Default, Clone)]
pub struct CaptureReporter {
    /// All events recorded so far, in emission order.
    pub events: Vec<ProgressEvent>,
}

impl ProgressReporter for CaptureReporter {
    fn on_stage(&mut self, event: ProgressEvent) {
        self.events.push(event);
    }
}

/// The detached-worker side of progress delivery: a pipeline running
/// under [`crate::cancellation::run_detach_on_cancel_with_events`] is
/// handed an `EventSink<ProgressEvent>` and can use it directly as a
/// reporter; the calling thread drains the events into the caller's
/// borrowed reporter.
impl ProgressReporter for EventSink<ProgressEvent> {
    fn on_stage(&mut self, event: ProgressEvent) {
        self.emit(event);
    }
}

/// Caller-side control bundle for a long-running operation: optional
/// cooperative cancellation and optional progress reporting, carried
/// as one parameter so operation signatures gain exactly one argument.
///
/// `OperationControl::default()` requests neither — a
/// `*_with_control` call with the default control behaves identically
/// to the plain entry point.
#[derive(Default)]
pub struct OperationControl<'a> {
    /// Cooperative cancellation token. When present, the operation
    /// runs detached (worker thread over private clones) so the
    /// < 1 s termination claim is structural.
    pub cancel: Option<&'a CancellationToken>,
    /// Per-stage progress reporter.
    pub reporter: Option<&'a mut dyn ProgressReporter>,
}

impl<'a> OperationControl<'a> {
    /// Control with cancellation only.
    pub fn with_cancel(cancel: &'a CancellationToken) -> Self {
        Self {
            cancel: Some(cancel),
            reporter: None,
        }
    }

    /// Control with progress reporting only.
    pub fn with_reporter(reporter: &'a mut dyn ProgressReporter) -> Self {
        Self {
            cancel: None,
            reporter: Some(reporter),
        }
    }

    /// Control with both cancellation and progress reporting.
    pub fn new(cancel: &'a CancellationToken, reporter: &'a mut dyn ProgressReporter) -> Self {
        Self {
            cancel: Some(cancel),
            reporter: Some(reporter),
        }
    }
}

impl std::fmt::Debug for OperationControl<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OperationControl")
            .field("cancel", &self.cancel.is_some())
            .field("reporter", &self.reporter.is_some())
            .finish()
    }
}

/// Threads cancellation polls and per-stage progress events through a
/// linear operation pipeline: `cancelled()` before each stage,
/// `finish_stage(label, started, count)` after it. Keeps the stage
/// cursor so pipelines never hand-number their events.
///
/// The reporter's trait-object pointee lifetime (`'r`) is independent
/// of the borrow (`'a`) so callers can reuse their reporter after the
/// tracker drops.
pub struct StageTracker<'a, 'r> {
    cancel: Option<&'a CancellationToken>,
    reporter: Option<&'a mut (dyn ProgressReporter + 'r)>,
    stage_index: usize,
    total_stages: usize,
}

impl<'a, 'r> StageTracker<'a, 'r> {
    /// New tracker over an optional token + reporter. `total_stages`
    /// is the honest per-run stage count (dynamic rounds included);
    /// an early-exiting pipeline may emit fewer events, never more.
    pub fn new(
        cancel: Option<&'a CancellationToken>,
        reporter: Option<&'a mut (dyn ProgressReporter + 'r)>,
        total_stages: usize,
    ) -> Self {
        Self {
            cancel,
            reporter,
            stage_index: 0,
            total_stages,
        }
    }

    /// A tracker that polls nothing and reports nowhere — for plain
    /// entry points delegating into a control-aware pipeline.
    pub fn inert() -> Self {
        Self::new(None, None, 0)
    }

    /// Cooperative poll: `true` once any holder of the token has
    /// cancelled. Pipelines call this at stage boundaries (and inside
    /// dominant loops) and return their `Cancelled` error.
    pub fn cancelled(&self) -> bool {
        self.cancel.is_some_and(CancellationToken::is_cancelled)
    }

    /// The tracked token, for pipelines that must poll inside code the
    /// tracker cannot reach (e.g. rayon-parallel per-item closures).
    pub fn token(&self) -> Option<&'a CancellationToken> {
        self.cancel
    }

    /// Emit the completed stage's event and advance the cursor.
    /// `started` is the instant the stage began.
    pub fn finish_stage(
        &mut self,
        stage: &'static str,
        started: std::time::Instant,
        entity_count: usize,
    ) {
        let event = ProgressEvent {
            stage,
            stage_index: self.stage_index,
            total_stages: self.total_stages,
            duration: started.elapsed(),
            entity_count,
        };
        self.stage_index += 1;
        if let Some(reporter) = self.reporter.as_mut() {
            reporter.on_stage(event);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn event(stage: &'static str, index: usize, total: usize) -> ProgressEvent {
        ProgressEvent {
            stage,
            stage_index: index,
            total_stages: total,
            duration: Duration::from_millis(1),
            entity_count: 3,
        }
    }

    #[test]
    fn percent_spans_zero_to_hundred() {
        assert_eq!(event("a", 0, 8).percent(), 12);
        assert_eq!(event("a", 7, 8).percent(), 100);
        assert_eq!(event("a", 0, 0).percent(), 100);
    }

    #[test]
    fn capture_reporter_records_in_order() {
        let mut reporter = CaptureReporter::default();
        reporter.on_stage(event("first", 0, 2));
        reporter.on_stage(event("second", 1, 2));
        assert_eq!(reporter.events.len(), 2);
        assert_eq!(reporter.events[0].stage, "first");
        assert_eq!(reporter.events[1].stage, "second");
    }

    #[test]
    fn null_reporter_ignores_events() {
        NullReporter.on_stage(event("ignored", 0, 1));
    }

    #[test]
    fn event_sink_reporter_delivers_through_detach_layer() {
        let token = CancellationToken::new();
        let mut received = CaptureReporter::default();
        let outcome = crate::cancellation::run_detach_on_cancel_with_events(
            &token,
            |mut sink| {
                sink.on_stage(ProgressEvent {
                    stage: "detached.stage",
                    stage_index: 0,
                    total_stages: 1,
                    duration: Duration::from_millis(2),
                    entity_count: 9,
                });
                "ok"
            },
            |event| received.on_stage(event),
        );
        assert!(matches!(
            outcome,
            crate::cancellation::DetachOutcome::Completed("ok")
        ));
        assert_eq!(received.events.len(), 1);
        assert_eq!(received.events[0].stage, "detached.stage");
        assert_eq!(received.events[0].entity_count, 9);
    }

    #[test]
    fn stage_tracker_polls_and_numbers_stages() {
        let token = CancellationToken::new();
        let mut reporter = CaptureReporter::default();
        {
            let mut tracker = StageTracker::new(Some(&token), Some(&mut reporter), 3);
            assert!(!tracker.cancelled());
            let t0 = std::time::Instant::now();
            tracker.finish_stage("op.first", t0, 5);
            tracker.finish_stage("op.second", t0, 6);
            token.cancel();
            assert!(tracker.cancelled());
        }
        assert_eq!(reporter.events.len(), 2);
        assert_eq!(reporter.events[0].stage, "op.first");
        assert_eq!(reporter.events[0].stage_index, 0);
        assert_eq!(reporter.events[1].stage_index, 1);
        assert_eq!(reporter.events[1].total_stages, 3);
        assert!(!StageTracker::inert().cancelled());
    }

    #[test]
    fn operation_control_constructors() {
        let token = CancellationToken::new();
        let control = OperationControl::with_cancel(&token);
        assert!(control.cancel.is_some());
        assert!(control.reporter.is_none());

        let mut reporter = CaptureReporter::default();
        let control = OperationControl::with_reporter(&mut reporter);
        assert!(control.cancel.is_none());
        assert!(control.reporter.is_some());

        let mut reporter = CaptureReporter::default();
        let control = OperationControl::new(&token, &mut reporter);
        assert!(control.cancel.is_some());
        assert!(control.reporter.is_some());

        let control = OperationControl::default();
        assert!(control.cancel.is_none());
        assert!(control.reporter.is_none());
        let printed = format!("{control:?}");
        assert!(printed.contains("OperationControl"));
    }
}
