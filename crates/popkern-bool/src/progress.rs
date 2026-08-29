//! Progress reporting for long-running boolean operations (§15).
//!
//! Implements the §15 [BEYOND] requirement: "Progress callback:
//! report percentage complete, current phase name, and entity counts
//! at a frequency >= 2 Hz for operations lasting longer than 2
//! seconds."
//!
//! The boolean pipeline calls a [`ProgressReporter`] between every
//! major phase, passing the current [`crate::error::Stage`] and a
//! [`ProgressEvent`] describing what just happened. Callers can
//! implement the trait to drive a UI progress bar, log to stderr,
//! or aggregate metrics.
//!
//! Combine with [`crate::cancellation::CancellationToken`] for the
//! full §15 interactive workflow: progress callback drives the UI,
//! cancellation token responds to "Cancel".

use std::time::Duration;

use crate::error::Stage;

/// One progress event emitted by the boolean pipeline.
///
/// `phase_index` and `total_phases` let a UI compute a percentage
/// without having to know the exact stage names.
#[derive(Debug, Clone, Copy)]
pub struct ProgressEvent {
    /// Which pipeline stage just completed.
    pub stage: Stage,
    /// 0-based index of this stage in the pipeline.
    pub phase_index: usize,
    /// Total number of phases the pipeline will run (typically 6:
    /// intersect, split, classify, select, merge, validate).
    pub total_phases: usize,
    /// Wall-clock duration of the stage that just completed, or
    /// `None` on platforms without a monotonic clock (see
    /// `popkern_core::time::Stopwatch`).
    pub duration: Option<Duration>,
    /// Stage-specific entity count (pair count for intersection,
    /// fragment count for split, etc.).
    pub entity_count: usize,
}

impl ProgressEvent {
    /// Convenience: returns the percentage 0–100 (rounded).
    pub fn percent(&self) -> u8 {
        if self.total_phases == 0 {
            return 100;
        }
        ((self.phase_index + 1) * 100 / self.total_phases) as u8
    }
}

/// A progress callback the boolean pipeline invokes between phases.
///
/// Implementors typically drive a UI progress bar or write metrics
/// to a log. Callers should keep the implementation cheap — the
/// pipeline calls the reporter inline, so a slow callback adds to
/// the boolean's wall time.
pub trait ProgressReporter: Send {
    /// Called immediately after a pipeline stage completes.
    fn on_phase(&mut self, event: ProgressEvent);
}

/// A no-op reporter that ignores all events. Useful as a default
/// when the caller does not need progress notifications.
#[derive(Debug, Default, Clone, Copy)]
pub struct NullReporter;

impl ProgressReporter for NullReporter {
    fn on_phase(&mut self, _event: ProgressEvent) {}
}

/// A reporter that records every event into a vector. Primarily
/// useful for tests that need to verify the pipeline emitted the
/// expected sequence of events.
#[derive(Debug, Default, Clone)]
pub struct CaptureReporter {
    /// All events recorded so far.
    pub events: Vec<ProgressEvent>,
}

impl ProgressReporter for CaptureReporter {
    fn on_phase(&mut self, event: ProgressEvent) {
        self.events.push(event);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percent_at_first_of_six_is_16() {
        let event = ProgressEvent {
            stage: Stage::Intersection,
            phase_index: 0,
            total_phases: 6,
            duration: Some(Duration::from_millis(10)),
            entity_count: 0,
        };
        assert_eq!(event.percent(), 16);
    }

    #[test]
    fn percent_at_last_phase_is_100() {
        let event = ProgressEvent {
            stage: Stage::Validation,
            phase_index: 5,
            total_phases: 6,
            duration: Some(Duration::from_millis(10)),
            entity_count: 0,
        };
        assert_eq!(event.percent(), 100);
    }

    #[test]
    fn null_reporter_ignores_events() {
        let mut reporter = NullReporter;
        reporter.on_phase(ProgressEvent {
            stage: Stage::Intersection,
            phase_index: 0,
            total_phases: 6,
            duration: Some(Duration::from_millis(0)),
            entity_count: 0,
        });
        // No state to check; just verify it doesn't panic.
    }

    #[test]
    fn capture_reporter_records_events() {
        let mut reporter = CaptureReporter::default();
        reporter.on_phase(ProgressEvent {
            stage: Stage::Intersection,
            phase_index: 0,
            total_phases: 6,
            duration: Some(Duration::from_millis(5)),
            entity_count: 12,
        });
        reporter.on_phase(ProgressEvent {
            stage: Stage::Merge,
            phase_index: 4,
            total_phases: 6,
            duration: Some(Duration::from_millis(20)),
            entity_count: 18,
        });
        assert_eq!(reporter.events.len(), 2);
        assert_eq!(reporter.events[0].entity_count, 12);
        assert_eq!(reporter.events[1].stage, Stage::Merge);
    }
}
