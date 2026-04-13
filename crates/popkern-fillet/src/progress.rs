//! Progress reporting for long-running fillet operations (§17).
//!
//! The fillet pipeline calls a [`ProgressReporter`] between every
//! major phase, passing a [`ProgressEvent`] describing what just
//! happened. Callers implement the trait to drive a UI progress
//! bar, log to stderr, or aggregate metrics.
//!
//! Combine with [`crate::cancellation::CancellationToken`] for the
//! full §17 interactive workflow: progress callback drives the UI,
//! cancellation token responds to "Cancel".

use std::time::Duration;

/// Named phases in the fillet pipeline — one value per stage that
/// emits a progress event. The full pipeline emits events in this
/// order: `Validate` → `Analyze` → `ComputeGeometry` →
/// `DetectOverlaps` → `ApplyTopology` → `RepairPcurves`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilletPhase {
    /// Parameter validation at the API boundary.
    Validate,
    /// Edge analysis (adjacent faces, convexity, seam detection).
    Analyze,
    /// Per-edge geometry computation (center curve, fillet surface, trim curves).
    ComputeGeometry,
    /// Overlap detection between adjacent fillet trim curves.
    DetectOverlaps,
    /// Topology modification (edge splitting, face splitting, surface assignment).
    ApplyTopology,
    /// Pcurve recomputation and repair after topology modifications.
    RepairPcurves,
}

impl FilletPhase {
    /// Returns the 0-based index of this phase in the pipeline.
    pub fn index(&self) -> usize {
        match self {
            FilletPhase::Validate => 0,
            FilletPhase::Analyze => 1,
            FilletPhase::ComputeGeometry => 2,
            FilletPhase::DetectOverlaps => 3,
            FilletPhase::ApplyTopology => 4,
            FilletPhase::RepairPcurves => 5,
        }
    }

    /// Returns the total number of phases in the pipeline (constant 6).
    pub fn total() -> usize {
        6
    }
}

/// One progress event emitted by the fillet pipeline.
///
/// `phase.index()` and `FilletPhase::total()` let a UI compute a
/// percentage without having to know every phase name.
#[derive(Debug, Clone, Copy)]
pub struct ProgressEvent {
    /// Which pipeline phase just completed.
    pub phase: FilletPhase,
    /// Wall-clock duration of the phase that just completed.
    pub duration: Duration,
    /// Phase-specific entity count (number of edges analysed,
    /// number of trim curves computed, etc.).
    pub entity_count: usize,
}

impl ProgressEvent {
    /// Returns the percentage `0..=100` (integer, rounded).
    pub fn percent(&self) -> u8 {
        ((self.phase.index() + 1) * 100 / FilletPhase::total()) as u8
    }
}

/// A progress callback the fillet pipeline invokes between phases.
///
/// Implementors typically drive a UI progress bar or write metrics
/// to a log. Callers should keep the implementation cheap — the
/// pipeline calls the reporter inline, so a slow callback adds
/// directly to the fillet's wall time.
pub trait ProgressReporter: Send {
    /// Called immediately after a pipeline phase completes.
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
    /// All events recorded so far, in phase order.
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
            phase: FilletPhase::Validate,
            duration: Duration::from_millis(10),
            entity_count: 0,
        };
        assert_eq!(event.percent(), 16);
    }

    #[test]
    fn percent_at_last_phase_is_100() {
        let event = ProgressEvent {
            phase: FilletPhase::RepairPcurves,
            duration: Duration::from_millis(10),
            entity_count: 0,
        };
        assert_eq!(event.percent(), 100);
    }

    #[test]
    fn null_reporter_ignores_events() {
        let mut reporter = NullReporter;
        reporter.on_phase(ProgressEvent {
            phase: FilletPhase::Analyze,
            duration: Duration::from_millis(0),
            entity_count: 0,
        });
    }

    #[test]
    fn capture_reporter_records_events_in_order() {
        let mut reporter = CaptureReporter::default();
        reporter.on_phase(ProgressEvent {
            phase: FilletPhase::Analyze,
            duration: Duration::from_millis(5),
            entity_count: 12,
        });
        reporter.on_phase(ProgressEvent {
            phase: FilletPhase::ApplyTopology,
            duration: Duration::from_millis(20),
            entity_count: 18,
        });
        assert_eq!(reporter.events.len(), 2);
        assert_eq!(reporter.events[0].phase, FilletPhase::Analyze);
        assert_eq!(reporter.events[1].phase, FilletPhase::ApplyTopology);
        assert_eq!(reporter.events[0].entity_count, 12);
    }

    #[test]
    fn phase_index_is_monotonic() {
        let phases = [
            FilletPhase::Validate,
            FilletPhase::Analyze,
            FilletPhase::ComputeGeometry,
            FilletPhase::DetectOverlaps,
            FilletPhase::ApplyTopology,
            FilletPhase::RepairPcurves,
        ];
        for (i, p) in phases.iter().enumerate() {
            assert_eq!(p.index(), i);
        }
    }

    #[test]
    fn total_is_six() {
        assert_eq!(FilletPhase::total(), 6);
    }
}
