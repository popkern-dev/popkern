//! Progress reporting and §18 telemetry for long-running fillet
//! operations.
//!
//! The fillet pipeline calls a [`ProgressReporter`] between every
//! major phase (phase-level events via [`ProgressEvent`]) AND after
//! every individual edge's geometry computation (edge-level events
//! via [`EdgeTelemetry`]). Callers implement the trait to drive a
//! UI progress bar, log to stderr, or aggregate metrics into a
//! structured [`FilletTelemetry`] report.
//!
//! Combine with [`crate::cancellation::CancellationToken`] for the
//! full §17 interactive workflow: progress callback drives the UI,
//! cancellation token responds to "Cancel".
//!
//! # §18 telemetry acceptance
//!
//! §18 of `docs/requirements/fillets_blends.md` requires the engine
//! to expose **per-edge timing, per-phase timing, Newton iteration
//! counts (mean, max, histogram), step rejection counts, memory
//! high-water, cache hit rates — all queryable**. This module ships
//! timing and iteration telemetry via [`CaptureReporter::telemetry`]
//! plus step rejection counts via [`EdgeTelemetry::step_rejections`]
//! and the `FilletTelemetry::*_step_rejections` aggregate accessors.
//! Memory high-water is measured by the `alloc_profile` integration
//! test's `CountingAlloc` at `crates/popkern-fillet/tests/alloc_profile.rs`.
//!
//! Cache hit rates are exposed by individual cache instances rather
//! than the fillet-wide telemetry struct. The fillet pipeline itself
//! currently maintains no caches — every Newton solve, trim curve,
//! and pcurve is built from scratch per edge — so there is nothing
//! for `FilletTelemetry` to aggregate on the fillet side. The
//! reference implementation of the cache-hit-rate pattern lives at
//! `popkern_core::healing::CacheStats` (and is wired into
//! `HealingCache`); any future cache added to the fillet pipeline
//! should reuse that type so the metric stays consistent across
//! the engine.

use popkern_core::topology::EdgeHandle;
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

/// One per-corner progress event emitted from inside the
/// `ApplyTopology` phase (PK-276). Vertex-blend corner construction
/// dominates that phase's wall clock, so phase-level events alone
/// leave the UI frozen through the longest stretch of a multi-edge
/// fillet; these events tick once per corner.
#[derive(Debug, Clone, Copy)]
pub struct CornerProgress {
    /// The vertex whose corner construction just completed (success
    /// or per-corner fallback — either way the corner is done).
    pub vertex: popkern_core::topology::VertexHandle,
    /// 0-based corner position in the deterministic (`BTreeMap`)
    /// vertex order.
    pub corner_index: usize,
    /// Number of corners this apply phase will construct.
    pub total_corners: usize,
    /// Wall-clock duration of this corner's construction.
    pub duration: Duration,
}

/// Per-edge §18 telemetry captured by the fillet pipeline's
/// compute loop. One record per edge in the spec list, emitted
/// immediately after `compute_fillet_geometry` returns for that
/// edge — whether the edge succeeded or failed.
#[derive(Debug, Clone, Copy)]
pub struct EdgeTelemetry {
    /// The edge handle the telemetry applies to.
    pub edge: EdgeHandle,
    /// Wall-clock duration of this edge's geometry computation
    /// (the time spent inside `compute_fillet_geometry` for this
    /// edge alone, excluding the pipeline-wide phases).
    pub duration: Duration,
    /// Total Newton iterations consumed by the rolling-ball
    /// contact solver across all walker samples for this edge.
    /// `None` on the legacy offset-intersect path (no walker, no
    /// Newton iteration telemetry). `Some(0)` would indicate the
    /// walker ran but produced no contact solutions — shouldn't
    /// happen in practice.
    pub newton_iterations: Option<usize>,
    /// Total backtracking-line-search rejections summed across
    /// every walker sample's contact solve. `None` on the legacy
    /// offset-intersect path. A high value relative to
    /// `newton_iterations` indicates ill-conditioned contact
    /// geometry — useful for flagging specific edges to a user.
    pub step_rejections: Option<usize>,
    /// `true` if the edge's geometry compute succeeded, `false`
    /// if it returned a `BlendFailure`.
    pub succeeded: bool,
}

/// Aggregated §18 telemetry summary. Produced by
/// [`CaptureReporter::telemetry`] after a fillet operation
/// completes. Gives queryable access to every metric §18 asks
/// for, without the caller having to iterate raw event vectors.
#[derive(Debug, Clone, Default)]
pub struct FilletTelemetry {
    /// Phase-level events in pipeline order.
    pub phase_events: Vec<ProgressEvent>,
    /// Per-edge events in the order edges were processed.
    pub edge_events: Vec<EdgeTelemetry>,
}

impl FilletTelemetry {
    /// Sum of all phase durations. Equivalent to the wall-clock
    /// time the reporter was called across, approximately equal
    /// to the total `fillet_edges` runtime.
    pub fn total_phase_duration(&self) -> Duration {
        self.phase_events.iter().map(|e| e.duration).sum()
    }

    /// Sum of all per-edge compute durations. This is *less than*
    /// [`Self::total_phase_duration`] because many phases (validate,
    /// analyze, apply-topology) run outside the per-edge loop.
    pub fn total_edge_duration(&self) -> Duration {
        self.edge_events.iter().map(|e| e.duration).sum()
    }

    /// Number of edges that succeeded (edge events with
    /// `succeeded = true`).
    pub fn edges_succeeded(&self) -> usize {
        self.edge_events.iter().filter(|e| e.succeeded).count()
    }

    /// Number of edges that failed.
    pub fn edges_failed(&self) -> usize {
        self.edge_events.iter().filter(|e| !e.succeeded).count()
    }

    /// Total Newton iterations summed across every walker edge
    /// that reported iteration counts. Edges on the legacy offset-
    /// intersect path (where `newton_iterations = None`) are
    /// excluded from the sum.
    pub fn total_newton_iterations(&self) -> usize {
        self.edge_events
            .iter()
            .filter_map(|e| e.newton_iterations)
            .sum()
    }

    /// Mean Newton iterations per walker edge, or `None` when no
    /// edge reported iteration counts.
    pub fn mean_newton_iterations(&self) -> Option<f64> {
        let counts: Vec<usize> = self
            .edge_events
            .iter()
            .filter_map(|e| e.newton_iterations)
            .collect();
        if counts.is_empty() {
            return None;
        }
        Some(counts.iter().sum::<usize>() as f64 / counts.len() as f64)
    }

    /// Maximum Newton iterations on any walker edge, or `None`
    /// when no edge reported iteration counts.
    pub fn max_newton_iterations(&self) -> Option<usize> {
        self.edge_events
            .iter()
            .filter_map(|e| e.newton_iterations)
            .max()
    }

    /// Newton-iteration histogram as a sorted `(count, frequency)`
    /// vector. Example: `[(3, 2), (5, 8), (7, 1)]` means 2 edges
    /// used 3 iterations, 8 edges used 5, 1 edge used 7.
    /// Useful for spotting outliers in the walker's convergence.
    pub fn newton_histogram(&self) -> Vec<(usize, usize)> {
        let mut hist: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
        for e in &self.edge_events {
            if let Some(count) = e.newton_iterations {
                *hist.entry(count).or_insert(0) += 1;
            }
        }
        hist.into_iter().collect()
    }

    /// Total backtracking line-search rejections summed across
    /// every walker edge. Edges on the legacy offset-intersect
    /// path (`step_rejections = None`) are excluded from the sum.
    pub fn total_step_rejections(&self) -> usize {
        self.edge_events
            .iter()
            .filter_map(|e| e.step_rejections)
            .sum()
    }

    /// Mean step rejections per walker edge, or `None` when no
    /// edge reported rejection counts.
    pub fn mean_step_rejections(&self) -> Option<f64> {
        let counts: Vec<usize> = self
            .edge_events
            .iter()
            .filter_map(|e| e.step_rejections)
            .collect();
        if counts.is_empty() {
            return None;
        }
        Some(counts.iter().sum::<usize>() as f64 / counts.len() as f64)
    }

    /// Maximum step rejections on any walker edge, or `None`
    /// when no edge reported rejection counts.
    pub fn max_step_rejections(&self) -> Option<usize> {
        self.edge_events
            .iter()
            .filter_map(|e| e.step_rejections)
            .max()
    }

    /// Ratio of step rejections to Newton iterations across every
    /// walker edge that reported both. Returns `None` when no walker
    /// edge ran. Intuition: 0.0 = line search always accepts on the
    /// first try (healthy); values above ~0.5 suggest the contact
    /// solver is burning half its iterations rejecting candidate
    /// steps — an ill-conditioned-surface flag.
    pub fn step_rejection_ratio(&self) -> Option<f64> {
        let iters: usize = self
            .edge_events
            .iter()
            .filter_map(|e| e.newton_iterations)
            .sum();
        if iters == 0 {
            return None;
        }
        let rejs: usize = self
            .edge_events
            .iter()
            .filter_map(|e| e.step_rejections)
            .sum();
        Some(rejs as f64 / iters as f64)
    }
}

/// A progress callback the fillet pipeline invokes between phases
/// and after each per-edge geometry computation.
///
/// Implementors typically drive a UI progress bar or write metrics
/// to a log. Callers should keep the implementation cheap — the
/// pipeline calls the reporter inline, so a slow callback adds
/// directly to the fillet's wall time.
pub trait ProgressReporter: Send {
    /// Called immediately after a pipeline phase completes.
    fn on_phase(&mut self, event: ProgressEvent);
    /// Called immediately after an edge's geometry computation
    /// completes (success or failure). Default implementation is
    /// a no-op so existing `ProgressReporter` impls don't need
    /// updating; override to capture per-edge telemetry.
    fn on_edge(&mut self, _telemetry: EdgeTelemetry) {}
    /// Called immediately after each vertex-blend corner construction
    /// inside the `ApplyTopology` phase (PK-276). Default is a no-op
    /// so existing impls don't need updating; override to drive
    /// fine-grained progress through the dominant phase.
    fn on_corner(&mut self, _event: CornerProgress) {}
}

/// A no-op reporter that ignores all events. Useful as a default
/// when the caller does not need progress notifications.
#[derive(Debug, Default, Clone, Copy)]
pub struct NullReporter;

impl ProgressReporter for NullReporter {
    fn on_phase(&mut self, _event: ProgressEvent) {}
}

/// A reporter that records every event into vectors. Primarily
/// useful for tests that need to verify the pipeline emitted the
/// expected sequence of events, and for callers who want a
/// structured §18 telemetry report after the fillet completes.
#[derive(Debug, Default, Clone)]
pub struct CaptureReporter {
    /// All phase events recorded so far, in pipeline order.
    pub events: Vec<ProgressEvent>,
    /// All edge events recorded so far, in the order edges were
    /// processed.
    pub edge_events: Vec<EdgeTelemetry>,
    /// All per-corner events recorded so far, in deterministic
    /// corner order (PK-276).
    pub corner_events: Vec<CornerProgress>,
}

impl CaptureReporter {
    /// Produce a structured [`FilletTelemetry`] report from the
    /// captured events. Clones the internal vectors — cheap for
    /// realistic batch sizes (hundreds of edges).
    pub fn telemetry(&self) -> FilletTelemetry {
        FilletTelemetry {
            phase_events: self.events.clone(),
            edge_events: self.edge_events.clone(),
        }
    }
}

impl ProgressReporter for CaptureReporter {
    fn on_phase(&mut self, event: ProgressEvent) {
        self.events.push(event);
    }
    fn on_edge(&mut self, telemetry: EdgeTelemetry) {
        self.edge_events.push(telemetry);
    }
    fn on_corner(&mut self, event: CornerProgress) {
        self.corner_events.push(event);
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
