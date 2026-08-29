//! Telemetry for healing operations: timing and operation counts.
//!
//! Every healing operation can optionally record timing data. This
//! enables performance profiling and regression detection.
//!
//! # Usage
//!
//! ```ignore
//! let mut telemetry = HealingTelemetry::new();
//! let start = HealingTelemetry::start_timer();
//! // ... run analysis ...
//! telemetry.record("topological_checks", start.elapsed(), 42);
//! println!("{telemetry}");
//! ```

use std::fmt;
use std::time::Duration;

use crate::time::Stopwatch;

/// Timing data for a single operation.
#[derive(Debug, Clone)]
pub struct OperationTiming {
    /// Name of the operation.
    pub name: String,
    /// How long the operation took.
    pub duration: Duration,
    /// Number of entities processed.
    pub entities_processed: usize,
}

/// Telemetry collector for a healing session.
///
/// Records timing and counts for all operations performed during
/// an analysis or repair session. The collector is owned by the
/// caller — no interior mutability or thread-safety overhead.
///
/// # Example
///
/// ```ignore
/// let mut t = HealingTelemetry::new();
/// let start = HealingTelemetry::start_timer();
/// // ... do work ...
/// t.record_elapsed("gap_repair", start, 12);
/// assert_eq!(t.operation_count(), 1);
/// ```
#[derive(Debug, Clone, Default)]
pub struct HealingTelemetry {
    timings: Vec<OperationTiming>,
    /// Operations whose duration could not be measured because the
    /// platform has no monotonic clock (clockless
    /// `wasm32-unknown-unknown`). Recorded as a count instead of
    /// fabricating zero-duration rows.
    unmeasured_operations: usize,
}

impl HealingTelemetry {
    /// Create an empty telemetry collector.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a completed operation.
    ///
    /// Call this after each discrete phase of analysis or repair
    /// to capture its timing and entity count.
    pub fn record(&mut self, name: impl Into<String>, duration: Duration, entities: usize) {
        self.timings.push(OperationTiming {
            name: name.into(),
            duration,
            entities_processed: entities,
        });
    }

    /// Start a timer — returns a [`Stopwatch`].
    ///
    /// Pass the stopwatch to [`record_elapsed`](Self::record_elapsed)
    /// when the operation completes.
    pub fn start_timer() -> Stopwatch {
        Stopwatch::start()
    }

    /// Record a completed operation timed by `started`. On platforms
    /// without a monotonic clock the row is counted in
    /// [`unmeasured_operations`](Self::unmeasured_operations) instead
    /// of being fabricated with a zero duration.
    pub fn record_elapsed(&mut self, name: impl Into<String>, started: Stopwatch, entities: usize) {
        match started.elapsed() {
            Some(duration) => self.record(name, duration, entities),
            None => self.unmeasured_operations += 1,
        }
    }

    /// Operations that completed without a measurable duration
    /// (clockless platform). Zero everywhere a monotonic clock
    /// exists.
    pub fn unmeasured_operations(&self) -> usize {
        self.unmeasured_operations
    }

    /// Total time across all recorded operations.
    pub fn total_duration(&self) -> Duration {
        self.timings.iter().map(|t| t.duration).sum()
    }

    /// Number of recorded operations.
    pub fn operation_count(&self) -> usize {
        self.timings.len()
    }

    /// Get all recorded timings.
    pub fn timings(&self) -> &[OperationTiming] {
        &self.timings
    }

    /// Find the slowest operation, or `None` if no operations recorded.
    pub fn slowest(&self) -> Option<&OperationTiming> {
        self.timings.iter().max_by_key(|t| t.duration)
    }

    /// Find the fastest operation, or `None` if no operations recorded.
    pub fn fastest(&self) -> Option<&OperationTiming> {
        self.timings.iter().min_by_key(|t| t.duration)
    }

    /// Total entities processed across all operations.
    pub fn total_entities(&self) -> usize {
        self.timings.iter().map(|t| t.entities_processed).sum()
    }
}

impl fmt::Display for HealingTelemetry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.timings.is_empty() {
            return write!(f, "HealingTelemetry: no operations recorded");
        }

        writeln!(f, "HealingTelemetry: {} operations", self.timings.len())?;
        writeln!(
            f,
            "{:<30} {:>12} {:>10}",
            "Operation", "Duration", "Entities"
        )?;
        writeln!(f, "{:-<30} {:->12} {:->10}", "", "", "")?;

        for timing in &self.timings {
            writeln!(
                f,
                "{:<30} {:>12.3?} {:>10}",
                timing.name, timing.duration, timing.entities_processed,
            )?;
        }

        writeln!(f, "{:-<30} {:->12} {:->10}", "", "", "")?;
        write!(
            f,
            "{:<30} {:>12.3?} {:>10}",
            "TOTAL",
            self.total_duration(),
            self.total_entities(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_telemetry_is_empty() {
        let t = HealingTelemetry::new();
        assert_eq!(t.operation_count(), 0);
        assert_eq!(t.total_duration(), Duration::ZERO);
        assert!(t.timings().is_empty());
        assert!(t.slowest().is_none());
        assert!(t.fastest().is_none());
        assert_eq!(t.total_entities(), 0);
    }

    #[test]
    fn record_and_retrieve() {
        let mut t = HealingTelemetry::new();
        t.record("op_a", Duration::from_millis(100), 10);
        t.record("op_b", Duration::from_millis(200), 20);

        assert_eq!(t.operation_count(), 2);
        assert_eq!(t.timings().len(), 2);
        assert_eq!(t.timings()[0].name, "op_a");
        assert_eq!(t.timings()[1].entities_processed, 20);
    }

    #[test]
    fn total_duration_sums_all() {
        let mut t = HealingTelemetry::new();
        t.record("a", Duration::from_millis(100), 5);
        t.record("b", Duration::from_millis(250), 10);
        t.record("c", Duration::from_millis(50), 3);

        assert_eq!(t.total_duration(), Duration::from_millis(400));
    }

    #[test]
    fn total_entities_sums_all() {
        let mut t = HealingTelemetry::new();
        t.record("a", Duration::from_millis(10), 5);
        t.record("b", Duration::from_millis(20), 15);

        assert_eq!(t.total_entities(), 20);
    }

    #[test]
    fn slowest_finds_max_duration() {
        let mut t = HealingTelemetry::new();
        t.record("fast", Duration::from_millis(10), 1);
        t.record("slow", Duration::from_millis(500), 1);
        t.record("medium", Duration::from_millis(100), 1);

        let slowest = t.slowest();
        assert!(slowest.is_some());
        assert_eq!(slowest.map(|o| &o.name), Some(&"slow".to_string()));
    }

    #[test]
    fn fastest_finds_min_duration() {
        let mut t = HealingTelemetry::new();
        t.record("fast", Duration::from_millis(10), 1);
        t.record("slow", Duration::from_millis(500), 1);
        t.record("medium", Duration::from_millis(100), 1);

        let fastest = t.fastest();
        assert!(fastest.is_some());
        assert_eq!(fastest.map(|o| &o.name), Some(&"fast".to_string()));
    }

    #[test]
    fn start_timer_measures_on_clocked_platforms() {
        let timer = HealingTelemetry::start_timer();
        let mut t = HealingTelemetry::new();
        t.record_elapsed("timed", timer, 2);
        assert_eq!(t.operation_count(), 1);
        assert_eq!(t.unmeasured_operations(), 0);
    }

    #[test]
    fn display_empty() {
        let t = HealingTelemetry::new();
        let s = t.to_string();
        assert!(s.contains("no operations recorded"));
    }

    #[test]
    fn display_with_operations() {
        let mut t = HealingTelemetry::new();
        t.record("topological_checks", Duration::from_millis(42), 12);
        t.record("geometric_checks", Duration::from_millis(108), 24);

        let s = t.to_string();
        assert!(s.contains("2 operations"));
        assert!(s.contains("topological_checks"));
        assert!(s.contains("geometric_checks"));
        assert!(s.contains("TOTAL"));
    }

    #[test]
    fn default_is_same_as_new() {
        let d = HealingTelemetry::default();
        let n = HealingTelemetry::new();
        assert_eq!(d.operation_count(), n.operation_count());
        assert_eq!(d.total_duration(), n.total_duration());
    }
}
