//! Progress reporting and cancellation for long-running healing operations.
//!
//! Provides a [`ProgressCallback`] trait for progress updates and a
//! [`CancellationToken`] for cooperative cancellation. These enable UI
//! integration without coupling the healing engine to any specific
//! framework.
//!
//! # Usage
//!
//! For simple cases where no progress reporting is needed, use [`NoProgress`]:
//!
//! ```ignore
//! let progress = NoProgress;
//! assert!(!progress.is_cancelled());
//! ```
//!
//! For cancellable operations, create a [`CancellationToken`] and share
//! it between the operation thread and the UI:
//!
//! ```ignore
//! let token = CancellationToken::new();
//! let token_clone = token.clone();
//! // In the UI thread:
//! token_clone.cancel();
//! // In the worker thread:
//! if token.is_cancelled() { return; }
//! ```

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

pub use crate::cancellation::CancellationToken;

/// Progress callback trait for healing operations.
///
/// Implementations receive periodic progress updates and can signal
/// cancellation. The trait is `Send + Sync` to support cross-thread use.
pub trait ProgressCallback: Send + Sync {
    /// Called periodically with a progress update.
    ///
    /// - `phase`: human-readable name of the current phase (e.g., "gap repair").
    /// - `current`: number of items processed so far in this phase.
    /// - `total`: total number of items expected in this phase.
    fn on_progress(&self, phase: &str, current: usize, total: usize);

    /// Check if the operation should be cancelled.
    ///
    /// Called before each major step. If this returns `true`, the
    /// operation should abort as soon as possible.
    fn is_cancelled(&self) -> bool;
}

// The healing `CancellationToken` is the kernel-wide token from
// `crate::cancellation` (PK-276 unification — this module carried its
// own near-identical definition before; re-exporting keeps every
// `popkern_core::healing::CancellationToken` caller compiling while
// letting one token cancel healing, booleans, fillets, and meshing).

/// No-op progress callback (for when progress reporting isn't needed).
///
/// Never reports cancellation. All progress updates are silently ignored.
pub struct NoProgress;

impl ProgressCallback for NoProgress {
    fn on_progress(&self, _phase: &str, _current: usize, _total: usize) {}
    fn is_cancelled(&self) -> bool {
        false
    }
}

/// A progress callback that stores the latest progress for polling.
///
/// Thread-safe: the UI thread can poll `snapshot()` while the worker
/// thread calls `on_progress()`. Cancellation is signalled via the
/// embedded [`CancellationToken`].
#[derive(Debug)]
pub struct ProgressTracker {
    phase: Mutex<String>,
    current: AtomicUsize,
    total: AtomicUsize,
    token: CancellationToken,
}

impl ProgressTracker {
    /// Create a new tracker with the given cancellation token.
    pub fn new(token: CancellationToken) -> Self {
        Self {
            phase: Mutex::new(String::new()),
            current: AtomicUsize::new(0),
            total: AtomicUsize::new(0),
            token,
        }
    }

    /// Take a snapshot of the current progress state.
    ///
    /// Returns `(phase, current, total)`.
    pub fn snapshot(&self) -> (String, usize, usize) {
        let phase = self.phase.lock().map(|g| g.clone()).unwrap_or_default();
        let current = self.current.load(Ordering::Relaxed);
        let total = self.total.load(Ordering::Relaxed);
        (phase, current, total)
    }

    /// Get a reference to the cancellation token.
    pub fn token(&self) -> &CancellationToken {
        &self.token
    }
}

impl ProgressCallback for ProgressTracker {
    fn on_progress(&self, phase: &str, current: usize, total: usize) {
        if let Ok(mut guard) = self.phase.lock() {
            guard.clear();
            guard.push_str(phase);
        }
        self.current.store(current, Ordering::Relaxed);
        self.total.store(total, Ordering::Relaxed);
    }

    fn is_cancelled(&self) -> bool {
        self.token.is_cancelled()
    }
}

/// A progress callback backed by a closure.
///
/// Useful for ad-hoc progress reporting without implementing the full trait.
pub struct ClosureProgress<F>
where
    F: Fn(&str, usize, usize) + Send + Sync,
{
    callback: F,
    token: CancellationToken,
}

impl<F> ClosureProgress<F>
where
    F: Fn(&str, usize, usize) + Send + Sync,
{
    /// Create a new closure-based progress callback.
    pub fn new(callback: F, token: CancellationToken) -> Self {
        Self { callback, token }
    }
}

impl<F> ProgressCallback for ClosureProgress<F>
where
    F: Fn(&str, usize, usize) + Send + Sync,
{
    fn on_progress(&self, phase: &str, current: usize, total: usize) {
        (self.callback)(phase, current, total);
    }

    fn is_cancelled(&self) -> bool {
        self.token.is_cancelled()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn cancellation_token_starts_uncancelled() {
        let token = CancellationToken::new();
        assert!(!token.is_cancelled());
    }

    #[test]
    fn cancellation_token_cancel_sets_flag() {
        let token = CancellationToken::new();
        token.cancel();
        assert!(token.is_cancelled());
    }

    #[test]
    fn cancellation_token_clone_shares_state() {
        let token = CancellationToken::new();
        let clone = token.clone();
        assert!(!clone.is_cancelled());
        token.cancel();
        assert!(clone.is_cancelled());
    }

    #[test]
    fn cancellation_token_default() {
        let token = CancellationToken::default();
        assert!(!token.is_cancelled());
    }

    #[test]
    fn no_progress_does_nothing() {
        let p = NoProgress;
        p.on_progress("test", 5, 10);
        assert!(!p.is_cancelled());
    }

    #[test]
    fn progress_tracker_stores_updates() {
        let token = CancellationToken::new();
        let tracker = ProgressTracker::new(token.clone());

        tracker.on_progress("gap repair", 3, 10);
        let (phase, current, total) = tracker.snapshot();
        assert_eq!(phase, "gap repair");
        assert_eq!(current, 3);
        assert_eq!(total, 10);

        tracker.on_progress("simplification", 7, 20);
        let (phase2, current2, total2) = tracker.snapshot();
        assert_eq!(phase2, "simplification");
        assert_eq!(current2, 7);
        assert_eq!(total2, 20);
    }

    #[test]
    fn progress_tracker_respects_cancellation() {
        let token = CancellationToken::new();
        let tracker = ProgressTracker::new(token.clone());

        assert!(!tracker.is_cancelled());
        token.cancel();
        assert!(tracker.is_cancelled());
    }

    #[test]
    fn progress_tracker_token_accessor() {
        let token = CancellationToken::new();
        let tracker = ProgressTracker::new(token);
        assert!(!tracker.token().is_cancelled());
        tracker.token().cancel();
        assert!(tracker.is_cancelled());
    }

    #[test]
    fn closure_progress_invokes_callback() {
        use std::sync::atomic::AtomicUsize;

        let count = Arc::new(AtomicUsize::new(0));
        let count_clone = count.clone();
        let token = CancellationToken::new();

        let progress = ClosureProgress::new(
            move |_phase, _cur, _tot| {
                count_clone.fetch_add(1, Ordering::Relaxed);
            },
            token.clone(),
        );

        progress.on_progress("test", 0, 5);
        progress.on_progress("test", 1, 5);
        assert_eq!(count.load(Ordering::Relaxed), 2);
        assert!(!progress.is_cancelled());

        token.cancel();
        assert!(progress.is_cancelled());
    }
}
