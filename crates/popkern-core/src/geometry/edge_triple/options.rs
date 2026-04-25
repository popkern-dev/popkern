//! Per-call options, branch policy, and cancellation token.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Per-call options. Defaults are tuned for production fillet
/// pipelines; UI callers customise `cancel`, `progress`, and
/// `prefer_branches`.
#[derive(Clone)]
pub struct EdgeTripleOptions {
    /// Cancellation token checked by the marching solver every 16
    /// iterations. Analytic paths are too fast to need cancellation.
    pub cancel: Option<CancellationToken>,
    /// Progress callback receiving 0.0..=1.0 from marching dispatch.
    /// Not invoked from analytic dispatch.
    pub progress: Option<Arc<dyn Fn(f32) + Send + Sync>>,
    /// If true, run [`crate::geometry::edge_triple::EdgeTriple::verify`]
    /// in release builds (always run in debug). Off by default for
    /// hot paths.
    pub verify_in_release: bool,
    /// Marching-solver iteration cap. Default: 250.
    pub max_marching_iterations: u32,
    /// Disambiguation policy when more than one branch matches.
    pub prefer_branches: BranchPolicy,
}

impl Default for EdgeTripleOptions {
    fn default() -> Self {
        Self {
            cancel: None,
            progress: None,
            verify_in_release: false,
            max_marching_iterations: 250,
            prefer_branches: BranchPolicy::ClosestToAnchors,
        }
    }
}

impl std::fmt::Debug for EdgeTripleOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EdgeTripleOptions")
            .field("cancel", &self.cancel.is_some())
            .field("progress", &self.progress.is_some())
            .field("verify_in_release", &self.verify_in_release)
            .field("max_marching_iterations", &self.max_marching_iterations)
            .field("prefer_branches", &self.prefer_branches)
            .finish()
    }
}

/// Branch disambiguation policy. See ADR §4 for the deterministic
/// 5-step rule.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum BranchPolicy {
    /// Default. Pick the branch whose midpoint is closest to the
    /// midpoint of the line through the anchors.
    ClosestToAnchors,
    /// Pick `branches[0]` after canonical enumeration. Useful when
    /// callers have already disambiguated by hint.
    First,
    /// Reject if more than one branch — caller must supply a hint.
    Reject,
}

/// Cooperative cancellation token. Cheap to clone (Arc-shared
/// AtomicBool). Setting `cancel()` from any thread aborts in-flight
/// marching solves at the next check.
#[derive(Clone, Debug, Default)]
pub struct CancellationToken {
    inner: Arc<AtomicBool>,
}

impl CancellationToken {
    /// Create a new, un-cancelled token.
    pub fn new() -> Self {
        Self::default()
    }

    /// Signal cancellation. Idempotent.
    pub fn cancel(&self) {
        self.inner.store(true, Ordering::Release);
    }

    /// True if `cancel()` has been called.
    pub fn is_cancelled(&self) -> bool {
        self.inner.load(Ordering::Acquire)
    }
}
