//! Cooperative cancellation for long-running boolean operations.
//!
//! Implements §15 of the boolean requirements: callers (typically a GUI
//! progress dialog or a CLI Ctrl-C handler) signal cancellation by
//! flipping an [`AtomicBool`]. The boolean pipeline checks the token
//! between major phases and returns [`BoolError::Cancelled`] if it
//! observes a request, leaving the input model unchanged.
//!
//! The token is cheap to clone (`Arc<AtomicBool>`) so the producer can
//! keep one handle and pass clones into worker threads or pipeline
//! stages without lifetime gymnastics.
//!
//! # Example
//!
//! ```ignore
//! use popkern_bool::cancellation::CancellationToken;
//! use std::sync::Arc;
//!
//! let token = CancellationToken::new();
//! let token_for_worker = token.clone();
//! std::thread::spawn(move || {
//!     // ... time passes, user clicks "Cancel"
//!     token_for_worker.cancel();
//! });
//!
//! // From the boolean call site:
//! let result = popkern_bool::boolean_cancellable(
//!     popkern_bool::BoolOperation::Union, a, b, tol, &token,
//! );
//! ```

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// A cooperative cancellation flag shared across threads.
///
/// The boolean pipeline polls [`Self::is_cancelled`] between phases.
/// Calling [`Self::cancel`] from any thread causes the next poll to
/// return `true`, which triggers a clean teardown — no inconsistent
/// state survives a cancellation, satisfying the §14 transactional
/// integrity requirement.
#[derive(Debug, Clone, Default)]
pub struct CancellationToken {
    flag: Arc<AtomicBool>,
}

impl CancellationToken {
    /// Create a fresh, non-cancelled token.
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Signal cancellation. All clones of this token observe the request.
    pub fn cancel(&self) {
        self.flag.store(true, Ordering::SeqCst);
    }

    /// Returns `true` if any holder of this token has called [`Self::cancel`].
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::SeqCst)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_token_is_not_cancelled() {
        let token = CancellationToken::new();
        assert!(!token.is_cancelled());
    }

    #[test]
    fn cancel_propagates_to_clones() {
        let token = CancellationToken::new();
        let clone = token.clone();
        token.cancel();
        assert!(clone.is_cancelled());
        assert!(token.is_cancelled());
    }

    #[test]
    fn default_token_is_not_cancelled() {
        let token = CancellationToken::default();
        assert!(!token.is_cancelled());
    }

    #[test]
    fn boolean_cancellable_returns_cancelled_when_pre_signaled() {
        use crate::{BoolOperation, boolean_cancellable};
        use popkern_core::geometry::Transform3;
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol).unwrap();
        let b = BRep::unit_cube(tol)
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap();
        let token = CancellationToken::new();
        token.cancel(); // Cancel before starting.
        let result = boolean_cancellable(BoolOperation::Union, &a, &b, tol, &token);
        assert!(matches!(result, Err(crate::BoolError::Cancelled)));
    }

    #[test]
    fn boolean_cancellable_succeeds_when_token_not_signaled() {
        use crate::{BoolOperation, boolean_cancellable};
        use popkern_core::geometry::Transform3;
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let tol = Tolerance::default();
        let a = BRep::unit_cube(tol).unwrap();
        let b = BRep::unit_cube(tol)
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap();
        let token = CancellationToken::new();
        let result = boolean_cancellable(BoolOperation::Union, &a, &b, tol, &token);
        assert!(result.is_ok(), "non-cancelled call should succeed");
    }
}
