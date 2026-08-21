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

// The boolean `CancellationToken` is the kernel-wide token from
// `popkern_core::cancellation` (PK-276 unification): one token held by
// a UI Cancel button cancels a boolean, a fillet, a healing pass, and
// a tessellation at once. The next poll after `cancel()` triggers a
// clean teardown — no inconsistent state survives a cancellation,
// satisfying the §14 transactional integrity requirement.
pub use popkern_core::cancellation::CancellationToken;

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
