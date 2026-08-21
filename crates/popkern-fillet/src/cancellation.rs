//! Cooperative cancellation for long-running fillet operations (§17).
//!
//! Callers (typically a GUI dialog or a CLI `Ctrl-C` handler) signal
//! cancellation by flipping an [`AtomicBool`]. The fillet pipeline
//! checks the token between major phases and returns
//! [`crate::BlendError::Cancelled`] if it observes a request, leaving
//! the input BRep unchanged.
//!
//! Tokens are cheap to clone (`Arc<AtomicBool>`) so the producer can
//! keep one handle and pass clones into worker threads or pipeline
//! stages without lifetime gymnastics. The pattern is lifted verbatim
//! from `popkern-bool::cancellation`.
//!
//! # Example
//!
//! ```ignore
//! use popkern_fillet::cancellation::CancellationToken;
//!
//! let token = CancellationToken::new();
//! let token_for_worker = token.clone();
//! std::thread::spawn(move || {
//!     // ... time passes, user clicks "Cancel"
//!     token_for_worker.cancel();
//! });
//!
//! // From the fillet call site:
//! let result = popkern_fillet::fillet_edges_cancellable(
//!     &brep, &specs, &options, tol, None, &token,
//! );
//! assert!(matches!(result, Err(popkern_fillet::BlendError::Cancelled)));
//! ```

// The fillet `CancellationToken` is the kernel-wide token from
// `popkern_core::cancellation` (PK-276 unification): one token held by
// a UI Cancel button cancels a fillet, a boolean, a healing pass, and
// a tessellation at once. The next poll after `cancel()` triggers a
// clean teardown — no inconsistent state survives a cancellation.
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
}
