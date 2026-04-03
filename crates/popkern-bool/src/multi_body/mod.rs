//! Multi-body boolean operations.
//!
//! Compute the boolean of N bodies simultaneously with globally consistent
//! intersection curves. Uses the cellular decomposition from Phase 6 with
//! N-bit membership to classify cells against all bodies at once.

pub mod multi_boolean;

#[cfg(test)]
#[path = "multi_body_tests.rs"]
mod tests;
