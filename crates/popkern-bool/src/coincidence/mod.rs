//! Coincidence engine for boolean operations.
//!
//! Detects and classifies coincident geometry between two BRep bodies.
//! Coincident faces (faces whose underlying surfaces are identical or
//! near-identical) require special handling during boolean operations —
//! the standard intersection algorithm produces garbage for these cases.
//!
//! The engine produces a [`CoincidenceMap`] that is consumed by the
//! classification and assembly phases to apply the correct ON-face rules
//! from the boolean spec (SS3.3).

pub mod types;

pub(crate) mod classify;
pub(crate) mod detect;

pub use types::CoincidenceMap;

#[cfg(test)]
#[path = "coincidence_tests.rs"]
mod tests;
