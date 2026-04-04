//! Entity provenance tracking for boolean operations.
//!
//! Maintains correspondence between input and output entities so that
//! parametric modeling workflows can trace features through booleans.

pub mod types;

pub use types::EntityMap;
