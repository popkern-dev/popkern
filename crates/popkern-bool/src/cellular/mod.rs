//! Cellular decomposition for boolean operations.
//!
//! The cellular approach partitions the space into disjoint volumetric
//! regions (cells) bounded by the refined face fragments. Each cell knows
//! which input bodies contain it. Boolean operations are reduced to
//! selecting cells by a predicate on their membership:
//!
//! - **Union** = cells in A OR B
//! - **Intersection** = cells in A AND B
//! - **Subtraction** = cells in A AND NOT B
//! - **Symmetric difference** = cells in A XOR B
//!
//! This architecture subsumes all boolean types and enables custom
//! predicates without modifying intersection, imprinting, or assembly.

pub(crate) mod assemble;
pub(crate) mod assemble_cleanup;
pub(crate) mod assemble_edge_geometry;
mod assemble_span_unify;
pub mod boolean;
pub mod cell;
pub(crate) mod classify_fragments;
pub mod decompose;
pub(crate) mod fragment_classify_one;
mod fragment_far;
pub(crate) mod fragment_interior;
pub(crate) mod fragment_propagation;
pub(crate) mod fragment_sample_helpers;
pub(crate) mod fragment_sampling;
pub(crate) mod fragment_torus_sample;
mod probe_clearance;
mod record_rotation;
pub mod select;

pub use cell::{Cell, CellDecomposition};

#[cfg(test)]
#[path = "cellular_tests.rs"]
mod tests;
