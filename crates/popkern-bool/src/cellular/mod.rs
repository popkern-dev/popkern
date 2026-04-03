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

pub mod boolean;
pub mod cell;
pub mod decompose;
pub mod select;

pub use cell::{Cell, CellDecomposition};

#[cfg(test)]
#[path = "cellular_tests.rs"]
mod tests;
