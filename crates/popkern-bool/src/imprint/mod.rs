//! Imprinting infrastructure for the boolean pipeline.
//!
//! Provides utilities for snap-to-existing (vertex merging), sliver detection,
//! and shared geometry validation. These are consumed by the graph-merge
//! pipeline and the assembly phase to produce watertight, sliver-free output.

pub mod micro_edge;
pub mod shared_geometry;
pub mod sliver;
pub mod snap;

#[cfg(test)]
#[path = "imprint_tests.rs"]
mod tests;
