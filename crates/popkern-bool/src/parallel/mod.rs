//! Parallel dispatch for face-face intersection (PK-277).
//!
//! Stage 1's per-pair intersection work is embarrassingly parallel:
//! each candidate face pair reads `&BRep` (plain `SlotMap` storage, no
//! interior mutability) and produces its own records. `ff_parallel`
//! holds the per-pair body and dispatches over rayon when the
//! `parallel` feature (default) is enabled, collecting results in
//! candidate order so output is identical at any thread count.

mod ff_junction_refine;
pub mod ff_parallel;

#[cfg(test)]
#[path = "parallel_tests.rs"]
mod parallel_tests;
