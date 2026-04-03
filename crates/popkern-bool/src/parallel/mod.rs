//! Parallel face-face intersection infrastructure.
//!
//! Provides a framework for parallelizing the most expensive part of the
//! boolean pipeline: face-face intersection computation. Face pairs are
//! independent and can be processed in parallel.
//!
//! The current implementation uses sequential processing. When the `rayon`
//! dependency is added, the parallel path can be enabled by switching
//! `process_pairs_sequential` to `process_pairs_parallel`.

pub mod ff_parallel;

#[cfg(test)]
#[path = "parallel_tests.rs"]
mod tests;
