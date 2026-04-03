//! Post-processing for boolean results: redundant edge removal and face simplification.

pub mod face_simplify;
pub mod redundant_edge;

#[cfg(test)]
#[path = "postprocess_tests.rs"]
mod tests;
