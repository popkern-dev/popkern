//! Extended boolean operations beyond union/intersection/subtraction.

pub mod face_level;
pub mod imprint;
pub mod interference;
pub mod section;
pub mod split_slice;

#[cfg(test)]
#[path = "ops_tests.rs"]
mod tests;
