//! Extended boolean operations beyond union/intersection/subtraction.

pub mod face_level;
pub(crate) mod half_space;
pub mod imprint;
pub mod interference;
pub mod plane_section;
pub mod section;
pub mod self_intersect;
pub mod self_trim;
mod self_trim_sharpen;
pub mod split_slice;

#[cfg(test)]
#[path = "ops_tests.rs"]
mod tests;
