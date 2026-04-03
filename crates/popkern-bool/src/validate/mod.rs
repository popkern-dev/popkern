//! Boolean-specific validation and repair.
//!
//! Pre-check validates inputs before the boolean, post-check validates
//! the result, and the repair pass closes small watertightness gaps.

pub mod post_check;
pub mod pre_check;
pub mod repair;

#[cfg(test)]
#[path = "validate_tests.rs"]
mod tests;
