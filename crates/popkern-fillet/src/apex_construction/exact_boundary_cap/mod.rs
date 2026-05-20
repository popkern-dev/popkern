//! Exact-boundary apex cap primitives.

pub(crate) mod boundary;
pub(crate) mod sector;

#[cfg(test)]
#[path = "boundary_tests.rs"]
mod boundary_tests;

#[cfg(test)]
#[path = "sector_tests.rs"]
mod sector_tests;
