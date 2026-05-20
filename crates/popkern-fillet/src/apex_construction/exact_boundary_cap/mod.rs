//! Exact-boundary apex cap primitives.

pub(crate) mod boundary;
pub(crate) mod corner;
pub(crate) mod ribbon;
pub(crate) mod sector;
pub(crate) mod solve;

#[cfg(test)]
#[path = "boundary_tests.rs"]
mod boundary_tests;

#[cfg(test)]
#[path = "corner_tests.rs"]
mod corner_tests;

#[cfg(test)]
#[path = "ribbon_tests.rs"]
mod ribbon_tests;

#[cfg(test)]
#[path = "sector_tests.rs"]
mod sector_tests;

#[cfg(test)]
#[path = "solve_tests.rs"]
mod solve_tests;
