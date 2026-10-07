//! Sweep operations — higher-level operations built from Euler operators.
//!
//! - [`extrude`](crate::topology::BRep::extrude): linear sweep along a direction vector.
//! - [`revolve`](crate::topology::BRep::revolve): rotational sweep around an axis.

mod extrude;
mod extrude_geometry;
#[cfg(test)]
#[path = "extrude_geometry_tests.rs"]
mod extrude_geometry_tests;
mod extrude_names;
mod revolve;

pub use extrude::ExtrudeResult;
pub use revolve::RevolveResult;
