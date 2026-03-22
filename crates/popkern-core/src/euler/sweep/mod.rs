//! Sweep operations — higher-level operations built from Euler operators.
//!
//! - [`extrude`](crate::topology::BRep::extrude): linear sweep along a direction vector.
//! - [`revolve`](crate::topology::BRep::revolve): rotational sweep around an axis.

mod extrude;
mod revolve;

pub use extrude::ExtrudeResult;
pub use revolve::RevolveResult;
