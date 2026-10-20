//! Topology construction from swept and lofted surfaces.
//!
//! Converts NURBS surfaces and boundary curves produced by the surface
//! generators into a complete BRep solid with proper half-edge wiring.

pub mod caps;
pub(crate) mod periodic_wall;
pub mod revolution_canonical;
pub mod revolution_faces;
pub mod side_faces;
pub(crate) mod side_faces_multi_wire;
pub mod sweep_naming;
