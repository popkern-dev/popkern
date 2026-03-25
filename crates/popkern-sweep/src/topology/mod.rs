//! Topology construction from swept and lofted surfaces.
//!
//! Converts NURBS surfaces and boundary curves produced by the surface
//! generators into a complete BRep solid with proper half-edge wiring.

pub mod caps;
pub mod side_faces;
