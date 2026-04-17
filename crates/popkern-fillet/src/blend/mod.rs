//! Vertex blending for multi-fillet junctions.
//!
//! When multiple fillets meet at a vertex, a blend surface fills the gap.
//! This module computes the blend region boundary, constructs the blend
//! surface (mitred or smooth), and handles fillet runout where a filleted
//! edge tapers to zero radius at an unfilleted endpoint.

pub mod boundary;
pub(crate) mod gregory;
pub mod runout;
pub mod surface;
mod surface_math;

#[cfg(test)]
#[path = "gregory_tests.rs"]
mod gregory_tests;

#[cfg(test)]
#[path = "gregory_diagnostic_tests.rs"]
mod gregory_diagnostic_tests;

pub use boundary::{BlendBoundary, BoundarySegment, compute_blend_boundary};
pub use runout::{RunoutEnd, runout_radius_fn};
pub use surface::{BlendSurfaceResult, build_blend_surface};
