//! Vertex blending for multi-fillet junctions.
//!
//! When multiple fillets meet at a vertex, a blend surface fills the gap.
//! This module computes the blend region boundary, constructs the blend
//! surface (mitred or smooth), and handles fillet runout where a filleted
//! edge tapers to zero radius at an unfilleted endpoint.

pub mod boundary;
pub(crate) mod boundary_hex;
pub(crate) mod boundary_ring_sense;
mod boundary_ring_stations;
pub(crate) mod boundary_sampling;
#[cfg(test)]
#[path = "boundary_tests.rs"]
mod boundary_tests;
pub(crate) mod gregory;
pub(crate) mod gregory_inner_cp_shrink;
pub(crate) mod gregory_parallelogram;
pub(crate) mod gregory_three_quad;
mod gregory_three_quad_geometry;
pub(crate) mod gregory_vertex_enclosure;
pub mod pns_adapter;
#[cfg(test)]
#[path = "pns_adapter_tests.rs"]
mod pns_adapter_tests;
pub(crate) mod pns_cap_authoring;
pub mod pns_promote;
#[cfg(test)]
#[path = "pns_promote_tests.rs"]
mod pns_promote_tests;
pub(crate) mod pns_rim_materialize;
pub mod runout;
pub(crate) mod sphere_fit_boundary;
pub(crate) mod sphere_fit_centre;
mod spoke_balancing;
pub mod surface;
pub(crate) mod surface_bridge;
mod surface_math;
mod surface_sphere;
pub(crate) mod tangent_corner;

#[cfg(test)]
#[path = "gregory_tests.rs"]
mod gregory_tests;

#[cfg(test)]
#[path = "gregory_diagnostic_metrics.rs"]
pub(crate) mod gregory_diagnostic_metrics;

#[cfg(test)]
#[path = "gregory_diagnostic_tests.rs"]
mod gregory_diagnostic_tests;

pub use boundary::{
    BlendBoundary, BlendBoundaryOptions, BoundarySegment, compute_blend_boundary,
    compute_blend_boundary_with_options,
};
pub use runout::{RunoutEnd, runout_radius_fn};
pub use surface::{BlendSurfaceResult, build_blend_surface};
