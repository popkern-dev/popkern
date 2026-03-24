//! Surface generation for sweep operations.
//!
//! Transforms a profile along a frame sequence to build skinned NURBS surfaces.

pub mod sweep_surface;

pub use sweep_surface::{SweepSurfaceResult, generate_sweep_surface};
