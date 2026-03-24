//! Surface generation for sweep operations.
//!
//! Transforms a profile along a frame sequence to build skinned NURBS surfaces.

pub mod loft_surface;
pub mod sweep_surface;

pub use loft_surface::{LoftSurfaceResult, generate_loft_surface};
pub use sweep_surface::{SweepSurfaceResult, generate_sweep_surface};
