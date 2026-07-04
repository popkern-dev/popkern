//! Surface generation for sweep operations.
//!
//! Transforms a profile along a frame sequence to build skinned NURBS surfaces.

pub(crate) mod canonical;
pub mod guided;
pub mod loft_continuity;
pub mod loft_surface;
pub mod sweep_surface;
pub mod swept_blend;

pub use guided::{GuideValidation, apply_guide_rails, validate_guide_rails};
pub use loft_surface::{LoftSurfaceResult, generate_loft_surface};
pub use sweep_surface::{SweepSurfaceResult, generate_sweep_surface};
pub use swept_blend::generate_swept_blend_surface;
