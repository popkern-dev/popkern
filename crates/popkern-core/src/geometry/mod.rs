//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod vector;

// Re-export Point3 from the v0 module so existing code still works
pub use crate::geometry_v0::Point3;
pub use vector::Vec3;
