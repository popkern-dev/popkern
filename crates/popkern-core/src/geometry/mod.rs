//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod analytic_curves;
pub mod analytic_surfaces;
pub mod curve;
pub mod point;
pub mod surface;
pub mod vector;

pub use curve::{CurveEval, CurveGeometry, CurveOnEdge};
pub use point::Point3;
pub use surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
pub use vector::Vec3;
