//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod analytic_curves;
pub mod analytic_surfaces;
pub mod curve;
pub mod nurbs_curve;
pub mod nurbs_ops;
pub mod nurbs_surface;
pub mod point;
pub mod surface;
pub mod vector;

pub use curve::{CurveEval, CurveGeometry, CurveOnEdge};
pub use nurbs_curve::NurbsCurve;
pub use nurbs_surface::NurbsSurface;
pub use point::Point3;
pub use surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
pub use vector::Vec3;
