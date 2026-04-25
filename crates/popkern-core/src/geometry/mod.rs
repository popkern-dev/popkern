//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod analytic_curves;
pub mod transform;
pub use transform::Transform3;
pub mod analytic_surfaces;
pub mod bernstein_poly;
pub mod curve;
pub mod curve_on_surface;
pub mod edge_triple;
pub mod fitting;
pub mod gregory;
pub mod intersection;
pub mod nurbs_curve;
pub mod nurbs_ops;
pub mod nurbs_surface;
pub mod point;
pub mod revolution;
pub mod surface;
pub mod uv_projection;
pub mod vector;

pub use curve::{CurveEval, CurveGeometry, CurveOnEdge};
pub use intersection::{CurveIntersection, CurveIntersectionResult, intersect_curves};
pub use nurbs_curve::NurbsCurve;
pub use nurbs_surface::NurbsSurface;
pub use point::Point3;
pub use surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
pub use vector::Vec3;

#[cfg(test)]
#[path = "surface_deriv2_tests.rs"]
mod surface_deriv2_tests;
