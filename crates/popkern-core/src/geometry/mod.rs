//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod analytic_curves;
pub mod analytic_surfaces;
pub mod curve;
pub mod curve_on_surface;
pub mod fitting;
pub mod intersection;
pub mod nurbs_curve;
pub mod nurbs_ops;
pub mod nurbs_surface;
pub mod point;
pub mod surface;
pub mod vector;

pub use curve::{CurveEval, CurveGeometry, CurveOnEdge};
pub use intersection::{CurveIntersection, CurveIntersectionResult, intersect_curves};
pub use nurbs_curve::NurbsCurve;
pub use nurbs_surface::NurbsSurface;
pub use point::Point3;
pub use surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
pub use vector::Vec3;
