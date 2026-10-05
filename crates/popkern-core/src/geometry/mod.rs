//! Geometric primitives: points, vectors, curves, and surfaces.

pub mod analytic_curves;
pub mod apex_blend;
pub(crate) mod dense;
pub mod transform;
pub use transform::Transform3;
pub mod analytic_surfaces;
pub mod bernstein_poly;
pub mod certify;
pub mod chart_pcurve;
pub mod chart_periods;
pub mod circle_recognition;
pub mod curve;
pub mod curve_inversion;
pub mod curve_nearest_param;
pub mod curve_on_surface;
pub mod curve_semantics;
mod curve_to_nurbs;
pub mod edge_triple;
pub mod fitting;
pub mod gregory;
pub mod intersection;
pub mod n_sided_fill;
pub mod nurbs_curve;
pub mod nurbs_ops;
pub mod nurbs_surface;
pub mod pcurve;
pub mod pcurve_split;
pub mod pcurve_trace;
pub mod point;
pub mod polyhedral_spline;
pub mod polyhedral_spline_patch;
pub mod revolution;
pub mod revolution_iso_pcurve;
pub mod surface;
pub mod surface_helpers;
pub mod trim_domain;
pub mod uv_projection;
pub(crate) mod uv_projection_face_box;
pub mod variational_cap;
pub mod vector;

pub use curve::{CurveEval, CurveGeometry, CurveOnEdge};
pub use curve_semantics::LineSegmentParams;
pub use intersection::{CurveIntersection, CurveIntersectionResult, intersect_curves};
pub use nurbs_curve::NurbsCurve;
pub use nurbs_surface::NurbsSurface;
pub use point::Point3;
/// Portable transcendentals (`x.psin()`, `y.patan2(x)`): the same bits
/// on every target. Kernel code calls these, never the inherent
/// `f64::sin` family (see `popkern_math`).
pub use popkern_math::PortableMath;
pub use surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
pub use vector::Vec3;

#[cfg(test)]
#[path = "surface_deriv2_tests.rs"]
mod surface_deriv2_tests;
pub use n_sided_fill::{FillBoundaryEdge, NSidedFill, NSidedFillError, fill_n_sided_boundary};
pub use polyhedral_spline_patch::{PolyhedralSplinePatch, PolyhedralSplineProvenance};
