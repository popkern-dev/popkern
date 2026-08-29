//! Geometry helpers for KP20 cap/support boundary solves.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::{BoundarySide, row_fit};

pub(super) fn cap_boundary_forward(
    cap_patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    start: f64,
    end: f64,
    half_edge_start: Point3,
    half_edge_end: Point3,
) -> bool {
    let p0 = row_fit::eval_boundary(cap_patch, side, start);
    let p1 = row_fit::eval_boundary(cap_patch, side, end);
    p0.distance_to(&half_edge_start) + p1.distance_to(&half_edge_end)
        <= p1.distance_to(&half_edge_start) + p0.distance_to(&half_edge_end)
}

pub(super) fn edge_curve_forward(
    curve: &CurveOnEdge,
    half_edge_start: Point3,
    half_edge_end: Point3,
) -> bool {
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    p0.distance_to(&half_edge_start) + p1.distance_to(&half_edge_end)
        <= p1.distance_to(&half_edge_start) + p0.distance_to(&half_edge_end)
}

pub(super) fn support_pcurve_forward(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    half_edge_start: Point3,
    half_edge_end: Point3,
) -> bool {
    let (t0, t1) = pcurve.parameter_range();
    let p0 = support_surface_point(surface, pcurve, t0);
    let p1 = support_surface_point(surface, pcurve, t1);
    p0.distance_to(&half_edge_start) + p1.distance_to(&half_edge_end)
        <= p1.distance_to(&half_edge_start) + p0.distance_to(&half_edge_end)
}

pub(super) fn support_surface_point(surface: &SurfaceGeometry, pcurve: &Pcurve, t: f64) -> Point3 {
    let uv = pcurve.point_at(t);
    surface.point_at(uv.x, uv.y)
}

pub(super) fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
