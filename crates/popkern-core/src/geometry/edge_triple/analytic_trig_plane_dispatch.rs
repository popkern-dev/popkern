use crate::geometry::edge_triple::analytic_trig_cone::try_plane_cone_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_cone_hyperbola::try_plane_cone_oblique_hyperbola;
use crate::geometry::edge_triple::analytic_trig_plane_cone_oblique::try_plane_cone_oblique_ellipse;
use crate::geometry::edge_triple::analytic_trig_plane_cone_parabola::try_plane_cone_oblique_parabola;
use crate::geometry::edge_triple::analytic_trig_plane_cyl::try_plane_cylinder_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_cyl_oblique::try_plane_cylinder_oblique;
use crate::geometry::edge_triple::analytic_trig_plane_sphere::try_plane_sphere_polar;
use crate::geometry::edge_triple::analytic_trig_plane_sphere_general::try_plane_sphere_general;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::EdgeTriple;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cylinder_pair(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    swap: bool,
    co: Point3,
    ca: Vec3,
    cr: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_plane_cylinder_perpendicular(po, pn, pu, pv, swap, co, ca, cr, anchors, tol)
        .or_else(|| try_plane_cylinder_oblique(po, pn, pu, pv, swap, co, ca, cr, anchors, tol))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_sphere_pair(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    swap: bool,
    center: Point3,
    radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_plane_sphere_polar(po, pn, pu, pv, swap, center, radius, anchors, tol)
        .or_else(|| try_plane_sphere_general(po, pn, pu, pv, swap, center, radius, anchors, tol))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cone_pair(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    swap: bool,
    apex: Point3,
    axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_plane_cone_perpendicular(po, pn, pu, pv, swap, apex, axis, half_angle, anchors, tol)
        .or_else(|| {
            try_plane_cone_oblique_ellipse(
                po, pn, pu, pv, swap, apex, axis, half_angle, anchors, tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_parabola(
                po, pn, pu, pv, swap, apex, axis, half_angle, anchors, tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_hyperbola(
                po, pn, pu, pv, swap, apex, axis, half_angle, anchors, tol,
            )
        })
}
