use crate::geometry::edge_triple::analytic_trig_cone_torus_axial::try_cone_torus_axial;
use crate::geometry::edge_triple::analytic_trig_cyl_torus_axial::try_cyl_torus_axial;
use crate::geometry::edge_triple::analytic_trig_plane_torus::try_plane_torus_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_torus_oblique::try_plane_torus_oblique;
use crate::geometry::edge_triple::analytic_trig_plane_torus_villarceau::try_plane_torus_villarceau;
use crate::geometry::edge_triple::analytic_trig_sphere_torus_axial::try_sphere_torus_axial;
use crate::geometry::edge_triple::analytic_trig_torus_torus_coaxial::try_torus_torus_coaxial;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::EdgeTriple;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_torus_pair(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    swap: bool,
    tc: Point3,
    ta: Vec3,
    t_big: f64,
    t_small: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_plane_torus_perpendicular(po, pn, pu, pv, swap, tc, ta, t_big, t_small, anchors, tol)
        .or_else(|| {
            try_plane_torus_villarceau(po, pn, pu, pv, swap, tc, ta, t_big, t_small, anchors, tol)
        })
        .or_else(|| {
            try_plane_torus_oblique(po, pn, pu, pv, swap, tc, ta, t_big, t_small, anchors, tol)
        })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_sphere_torus_pair(
    center: Point3,
    radius: f64,
    swap: bool,
    tc: Point3,
    ta: Vec3,
    tor_big: f64,
    tor_small: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_sphere_torus_axial(
        center, radius, swap, tc, ta, tor_big, tor_small, anchors, tol,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cone_torus_pair(
    apex: Point3,
    axis: Vec3,
    half_angle: f64,
    swap: bool,
    tc: Point3,
    ta: Vec3,
    tor_big: f64,
    tor_small: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_cone_torus_axial(
        apex, axis, half_angle, swap, tc, ta, tor_big, tor_small, anchors, tol,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cyl_torus_pair(
    co: Point3,
    ca: Vec3,
    cr: f64,
    swap: bool,
    tc: Point3,
    ta: Vec3,
    tor_big: f64,
    tor_small: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_cyl_torus_axial(co, ca, cr, swap, tc, ta, tor_big, tor_small, anchors, tol)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn try_torus_torus_pair(
    ca_t: Point3,
    aa_t: Vec3,
    big_a: f64,
    small_a: f64,
    cb_t: Point3,
    ab_t: Vec3,
    big_b: f64,
    small_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    try_torus_torus_coaxial(
        ca_t, aa_t, big_a, small_a, cb_t, ab_t, big_b, small_b, anchors, tol,
    )
}
