//! Private helpers shared by the `sheet_*` constructors in
//! `sheet.rs`. Split out so `sheet.rs` stays under the 500-line
//! structural cap while still hosting the public constructor
//! surface and its doc-comments.

use crate::error::CoreError;
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;
use popkern_math::PortableMath;

/// Locate the edge whose canonical half-edge runs from `p_from`
/// to `p_to` (within linear tolerance) and replace its `Line`
/// curve with an `Arc` on the specified circle.
///
/// Preserves the canonical half-edge's direction: the arc
/// parameterisation `(t_start, t_end)` matches the canonical
/// endpoint order.
pub(super) fn rewrite_edge_as_arc(
    brep: &mut BRep,
    p_from: Point3,
    p_to: Point3,
    center: Point3,
    normal: Vec3,
    radius: f64,
) -> Result<(), CoreError> {
    let edge_h = find_edge_between_points(brep, p_from, p_to, 1.0e-9).ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "arched_triangle: cannot find edge between {p_from:?} and {p_to:?}"
        ))
    })?;

    // The canonical half-edge's origin gives the "first" vertex
    // of the edge in parameter-range order. We use that to
    // compute the angle at t_start.
    let canonical_he = brep.edges[edge_h].half_edge;
    let canonical_origin_vh = brep.half_edges[canonical_he].origin;
    let canonical_twin_h = brep.half_edges[canonical_he].twin;
    let canonical_twin_origin_vh = brep.half_edges[canonical_twin_h].origin;
    let canonical_origin_pos = brep.vertices[canonical_origin_vh].position;
    let canonical_twin_origin_pos = brep.vertices[canonical_twin_origin_vh].position;

    let (u_basis, v_basis) = circle_basis_local(normal);
    let t_start = angle_in_basis(canonical_origin_pos - center, u_basis, v_basis);
    let t_end_raw = angle_in_basis(canonical_twin_origin_pos - center, u_basis, v_basis);
    // Choose the "short" arc by ensuring `t_end` is within `π`
    // of `t_start` (add or subtract `2π` as needed).
    let t_end = wrap_nearest(t_start, t_end_raw);

    let arc = crate::geometry::curve::CurveGeometry::Arc {
        center,
        normal,
        radius,
        start_angle: t_start,
        end_angle: t_end,
    };
    brep.edges[edge_h].curve = Some(crate::geometry::curve::CurveOnEdge {
        geometry: arc,
        t_start,
        t_end,
    });
    Ok(())
}

fn find_edge_between_points(
    brep: &BRep,
    p_from: Point3,
    p_to: Point3,
    tol_eps: f64,
) -> Option<crate::topology::EdgeHandle> {
    for (eh, edge) in &brep.edges {
        let he_h = edge.half_edge;
        let twin_h = brep.half_edges[he_h].twin;
        let a = brep.vertices[brep.half_edges[he_h].origin].position;
        let b = brep.vertices[brep.half_edges[twin_h].origin].position;
        let matches = (a.distance_to(&p_from) < tol_eps && b.distance_to(&p_to) < tol_eps)
            || (a.distance_to(&p_to) < tol_eps && b.distance_to(&p_from) < tol_eps);
        if matches {
            return Some(eh);
        }
    }
    None
}

fn angle_in_basis(rel: Vec3, u: Vec3, v: Vec3) -> f64 {
    rel.dot(&v).patan2(rel.dot(&u))
}

/// Wrap `t` so it's within `±π` of `anchor`. Used to pick the
/// short-arc parameter for an otherwise-ambiguous endpoint angle.
fn wrap_nearest(anchor: f64, t: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let mut result = t;
    while result - anchor > PI {
        result -= TAU;
    }
    while result - anchor < -PI {
        result += TAU;
    }
    result
}

/// Local replica of the (private) `circle_basis` helper in
/// `analytic_curves`, used here so the arc parameterisation we
/// compute matches the frame `CurveGeometry::Arc::point_at`
/// evaluates in at runtime.
fn circle_basis_local(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = normal.cross(&axis);
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}
