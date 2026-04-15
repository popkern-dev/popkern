//! Shared helpers for 2D corner operations (`fillet_2d`,
//! `chamfer_2d`).
//!
//! All the curve-agnostic topology-surgery primitives that both
//! `fillet_corner_2d` and `chamfer_corner_2d` need live here so
//! each consuming module stays under the 500-line structural cap
//! and the curved-edge support (Arc adjacent edges in the
//! chamfer path) is easy to extend with other curve kinds later.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::error::{BlendError, BlendFailure};

/// Rebuild `edge`'s curve as a straight `Line` whose direction
/// runs from the canonical half-edge's origin vertex to its
/// twin's origin. The origin-vertex positions are read from the
/// current BRep, so this picks up any vertex moves that happened
/// earlier in the surgery.
pub(crate) fn rebuild_line_curve(brep: &mut BRep, edge: popkern_core::topology::EdgeHandle) {
    let canonical = brep.edges[edge].half_edge;
    let twin = brep.half_edges[canonical].twin;
    let v_from = brep.half_edges[canonical].origin;
    let v_to = brep.half_edges[twin].origin;
    let p_from = brep.vertices[v_from].position;
    let p_to = brep.vertices[v_to].position;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p_from,
            direction: p_to - p_from,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// Unit tangent direction at the corner vertex, walking outward
/// along `he_out` (the half-edge whose `origin` is the corner).
/// Handles both `Line` and `Arc` curves by picking up the
/// appropriate direction from the underlying `CurveEval::tangent_at`.
/// Returns an error for `Nurbs` curves (not yet supported) or
/// missing geometry.
pub(crate) fn tangent_at_corner_outgoing(
    brep: &BRep,
    he_out: HalfEdgeHandle,
) -> Result<Vec3, String> {
    use popkern_core::geometry::curve::CurveEval;
    let edge_h = brep.half_edges[he_out].edge;
    let curve = brep.edges[edge_h]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_h:?} has no curve geometry"))?;
    let corner_vh = brep.half_edges[he_out].origin;
    let corner_pos = brep.vertices[corner_vh].position;

    // Which parameter endpoint sits at the corner? The curve
    // parameterises `point_at(t_start) → point_at(t_end)`. Pick
    // whichever endpoint is closer to the corner vertex.
    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let d_start = corner_pos.distance_to(&p_start);
    let d_end = corner_pos.distance_to(&p_end);

    let (t_at_corner, sign) = if d_start <= d_end {
        (curve.t_start, 1.0)
    } else {
        (curve.t_end, -1.0)
    };

    match &curve.geometry {
        CurveGeometry::Line { .. } | CurveGeometry::Arc { .. } => {}
        other => {
            return Err(format!(
                "unsupported curve type for 2D corner tangent: {other:?}"
            ));
        }
    }

    let tan = curve.geometry.tangent_at(t_at_corner) * sign;
    let unit = tan
        .normalized()
        .ok_or_else(|| "tangent degenerate (zero length)".to_string())?;
    Ok(unit)
}

/// Point on the edge at arc-length distance `dist` from the
/// corner vertex, walking outward along `he_out`. For a `Line`
/// edge this is `corner + tangent * dist`; for an `Arc` it's
/// `arc.point_at(t_at_corner + sign · dist / radius)`. Returns
/// an error for unsupported curve types.
pub(crate) fn point_at_arc_length_from_corner(
    brep: &BRep,
    he_out: HalfEdgeHandle,
    dist: f64,
) -> Result<popkern_core::geometry::Point3, String> {
    use popkern_core::geometry::curve::CurveEval;
    let edge_h = brep.half_edges[he_out].edge;
    let curve = brep.edges[edge_h]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_h:?} has no curve geometry"))?;
    let corner_vh = brep.half_edges[he_out].origin;
    let corner_pos = brep.vertices[corner_vh].position;

    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let d_start = corner_pos.distance_to(&p_start);
    let d_end = corner_pos.distance_to(&p_end);
    let (t_at_corner, sign) = if d_start <= d_end {
        (curve.t_start, 1.0)
    } else {
        (curve.t_end, -1.0)
    };

    match &curve.geometry {
        CurveGeometry::Line { .. } => {
            let tan = tangent_at_corner_outgoing(brep, he_out)?;
            Ok(corner_pos + tan * dist)
        }
        CurveGeometry::Arc { radius, .. } => {
            let dtheta = dist / *radius;
            let new_t = t_at_corner + sign * dtheta;
            Ok(curve.geometry.point_at(new_t))
        }
        other => Err(format!(
            "unsupported curve type for 2D corner arc-length walk: {other:?}"
        )),
    }
}

/// Arc-length of an edge, used by 2D corner operations to check
/// that a requested setback or chamfer distance doesn't exceed
/// the available budget on the adjacent edge.
pub(crate) fn edge_arc_length(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Result<f64, String> {
    let curve = brep.edges[edge]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge:?} has no curve geometry"))?;
    match &curve.geometry {
        CurveGeometry::Line { direction, .. } => {
            let scale = curve.t_end - curve.t_start;
            Ok(direction.length() * scale.abs())
        }
        CurveGeometry::Arc { radius, .. } => {
            let dtheta = (curve.t_end - curve.t_start).abs();
            Ok(radius * dtheta)
        }
        other => Err(format!(
            "unsupported curve type for arc-length query: {other:?}"
        )),
    }
}

/// After a 2D corner operation moves one of an edge's endpoint
/// vertices onto a new tangent-point position, rebuild the
/// edge's curve to span the new endpoint pair. Dispatches on
/// curve kind:
///
/// - `Line` — rebuild as a new `Line` between the canonical
///   half-edge's current endpoint positions.
/// - `Arc` — preserve the original circle (`center`, `normal`,
///   `radius`) and clip `start_angle`/`end_angle` to match the
///   new endpoint angles. Uses the core's `circle_basis` frame
///   (replicated locally via [`circle_basis_like_core`]) so the
///   clipped parameters continue to evaluate correctly.
/// - Any other curve type — returns an error.
///
/// New angles are chosen via [`wrap_nearest_to`] against the old
/// `t_start`/`t_end` so the sweep direction is preserved.
pub(crate) fn rebuild_edge_curve_clipped(
    brep: &mut BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Result<(), String> {
    let canonical = brep.edges[edge].half_edge;
    let twin = brep.half_edges[canonical].twin;
    let v_from = brep.half_edges[canonical].origin;
    let v_to = brep.half_edges[twin].origin;
    let p_from = brep.vertices[v_from].position;
    let p_to = brep.vertices[v_to].position;

    let existing = brep.edges[edge]
        .curve
        .clone()
        .ok_or_else(|| format!("edge {edge:?} has no curve geometry"))?;

    match &existing.geometry {
        CurveGeometry::Line { .. } => {
            brep.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: p_from,
                    direction: p_to - p_from,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
            Ok(())
        }
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let (u_basis, v_basis) = circle_basis_like_core(*normal);
            let rel_from = p_from - *center;
            let rel_to = p_to - *center;
            let t_from_raw = rel_from.dot(&v_basis).atan2(rel_from.dot(&u_basis));
            let t_to_raw = rel_to.dot(&v_basis).atan2(rel_to.dot(&u_basis));
            let t_from_new = wrap_nearest_to(existing.t_start, t_from_raw);
            let t_to_new = wrap_nearest_to(existing.t_end, t_to_raw);
            brep.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center: *center,
                    normal: *normal,
                    radius: *radius,
                    start_angle: t_from_new,
                    end_angle: t_to_new,
                },
                t_start: t_from_new,
                t_end: t_to_new,
            });
            Ok(())
        }
        other => Err(format!(
            "cannot rebuild curve after 2D corner surgery: unsupported type {other:?}"
        )),
    }
}

/// Wrap `t` so it's within `±π` of `anchor` — used when clipping
/// arc angles to preserve the previous sweep direction.
fn wrap_nearest_to(anchor: f64, t: f64) -> f64 {
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

/// Replica of `popkern_core`'s private
/// `analytic_curves::circle_basis`. Used so the angles we compute
/// for clipping an `Arc` match the frame that
/// `CurveGeometry::Arc::point_at` evaluates in at runtime.
pub(crate) fn circle_basis_like_core(normal: Vec3) -> (Vec3, Vec3) {
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

/// Validate that the given edge's curve is a straight `Line`.
/// Used by `fillet_corner_2d` and `chamfer_corner_2d_by_angle`
/// which don't yet support curved adjacent edges (the curved-
/// edge path for basic chamfer goes through
/// [`tangent_at_corner_outgoing`] / [`point_at_arc_length_from_corner`]
/// instead).
pub(crate) fn require_line(
    brep: &BRep,
    edge_h: popkern_core::topology::EdgeHandle,
) -> Result<(), String> {
    let edge = brep
        .edges
        .get(edge_h)
        .ok_or_else(|| format!("edge {edge_h:?} not found"))?;
    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_h:?} has no curve geometry"))?;
    match &curve.geometry {
        CurveGeometry::Line { .. } => Ok(()),
        other => Err(format!(
            "edge {edge_h:?} curve must be a Line, got {other:?}"
        )),
    }
}

/// Build a uniform `PartialFailure` wrapping an `InvalidTopology`
/// blend failure for a 2D corner operation. Shared between
/// fillet and chamfer so both emit the same error shape.
pub(crate) fn topology_failure(vertex: VertexHandle, reason: String) -> BlendError {
    BlendError::PartialFailure {
        failed: 1,
        total: 1,
        partial_result: Box::new(BRep::new(popkern_core::tolerance::Tolerance::default())),
        failures: vec![(
            popkern_core::topology::EdgeHandle::default(),
            BlendFailure::InvalidTopology {
                reason: format!("2D corner at {vertex:?}: {reason}"),
                location: None,
            },
        )],
    }
}
