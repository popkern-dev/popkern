//! Shared helpers for 2D corner operations (`fillet_2d`,
//! `chamfer_2d`).
//!
//! All the curve-agnostic topology-surgery primitives that both
//! `fillet_corner_2d` and `chamfer_corner_2d` need live here so
//! each consuming module stays under the 500-line structural cap
//! and the curved-edge support (Arc adjacent edges in the
//! chamfer path) is easy to extend with other curve kinds later.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::error::{BlendError, BlendFailure};

/// Find the parameter `t` on a NURBS curve such that the arc
/// length from `t_start_walk` to `t` equals `target_dist`,
/// walking monotonically toward `t_far_walk`. Uses Newton's
/// method with a bisection fallback — the derivative of arc
/// length w.r.t. the parameter is the tangent's magnitude
/// (speed), which we have via `CurveEval::tangent_at`.
///
/// `t_far_walk` may be less than `t_start_walk` (the walk runs
/// in the negative parameter direction); the helper normalises
/// via a `[0, 1]` walk-progress `s`.
pub(crate) fn invert_nurbs_arc_length(
    nurbs: &NurbsCurve,
    t_start_walk: f64,
    t_far_walk: f64,
    target_dist: f64,
) -> Result<f64, String> {
    use popkern_core::geometry::curve::CurveEval;

    if target_dist < 0.0 {
        return Err(format!(
            "invert_nurbs_arc_length: target_dist {target_dist} must be non-negative"
        ));
    }
    // `length_between(a, b)` is a signed definite integral; the
    // Gauss-Legendre quadrature returns `-|length|` when a > b.
    // We want absolute arc-length here.
    let total_dist = nurbs.length_between(t_start_walk, t_far_walk).abs();
    if target_dist > total_dist + 1.0e-9 {
        return Err(format!(
            "target_dist {target_dist:.9} exceeds walk length {total_dist:.9}"
        ));
    }
    if target_dist <= 1.0e-12 {
        return Ok(t_start_walk);
    }
    if (t_far_walk - t_start_walk).abs() < 1.0e-14 {
        return Err("invert_nurbs_arc_length: walk span is degenerate".into());
    }

    let span = t_far_walk - t_start_walk;
    let mut lo_s = 0.0_f64;
    let mut hi_s = 1.0_f64;
    let mut s = target_dist / total_dist;

    for _ in 0..50 {
        let t = t_start_walk + s * span;
        let acc = nurbs.length_between(t_start_walk, t).abs();
        let error = acc - target_dist;
        if error.abs() < 1.0e-10 {
            return Ok(t);
        }
        if error < 0.0 {
            lo_s = s;
        } else {
            hi_s = s;
        }
        let speed = nurbs.tangent_at(t).length();
        let s_newton = if speed > 1.0e-12 {
            s - error / (speed * span.abs())
        } else {
            f64::NAN
        };
        s = if s_newton.is_finite() && s_newton > lo_s && s_newton < hi_s {
            s_newton
        } else {
            0.5 * (lo_s + hi_s)
        };
        if hi_s - lo_s < 1.0e-14 {
            return Ok(t_start_walk + s * span);
        }
    }
    Err(format!(
        "invert_nurbs_arc_length failed to converge for target={target_dist:.6}"
    ))
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
        CurveGeometry::Line { .. } | CurveGeometry::Arc { .. } | CurveGeometry::Nurbs(_) => {}
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
    let (t_at_corner, t_far, sign) = if d_start <= d_end {
        (curve.t_start, curve.t_end, 1.0)
    } else {
        (curve.t_end, curve.t_start, -1.0)
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
        CurveGeometry::Nurbs(nurbs) => {
            let new_t = invert_nurbs_arc_length(nurbs, t_at_corner, t_far, dist)?;
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
        CurveGeometry::Nurbs(nurbs) => {
            use popkern_core::geometry::curve::CurveEval;
            Ok(nurbs.length_between(curve.t_start, curve.t_end).abs())
        }
        other => Err(format!(
            "unsupported curve type for arc-length query: {other:?}"
        )),
    }
}

/// Signed arc-length distance from the corner (origin of `he_out`)
/// to `target`, measured along the edge's forward direction.
///
/// - **Line edge** — `(target − corner_pos) · tangent`. Positive
///   means `target` lies ahead of the corner along the edge;
///   negative means behind.
/// - **Arc edge** — the angular difference between the corner's
///   parameter and `target`'s angle on the same circle (wrapped
///   into `±π` via [`wrap_nearest_to`]) scaled by the radius,
///   with the sign flipped so that "toward the far end of the
///   edge" is always positive regardless of whether the corner
///   sits at `t_start` or `t_end` in the curve parameterisation.
///
/// Used by `fillet_corner_2d` to validate that a computed
/// tangent point falls within the supporting edge's arc-length
/// budget. The caller compares the return value against
/// [`edge_arc_length`] for the same edge.
pub(crate) fn setback_from_corner(
    brep: &BRep,
    he_out: HalfEdgeHandle,
    target: Point3,
) -> Result<f64, String> {
    use popkern_core::geometry::curve::CurveEval;
    let edge_h = brep.half_edges[he_out].edge;
    let curve = brep.edges[edge_h]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_h:?} has no curve geometry"))?;
    let corner_vh = brep.half_edges[he_out].origin;
    let corner_pos = brep.vertices[corner_vh].position;

    match &curve.geometry {
        CurveGeometry::Line { .. } => {
            let tangent = tangent_at_corner_outgoing(brep, he_out)?;
            Ok((target - corner_pos).dot(&tangent))
        }
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let p_start = curve.geometry.point_at(curve.t_start);
            let p_end = curve.geometry.point_at(curve.t_end);
            let d_start = corner_pos.distance_to(&p_start);
            let d_end = corner_pos.distance_to(&p_end);
            let (t_corner, t_far) = if d_start <= d_end {
                (curve.t_start, curve.t_end)
            } else {
                (curve.t_end, curve.t_start)
            };

            let (u, v) = circle_basis_like_core(*normal);
            let target_rel = target - *center;
            let target_angle_raw = target_rel.dot(&v).atan2(target_rel.dot(&u));
            let target_angle = wrap_nearest_to(t_corner, target_angle_raw);

            let sweep_sign = (t_far - t_corner).signum();
            let progress = (target_angle - t_corner) * sweep_sign;
            Ok(progress * radius)
        }
        CurveGeometry::Nurbs(nurbs) => {
            let p_start = curve.geometry.point_at(curve.t_start);
            let p_end = curve.geometry.point_at(curve.t_end);
            let d_start = corner_pos.distance_to(&p_start);
            let d_end = corner_pos.distance_to(&p_end);
            let (t_corner, t_far) = if d_start <= d_end {
                (curve.t_start, curve.t_end)
            } else {
                (curve.t_end, curve.t_start)
            };
            // Locate the target's parameter via projection.
            let target_t = nurbs
                .closest_point(&target)
                .ok_or_else(|| "closest_point on NURBS returned None".to_string())?;
            let abs_dist = nurbs.length_between(t_corner, target_t).abs();
            // Sign = +1 if target_t lies on the "toward t_far" side of t_corner.
            let walk_dir = (t_far - t_corner).signum();
            let target_dir = (target_t - t_corner).signum();
            let signed = if walk_dir * target_dir >= 0.0 {
                abs_dist
            } else {
                -abs_dist
            };
            Ok(signed)
        }
        other => Err(format!(
            "unsupported curve type for setback query: {other:?}"
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
        CurveGeometry::Nurbs(nurbs) => {
            use popkern_core::geometry::curve::CurveEval;
            // Parameters of the new canonical endpoints on the
            // original NURBS, located via projection.
            let t_from_new = nurbs
                .closest_point(&p_from)
                .ok_or_else(|| "closest_point failed for p_from".to_string())?;
            let t_to_new = nurbs
                .closest_point(&p_to)
                .ok_or_else(|| "closest_point failed for p_to".to_string())?;
            let (tmin, tmax) = if t_from_new < t_to_new {
                (t_from_new, t_to_new)
            } else {
                (t_to_new, t_from_new)
            };
            if (tmax - tmin) < 1.0e-12 {
                return Err(format!(
                    "rebuild_edge_curve_clipped: degenerate NURBS trim range [{tmin}, {tmax}]"
                ));
            }
            // `split_at` rejects splits at the curve's existing
            // parameter boundaries, so we skip redundant splits
            // when a new endpoint coincides with an original one
            // (common when only ONE of the two canonical
            // endpoints was moved by the corner op).
            let (nurbs_start, nurbs_end) = nurbs.parameter_range();
            let tol = 1.0e-12_f64;
            let trimmed = if (tmin - nurbs_start).abs() < tol && (tmax - nurbs_end).abs() < tol {
                nurbs.clone()
            } else if (tmin - nurbs_start).abs() < tol {
                let (left, _) = nurbs
                    .split_at(tmax)
                    .map_err(|e| format!("NURBS split_at({tmax}): {e:?}"))?;
                left
            } else if (tmax - nurbs_end).abs() < tol {
                let (_, right) = nurbs
                    .split_at(tmin)
                    .map_err(|e| format!("NURBS split_at({tmin}): {e:?}"))?;
                right
            } else {
                let (_, right) = nurbs
                    .split_at(tmin)
                    .map_err(|e| format!("NURBS split_at({tmin}): {e:?}"))?;
                let (trimmed, _) = right
                    .split_at(tmax)
                    .map_err(|e| format!("NURBS split_at({tmax}): {e:?}"))?;
                trimmed
            };
            brep.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(trimmed),
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

/// Build a uniform `PartialFailure` wrapping an `InvalidTopology`
/// blend failure for a 2D corner operation. Shared between
/// fillet and chamfer so both emit the same error shape.
///
/// The location is looked up from `brep` when `vertex` resolves;
/// callers that fail outside the context of a specific corner
/// (e.g. chain-walk failures in `edges_2d`) pass
/// `VertexHandle::default()`, which yields `location: None`.
pub(crate) fn topology_failure(brep: &BRep, vertex: VertexHandle, reason: String) -> BlendError {
    let location = brep.vertices.get(vertex).map(|v| v.position);
    BlendError::PartialFailure {
        failed: 1,
        total: 1,
        partial_result: Box::new(BRep::new(popkern_core::tolerance::Tolerance::default())),
        failures: vec![(
            popkern_core::topology::EdgeHandle::default(),
            BlendFailure::InvalidTopology {
                reason: format!("2D corner at {vertex:?}: {reason}"),
                location,
            },
        )],
    }
}
