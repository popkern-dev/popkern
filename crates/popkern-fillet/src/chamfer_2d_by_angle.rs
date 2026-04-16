//! §5.7 2D corner **distance + angle** chamfer mode —
//! the distance-offset-angle path
//! style. Specifies the chamfer as `dist_a` along the first
//! adjacent edge plus `angle_rad`, the angle at the first tangent
//! point between the *remaining* edge direction and the chamfer
//! line direction toward the second tangent point.
//!
//! Split out of `chamfer_2d.rs` so the curved-support
//! ray-intersection math — which needs its own trio of in-plane
//! geometry helpers — doesn't push the parent file over the
//! 500-line structural cap.
//!
//! ## Curved support behaviour
//!
//! For **Line** adjacent edges the math reduces to the classic
//! law of sines in the triangle `(corner, t1, t2)`. For **Arc**
//! adjacent edges the law of sines no longer applies — the
//! triangle has curved sides from the corner to the two tangent
//! points — so instead:
//!
//! 1. Walk arc-length `dist_a` along `edge_a` from the corner to
//!    get `t1`.
//! 2. Compute the tangent direction at `t1` pointing outward
//!    along `edge_a` (same sense as the outgoing tangent at the
//!    corner, but rotated along the arc).
//! 3. Rotate that tangent by `angle_rad` around the face normal
//!    to get the chamfer ray direction from `t1`.
//! 4. Intersect the ray with `edge_b` in closed form: a 2D
//!    line-line Cramer solve for Line edges, a ray-circle
//!    quadratic for Arc edges.
//! 5. The intersection point is `t2`; verify its arc-length
//!    setback on `edge_b` fits within that edge's budget.
//!
//! The (Line, Line) case evaluates the unified path and still
//! gives the same answer as the law-of-sines formula — a
//! regression-guarded invariant in the tests.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::chamfer_2d::{ChamferPlan, apply_chamfer_plan};
use crate::chamfer_2d_ray::{
    ray_arc_edge_intersection, ray_line_edge_intersection, ray_nurbs_edge_intersection,
};
use crate::corner_2d_helpers::{
    edge_arc_length, invert_nurbs_arc_length, point_at_arc_length_from_corner, setback_from_corner,
    tangent_at_corner_outgoing, topology_failure,
};
use crate::error::BlendError;

/// §2.3 / §5.7 **distance + angle** chamfer mode. Curved
/// supports (Line or Arc) are handled uniformly — see the
/// module-level docs.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] for non-finite /
///   non-positive `dist_a` or `angle_rad`, or an unknown vertex.
/// - [`BlendError::PartialFailure`] wrapping
///   `BlendFailure::InvalidTopology` for geometric / topological
///   failures: corner has ≠ 2 incident edges, edges are
///   colinear, `angle_rad` outside `(corner_angle, π)`, the
///   chamfer ray doesn't hit `edge_b`, or the derived `dist_b`
///   exceeds the edge-b budget.
pub fn chamfer_corner_2d_by_angle(
    brep: &BRep,
    vertex: VertexHandle,
    dist_a: f64,
    angle_rad: f64,
) -> Result<BRep, BlendError> {
    if !dist_a.is_finite() || dist_a <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "dist_a".into(),
            value: format!("{dist_a}"),
            reason: "chamfer distance must be a positive finite value".into(),
        });
    }
    if !angle_rad.is_finite() || angle_rad <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "angle_rad".into(),
            value: format!("{angle_rad}"),
            reason: "chamfer angle must be a positive finite value in radians".into(),
        });
    }
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }

    let plan = chamfer_plan_by_angle(brep, vertex, dist_a, angle_rad)
        .map_err(|reason| topology_failure(brep, vertex, reason))?;

    let mut out = brep.clone();
    apply_chamfer_plan(&mut out, vertex, &plan)
        .map_err(|reason| topology_failure(brep, vertex, reason))?;

    out.attach_planar_pcurves()
        .map_err(|e| topology_failure(brep, vertex, format!("pcurve refresh failed: {e:?}")))?;

    out.validate()
        .map_err(|e| topology_failure(brep, vertex, format!("validation failed: {e:?}")))?;

    Ok(out)
}

fn chamfer_plan_by_angle(
    brep: &BRep,
    vertex: VertexHandle,
    dist_a: f64,
    angle_rad: f64,
) -> Result<ChamferPlan, String> {
    let mut outgoing: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect();
    outgoing.sort();

    if outgoing.len() != 2 {
        return Err(format!(
            "corner vertex must have exactly 2 outgoing half-edges, got {}",
            outgoing.len()
        ));
    }

    let edge_a = brep.half_edges[outgoing[0]].edge;
    let edge_b = brep.half_edges[outgoing[1]].edge;
    if edge_a == edge_b {
        return Err("both outgoing half-edges belong to the same edge (degenerate loop)".into());
    }

    // Tangent directions at the corner, outward along each
    // adjacent edge (Line or Arc both supported).
    let dir_a = tangent_at_corner_outgoing(brep, outgoing[0])?;
    let dir_b = tangent_at_corner_outgoing(brep, outgoing[1])?;

    let len_a = edge_arc_length(brep, edge_a)?;
    let len_b = edge_arc_length(brep, edge_b)?;
    if len_a < 1e-12 || len_b < 1e-12 {
        return Err("adjacent edge has near-zero length".into());
    }

    let cos_corner = dir_a.dot(&dir_b).clamp(-1.0, 1.0);
    let corner_angle = cos_corner.acos();
    if corner_angle < 1.0e-4 {
        return Err("edges are nearly colinear (no meaningful corner)".into());
    }
    if angle_rad <= corner_angle + 1.0e-6 {
        return Err(format!(
            "angle_rad {angle_rad:.6} must exceed corner angle {corner_angle:.6}"
        ));
    }
    if angle_rad >= std::f64::consts::PI - 1.0e-6 {
        return Err(format!("angle_rad {angle_rad:.6} must be < π"));
    }
    if dist_a >= len_a {
        return Err(format!(
            "dist_a {dist_a:.6} ≥ adjacent edge length {len_a:.6}"
        ));
    }

    // First tangent point and the tangent direction at that point.
    let t1 = point_at_arc_length_from_corner(brep, outgoing[0], dist_a)?;
    let tangent_at_t1 = tangent_at_forward_distance(brep, outgoing[0], dist_a)?;

    // Face-plane normal from the corner tangent frame. arc_normal
    // = dir_a × dir_b is oriented so that rotation by `+corner_angle`
    // from `dir_a` around `arc_normal` reaches `dir_b` — the same
    // rotation sense we use for the chamfer below.
    let arc_normal = dir_a
        .cross(&dir_b)
        .normalized()
        .ok_or_else(|| "face normal degenerate".to_string())?;

    // Rotate `tangent_at_t1` by `+angle_rad` around `arc_normal`
    // to get the chamfer ray direction. Rodrigues formula for a
    // vector in the plane (so the parallel component is zero):
    //    v_rot = v · cos(α) + (n × v) · sin(α)
    let cos_a = angle_rad.cos();
    let sin_a = angle_rad.sin();
    let chamfer_dir = tangent_at_t1 * cos_a + arc_normal.cross(&tangent_at_t1) * sin_a;
    let chamfer_dir = chamfer_dir
        .normalized()
        .ok_or_else(|| "chamfer direction degenerate".to_string())?;

    // Ray-edge_b intersection in the face plane.
    let curve_b = brep.edges[edge_b]
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_b:?} has no curve geometry"))?;
    let r = match &curve_b.geometry {
        CurveGeometry::Line { .. } => {
            ray_line_edge_intersection(brep, edge_b, t1, chamfer_dir, arc_normal)?
        }
        CurveGeometry::Arc { .. } => ray_arc_edge_intersection(brep, edge_b, t1, chamfer_dir)?,
        CurveGeometry::Nurbs(_) => {
            ray_nurbs_edge_intersection(brep, edge_b, t1, chamfer_dir, arc_normal)?
        }
        other => {
            return Err(format!(
                "chamfer_by_angle: unsupported curve type on edge_b: {other:?}"
            ));
        }
    };
    let t2 = t1 + chamfer_dir * r;

    // Verify the derived tangent point lies ahead of the corner
    // on edge_b and within its arc-length budget.
    let dist_b = setback_from_corner(brep, outgoing[1], t2)?;
    if dist_b <= 1.0e-9 {
        return Err(format!(
            "derived dist_b {dist_b:.6} is non-positive (chamfer ray doesn't hit edge_b forward)"
        ));
    }
    if dist_b >= len_b {
        return Err(format!(
            "derived dist_b {dist_b:.6} ≥ adjacent edge length {len_b:.6}"
        ));
    }

    Ok(ChamferPlan {
        he_out_a: outgoing[0],
        he_out_b: outgoing[1],
        t1,
        t2,
    })
}

/// Tangent direction at arc-length `dist` along the edge from
/// the corner, pointing outward (same sense as the outgoing
/// tangent at the corner). For a **Line** edge this is constant
/// and equals [`tangent_at_corner_outgoing`]; for an **Arc**
/// edge the tangent rotates by `dist / radius` along the sweep.
fn tangent_at_forward_distance(
    brep: &BRep,
    he_out: HalfEdgeHandle,
    dist: f64,
) -> Result<Vec3, String> {
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
        CurveGeometry::Line { .. } => tangent_at_corner_outgoing(brep, he_out),
        CurveGeometry::Arc { radius, .. } => {
            let dtheta = dist / *radius;
            let new_t = t_at_corner + sign * dtheta;
            let tan = curve.geometry.tangent_at(new_t) * sign;
            tan.normalized()
                .ok_or_else(|| "forward tangent degenerate".to_string())
        }
        CurveGeometry::Nurbs(nurbs) => {
            let new_t = invert_nurbs_arc_length(nurbs, t_at_corner, t_far, dist)?;
            let tan = curve.geometry.tangent_at(new_t) * sign;
            tan.normalized()
                .ok_or_else(|| "forward tangent degenerate".to_string())
        }
        other => Err(format!(
            "unsupported curve type for forward tangent: {other:?}"
        )),
    }
}

#[cfg(test)]
#[path = "chamfer_2d_by_angle_tests.rs"]
mod tests;
