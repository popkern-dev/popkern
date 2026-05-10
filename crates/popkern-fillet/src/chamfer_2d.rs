//! §5.7 2D corner chamfers for sheet-body planar wires.
//!
//! A **2D corner chamfer** replaces a sharp convex corner of a
//! planar polygon with a straight-line segment ("bevel") between
//! two tangent points on the adjacent edges. Given a sheet body
//! (produced by e.g. `BRep::sheet_quad`) and a corner vertex,
//! [`chamfer_corner_2d`] removes the vertex and inserts a new
//! line edge whose endpoints sit at caller-specified distances
//! along each adjacent edge.
//!
//! This is the companion to `fillet_2d::fillet_corner_2d`; the two
//! share tangent-point math and topology surgery via the helpers
//! in [`crate::fillet_2d`]. The only real differences are:
//!
//! - **Tangent placement.** A fillet sets `setback = r / tan(θ/2)`
//!   to make the arc tangent to both edges; a chamfer just uses
//!   caller-provided distances `dist_a` and `dist_b` along each
//!   edge from the corner.
//! - **New-edge curve.** A fillet's new edge carries a
//!   [`CurveGeometry::Arc`]; a chamfer's new edge is a
//!   [`CurveGeometry::Line`] from one tangent point to the other.
//!
//! With equal distances (`dist_a == dist_b`) on a 90° corner, a
//! chamfer produces the standard 45° bevel of the 2D CAD idiom.
//!
//! # Example
//!
//! ```no_run
//! use popkern_core::tolerance::Tolerance;
//! use popkern_core::topology::BRep;
//! use popkern_fillet::chamfer_corner_2d;
//!
//! let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, Tolerance::default()).unwrap();
//! let corner = sheet.vertices.keys().next().unwrap();
//! // 0.1 × 0.1 equal chamfer.
//! let _chamfered = chamfer_corner_2d(&sheet, corner, 0.1, 0.1).unwrap();
//! ```

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, Edge, HalfEdge, HalfEdgeHandle, VertexHandle};

use crate::corner_2d_helpers::{
    edge_arc_length, point_at_arc_length_from_corner, rebuild_edge_curve_clipped,
    tangent_at_corner_outgoing, topology_failure,
};
use crate::error::BlendError;

/// Re-export of [`crate::chamfer_2d_by_angle::chamfer_corner_2d_by_angle`]
/// — the distance+angle chamfer mode lives in a sibling module
/// to keep `chamfer_2d.rs` under the 500-line cap, but we expose
/// it through this module so the `popkern_fillet::chamfer_2d`
/// namespace remains the single entry point for 2D chamfer ops.
pub use crate::chamfer_2d_by_angle::chamfer_corner_2d_by_angle;

/// Chamfer a single sharp corner of a sheet body by distances
/// `dist_a` along the first adjacent edge and `dist_b` along the
/// second.
///
/// Returns a new BRep with the corner vertex replaced by two
/// tangent points joined by a straight-line edge. The input
/// `brep` is unchanged.
///
/// **Distance ordering.** The "first" and "second" adjacent edges
/// are picked by the half-edge slot-map ordering, which is stable
/// but not otherwise meaningful. For symmetric chamfers
/// (`dist_a == dist_b`) the ordering is immaterial. For asymmetric
/// chamfers, callers that need a specific orientation should
/// either use the symmetric form or swap distances if the result
/// comes out mirrored.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] — either distance is NaN,
///   infinite, or non-positive; or `vertex` isn't a handle in `brep`.
/// - [`BlendError::PartialFailure`] wrapping
///   `BlendFailure::InvalidTopology` for:
///   - corner has ≠ 2 incident edges
///   - an incident edge is not straight (not a `Line` curve)
///   - a distance exceeds the corresponding edge length
///   - post-surgery validation fails
pub fn chamfer_corner_2d(
    brep: &BRep,
    vertex: VertexHandle,
    dist_a: f64,
    dist_b: f64,
) -> Result<BRep, BlendError> {
    if !dist_a.is_finite() || dist_a <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "dist_a".into(),
            value: format!("{dist_a}"),
            reason: "chamfer distance must be a positive finite value".into(),
        });
    }
    if !dist_b.is_finite() || dist_b <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "dist_b".into(),
            value: format!("{dist_b}"),
            reason: "chamfer distance must be a positive finite value".into(),
        });
    }
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }

    let plan = chamfer_plan(brep, vertex, dist_a, dist_b)
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

/// Batch variant of [`chamfer_corner_2d`] — chamfer every
/// `(vertex, dist_a, dist_b)` spec in sequence, chaining the BRep
/// through each call.
///
/// Collision detection is the same implicit
/// "setback-beats-remaining-edge-length" check as in the single-
/// corner path. All-or-nothing semantics: any failure aborts the
/// batch and discards the partial result.
pub fn chamfer_corners_2d(
    brep: &BRep,
    specs: &[(VertexHandle, f64, f64)],
) -> Result<BRep, BlendError> {
    let mut current = brep.clone();
    for &(vertex, dist_a, dist_b) in specs {
        current = chamfer_corner_2d(&current, vertex, dist_a, dist_b)?;
    }
    Ok(current)
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ChamferPlan {
    pub(crate) he_out_a: HalfEdgeHandle,
    pub(crate) he_out_b: HalfEdgeHandle,
    pub(crate) t1: Point3,
    pub(crate) t2: Point3,
}

fn chamfer_plan(
    brep: &BRep,
    vertex: VertexHandle,
    dist_a: f64,
    dist_b: f64,
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

    // Curve-type gate: `chamfer_corner_2d` supports `Line` and
    // `Arc` adjacent edges. `tangent_at_corner_outgoing` and
    // `point_at_arc_length_from_corner` reject unsupported types
    // (Circle, Nurbs) via their own error paths.
    let dir_a = tangent_at_corner_outgoing(brep, outgoing[0])?;
    let dir_b = tangent_at_corner_outgoing(brep, outgoing[1])?;

    let len_a = edge_arc_length(brep, edge_a)?;
    let len_b = edge_arc_length(brep, edge_b)?;
    if len_a < 1e-12 || len_b < 1e-12 {
        return Err("adjacent edge has near-zero length".into());
    }

    if dist_a >= len_a {
        return Err(format!(
            "dist_a {dist_a:.6} ≥ adjacent edge arc-length {len_a:.6}"
        ));
    }
    if dist_b >= len_b {
        return Err(format!(
            "dist_b {dist_b:.6} ≥ adjacent edge arc-length {len_b:.6}"
        ));
    }

    // Check the tangent directions at the corner aren't colinear
    // — a chamfer on a colinear corner is pointless and the new
    // line would overlap the support edges. For Line edges this
    // is the straight-line direction; for Arc edges it's the
    // tangent at the corner endpoint.
    let cos_angle = dir_a.dot(&dir_b).clamp(-1.0, 1.0);
    let angle = cos_angle.acos();
    if angle < 1.0e-4 {
        return Err("edges are nearly colinear (no meaningful corner)".into());
    }

    let t1 = point_at_arc_length_from_corner(brep, outgoing[0], dist_a)?;
    let t2 = point_at_arc_length_from_corner(brep, outgoing[1], dist_b)?;

    Ok(ChamferPlan {
        he_out_a: outgoing[0],
        he_out_b: outgoing[1],
        t1,
        t2,
    })
}

pub(crate) fn apply_chamfer_plan(
    brep: &mut BRep,
    corner: VertexHandle,
    plan: &ChamferPlan,
) -> Result<(), String> {
    // 1. Insert the two new tangent-point vertices.
    let v_t1 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: plan.t1,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v_t2 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: plan.t2,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    // 2. Reparent the outgoing half-edges at the corner onto the
    //    new tangent vertices and rebuild their edge curves. For
    //    Line edges this is a straight Line rebuild; for Arc
    //    edges it clips the arc's start/end angles to the new
    //    tangent-point angles on the original circle.
    brep.half_edges[plan.he_out_a].origin = v_t1;
    brep.half_edges[plan.he_out_b].origin = v_t2;
    let edge_a_h = brep.half_edges[plan.he_out_a].edge;
    let edge_b_h = brep.half_edges[plan.he_out_b].edge;
    rebuild_edge_curve_clipped(brep, edge_a_h)?;
    rebuild_edge_curve_clipped(brep, edge_b_h)?;

    // 3. Create the new chamfer edge with a straight-line curve.
    let line = CurveGeometry::Line {
        origin: plan.t1,
        direction: plan.t2 - plan.t1,
    };
    let chamfer_edge = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(CurveOnEdge {
            geometry: line,
            t_start: 0.0,
            t_end: 1.0,
        }),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });

    // 4. Create the two new half-edges for the chamfer edge and
    //    splice them into the two face loops. The canonical
    //    half-edge of the new edge has origin `v_t1` so the Line
    //    curve's natural direction matches.
    let loop_front = brep.half_edges[plan.he_out_a].loop_;
    let loop_back = brep.half_edges[plan.he_out_b].loop_;
    let prev_front = brep.half_edges[plan.he_out_a].prev;
    let prev_back = brep.half_edges[plan.he_out_b].prev;

    let new_front = brep.half_edges.insert(HalfEdge {
        origin: v_t2,
        twin: Default::default(),
        next: plan.he_out_a,
        prev: prev_front,
        edge: chamfer_edge,
        loop_: loop_front,
        pcurve: None,
    });
    let new_back = brep.half_edges.insert(HalfEdge {
        origin: v_t1,
        twin: new_front,
        next: plan.he_out_b,
        prev: prev_back,
        edge: chamfer_edge,
        loop_: loop_back,
        pcurve: None,
    });
    brep.half_edges[new_front].twin = new_back;
    // Canonical half-edge matches the Line curve direction
    // `t1 → t2` (origin = `v_t1` ⇒ `new_back`).
    brep.edges[chamfer_edge].half_edge = new_back;

    // 5. Rewire next/prev for both splices.
    brep.half_edges[prev_front].next = new_front;
    brep.half_edges[plan.he_out_a].prev = new_front;
    brep.half_edges[prev_back].next = new_back;
    brep.half_edges[plan.he_out_b].prev = new_back;

    // 6. Update the new tangent vertices' outgoing handles.
    brep.vertices[v_t1].outgoing = plan.he_out_a;
    brep.vertices[v_t2].outgoing = plan.he_out_b;

    // 7. Delete the corner vertex.
    brep.vertices.remove(corner);

    Ok(())
}

#[cfg(test)]
#[path = "chamfer_2d_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "chamfer_2d_curved_tests.rs"]
mod curved_tests;

#[cfg(test)]
#[path = "chamfer_2d_nurbs_tests.rs"]
mod nurbs_tests;
