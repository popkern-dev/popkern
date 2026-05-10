//! §5.7 2D corner fillets for sheet-body planar wires.
//!
//! A **2D corner fillet** replaces a sharp convex corner of a
//! planar polygon with a circular arc of the requested radius.
//! Given a sheet body (produced by e.g. `BRep::sheet_quad`) and a
//! vertex where two straight-line edges meet, [`fillet_corner_2d`]
//! builds two tangent points (one on each adjacent edge) and a
//! new arc edge that spans between them, then deletes the original
//! corner vertex.
//!
//! This is the "sheet body" half of §5.7 Planar Wire Fillets:
//! A commercial kernel's 2D profile fillet does the analogous
//! operation for 2D wires. Since popkern represents sheet bodies
//! as two-face paper-thin shells (front + backside of the same
//! plane, opposite orientations), the surgery runs **twice** —
//! once on the front face's outer loop and once on the backside's
//! — keeping the manifold pairing intact.
//!
//! ## Algorithm
//!
//! 1. Find the two half-edges whose origin is the corner vertex.
//!    They live on distinct edges `edge_a` and `edge_b`.
//! 2. Verify both edges carry `CurveGeometry::Line`. 2D corner
//!    fillets only operate on straight-line edges today; curved
//!    boundary edges are a future extension.
//! 3. Compute unit directions `dir_a`, `dir_b` from the corner
//!    toward each edge's other endpoint.
//! 4. Compute the interior half-angle `θ/2` from
//!    `cos(θ) = dir_a · dir_b`. Reject if `θ ≈ 0` (colinear) or
//!    `θ ≈ π` (near-straight, no meaningful corner).
//! 5. Tangent setback = `r / tan(θ/2)`; arc centre offset = `r /
//!    sin(θ/2)` along the interior bisector.
//! 6. Reject if the setback exceeds either adjacent edge's length
//!    (`RadiusTooLarge`).
//! 7. Create new vertices `T1` (on edge_a) and `T2` (on edge_b)
//!    at the tangent points.
//! 8. Update the two half-edges with origin = corner: their
//!    origins now become `T1` and `T2` respectively, and the
//!    underlying edge curves get their origin-side endpoint
//!    shortened.
//! 9. Create a new arc edge with a `CurveGeometry::Arc` and two
//!    fresh half-edges (one per face loop). Splice both into
//!    their loops between the predecessor of the old corner
//!    half-edge and the updated half-edge.
//! 10. Delete the corner vertex (no longer referenced).
//! 11. Re-attach planar pcurves on the affected faces.
//! 12. Validate the result; on failure return `InvalidTopology`.
//!
//! The caller receives a fresh [`BRep`] on success (the input is
//! never mutated).

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, Edge, HalfEdge, VertexHandle};

use crate::corner_2d_helpers::{rebuild_edge_curve_clipped, topology_failure};
use crate::error::BlendError;
use crate::fillet_2d_plan::{CornerPlan, corner_plan};

/// Fillet a single sharp corner of a sheet body.
///
/// Returns a new BRep with the corner replaced by a circular arc
/// of the requested radius. The input `brep` is unchanged.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] — `radius ≤ 0`, `radius` is
///   NaN or infinite, or `vertex` isn't a handle in `brep`.
/// - [`BlendError::PartialFailure`] with a `BlendFailure::InvalidTopology`
///   wrapping one of:
///   - corner has ≠ 2 incident edges
///   - an incident edge is not straight (not a `Line` curve)
///   - edges are colinear / near-parallel (no meaningful corner)
///   - radius is too large for the shorter adjacent edge
///   - post-surgery validation fails
pub fn fillet_corner_2d(
    brep: &BRep,
    vertex: VertexHandle,
    radius: f64,
) -> Result<BRep, BlendError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{radius}"),
            reason: "radius must be a positive finite value".into(),
        });
    }
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }

    let plan = corner_plan(brep, vertex, radius)
        .map_err(|reason| topology_failure(brep, vertex, reason))?;

    let mut out = brep.clone();
    apply_corner_plan(&mut out, vertex, &plan)
        .map_err(|reason| topology_failure(brep, vertex, reason))?;

    out.attach_planar_pcurves()
        .map_err(|e| topology_failure(brep, vertex, format!("pcurve refresh failed: {e:?}")))?;

    out.validate()
        .map_err(|e| topology_failure(brep, vertex, format!("validation failed: {e:?}")))?;

    Ok(out)
}

/// Batch variant of [`fillet_corner_2d`] — fillet every
/// `(vertex, radius)` spec in sequence, chaining the BRep through
/// each call so later specs see the topology produced by earlier
/// ones.
///
/// **Collision detection** is implicit: after filleting corner A,
/// its adjacent edges are shortened; filleting a neighbouring
/// corner B then sees the reduced budget via the existing
/// "setback ≥ edge length" check in [`fillet_corner_2d`]. Callers
/// that want predictable behaviour across sibling corners should
/// order the specs deterministically and pick radii that leave
/// the shared edges with strictly positive length after both
/// operations.
///
/// The caller's input `brep` is unchanged on both success and
/// failure.
///
/// # Errors
///
/// Returns the first per-corner error encountered, wrapped
/// unchanged. On any failure the partial BRep is discarded —
/// this is an all-or-nothing batch, not a best-effort one.
pub fn fillet_corners_2d(brep: &BRep, specs: &[(VertexHandle, f64)]) -> Result<BRep, BlendError> {
    let mut current = brep.clone();
    for &(vertex, radius) in specs {
        current = fillet_corner_2d(&current, vertex, radius)?;
    }
    Ok(current)
}

fn apply_corner_plan(
    brep: &mut BRep,
    corner: VertexHandle,
    plan: &CornerPlan,
) -> Result<(), String> {
    // 1. Insert the new tangent-point vertices.
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

    // 2. Move the origin of the outgoing half-edges onto the new
    //    tangent vertices. `he_out_a` now starts at `t1`, and its
    //    underlying edge's curve is a line from `t1` to `v_a_far`.
    //    Same for `he_out_b` with `t2` and `v_b_far`.
    brep.half_edges[plan.he_out_a].origin = v_t1;
    brep.half_edges[plan.he_out_b].origin = v_t2;

    // Rebuild the two support edges' curves after their endpoints
    // moved onto the tangent points. For `Line` edges this writes
    // a fresh straight-line curve; for `Arc` edges it clips the
    // arc's start/end angles to match the new endpoint on the
    // same circle, preserving the original sweep direction.
    let edge_a_h = brep.half_edges[plan.he_out_a].edge;
    let edge_b_h = brep.half_edges[plan.he_out_b].edge;
    rebuild_edge_curve_clipped(brep, edge_a_h)?;
    rebuild_edge_curve_clipped(brep, edge_b_h)?;

    // 3. Create the arc edge and its two half-edges (one per face
    //    loop). `new_front` splices into the loop containing
    //    `he_out_a` (just before it) and points from `t2` → `t1`.
    //    `new_back` splices into the loop containing `he_out_b`
    //    and points from `t1` → `t2`.
    let arc_geometry = CurveGeometry::Arc {
        center: plan.arc_center,
        normal: plan.arc_normal,
        radius: plan.radius,
        start_angle: plan.start_angle,
        end_angle: plan.end_angle,
    };
    let arc_edge = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(CurveOnEdge {
            geometry: arc_geometry,
            t_start: plan.start_angle,
            t_end: plan.end_angle,
        }),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });

    // Gather the loop and predecessor data we need BEFORE inserting
    // the new half-edges, because both belong to the same half-edge
    // slot-map and we mustn't borrow while inserting.
    let loop_front = brep.half_edges[plan.he_out_a].loop_;
    let loop_back = brep.half_edges[plan.he_out_b].loop_;
    let prev_front = brep.half_edges[plan.he_out_a].prev;
    let prev_back = brep.half_edges[plan.he_out_b].prev;

    // Front half-edge: origin `t2`, destination `t1`. Inserted
    // before `he_out_a` in the front loop.
    let new_front = brep.half_edges.insert(HalfEdge {
        origin: v_t2,
        twin: Default::default(),
        next: plan.he_out_a,
        prev: prev_front,
        edge: arc_edge,
        loop_: loop_front,
        pcurve: None,
    });
    // Back half-edge: origin `t1`, destination `t2`. Inserted
    // before `he_out_b` in the back loop.
    let new_back = brep.half_edges.insert(HalfEdge {
        origin: v_t1,
        twin: new_front,
        next: plan.he_out_b,
        prev: prev_back,
        edge: arc_edge,
        loop_: loop_back,
        pcurve: None,
    });
    brep.half_edges[new_front].twin = new_back;
    // Canonical half-edge direction matches the arc curve's
    // parameterisation: start_angle → t1, end_angle → t2. The
    // half-edge that starts at t1 is `new_back` (origin = v_t1).
    brep.edges[arc_edge].half_edge = new_back;

    // 4. Rewire next/prev across both splices.
    brep.half_edges[prev_front].next = new_front;
    brep.half_edges[plan.he_out_a].prev = new_front;
    brep.half_edges[prev_back].next = new_back;
    brep.half_edges[plan.he_out_b].prev = new_back;

    // 5. Make sure each vertex's `outgoing` handle points at a
    //    valid half-edge. The new tangent vertices emit the
    //    modified outgoing half-edges; the old corner vertex loses
    //    its references.
    brep.vertices[v_t1].outgoing = plan.he_out_a;
    brep.vertices[v_t2].outgoing = plan.he_out_b;

    // 6. Any loop anchor that pointed at `he_out_a` / `he_out_b`
    //    is still valid (the handles are unchanged), but anchors
    //    that happened to point at the front/back loop via the
    //    predecessor are also valid. Nothing to do.

    // 7. Delete the corner vertex. It is no longer referenced —
    //    any surviving reference would indicate a wiring bug and
    //    validation will catch it.
    brep.vertices.remove(corner);

    Ok(())
}

#[cfg(test)]
#[path = "fillet_2d_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "fillet_2d_curved_tests.rs"]
mod curved_tests;

#[cfg(test)]
#[path = "fillet_2d_nurbs_tests.rs"]
mod nurbs_tests;
