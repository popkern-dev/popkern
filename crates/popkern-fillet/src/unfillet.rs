//! Unfillet operations (§14.2).
//!
//! popkern-fillet exposes three flavours of "remove a fillet":
//!
//! 1. **Snapshot undo** ([`unfillet_via_record`]) — when the caller
//!    has a [`crate::FilletOperationRecord`] from a previous
//!    `fillet_recorded` call, the original input BRep is already
//!    cached and can be restored in O(1) with no geometry work.
//!    This is the right path for *interactive* in-session undo.
//!
//! 2. **Algorithmic plan** ([`compute_unfillet_plan`]) — pure
//!    non-destructive geometry query that returns an
//!    [`UnfilletPlan`] describing what an algorithmic unfillet
//!    *would* produce: which face is the blend, which two faces
//!    are its supports, where the restored sharp edge would lie,
//!    and which face boundaries would need to be re-stitched. This
//!    is the right path for *previewing* an unfillet (e.g. a UI
//!    that wants to highlight the restored edge before the user
//!    commits) and for the future Round 5c full execution.
//!
//! 3. **Algorithmic execution** ([`unfillet_face`]) — actually
//!    runs the topology surgery for plane-plane fillets via the
//!    `kef` / `mef` Euler operators plus targeted vertex
//!    repositioning. **Round 5b ships only the snapshot path and
//!    the plan; the surgery is deferred to Round 5c** because
//!    geometrically-exact unfillet requires vertex welding
//!    primitives that popkern-core does not yet expose, and a
//!    purely topology-correct (but geometrically-approximate)
//!    fallback is worse than no implementation. The function
//!    therefore returns a structured `Err(BlendError::Cancelled)`
//!    today with a clear remediation hint pointing at the two
//!    working alternatives.
//!
//! ## Why no in-place modification yet
//!
//! In-place radius modification (§14.7) is the natural follow-up
//! to algorithmic unfillet: change the radius by unfillet → fillet
//! at the new radius. Until [`unfillet_face`] works for a general
//! BRep, in-place modification is a thin wrapper around
//! "unfillet_via_record + new fillet" and that path already exists
//! through [`crate::FilletOperationRecord::re_execute`].

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::error::BlendError;
use crate::history::FilletOperationRecord;
use crate::identify::{BlendKind, identify_blend};
use crate::unfillet_helpers::{
    collect_end_face_planes, extract_plane, find_torus_support_faces, intersect_line_plane,
    intersect_planes, pick_extreme_endpoints, walk_blend_boundary,
};
use crate::unfillet_surgery::{
    corner_vertex_pair_near, find_he_in_face_with_origin, rebuild_curves_touching,
    remove_degenerate_self_loop_edges,
};

/// A computed plan for an unfillet operation.
///
/// Returned by [`compute_unfillet_plan`]. Pure data — no references
/// to the input BRep are retained — so the plan outlives any
/// borrows of the source model.
#[derive(Debug, Clone)]
pub struct UnfilletPlan {
    /// The blend face that would be removed.
    pub blend_face: FaceHandle,
    /// First support face (one of the two G1-tangent neighbours).
    pub support_a: FaceHandle,
    /// Second support face.
    pub support_b: FaceHandle,
    /// The 3D start position of the restored sharp edge.
    pub restored_edge_start: Point3,
    /// The 3D end position of the restored sharp edge.
    pub restored_edge_end: Point3,
    /// Direction vector of the restored sharp edge (unit length).
    pub restored_edge_direction: Vec3,
    /// The two long G1-tangent boundary half-edges of the blend
    /// face. These are the contact-curve edges that the surgery
    /// would dissolve via `kef`.
    pub long_he_to_a: Option<HalfEdgeHandle>,
    /// Long boundary to support B.
    pub long_he_to_b: Option<HalfEdgeHandle>,
    /// Edges that would need to be removed by the surgery (the
    /// fillet's contact-curve edges plus any short end-cap edges).
    pub edges_to_remove: Vec<EdgeHandle>,
}

/// Compute an unfillet plan without modifying the input BRep.
///
/// Pure function. Plane-plane fillets only — returns
/// `Err(BlendError::InvalidParameter)` for blends on curved support
/// surfaces, which require Round 5c.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] when `blend_face` is not
///   recognised as a blend, when its kind is not
///   [`BlendKind::Fillet`], or when one of its supports is not a
///   plane.
pub fn compute_unfillet_plan(
    brep: &BRep,
    blend_face: FaceHandle,
) -> Result<UnfilletPlan, BlendError> {
    let info = identify_blend(brep, blend_face).ok_or_else(|| BlendError::InvalidParameter {
        parameter: "blend_face".to_string(),
        value: format!("{blend_face:?}"),
        reason: "face is not recognised as a blend".to_string(),
    })?;

    if info.kind != BlendKind::Fillet {
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face.kind".to_string(),
            value: format!("{:?}", info.kind),
            reason: "unfillet currently supports BlendKind::Fillet only".to_string(),
        });
    }

    // For a Torus-surface fillet the canonical support detection is
    // "plane neighbours whose normals are perpendicular to the torus
    // axis". This is more reliable than reusing the
    // tangency-sample-based detection from identify_blend, which can
    // miss G1 boundaries when the contact-curve midpoint doesn't
    // land precisely on the tangent line. Falls back to
    // identify_blend's pair if the geometric search fails.
    let (support_a, support_b) = find_torus_support_faces(brep, blend_face)
        .or(info.support_faces)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "could not determine the support face pair".to_string(),
        })?;

    let plane_a = extract_plane(brep, support_a)?;
    let plane_b = extract_plane(brep, support_b)?;

    // Compute the line of intersection between the two support
    // planes — that line is the locus of the restored sharp edge.
    let (line_origin, line_direction) = intersect_planes(&plane_a, &plane_b)?;

    // Bound the restored edge to the segment that the original
    // fillet covered. The two endpoints are determined by the
    // blend's "end" boundaries — short edges that connect the
    // blend face to the perpendicular end faces (faces C / D in
    // the unit-cube example). Each end face's plane intersected
    // with the restored line gives one endpoint.
    let end_face_planes = collect_end_face_planes(brep, blend_face, support_a, support_b)?;

    let mut endpoints = Vec::new();
    for end_plane in &end_face_planes {
        if let Some(p) = intersect_line_plane(&line_origin, &line_direction, end_plane) {
            endpoints.push(p);
        }
    }

    if endpoints.len() < 2 {
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: format!(
                "found {} valid restored-edge endpoints; need exactly 2",
                endpoints.len()
            ),
        });
    }

    // Pick the two endpoints that are furthest apart along the
    // intersection line (covers the full extent of the original
    // fillet edge).
    let (p_start, p_end) = pick_extreme_endpoints(&endpoints, &line_direction);

    // Walk the blend's outer loop to identify which boundary
    // half-edges are the "long" G1-tangent ones (contact curves
    // to support A and support B respectively) and collect the
    // edges that would need to be removed by the surgery.
    let (long_he_to_a, long_he_to_b, edges_to_remove) =
        walk_blend_boundary(brep, blend_face, support_a, support_b);

    Ok(UnfilletPlan {
        blend_face,
        support_a,
        support_b,
        restored_edge_start: p_start,
        restored_edge_end: p_end,
        restored_edge_direction: line_direction,
        long_he_to_a,
        long_he_to_b,
        edges_to_remove,
    })
}

/// Run an algorithmic unfillet on `blend_face`.
///
/// **Round 5b status:** this entry point currently returns a
/// structured error directing callers to one of the two working
/// paths (snapshot undo via [`unfillet_via_record`], or the plan
/// query via [`compute_unfillet_plan`]). The full topology surgery
/// is Round 5c and requires vertex-welding primitives that
/// popkern-core does not yet expose.
///
/// The function is shipped today so callers can wire their UI
/// against the eventual public shape — the signature is stable
/// and Round 5c will replace the body without changing the
/// signature.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] for the same reasons as
///   [`compute_unfillet_plan`].
/// - A structured "not yet supported" `BlendError::Cancelled`
///   variant when the plan was computable but the surgery isn't
///   implemented yet. The remediation hint points to the working
///   alternatives.
pub fn unfillet_face(brep: &BRep, blend_face: FaceHandle) -> Result<BRep, BlendError> {
    // 1. Plan validates that the face is a blend and gives us the
    //    geometry we need (restored edge endpoints, support pair,
    //    long boundary half-edges).
    let plan = compute_unfillet_plan(brep, blend_face)?;

    let long_a = plan
        .long_he_to_a
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "blend has no detected long boundary to support A".to_string(),
        })?;
    let long_b = plan
        .long_he_to_b
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "blend has no detected long boundary to support B".to_string(),
        })?;

    // 2. Snapshot the second support's plane geometry — we'll
    //    re-apply it after the topology surgery dissolves and re-
    //    splits the merged face.
    let plane_b_geometry = brep
        .faces
        .get(plan.support_b)
        .and_then(|f| f.surface.as_ref())
        .map(|s| (s.geometry.clone(), s.orientation))
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "support_b.surface".to_string(),
            value: format!("{:?}", plan.support_b),
            reason: "support B has no attached surface geometry".to_string(),
        })?;

    // 3. Identify the four "fillet corner vertices" — the vertices
    //    on the blend face's outer loop that are closest to each
    //    restored sharp-edge endpoint. After the KEFs we'll move
    //    each pair to the exact restored corner position and weld
    //    them, leaving one vertex per corner.
    let corner_pair_at_start =
        corner_vertex_pair_near(brep, blend_face, &plan.restored_edge_start)?;
    let corner_pair_at_end = corner_vertex_pair_near(brep, blend_face, &plan.restored_edge_end)?;

    // 4. Topology surgery. We do all mutation on a clone so the
    //    input BRep is unchanged on failure.
    let mut out = brep.clone();

    // 4a. KEF the long boundary on F-A so that F is killed and A's
    //     loop wraps around F's old area. KEF kills the *twin's*
    //     face, so we pass the twin of long_a (which lives on F):
    //     that makes long_a the half-edge on A and its twin's loop
    //     (F) gets killed.
    let twin_of_long_a = out
        .half_edges
        .get(long_a)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "long_he_to_a".to_string(),
            value: format!("{long_a:?}"),
            reason: "half-edge handle does not exist".to_string(),
        })?
        .twin;
    out.kef(twin_of_long_a).map_err(BlendError::Topology)?;

    // 4b. KEF the long boundary on F-B. After step 4a, long_b is
    //     in A's surviving loop (its `loop_` field was rewritten
    //     by the kef). To kill B, we kef long_b directly — its
    //     twin is on B's loop, which gets killed.
    out.kef(long_b).map_err(BlendError::Topology)?;

    // 5. Move the corner vertices to their exact restored positions,
    //    then weld coincident pairs into single representatives.
    //    We move BOTH vertices of each pair onto the same point,
    //    then weld_vertex_pair collapses them.
    if let [v0, v1] = corner_pair_at_start[..] {
        out.vertices[v0].position = plan.restored_edge_start;
        out.vertices[v1].position = plan.restored_edge_start;
        out.weld_vertex_pair(v0, v1).map_err(BlendError::Topology)?;
    }
    if let [v0, v1] = corner_pair_at_end[..] {
        out.vertices[v0].position = plan.restored_edge_end;
        out.vertices[v1].position = plan.restored_edge_end;
        out.weld_vertex_pair(v0, v1).map_err(BlendError::Topology)?;
    }

    // 5b. Welding leaves behind "degenerate" edges where both
    //     endpoints collapsed to the same canonical vertex (the
    //     short end-cap arcs of the fillet, whose two endpoints
    //     just got welded together). These contribute 1 to the
    //     edge count without representing real topology, breaking
    //     the Euler-Poincaré invariant. Remove them.
    remove_degenerate_self_loop_edges(&mut out);

    // 5c. Welding moved vertices but did NOT update the curve
    //     geometry of edges that touch them. The result is a
    //     loop where each half-edge's curve start/end points still
    //     reference the OLD vertex positions, producing a
    //     `gap = fillet_radius` validation failure. Walk the
    //     surviving edges and rebuild any curve to fit the new
    //     endpoint positions as a straight line (correct for the
    //     plane-plane scope).
    let welded_vertices = [corner_pair_at_start[0], corner_pair_at_end[0]];
    rebuild_curves_touching(&mut out, &welded_vertices);

    // 6. After the welds, the two surviving canonical vertices are
    //    `corner_pair_at_start[0]` and `corner_pair_at_end[0]`
    //    (weld_vertex_pair keeps the primary, removes the secondary).
    //    Find half-edges in the merged face whose origin is each
    //    surviving vertex — these are the MEF endpoints.
    //
    //    The merged face is support_a — the kef chain killed F and
    //    B but kept A as the surviving face holding the merged loop.
    let merged_face = plan.support_a;
    let surviving_start = corner_pair_at_start[0];
    let surviving_end = corner_pair_at_end[0];

    let he_at_start =
        find_he_in_face_with_origin(&out, merged_face, surviving_start).ok_or_else(|| {
            BlendError::InvalidParameter {
                parameter: "merged_loop".to_string(),
                value: format!("{merged_face:?}"),
                reason: "could not find a half-edge in the merged loop with the welded \
                     start vertex as its origin"
                    .to_string(),
            }
        })?;
    let he_at_end =
        find_he_in_face_with_origin(&out, merged_face, surviving_end).ok_or_else(|| {
            BlendError::InvalidParameter {
                parameter: "merged_loop".to_string(),
                value: format!("{merged_face:?}"),
                reason: "could not find a half-edge in the merged loop with the welded \
                     end vertex as its origin"
                    .to_string(),
            }
        })?;

    // 7. MEF splits the merged face along the new edge. The new
    //    face inherits the parent face's surface — A's plane.
    let mef_result = out
        .mef(he_at_start, he_at_end)
        .map_err(BlendError::Topology)?;

    // 8. The MEF result has TWO faces sharing the new edge: the
    //    parent (still A) and `mef_result.face` (also A's plane,
    //    inherited). Reassign B's plane to the new face.
    use popkern_core::geometry::surface::SurfaceOnFace;
    out.faces[mef_result.face].surface = Some(SurfaceOnFace {
        geometry: plane_b_geometry.0,
        orientation: plane_b_geometry.1,
    });

    Ok(out)
}

/// Snapshot-based unfillet — undo a fillet by reading the input
/// BRep cached inside the [`FilletOperationRecord`].
///
/// This is the right path for *interactive* in-session undo. It is
/// O(1), exact, and works regardless of surface type — the
/// algorithmic limitations of [`unfillet_face`] do not apply.
///
/// Returns a reference into the record. Clone the result if the
/// caller needs an owned BRep.
pub fn unfillet_via_record(record: &FilletOperationRecord) -> &BRep {
    record.undo_state()
}

// Internal helpers live in `unfillet_helpers` — kept separate so
// each source file stays under the 500-line structural cap.

#[cfg(test)]
#[path = "unfillet_tests.rs"]
mod tests;
