//! Closed-fillet unfillet — the inverse of
//! [`crate::construct::closed_fillet::apply_closed_fillet`].
//!
//! Given an annulus fillet face created by the Round 5c-3
//! closed-fillet construction (one outer loop + one inner loop,
//! both single-element self-closed half-edges), this module
//! reverses the surgery:
//!
//! 1. Identify the two contact-circle edges (one per boundary
//!    loop) and the two support faces they bridge.
//! 2. Pick the seam-bearing support — its rim vertex stays put
//!    and gets moved back to the original sharp-rim position
//!    (the plane-cylinder intersection circle).
//! 3. The other support's contact-circle vertex gets deleted
//!    (it was a fresh vertex created by the construction path).
//! 4. The two contact-circle edges merge back into a single
//!    closed rim edge — its half-edges live in the two support
//!    face loops.
//! 5. The seam edge gets its curve geometry rebuilt to span back
//!    to the new rim position.
//! 6. The fillet face and its loops are deleted.
//!
//! Round 5c-4 scope: handles plane-cylinder closed-fillet annuli
//! (the case the Round 5c-3 path produces). The rim is computed
//! analytically from the support surfaces — Plane and Cylinder
//! geometries are required.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BlendError;
use crate::unfillet_closed_helpers::{
    Side, build_rim_circle_nurbs, classify_support_pair, compute_rim_vertex_position,
    extract_cylinder_geometry, extract_plane_geometry, find_seam_edge_at_vertex_excluding,
    plane_cylinder_intersection_circle, rebuild_line_curve,
};

/// Reverse the closed-fillet construction for `fillet_face`.
///
/// Returns a fresh BRep (clone of `brep`) with the annulus face
/// removed, the original sharp-rim edge restored, the seam edge
/// extended back to the new rim height, and the redundant
/// contact-circle vertex deleted.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if `fillet_face` is not an
///   annulus, has no surface, or its supports aren't the
///   plane-cylinder pair the closed-fillet path produces.
pub fn unfillet_closed_annulus(brep: &BRep, fillet_face: FaceHandle) -> Result<BRep, BlendError> {
    // 1. Validate the fillet face is an annulus.
    let face_data = brep
        .faces
        .get(fillet_face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "face handle does not exist".to_string(),
        })?;

    if face_data.inner_loops.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "face is not an annulus (has no inner loops)".to_string(),
        });
    }

    // The closed-fillet construction produces exactly one inner loop.
    let inner_loop_h = face_data.inner_loops[0];
    let outer_loop_h = face_data.outer_loop;

    // 2. Walk each loop to find its (single) half-edge → edge → twin → support face.
    let outer_he_h = brep
        .loops
        .get(outer_loop_h)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "fillet_face.outer_loop".to_string(),
            value: format!("{outer_loop_h:?}"),
            reason: "outer loop missing".to_string(),
        })?
        .half_edge;
    let inner_he_h = brep
        .loops
        .get(inner_loop_h)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "fillet_face.inner_loop".to_string(),
            value: format!("{inner_loop_h:?}"),
            reason: "inner loop missing".to_string(),
        })?
        .half_edge;

    let outer_he =
        brep.half_edges
            .get(outer_he_h)
            .cloned()
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "outer_loop.half_edge".to_string(),
                value: format!("{outer_he_h:?}"),
                reason: "outer loop's half-edge missing".to_string(),
            })?;
    let inner_he =
        brep.half_edges
            .get(inner_he_h)
            .cloned()
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "inner_loop.half_edge".to_string(),
                value: format!("{inner_he_h:?}"),
                reason: "inner loop's half-edge missing".to_string(),
            })?;

    let outer_twin_h = outer_he.twin;
    let inner_twin_h = inner_he.twin;
    let outer_support_face = brep
        .loops
        .get(brep.half_edges[outer_twin_h].loop_)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "outer_twin.loop".to_string(),
            value: format!("{outer_twin_h:?}"),
            reason: "outer twin's loop missing".to_string(),
        })?
        .face;
    let inner_support_face = brep
        .loops
        .get(brep.half_edges[inner_twin_h].loop_)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "inner_twin.loop".to_string(),
            value: format!("{inner_twin_h:?}"),
            reason: "inner twin's loop missing".to_string(),
        })?
        .face;

    // 3. Identify which support is the cylinder (has a seam edge,
    //    multi-half-edge loop) and which is the plane (single-
    //    half-edge inner loop on a face with a hole).
    let (
        plane_support,
        cyl_support,
        plane_twin_h,
        cyl_twin_h,
        plane_edge,
        cyl_edge,
        v_plane_side,
        v_cyl_side,
    ) = match classify_support_pair(brep, outer_support_face, inner_support_face) {
        Side::Outer => (
            outer_support_face,
            inner_support_face,
            outer_twin_h,
            inner_twin_h,
            outer_he.edge,
            inner_he.edge,
            outer_he.origin,
            inner_he.origin,
        ),
        Side::Inner => (
            inner_support_face,
            outer_support_face,
            inner_twin_h,
            outer_twin_h,
            inner_he.edge,
            outer_he.edge,
            inner_he.origin,
            outer_he.origin,
        ),
        Side::Ambiguous => {
            return Err(BlendError::InvalidParameter {
                parameter: "fillet_face.supports".to_string(),
                value: format!("{outer_support_face:?}, {inner_support_face:?}"),
                reason: "could not classify which support is the plane and which is the cylinder"
                    .to_string(),
            });
        }
    };

    // 4. Compute the original rim position from the plane and
    //    cylinder geometries. The rim is a circle on the plane
    //    at the cylinder radius, centered on the cylinder axis
    //    projected onto the plane.
    let plane_geom = extract_plane_geometry(brep, plane_support)?;
    let cyl_geom = extract_cylinder_geometry(brep, cyl_support)?;
    let (rim_center, rim_radius) = plane_cylinder_intersection_circle(&plane_geom, &cyl_geom)?;

    // The vertex that stays is the cylinder-side vertex — it's
    // attached to the seam edge whose other endpoint is on the
    // cylinder's bottom rim. Move it to a point on the new rim
    // circle (preserve its existing angular position around the
    // cylinder axis).
    let v_keep = v_cyl_side;
    let v_delete = v_plane_side;

    let v_keep_pos_old = brep.vertices[v_keep].position;
    let new_rim_vertex_pos =
        compute_rim_vertex_position(&v_keep_pos_old, &rim_center, &cyl_geom, rim_radius);

    let mut out = brep.clone();

    // 5. Move v_keep to the new rim position.
    out.vertices[v_keep].position = new_rim_vertex_pos;

    // 6. Find and rebuild the seam edge curve so its endpoint at
    //    v_keep is at the new rim height.
    if let Some(seam_edge) = find_seam_edge_at_vertex_excluding(
        &out,
        out.loops[brep.half_edges[cyl_twin_h].loop_].half_edge,
        cyl_twin_h,
        v_keep,
    ) {
        rebuild_line_curve(&mut out, seam_edge);
    }

    // 7. Update the cyl_edge curve geometry to be the new rim
    //    circle (NURBS approximation centered on cylinder axis).
    let rim_curve = build_rim_circle_nurbs(&rim_center, rim_radius, &cyl_geom.axis);
    out.edges[cyl_edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(rim_curve),
        t_start: 0.0,
        t_end: 1.0,
    });

    // 8. Re-twin the support half-edges to share `cyl_edge` as the
    //    canonical rim edge. The plane-side half-edge currently
    //    references `plane_edge` and is twinned with the inner
    //    fillet-face half-edge (which we're about to delete).
    //    Rewire it to reference cyl_edge and twin with the
    //    cyl-side half-edge.
    out.half_edges[plane_twin_h].twin = cyl_twin_h;
    out.half_edges[plane_twin_h].edge = cyl_edge;
    out.half_edges[plane_twin_h].origin = v_keep;
    out.half_edges[cyl_twin_h].twin = plane_twin_h;

    // 9. Delete the fillet face's half-edges (outer_he, inner_he),
    //    the obsolete plane_edge, the obsolete v_delete, and the
    //    fillet face + its loops.
    out.half_edges.remove(outer_he_h);
    out.half_edges.remove(inner_he_h);
    out.edges.remove(plane_edge);
    out.vertices.remove(v_delete);
    out.loops.remove(outer_loop_h);
    out.loops.remove(inner_loop_h);
    out.faces.remove(fillet_face);

    // 10. Update v_keep's outgoing pointer if needed.
    if !out.half_edges.contains_key(out.vertices[v_keep].outgoing) {
        out.vertices[v_keep].outgoing = cyl_twin_h;
    }

    Ok(out)
}

#[cfg(test)]
#[path = "unfillet_closed_tests.rs"]
mod tests;
