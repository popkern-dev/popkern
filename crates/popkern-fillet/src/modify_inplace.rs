//! True in-place geometry update for §14.7 — modifies a fillet's
//! radius without any topology surgery. This is the ≥2× speedup
//! path: instead of unfilleting and re-filleting (which rebuilds
//! topology, recomputes pcurves, and re-runs validation), we just
//! move the contact-circle vertices, rebuild the contact-circle
//! edge curves, rebuild the seam edge, and update the fillet face's
//! surface.
//!
//! Round 5d-2 scope: closed plane-cylinder annulus only (the case
//! Round 5c-3 produces). Open plane-plane is a future round.
//!
//! # Why "in place" is faster
//!
//! The remove+recreate path through `modify_fillet_radius` runs:
//!   1. Validate (~µs)
//!   2. `compute_unfillet_plan` or `unfillet_closed_annulus` (~ms — surface intersections, vertex welding)
//!   3. `fillet_edge` end-to-end (~ms — analyze, geometry compute, topology surgery, pcurve recompute, validation)
//!
//! In-place skips steps 2 and most of step 3 — no topology surgery,
//! no pcurve recompute, no second validation. The geometry compute
//! itself is short-circuited by reusing the fact that for a known
//! plane+cylinder pair, the new contact circles are determined
//! analytically by the new radius.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;
use crate::modify_inplace_geometry::{
    determine_fillet_signs, move_cyl_side_vertex, move_plane_side_vertex,
    project_to_cylinder_axis_at_height, project_to_plane_along_axis, rebuild_circle_edge,
};
use crate::unfillet_closed_helpers::{
    Side, classify_support_pair, extract_cylinder_geometry, extract_plane_geometry,
    find_seam_edge_at_vertex_excluding, rebuild_line_curve,
};

/// Modify a closed-fillet annulus's radius in place. Returns a
/// fresh BRep (clone of `brep`) with the contact-circle vertices,
/// edge curves, seam edge, and fillet surface updated for the new
/// radius. No topology changes.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if `fillet_face` isn't a
///   closed-fillet annulus, if the supports aren't a plane-cylinder
///   pair, or if the cylinder axis is parallel to the plane.
pub(crate) fn modify_closed_in_place(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: f64,
) -> Result<BRep, BlendError> {
    let layout = collect_closed_layout(brep, fillet_face)?;
    let plane_geom = extract_plane_geometry(brep, layout.plane_support)?;
    let cyl_geom = extract_cylinder_geometry(brep, layout.cyl_support)?;

    // The fillet sits in the corner where the plane meets the
    // cylinder. Its sign conventions are determined by where the
    // existing contact vertices lie: s_plane = ±1 along n_p (which
    // side of the plane the cylinder vertex is on), s_radial = ±1
    // along the radial direction (whether the plane vertex is
    // outside or inside the cylinder).
    let v_plane_pos_old = brep.vertices[layout.v_plane_side].position;
    let v_cyl_pos_old = brep.vertices[layout.v_cyl_side].position;

    let signs = determine_fillet_signs(&plane_geom, &cyl_geom, &v_plane_pos_old, &v_cyl_pos_old)?;

    // Compute new vertex positions, preserving angular position
    // around the cylinder axis.
    let new_plane_pos = move_plane_side_vertex(
        &v_plane_pos_old,
        &plane_geom,
        &cyl_geom,
        signs.s_radial,
        new_radius,
    );
    let new_cyl_pos = move_cyl_side_vertex(
        &v_cyl_pos_old,
        &plane_geom,
        &cyl_geom,
        signs.s_plane,
        new_radius,
    );

    let mut out = brep.clone();
    out.vertices[layout.v_plane_side].position = new_plane_pos;
    out.vertices[layout.v_cyl_side].position = new_cyl_pos;

    // Rebuild the two contact-circle edges. Each is a closed
    // circular edge centered on the cylinder axis.
    let plane_circle_centre = project_to_plane_along_axis(&plane_geom, &cyl_geom);
    let cyl_circle_centre =
        project_to_cylinder_axis_at_height(&cyl_geom, &plane_geom, signs.s_plane, new_radius);

    rebuild_circle_edge(
        &mut out,
        layout.plane_edge,
        &plane_circle_centre,
        cyl_geom.radius + signs.s_radial * new_radius,
        &cyl_geom.axis,
    );
    rebuild_circle_edge(
        &mut out,
        layout.cyl_edge,
        &cyl_circle_centre,
        cyl_geom.radius,
        &cyl_geom.axis,
    );

    // Rebuild the seam edge — its v_cyl_side endpoint just moved, so
    // its straight-line curve needs to update too.
    let cyl_loop_start = out.loops[brep.half_edges[layout.cyl_twin_h].loop_].half_edge;
    if let Some(seam_edge) = find_seam_edge_at_vertex_excluding(
        &out,
        cyl_loop_start,
        layout.cyl_twin_h,
        layout.v_cyl_side,
    ) {
        rebuild_line_curve(&mut out, seam_edge);
    }

    // Update the fillet face's surface to a new torus that matches
    // the new radius. The torus axis is the cylinder axis; the
    // major radius traces the rolling-ball center curve; the minor
    // radius is the new fillet radius.
    let torus_centre = Point3::new(
        plane_circle_centre.x + signs.s_plane * new_radius * cyl_geom.axis.x,
        plane_circle_centre.y + signs.s_plane * new_radius * cyl_geom.axis.y,
        plane_circle_centre.z + signs.s_plane * new_radius * cyl_geom.axis.z,
    );
    let torus = SurfaceGeometry::Torus {
        center: torus_centre,
        axis: cyl_geom.axis,
        major_radius: cyl_geom.radius + signs.s_radial * new_radius,
        minor_radius: new_radius,
    };
    out.faces[fillet_face].surface = Some(SurfaceOnFace {
        geometry: torus,
        orientation: true,
    });

    Ok(out)
}

/// Topology layout of a closed-fillet annulus, captured by walking
/// the outer + inner loops. All handles refer to the input BRep.
struct ClosedLayout {
    plane_support: FaceHandle,
    cyl_support: FaceHandle,
    cyl_twin_h: HalfEdgeHandle,
    plane_edge: EdgeHandle,
    cyl_edge: EdgeHandle,
    v_plane_side: VertexHandle,
    v_cyl_side: VertexHandle,
}

fn collect_closed_layout(brep: &BRep, fillet_face: FaceHandle) -> Result<ClosedLayout, BlendError> {
    let face_data = brep
        .faces
        .get(fillet_face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "face handle does not exist".to_string(),
        })?;

    if face_data.inner_loops.len() != 1 {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "in-place modify requires exactly one inner loop (closed-fillet annulus)"
                .to_string(),
        });
    }

    let outer_loop_h = face_data.outer_loop;
    let inner_loop_h = face_data.inner_loops[0];

    let outer_he_h = brep.loops[outer_loop_h].half_edge;
    let inner_he_h = brep.loops[inner_loop_h].half_edge;

    let outer_he = brep.half_edges[outer_he_h].clone();
    let inner_he = brep.half_edges[inner_he_h].clone();

    let outer_twin_h = outer_he.twin;
    let inner_twin_h = inner_he.twin;
    let outer_support = brep.loops[brep.half_edges[outer_twin_h].loop_].face;
    let inner_support = brep.loops[brep.half_edges[inner_twin_h].loop_].face;

    let (
        plane_support,
        cyl_support,
        _plane_twin_h,
        cyl_twin_h,
        plane_edge,
        cyl_edge,
        v_plane_side,
        v_cyl_side,
    ) = match classify_support_pair(brep, outer_support, inner_support) {
        Side::Outer => (
            outer_support,
            inner_support,
            outer_twin_h,
            inner_twin_h,
            outer_he.edge,
            inner_he.edge,
            outer_he.origin,
            inner_he.origin,
        ),
        Side::Inner => (
            inner_support,
            outer_support,
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
                value: format!("{outer_support:?}, {inner_support:?}"),
                reason: "could not classify supports — in-place modify needs a plane-cylinder pair"
                    .to_string(),
            });
        }
    };

    Ok(ClosedLayout {
        plane_support,
        cyl_support,
        cyl_twin_h,
        plane_edge,
        cyl_edge,
        v_plane_side,
        v_cyl_side,
    })
}
