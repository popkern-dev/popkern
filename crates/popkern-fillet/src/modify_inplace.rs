//! True in-place geometry update for §14.7 — modifies a fillet's
//! radius without any topology surgery. The ≥2× speedup path:
//! instead of unfillet+refillet, just move the contact vertices,
//! rebuild touching edges, and update the fillet face's surface.
//!
//! - Round 5d-2: closed plane-cylinder annulus
//!   ([`modify_closed_in_place`])
//! - Round 5d-3: open plane-plane fillet
//!   ([`modify_open_in_place`])

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;
use crate::modify_inplace_geometry::{
    determine_fillet_signs, move_cyl_side_vertex, move_plane_side_vertex,
    project_to_cylinder_axis_at_height, project_to_plane_along_axis, rebuild_circle_edge,
};
use crate::unfillet::{UnfilletPlan, compute_unfillet_plan};
use crate::unfillet_closed_helpers::{
    Side, classify_support_pair, extract_cylinder_geometry, extract_plane_geometry,
    find_seam_edge_at_vertex_excluding, rebuild_line_curve,
};
use crate::unfillet_helpers::extract_plane;
use crate::unfillet_surgery::rebuild_curves_touching;

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

/// Modify an open plane-plane fillet's radius in place. The fillet
/// face must have a single outer loop (no inner loops) bounding two
/// long contact lines (parallel to the original sharp edge) and two
/// short end-cap lines.
///
/// The implementation:
///   1. Uses [`compute_unfillet_plan`] to capture the support faces
///      and the original sharp-edge endpoints.
///   2. Calls [`corner_vertex_pair_near`] at each endpoint to find
///      the two corner vertices on the fillet face boundary.
///   3. Each corner is at distance `R` from the sharp-edge endpoint
///      along a perpendicular direction in its support plane —
///      moves it to distance `new_radius` along the same direction.
///   4. Calls [`rebuild_curves_touching`] to refit every edge whose
///      endpoint just moved (the long contact edges, the short
///      end-cap edges, and the perpendicular cube-edge segments).
///   5. Reassigns the fillet face's surface to a new analytic torus
///      with `minor_radius = new_radius`.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if `fillet_face` isn't an
///   identifiable plane-plane fillet, the corners can't be found,
///   or the support planes can't be extracted.
pub(crate) fn modify_open_in_place(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: f64,
) -> Result<BRep, BlendError> {
    modify_open_in_place_general(brep, fillet_face, new_radius, new_radius)
}

/// Variable-radius variant of [`modify_open_in_place`] (Round 5h).
/// Moves the start-end corners to distance `r_start` from the
/// start sharp-edge endpoint and the end-end corners to distance
/// `r_end` from the end sharp-edge endpoint. The fillet face's
/// surface is approximated as a Torus with the average minor radius;
/// callers needing exact NURBS-sweep surface fidelity should use
/// the remove+recreate path via [`crate::modify_fillet_radius`]
/// with a `FilletRadius::Variable` spec.
pub(crate) fn modify_open_in_place_variable(
    brep: &BRep,
    fillet_face: FaceHandle,
    r_start: f64,
    r_end: f64,
) -> Result<BRep, BlendError> {
    modify_open_in_place_general(brep, fillet_face, r_start, r_end)
}

fn modify_open_in_place_general(
    brep: &BRep,
    fillet_face: FaceHandle,
    r_start: f64,
    r_end: f64,
) -> Result<BRep, BlendError> {
    let plan = compute_unfillet_plan(brep, fillet_face)?;
    let plane_a = extract_plane(brep, plan.support_a)?;
    let plane_b = extract_plane(brep, plan.support_b)?;

    // Identify the 4 fillet corners by analytic direction matching.
    // For each end of the sharp edge we use the local radius
    // (`r_start` at the start, `r_end` at the end). The constant
    // case is the special case `r_start == r_end`.
    let edge_dir = {
        let d = plan.restored_edge_direction;
        let len = d.length();
        if len > 1e-12 {
            Vec3::new(d.x / len, d.y / len, d.z / len)
        } else {
            return Err(BlendError::InvalidParameter {
                parameter: "fillet_face".to_string(),
                value: format!("{fillet_face:?}"),
                reason: "restored edge direction is degenerate".to_string(),
            });
        }
    };

    let perp_a = perpendicular_in_plane(&plane_a.normal, &edge_dir);
    let perp_b = perpendicular_in_plane(&plane_b.normal, &edge_dir);

    let loop_vertices = collect_outer_loop_vertices(brep, fillet_face);
    let signed_perp_a = pick_signed_perpendicular(&perp_a, &plan, &loop_vertices, brep);
    let signed_perp_b = pick_signed_perpendicular(&perp_b, &plan, &loop_vertices, brep);

    let mut out = brep.clone();
    let mut welded: Vec<VertexHandle> = Vec::with_capacity(4);

    for (endpoint, local_radius) in [
        (&plan.restored_edge_start, r_start),
        (&plan.restored_edge_end, r_end),
    ] {
        for &perp in [&signed_perp_a, &signed_perp_b].iter() {
            let Some(v) = find_corner_vertex(brep, &loop_vertices, endpoint, perp) else {
                return Err(BlendError::InvalidParameter {
                    parameter: "fillet_face".to_string(),
                    value: format!("{fillet_face:?}"),
                    reason: "could not locate a fillet corner along an expected perpendicular"
                        .to_string(),
                });
            };
            out.vertices[v].position = Point3::new(
                endpoint.x + perp.x * local_radius,
                endpoint.y + perp.y * local_radius,
                endpoint.z + perp.z * local_radius,
            );
            welded.push(v);
        }
    }

    rebuild_curves_touching(&mut out, &welded);

    // Reassign the fillet face surface. For constant radius this is
    // an exact analytic torus. For variable radius the surface is
    // an approximation — a Torus with the *average* minor radius —
    // which renders correctly enough for visualisation but is not
    // geometrically exact. Use remove+recreate for exact fidelity.
    let avg_radius = (r_start + r_end) * 0.5;
    let cos_angle = plane_a.normal.dot(&plane_b.normal).clamp(-1.0, 1.0);
    let dihedral = cos_angle.acos();
    let half_angle = dihedral * 0.5;
    let sin_half = half_angle.sin();
    if sin_half.abs() > 1e-10 {
        let axis_raw = plane_a.normal.cross(&plane_b.normal);
        let axis_len = axis_raw.length();
        let axis = if axis_len > 1e-12 {
            Vec3::new(
                axis_raw.x / axis_len,
                axis_raw.y / axis_len,
                axis_raw.z / axis_len,
            )
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let major = avg_radius / sin_half;
        let centre = Point3::new(
            (plan.restored_edge_start.x + plan.restored_edge_end.x) * 0.5,
            (plan.restored_edge_start.y + plan.restored_edge_end.y) * 0.5,
            (plan.restored_edge_start.z + plan.restored_edge_end.z) * 0.5,
        );
        out.faces[fillet_face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Torus {
                center: centre,
                axis,
                major_radius: major,
                minor_radius: avg_radius,
            },
            orientation: true,
        });
    }

    Ok(out)
}

/// Compute a unit perpendicular vector that lies in `plane_normal`'s
/// plane and is perpendicular to `edge_dir`. Sign is arbitrary.
fn perpendicular_in_plane(plane_normal: &Vec3, edge_dir: &Vec3) -> Vec3 {
    let raw = plane_normal.cross(edge_dir);
    let len = raw.length();
    if len > 1e-12 {
        Vec3::new(raw.x / len, raw.y / len, raw.z / len)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

/// Collect the unique vertices in `fillet_face`'s outer loop.
fn collect_outer_loop_vertices(brep: &BRep, fillet_face: FaceHandle) -> Vec<VertexHandle> {
    let mut out = Vec::new();
    let lh = brep.faces[fillet_face].outer_loop;
    let start_he = brep.loops[lh].half_edge;
    let mut he = start_he;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let h = &brep.half_edges[he];
        if !out.contains(&h.origin) {
            out.push(h.origin);
        }
        he = h.next;
        if he == start_he {
            break;
        }
    }
    out
}

/// Resolve the sign of `perp` so that it points toward the side of
/// the sharp edge where the fillet sits. Picks whichever sign has
/// more loop vertices' offsets aligning positively with it.
fn pick_signed_perpendicular(
    perp: &Vec3,
    plan: &UnfilletPlan,
    loop_vertices: &[VertexHandle],
    brep: &BRep,
) -> Vec3 {
    let mut pos_count = 0;
    let mut neg_count = 0;
    for &v in loop_vertices {
        let pos = brep.vertices[v].position;
        let d_start = (pos.x - plan.restored_edge_start.x).powi(2)
            + (pos.y - plan.restored_edge_start.y).powi(2)
            + (pos.z - plan.restored_edge_start.z).powi(2);
        let d_end = (pos.x - plan.restored_edge_end.x).powi(2)
            + (pos.y - plan.restored_edge_end.y).powi(2)
            + (pos.z - plan.restored_edge_end.z).powi(2);
        let endpoint = if d_start <= d_end {
            plan.restored_edge_start
        } else {
            plan.restored_edge_end
        };
        let offset = Vec3::new(pos.x - endpoint.x, pos.y - endpoint.y, pos.z - endpoint.z);
        if offset.length() < 1e-9 {
            continue;
        }
        let aligned = offset.dot(perp);
        if aligned > 1e-6 {
            pos_count += 1;
        } else if aligned < -1e-6 {
            neg_count += 1;
        }
    }
    if pos_count >= neg_count {
        *perp
    } else {
        Vec3::new(-perp.x, -perp.y, -perp.z)
    }
}

/// Find the loop vertex whose direction from `endpoint` best aligns
/// with `perp`. Cosine-similarity scored.
fn find_corner_vertex(
    brep: &BRep,
    loop_vertices: &[VertexHandle],
    endpoint: &Point3,
    perp: &Vec3,
) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    for &v in loop_vertices {
        let pos = brep.vertices[v].position;
        let offset = Vec3::new(pos.x - endpoint.x, pos.y - endpoint.y, pos.z - endpoint.z);
        let len = offset.length();
        if len < 1e-9 {
            continue;
        }
        let cos_sim = offset.dot(perp) / len;
        if cos_sim < 0.5 {
            continue;
        }
        match &best {
            Some((_, prev)) if cos_sim <= *prev => {}
            _ => best = Some((v, cos_sim)),
        }
    }
    best.map(|(v, _)| v)
}

#[cfg(test)]
#[path = "modify_inplace_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "modify_equivalence_tests.rs"]
mod equivalence_tests;
