//! True in-place geometry update for §14.7 — modifies a fillet's
//! radius without any topology surgery. The ≥2× speedup path:
//! instead of unfillet+refillet, just move the contact vertices,
//! rebuild touching edges, and update the fillet face's surface.
//!
//! - closed plane-cylinder annulus
//!   ([`modify_closed_in_place`])
//! - open plane-plane fillet
//!   ([`modify_open_in_place`])

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use popkern_math::PortableMath;

use crate::error::BlendError;
use crate::modify_inplace_geometry::{
    collect_outer_loop_vertices, determine_fillet_signs, find_corner_vertex, move_cyl_side_vertex,
    move_plane_side_vertex, perpendicular_in_plane, pick_signed_perpendicular,
    project_to_cylinder_axis_at_height, project_to_plane_along_axis, rebuild_circle_edge,
};
use crate::modify_inplace_sweep::{
    VariableSweepInputs, build_parametric_plane_plane_sweep, build_variable_plane_plane_sweep,
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
pub fn modify_closed_in_place(
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
    crate::construct::strip_orientation::reconcile_strip_orientation_boundary(
        &mut out,
        fillet_face,
    );

    Ok(out)
}

pub use crate::modify_inplace_closed_sweep::modify_closed_in_place_bypoints;

/// Topology layout of a closed-fillet annulus, captured by walking
/// the outer + inner loops. All handles refer to the input BRep.
pub(crate) struct ClosedLayout {
    pub(crate) plane_support: FaceHandle,
    pub(crate) cyl_support: FaceHandle,
    pub(crate) cyl_twin_h: HalfEdgeHandle,
    pub(crate) plane_edge: EdgeHandle,
    pub(crate) cyl_edge: EdgeHandle,
    pub(crate) v_plane_side: VertexHandle,
    pub(crate) v_cyl_side: VertexHandle,
}

pub(crate) fn collect_closed_layout(
    brep: &BRep,
    fillet_face: FaceHandle,
) -> Result<ClosedLayout, BlendError> {
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
pub fn modify_open_in_place(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: f64,
) -> Result<BRep, BlendError> {
    modify_open_in_place_general(brep, fillet_face, new_radius, new_radius, None)
}

/// Variable-radius variant of [`modify_open_in_place`] .
/// Moves the start-end corners to distance `r_start` from the
/// start sharp-edge endpoint and the end-end corners to distance
/// `r_end` from the end sharp-edge endpoint. The fillet face's
/// surface is approximated as a Torus with the average minor radius;
/// callers needing exact NURBS-sweep surface fidelity should use
/// the remove+recreate path via [`crate::modify_fillet_radius`]
/// with a `FilletRadius::Variable` spec.
pub fn modify_open_in_place_variable(
    brep: &BRep,
    fillet_face: FaceHandle,
    r_start: f64,
    r_end: f64,
) -> Result<BRep, BlendError> {
    modify_open_in_place_general(brep, fillet_face, r_start, r_end, None)
}

/// ByPoints variant of the in-place open path. The
/// corners are positioned by the radii at parameters 0 and 1
/// (matching the perp-split convention used for `Constant` and
/// `Variable`), and the fillet surface is built as a NURBS sweep
/// driven by a piecewise-linear closure over the *full* points
/// list — so intermediate radius points shape the surface.
pub fn modify_open_in_place_bypoints(
    brep: &BRep,
    fillet_face: FaceHandle,
    points: &[(f64, f64)],
) -> Result<BRep, BlendError> {
    if points.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "points".to_string(),
            value: "empty".to_string(),
            reason: "ByPoints radius profile must have at least one point".to_string(),
        });
    }
    let r_start = crate::center_curve::bypoints_local_radius(points, 0.0);
    let r_end = crate::center_curve::bypoints_local_radius(points, 1.0);
    let pts: Vec<(f64, f64)> = points.to_vec();
    let closure: Box<dyn Fn(f64) -> f64> =
        Box::new(move |t| crate::center_curve::bypoints_local_radius(&pts, t));
    modify_open_in_place_general(brep, fillet_face, r_start, r_end, Some(closure))
}

fn modify_open_in_place_general(
    brep: &BRep,
    fillet_face: FaceHandle,
    r_start: f64,
    r_end: f64,
    radius_closure: Option<Box<dyn Fn(f64) -> f64>>,
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

    // Reassign the fillet face surface.
    //   - ByPoints (closure provided) → tensor-product NURBS sweep
    //     driven by the piecewise-linear closure over the full
    //     points list. (ByPoints in-place path.)
    //   - Variable (no closure, r_start ≠ r_end) → linear-interp
    //     NURBS sweep. (variable in-place path.)
    //   - Constant (no closure, r_start == r_end) → analytic torus
    //     matching `try_plane_plane_analytic`. (open plane-plane in-place path.)
    let sweep_inputs = VariableSweepInputs {
        p_start: plan.restored_edge_start,
        p_end: plan.restored_edge_end,
        plane_a_normal: plane_a.normal,
        plane_b_normal: plane_b.normal,
        signed_perp_a,
        signed_perp_b,
        r_start,
        r_end,
        n_sections: 8,
        n_arc_samples: 8,
    };
    let new_surface = if let Some(radius_at) = radius_closure.as_deref() {
        build_parametric_plane_plane_sweep(&sweep_inputs, radius_at)
    } else if (r_start - r_end).abs() > 1e-12 {
        build_variable_plane_plane_sweep(&sweep_inputs)
    } else {
        constant_plane_plane_cylinder(
            &plane_a.normal,
            &plane_b.normal,
            &signed_perp_a,
            &signed_perp_b,
            &plan,
            r_start,
        )
    };
    out.faces[fillet_face].surface = Some(SurfaceOnFace {
        geometry: new_surface,
        orientation: true,
    });
    crate::construct::strip_orientation::reconcile_strip_orientation_boundary(
        &mut out,
        fillet_face,
    );

    // The swap retires the old strip chart — re-author the face's
    // pcurves against the NEW surface before any settlement reads
    // them (PK-179; see `modify_inplace_chart`).
    crate::modify_inplace_chart::reauthor_strip_chart_pcurves(&mut out, fillet_face)?;

    // PK-148: the chain-end mid-arc corner vertices (settled onto
    // the OLD terminal arc by the construction pipeline) must land
    // on the NEW radius' arc, with stub curves and cap pcurves
    // re-authored to match. Same settlement as the fresh-blend
    // pipeline, scoped to the modified strip. Uses the body's own
    // tolerance — the in-place API carries no caller tolerance.
    let body_tolerance = out.tolerance;
    let _ = crate::construct::terminal_arc_trim::settle_terminal_arc_faces(
        &mut out,
        &[fillet_face],
        &body_tolerance,
    );

    Ok(out)
}

/// Build the analytic Cylinder surface for a constant-radius
/// plane-plane fillet, mirroring `try_plane_plane_analytic` in
/// `fillet_surface_analytic.rs`. The center curve for a straight
/// dihedral edge is a line, so the fillet envelope is a cylinder
/// whose axis is the edge direction and radius equals the fillet
/// radius. The rolling-ball centre sits on the MATERIAL-side
/// bisector (the signed perpendiculars, mirroring the variable
/// sweep) at distance `r / sin(half_angle)` from the sharp edge —
/// distance `r` from EACH support plane, not `r` along the
/// bisector (the old `bisector · r` placement put the axis
/// `r·(1−1/√2)` too close to a perpendicular dihedral's edge, so
/// the authored cylinder cut the supports instead of touching
/// them).
fn constant_plane_plane_cylinder(
    n_a: &Vec3,
    n_b: &Vec3,
    signed_perp_a: &Vec3,
    signed_perp_b: &Vec3,
    plan: &UnfilletPlan,
    radius: f64,
) -> SurfaceGeometry {
    let cos_angle = n_a.dot(n_b).clamp(-1.0, 1.0);
    let half_angle = cos_angle.pacos() * 0.5;
    let sin_half = half_angle.psin().max(1e-12);
    let bisector = (*signed_perp_a + *signed_perp_b)
        .normalized()
        .unwrap_or(*n_a);
    let origin = plan.restored_edge_start + bisector * (radius / sin_half);
    let chord = plan.restored_edge_end - plan.restored_edge_start;
    let axis = chord
        .normalized()
        .unwrap_or_else(|| n_a.cross(n_b).normalized().unwrap_or(*n_a));
    SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }
}

#[cfg(test)]
#[path = "modify_inplace_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "modify_inplace_variable_tests.rs"]
mod variable_tests;

#[cfg(test)]
#[path = "modify_equivalence_tests.rs"]
mod equivalence_tests;
