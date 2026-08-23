//! In-place chamfer distance modification (PK-280) — the chamfer
//! sibling of `modify_inplace.rs`'s §14.7 fillet path.
//!
//! A plane-plane chamfer strip is fully analytic: its corners sit at
//! `restored_edge_endpoint + perp_side * distance_side` exactly like
//! a fillet's contact corners (with the chamfer distance in place of
//! the radius), and the strip surface is the PLANE through the moved
//! corners. Changing the distance therefore needs no topology
//! surgery: move the four corner vertices, rebuild the touched
//! curves, swap the strip plane, and re-author the strip chart —
//! the same skeleton `modify_open_in_place` uses, measured there at
//! 35–2693× over remove+recreate.
//!
//! Scope (stage 1 of PK-280): OPEN plane-plane chamfer strips with
//! `Equal` or `TwoDistance` specs. Closed rim chamfers (cone
//! frustums) and `DistanceAngle` decline with a structured reason
//! naming the remove+recreate route — recorded follow-up in the
//! PK-280 ledger.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::error::BlendError;
use crate::identify::BlendKind;
use crate::modify_inplace_geometry::{
    collect_outer_loop_vertices, find_corner_vertex, perpendicular_in_plane,
    pick_signed_perpendicular,
};
use crate::types::ChamferKind;
use crate::unfillet::compute_blend_strip_plan_with_supports;
use crate::unfillet_helpers::extract_plane;
use crate::unfillet_surgery::rebuild_curves_touching;

/// Modify a chamfer strip's distances in place.
///
/// `chamfer_face` must be an OPEN plane-plane chamfer strip (the
/// common box-edge bevel). The new spec is interpreted exactly as
/// [`crate::chamfer_edges`] interprets it:
/// [`ChamferKind::Equal`] applies one distance on both supports;
/// [`ChamferKind::TwoDistance`] applies `dist_1` on `face_1` (which
/// must be one of the strip's two support faces) and `dist_2` on the
/// other support.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] when the face is not an open
///   plane-plane chamfer strip, when a `TwoDistance` reference face
///   is not one of the supports, when a distance is not strictly
///   positive, or for the (declined) `DistanceAngle` spec.
pub fn modify_chamfer_distance(
    brep: &BRep,
    chamfer_face: FaceHandle,
    new_kind: ChamferKind,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    // Kind gate FIRST: support derivation reads the strip's plane and
    // would report a confusing surface error on e.g. a cylinder
    // fillet strip.
    let info = crate::identify::identify_blend(brep, chamfer_face).ok_or_else(|| {
        BlendError::InvalidParameter {
            parameter: "chamfer_face".to_string(),
            value: format!("{chamfer_face:?}"),
            reason: "face is not recognised as a blend".to_string(),
        }
    })?;
    if info.kind != BlendKind::Chamfer {
        return Err(BlendError::InvalidParameter {
            parameter: "chamfer_face.kind".to_string(),
            value: format!("{:?}", info.kind),
            reason: "this operation supports Chamfer blend faces only".to_string(),
        });
    }
    let supports = chamfer_supports_by_dihedral(brep, chamfer_face)?;
    let plan = compute_blend_strip_plan_with_supports(
        brep,
        chamfer_face,
        &[BlendKind::Chamfer],
        Some(supports),
    )?;
    let plane_a = extract_plane(brep, plan.support_a)?;
    let plane_b = extract_plane(brep, plan.support_b)?;

    let (dist_a, dist_b) = match &new_kind {
        ChamferKind::Equal(d) => (*d, *d),
        ChamferKind::TwoDistance {
            face_1,
            dist_1,
            dist_2,
        } => {
            if *face_1 == plan.support_a {
                (*dist_1, *dist_2)
            } else if *face_1 == plan.support_b {
                (*dist_2, *dist_1)
            } else {
                return Err(BlendError::InvalidParameter {
                    parameter: "new_kind.face_1".to_string(),
                    value: format!("{face_1:?}"),
                    reason: format!(
                        "TwoDistance reference face must be one of the strip's supports \
                         ({:?} / {:?})",
                        plan.support_a, plan.support_b
                    ),
                });
            }
        }
        ChamferKind::DistanceAngle { .. } => {
            return Err(BlendError::InvalidParameter {
                parameter: "new_kind".to_string(),
                value: "DistanceAngle".to_string(),
                reason: "DistanceAngle in-place modification is not derivable from the strip \
                         alone (the angle re-orients the bevel against one reference face); \
                         remove and re-chamfer via chamfer_edges"
                    .to_string(),
            });
        }
    };
    for (label, d) in [("dist_a", dist_a), ("dist_b", dist_b)] {
        if !(d.is_finite() && d > 0.0) {
            return Err(BlendError::InvalidParameter {
                parameter: label.to_string(),
                value: format!("{d}"),
                reason: "chamfer distance must be strictly positive and finite".to_string(),
            });
        }
    }

    let edge_dir =
        plan.restored_edge_direction
            .normalized()
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "chamfer_face".to_string(),
                value: format!("{chamfer_face:?}"),
                reason: "restored edge direction is degenerate".to_string(),
            })?;

    let perp_a = perpendicular_in_plane(&plane_a.normal, &edge_dir);
    let perp_b = perpendicular_in_plane(&plane_b.normal, &edge_dir);
    let loop_vertices = collect_outer_loop_vertices(brep, chamfer_face);
    let signed_perp_a = pick_signed_perpendicular(&perp_a, &plan, &loop_vertices, brep);
    let signed_perp_b = pick_signed_perpendicular(&perp_b, &plan, &loop_vertices, brep);

    // OLD per-side distances, measured from the existing contact
    // corners (needed to parameterise every boundary vertex across
    // the old cross-section, mid-bevel end vertices included — the
    // measured cube bevel carries a mid vertex on each end edge, and
    // moving only the four contact corners left them off the new
    // plane).
    let corner_a = find_corner_vertex(
        brep,
        &loop_vertices,
        &plan.restored_edge_start,
        &signed_perp_a,
    );
    let corner_b = find_corner_vertex(
        brep,
        &loop_vertices,
        &plan.restored_edge_start,
        &signed_perp_b,
    );
    let (Some(ca), Some(cb)) = (corner_a, corner_b) else {
        return Err(BlendError::InvalidParameter {
            parameter: "chamfer_face".to_string(),
            value: format!("{chamfer_face:?}"),
            reason: "could not locate the strip's contact corners".to_string(),
        });
    };
    let old_dist = |v: VertexHandle, perp: &Vec3| -> f64 {
        (brep.vertices[v].position - plan.restored_edge_start).dot(perp)
    };
    let d_a_old = old_dist(ca, &signed_perp_a);
    let d_b_old = old_dist(cb, &signed_perp_b);
    for (label, d) in [("old dist a", d_a_old), ("old dist b", d_b_old)] {
        if d <= 0.0 {
            return Err(BlendError::InvalidParameter {
                parameter: "chamfer_face".to_string(),
                value: format!("{chamfer_face:?}"),
                reason: format!("{label} measured non-positive — unexpected strip layout"),
            });
        }
    }

    // Every boundary vertex maps across the cross-section: express
    // its offset from the restored edge in the (perp_a, perp_b)
    // frame, recover its old span fraction t (t = 0 at the A
    // contact, t = 1 at the B contact), and re-emit it on the new
    // segment. Contact corners land exactly at the new distances;
    // mid-bevel vertices interpolate.
    let gram = signed_perp_a.dot(&signed_perp_b);
    let det = 1.0 - gram * gram;
    if det.abs() < 1.0e-12 {
        return Err(BlendError::InvalidParameter {
            parameter: "chamfer_face".to_string(),
            value: format!("{chamfer_face:?}"),
            reason: "support perpendiculars are parallel — degenerate dihedral".to_string(),
        });
    }
    let mut out = brep.clone();
    let mut welded: Vec<VertexHandle> = Vec::with_capacity(loop_vertices.len());
    for &v in &loop_vertices {
        let p_old = brep.vertices[v].position;
        let along = (p_old - plan.restored_edge_start).dot(&edge_dir);
        let base = plan.restored_edge_start + edge_dir * along;
        let delta = p_old - base;
        let da = delta.dot(&signed_perp_a);
        let db = delta.dot(&signed_perp_b);
        let alpha = (da - gram * db) / det;
        let beta = (db - gram * da) / det;
        let t = (beta / d_b_old).clamp(0.0, 1.0);
        // Cross-section honesty: on a straight bevel every boundary
        // vertex satisfies alpha/d_a + beta/d_b = 1. A vertex off
        // that line means the strip is NOT the plane-plane bevel
        // this path models — decline rather than reshape it.
        let span_residual = (alpha / d_a_old + beta / d_b_old - 1.0).abs();
        if span_residual > 0.05 {
            return Err(BlendError::InvalidParameter {
                parameter: "chamfer_face".to_string(),
                value: format!("{chamfer_face:?}"),
                reason: format!(
                    "boundary vertex off the straight bevel cross-section                      (span residual {span_residual:.3}) — not a plane-plane chamfer strip"
                ),
            });
        }
        let p_new = base + signed_perp_a * ((1.0 - t) * dist_a) + signed_perp_b * (t * dist_b);
        out.vertices[v].position = p_new;
        welded.push(v);
    }
    let corners = [
        plan.restored_edge_start + signed_perp_a * dist_a,
        plan.restored_edge_start + signed_perp_b * dist_b,
        plan.restored_edge_end + signed_perp_a * dist_a,
        plan.restored_edge_end + signed_perp_b * dist_b,
    ];

    rebuild_curves_touching(&mut out, &welded);

    // The strip plane through the moved corners: u along the restored
    // edge, v across the bevel from the A-side contact toward the
    // B-side contact.
    let across =
        (corners[1] - corners[0])
            .normalized()
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "chamfer_face".to_string(),
                value: format!("{chamfer_face:?}"),
                reason: "degenerate bevel cross direction (coincident contact corners)".to_string(),
            })?;
    let normal =
        edge_dir
            .cross(&across)
            .normalized()
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "chamfer_face".to_string(),
                value: format!("{chamfer_face:?}"),
                reason: "degenerate bevel plane (edge direction parallel to cross direction)"
                    .to_string(),
            })?;
    let new_surface = SurfaceGeometry::Plane {
        origin: corners[0],
        normal,
        u_axis: edge_dir,
        v_axis: across,
    };
    out.faces[chamfer_face].surface = Some(SurfaceOnFace {
        geometry: new_surface,
        orientation: true,
    });
    crate::construct::strip_orientation::reconcile_strip_orientation_boundary(
        &mut out,
        chamfer_face,
    );
    crate::modify_inplace_chart::reauthor_strip_chart_pcurves(&mut out, chamfer_face)?;

    // Same acceptance the fillet path relies on: the update must
    // leave a structurally valid body — a failed rebuild here must
    // never serve (the caller's fallback is remove+recreate).
    out.validate().map_err(|e| BlendError::InvalidParameter {
        parameter: "chamfer_face".to_string(),
        value: format!("{chamfer_face:?}"),
        reason: format!("in-place chamfer update failed validation: {e}"),
    })?;
    let _ = tolerance;
    Ok(out)
}

/// The bevel's true supports are its two neighbours whose planes
/// meet the strip plane at a genuinely oblique dihedral (45° for an
/// equal chamfer; anywhere clearly away from 0° and 90° for a
/// two-distance one). End faces meet the strip at ~90° and coplanar
/// continuations at ~0°, so the oblique pair is unambiguous.
pub(crate) fn chamfer_supports_by_dihedral(
    brep: &BRep,
    chamfer_face: FaceHandle,
) -> Result<(FaceHandle, FaceHandle), BlendError> {
    let strip_plane = extract_plane(brep, chamfer_face)?;
    let mut oblique: Vec<FaceHandle> = Vec::new();
    for n in crate::identify::blend_neighbours(brep, chamfer_face) {
        let Ok(p) = extract_plane(brep, n) else {
            continue;
        };
        let cosine = p.normal.dot(&strip_plane.normal).abs();
        if (0.05..=0.95).contains(&cosine) {
            oblique.push(n);
        }
    }
    if oblique.len() != 2 {
        return Err(BlendError::InvalidParameter {
            parameter: "chamfer_face".to_string(),
            value: format!("{chamfer_face:?}"),
            reason: format!(
                "expected exactly 2 oblique support neighbours, found {} — not an open                  plane-plane chamfer strip",
                oblique.len()
            ),
        });
    }
    Ok((oblique[0], oblique[1]))
}

#[cfg(test)]
#[path = "modify_chamfer_inplace_tests.rs"]
mod tests;
