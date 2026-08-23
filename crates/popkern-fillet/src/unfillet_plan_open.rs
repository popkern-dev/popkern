//! Open-blend strip plan derivation (split from `unfillet.rs` for
//! the 500-line structural cap): [`compute_unfillet_plan`] and the
//! kind-parameterised [`compute_blend_strip_plan`] family shared by
//! the fillet and chamfer in-place modification paths (PK-280).

use popkern_core::topology::{BRep, FaceHandle};

use super::UnfilletPlan;
use crate::error::BlendError;
use crate::identify::{BlendKind, identify_blend};
use crate::unfillet_helpers::{
    collect_end_face_planes, extract_plane, find_torus_support_faces, intersect_line_plane,
    intersect_planes, is_closed_fillet_annulus, pick_extreme_endpoints, walk_blend_boundary,
};

/// Compute an unfillet plan without modifying the input BRep.
///
/// Pure function. Open plane-plane blends plan a restored
/// [`RestoredEdgeGeometry::Line`]; closed plane-cylinder rims
/// (annulus or seamed band) plan the restored intersection
/// [`RestoredEdgeGeometry::Circle`]. Open blends on curved supports
/// decline with the configuration named — the operation record
/// (`unfillet_via_record`) is the route for those.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] when `blend_face` is not
///   recognised as a blend, when its kind is not
///   [`BlendKind::Fillet`], or when the support configuration has
///   no implemented plan class.
pub fn compute_unfillet_plan(
    brep: &BRep,
    blend_face: FaceHandle,
) -> Result<UnfilletPlan, BlendError> {
    compute_blend_strip_plan(brep, blend_face, &[BlendKind::Fillet])
}

/// Kind-parameterised strip plan (PK-280): everything downstream of
/// the kind gate — support-plane intersection, end-plane endpoint
/// derivation, boundary walk — is identical for fillets and chamfers,
/// so the chamfer in-place path shares this plan instead of
/// duplicating it. `unfillet_face` keeps its fillet-only contract via
/// [`compute_unfillet_plan`].
pub(crate) fn compute_blend_strip_plan(
    brep: &BRep,
    blend_face: FaceHandle,
    accepted_kinds: &[BlendKind],
) -> Result<UnfilletPlan, BlendError> {
    compute_blend_strip_plan_with_supports(brep, blend_face, accepted_kinds, None)
}

/// [`compute_blend_strip_plan`] with an optional caller-supplied
/// support pair. Curvature-analysis identification cannot reliably
/// pick a CHAMFER's supports (every plane flanked by non-tangent
/// planes pattern-matches; the measured cube bevel got one support
/// and one END face, leaving a single restored-edge endpoint) — the
/// chamfer in-place path derives supports geometrically (the
/// neighbours meeting the bevel at 45°) and passes them here.
pub(crate) fn compute_blend_strip_plan_with_supports(
    brep: &BRep,
    blend_face: FaceHandle,
    accepted_kinds: &[BlendKind],
    support_override: Option<(FaceHandle, FaceHandle)>,
) -> Result<UnfilletPlan, BlendError> {
    // Closed rims first (same dispatch `unfillet_face` uses): their
    // plane+cylinder supports would misread as open-blend declines.
    if is_closed_fillet_annulus(brep, blend_face)
        || crate::unfillet_closed::seamed_band_bridge(brep, blend_face).is_some()
    {
        return crate::unfillet_plan_closed::compute_closed_rim_plan(brep, blend_face);
    }
    let info = identify_blend(brep, blend_face).ok_or_else(|| BlendError::InvalidParameter {
        parameter: "blend_face".to_string(),
        value: format!("{blend_face:?}"),
        reason: "face is not recognised as a blend".to_string(),
    })?;

    if !accepted_kinds.contains(&info.kind) {
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face.kind".to_string(),
            value: format!("{:?}", info.kind),
            reason: format!("this operation supports {accepted_kinds:?} blend faces only"),
        });
    }

    // For a Torus- or Cylinder-surface fillet the canonical support
    // detection is "plane neighbours whose normals are perpendicular
    // to the fillet axis". This is more reliable than reusing the
    // tangency-sample-based detection from identify_blend, which can
    // miss G1 boundaries when the contact-curve midpoint doesn't
    // land precisely on the tangent line. Falls back to
    // identify_blend's pair if the geometric search fails.
    let (support_a, support_b) = support_override
        .or_else(|| find_torus_support_faces(brep, blend_face))
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
        // PK-143 explicit fallback boundary: a strip end that meets
        // a corner CAP has no planar end face to intersect — the
        // corner's radii are coupled and the group must re-dispatch
        // together. Name the record path so callers route there.
        // A failed endpoint derivation means the recordless restore
        // cannot serve this strip's class (corner caps, setback
        // sectors, constructively-reorganized neighborhoods) — the
        // operation record is the route, so name it unconditionally.
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: format!(
                "found {} valid restored-edge endpoints; need exactly 2 — the recordless \
                 restore cannot derive this class; unfillet/modify via the operation record \
                 (unfillet_via_record / modify_fillet_radius_via_record; coupled corner \
                 groups share one record)",
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
        restored_edge: crate::unfillet_plan_closed::RestoredEdgeGeometry::Line {
            start: p_start,
            end: p_end,
            direction: line_direction,
        },
        long_he_to_a,
        long_he_to_b,
        edges_to_remove,
    })
}
