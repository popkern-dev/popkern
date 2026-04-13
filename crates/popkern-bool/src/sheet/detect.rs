//! Sheet-body detection and geometric extraction helpers.
//!
//! [`is_sheet_body`] is the entry point used by the boolean engine's
//! pre-check to reject open shells upfront. [`is_planar_sheet`] is the
//! fast-path detector used by the multi-piece dispatcher to identify
//! planar sheets that can be split via half-space boolean operations.

use popkern_core::topology::BRep;

/// Returns `true` if `brep` is a sheet body (open shell).
///
/// Detection: an edge whose two half-edges share the same face is
/// "free" — but only if the underlying surface is NOT periodic in
/// the direction of that edge. Periodic surfaces (cylinder seam,
/// sphere meridian, torus seams) legitimately have face_a == face_b
/// at their seam edges, and the body is still closed.
///
/// The check therefore: walks every edge, identifies free edges by
/// the same-face criterion, and skips those whose hosting surface
/// reports periodicity. Anything left over is a real open boundary.
///
/// Empty BReps (no edges, no faces) are NOT considered sheet bodies
/// — they are simply degenerate, and the existing pre-check rejects
/// them via the "no faces" finding.
pub fn is_sheet_body(brep: &BRep) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;

    if brep.edges.is_empty() {
        return false;
    }

    for (_, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            return true;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            return true;
        };
        let face_a = brep.loops.get(he.loop_).map(|l| l.face);
        let face_b = brep.loops.get(twin.loop_).map(|l| l.face);
        let (face_a, face_b) = match (face_a, face_b) {
            (Some(fa), Some(fb)) => (fa, fb),
            _ => return true,
        };

        if face_a != face_b {
            continue;
        }

        let Some(face) = brep.faces.get(face_a) else {
            return true;
        };
        let Some(surf_on_face) = face.surface.as_ref() else {
            return true;
        };
        let is_periodic = matches!(
            surf_on_face.geometry,
            SurfaceGeometry::Cylinder { .. }
                | SurfaceGeometry::Cone { .. }
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
        );
        if !is_periodic {
            return true;
        }
    }
    false
}

/// Returns `true` if `brep` is a closed solid (no free edges).
///
/// Inverse of [`is_sheet_body`] for non-empty BReps. An empty BRep
/// (no edges) is reported as closed because there is nothing to
/// breach the closure condition.
pub fn is_closed_solid(brep: &BRep) -> bool {
    !is_sheet_body(brep)
}

/// Detect whether `sheet` is "effectively planar" and extract its plane.
///
/// Returns `Some((origin, outward_normal))` when every face of the
/// sheet lies on the same supporting plane within `tolerance.linear()`
/// — either because the sheet has a single planar face, or because it
/// has multiple planar faces that all share the same normal and offset
/// along that normal. The returned normal is the face's outward-facing
/// normal (orientation flag applied). Returns `None` for sheets with
/// non-planar surfaces, mixed normals, missing surface geometry, or
/// zero faces.
///
/// Used by [`crate::sheet::boolean_with_sheets`] to dispatch the
/// planar fast path for `SolidMinusSheet` — the resulting plane is
/// passed to [`crate::sheet::split_solid_by_plane`] to produce
/// multiple closed sub-solids via half-space split.
pub fn is_planar_sheet(
    sheet: &BRep,
) -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    if sheet.faces.is_empty() {
        return None;
    }

    // Tolerance floor for "same plane" checks. We compare the signed
    // offset of each face's own plane origin along the reference
    // normal — if every face agrees within ~1e-6, treat them as a
    // single plane.
    let offset_tol = sheet.tolerance.linear();
    // Cosine threshold for "same normal" (normals within ~1e-4 rad).
    let angle_tol = 1.0e-4;

    let mut ref_origin: Option<popkern_core::geometry::Point3> = None;
    let mut ref_normal: Option<popkern_core::geometry::Vec3> = None;

    for face in sheet.faces.values() {
        let surf_on_face = face.surface.as_ref()?;
        let SurfaceGeometry::Plane { origin, normal, .. } = &surf_on_face.geometry else {
            return None;
        };
        let face_normal = if surf_on_face.orientation {
            *normal
        } else {
            *normal * -1.0
        };

        match (ref_origin, ref_normal) {
            (None, None) => {
                ref_origin = Some(*origin);
                ref_normal = Some(face_normal);
            }
            (Some(_), Some(rn)) => {
                // Compare normals (absolute value — opposite-facing
                // planar faces are still on the same plane; we keep
                // the first-seen direction as the reference normal).
                let cos = rn.dot(&face_normal);
                if (cos.abs() - 1.0).abs() > angle_tol {
                    return None;
                }
                // Project this face's origin onto the reference
                // plane's signed-distance axis. If the offset is
                // larger than `offset_tol`, the face sits on a
                // different plane and we bail out.
                let delta = *origin - ref_origin.unwrap();
                let signed = delta.dot(&rn);
                if signed.abs() > offset_tol {
                    return None;
                }
            }
            _ => unreachable!(),
        }
    }

    Some((ref_origin?, ref_normal?))
}
