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

/// Detect whether `sheet` is a single planar face and extract its plane.
///
/// Returns `Some((origin, outward_normal))` if `sheet` has exactly one
/// face whose surface is a [`popkern_core::geometry::surface::SurfaceGeometry::Plane`],
/// with the outward normal already adjusted by the face's `orientation`
/// flag. Returns `None` for multi-face sheets, non-planar surfaces, or
/// missing surface geometry.
///
/// Used by [`crate::sheet::boolean_with_sheets`] to dispatch the planar
/// fast path for `SolidMinusSheet` to a half-space split that produces
/// multiple closed sub-solids.
pub fn is_planar_sheet(
    sheet: &BRep,
) -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    if sheet.faces.len() != 1 {
        return None;
    }
    let face = sheet.faces.values().next()?;
    let surf_on_face = face.surface.as_ref()?;
    let SurfaceGeometry::Plane { origin, normal, .. } = &surf_on_face.geometry else {
        return None;
    };
    let n = if surf_on_face.orientation {
        *normal
    } else {
        *normal * -1.0
    };
    Some((*origin, n))
}
