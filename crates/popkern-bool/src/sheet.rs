//! Sheet body / open shell detection (§2.3).
//!
//! A sheet body is a BRep whose boundary is not closed: at least one
//! edge is "free" — used by only one face. The boolean engine's
//! classification logic assumes a closed volume (point-set membership)
//! and produces undefined results on sheet bodies. This module
//! provides the detector that the public entry points use to reject
//! sheet inputs with a structured [`crate::BoolError::OpenShellNotSupported`]
//! error rather than silently corrupting the result.
//!
//! Future work will add boundary-relative classification semantics
//! and a sheet-aware operation matrix per §2.3 / §6.x. Until then,
//! this module's role is purely diagnostic.

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
            // Two distinct face uses → manifold edge, not free.
            continue;
        }

        // Same face on both sides. This is either a seam edge on a
        // periodic surface (closed) or a real free edge (open).
        // Inspect the face's surface to disambiguate.
        let Some(face) = brep.faces.get(face_a) else {
            return true;
        };
        let Some(surf_on_face) = face.surface.as_ref() else {
            // No surface attached: can't determine periodicity, so
            // treat as a free edge for safety.
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
            // Same face on both sides AND surface not periodic =
            // genuine free boundary.
            return true;
        }
        // Otherwise it is a seam edge — does not break closure.
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn unit_cube_is_closed() {
        let cube = BRep::unit_cube(tol()).unwrap();
        assert!(is_closed_solid(&cube));
        assert!(!is_sheet_body(&cube));
    }

    #[test]
    fn cylinder_is_closed() {
        let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
        assert!(is_closed_solid(&cyl));
        assert!(!is_sheet_body(&cyl));
    }

    #[test]
    fn empty_brep_is_not_a_sheet() {
        let empty = BRep::new(tol());
        assert!(!is_sheet_body(&empty));
        assert!(is_closed_solid(&empty));
    }

    #[test]
    fn dangling_edge_brep_is_a_sheet() {
        // mvfs + mev creates a single face holding a dangling edge.
        // The edge's two half-edges both belong to the same face;
        // the face has no surface attached, so the periodic-skip
        // branch does not fire and the detector reports a sheet.
        use popkern_core::geometry::Point3;
        let mut brep = BRep::new(tol());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        assert!(
            is_sheet_body(&brep),
            "dangling-edge BRep should be a sheet body"
        );
    }

    #[test]
    fn boolean_on_sheet_body_returns_open_shell_error() {
        use popkern_core::geometry::Point3;
        let mut sheet = BRep::new(tol());
        let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = sheet
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        let cube = BRep::unit_cube(tol()).unwrap();
        let result = crate::union(sheet, cube, tol(), None);
        assert!(
            matches!(
                result,
                Err(crate::BoolError::OpenShellNotSupported { a: true, b: false })
            ),
            "expected OpenShellNotSupported, got {result:?}"
        );
    }
}
