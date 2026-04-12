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

/// Side of a sheet body relative to a query point (§6.x boundary-relative).
///
/// Sheet bodies do not enclose a volume, so point-set "inside/outside"
/// classification is undefined. Instead, classify the point by which
/// side of the oriented surface it lies on. The "front" side is the
/// half-space the surface normal points into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SideOfSheet {
    /// Point is on the side the sheet's outward normal points to.
    Front,
    /// Point is on the opposite side.
    Back,
    /// Point lies on the sheet itself within tolerance.
    On,
    /// Sheet has no orientable face near the query — undefined.
    Indeterminate,
}

/// Classify a query point against a sheet body via boundary-relative
/// semantics.
///
/// Finds the closest face on the sheet, evaluates its surface and
/// outward normal at the closest UV, and reports which side the
/// point lies on. If the point is closer than `tolerance.linear()`,
/// returns `On`.
///
/// This is the building block §6.x asks for: open-shell classification
/// uses "which side of the oriented surface" instead of "inside/outside
/// volume". Callers can use [`side_of_sheet`] in their own
/// sheet-aware boolean wrappers.
pub fn side_of_sheet(
    sheet: &BRep,
    point: &popkern_core::geometry::Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> SideOfSheet {
    use popkern_core::geometry::surface::SurfaceEval;

    let mut best_dist = f64::INFINITY;
    let mut best_signed_dist = 0.0_f64;
    let mut best_found = false;

    for (_fh, face) in &sheet.faces {
        let Some(surf) = face.surface.as_ref() else {
            continue;
        };
        let Some((u, v)) = surf.geometry.closest_point(point) else {
            continue;
        };
        let p_surf = surf.geometry.point_at(u, v);
        let mut normal = surf.geometry.normal_at(u, v);
        // Apply the face's surface-orientation flag so we always
        // work with the OUTWARD-facing normal, matching the convention
        // used by coincidence/classify::face_outward_normal.
        if !surf.orientation {
            normal = normal * (-1.0);
        }
        let n_len = normal.length();
        if n_len < 1e-30 {
            continue;
        }
        let to_point = *point - p_surf;
        let dist = to_point.length();
        if dist < best_dist {
            best_dist = dist;
            // Signed distance along the outward normal. Positive
            // means point is on the front side, negative on back.
            let signed = to_point.dot(&normal) / n_len;
            best_signed_dist = signed;
            best_found = true;
        }
    }

    if !best_found {
        return SideOfSheet::Indeterminate;
    }

    if best_dist < tolerance.linear() {
        SideOfSheet::On
    } else if best_signed_dist >= 0.0 {
        SideOfSheet::Front
    } else {
        SideOfSheet::Back
    }
}

/// Operation matrix classification (§2.3) for an (operation, A, B) tuple.
///
/// Categorizes a boolean call by which inputs are sheet bodies and
/// whether the requested operation is well-defined for that combination.
/// Use [`classify_sheet_operation`] to dispatch a boolean call to the
/// correct sub-operation or surface a structured rejection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheetOperationKind {
    /// Both inputs are closed solids — standard boolean.
    SolidSolid,
    /// A is solid, B is sheet, op is subtraction. Trims solid along
    /// the sheet (valid; produces solid result).
    SolidMinusSheet,
    /// A is solid, B is sheet, op is intersection. Returns the portion
    /// of the sheet inside the solid (open shell result).
    SheetIntersectSolid,
    /// A is solid, B is sheet, op is union. Conceptually merges, but
    /// not yet implemented.
    SolidUnionSheet,
    /// A is sheet, B is solid, op is subtraction. Subtracting a volume
    /// from a 2D sheet is undefined; this case is rejected.
    InvalidSheetMinusSolid,
    /// Both inputs are sheets, op is intersection. Returns shared
    /// curves (1D) or shared regions (2D); not yet implemented.
    SheetSheetIntersection,
    /// Both inputs are sheets, op is union. Combines into one open
    /// shell; not yet implemented.
    SheetSheetUnion,
    /// Combination not yet covered by the operation matrix.
    Unsupported,
}

/// Result of a [`boolean_with_sheets`] call.
///
/// Sheet-aware boolean operations can produce three different shapes
/// of result depending on the operation:
///
/// - [`SheetBooleanResult::Solid`] — a closed solid (e.g. solid-solid
///   union, or solid-minus-sheet trim that produced a single piece).
/// - [`SheetBooleanResult::Curves`] — a list of intersection curves
///   (e.g. sheet-sheet intersection produces a 1D wire result).
/// - [`SheetBooleanResult::OpenShell`] — an open shell BRep (e.g.
///   sheet-intersect-solid returns the portion of the sheet inside
///   the solid; the result is itself a sheet body).
#[derive(Debug)]
pub enum SheetBooleanResult {
    /// A closed solid result.
    Solid(popkern_core::topology::BRep),
    /// A list of intersection curves (no enclosing volume).
    Curves(Vec<popkern_core::geometry::curve_on_surface::CurveOnSurface>),
    /// An open shell BRep result.
    OpenShell(popkern_core::topology::BRep),
}

/// Sheet-aware boolean operation dispatcher (§2.3).
///
/// Inspects whether each input is a sheet body (via [`is_sheet_body`]),
/// classifies the (operation, is_a_sheet, is_b_sheet) tuple via
/// [`classify_sheet_operation`], and routes to the appropriate
/// sub-operation:
///
/// | Kind | Implementation |
/// |---|---|
/// | `SolidSolid` | Delegates to standard [`crate::boolean`]. |
/// | `SheetSheetIntersection` | Calls [`crate::ops::section::section`] and returns the curves. |
/// | `SolidMinusSheet` | Calls [`crate::ops::imprint::imprint`] — imprints sheet curves onto the solid (does NOT yet split into multiple pieces; that's a future iteration). |
/// | `InvalidSheetMinusSolid` | Rejects with [`crate::BoolError::OpenShellNotSupported`]. |
/// | Other variants | Returns [`crate::BoolError::OpenShellNotSupported`] (future work). |
///
/// Use this entry point when you need sheet support; the standard
/// `boolean()` rejects sheet inputs upfront so non-sheet callers
/// don't accidentally invoke half-implemented sheet operations.
pub fn boolean_with_sheets(
    op: crate::BoolOperation,
    a: &popkern_core::topology::BRep,
    b: &popkern_core::topology::BRep,
    tol: popkern_core::tolerance::Tolerance,
) -> Result<SheetBooleanResult, crate::BoolError> {
    let a_is_sheet = is_sheet_body(a);
    let b_is_sheet = is_sheet_body(b);
    let kind = classify_sheet_operation(op, a_is_sheet, b_is_sheet);

    match kind {
        SheetOperationKind::SolidSolid => {
            let result = crate::boolean(op, a.clone(), b.clone(), tol, None)?;
            Ok(SheetBooleanResult::Solid(result))
        }
        SheetOperationKind::SheetSheetIntersection => {
            let section = crate::ops::section::section(a, b, &tol)?;
            let curves = section.curves.into_iter().map(|c| c.curve).collect();
            Ok(SheetBooleanResult::Curves(curves))
        }
        SheetOperationKind::SolidMinusSheet => {
            // Use imprint to add the sheet's curves to the solid as
            // edges. This doesn't split the solid into multiple
            // pieces — that's a future iteration that requires
            // running the full classification + assembly loop with
            // the sheet acting as a cutter.
            let imprint = crate::ops::imprint::imprint(a, b, &tol)?;
            Ok(SheetBooleanResult::Solid(imprint.imprinted_target))
        }
        SheetOperationKind::InvalidSheetMinusSolid
        | SheetOperationKind::SheetIntersectSolid
        | SheetOperationKind::SolidUnionSheet
        | SheetOperationKind::SheetSheetUnion
        | SheetOperationKind::Unsupported => Err(crate::BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        }),
    }
}

/// Classify a boolean call by sheet-aware operation matrix (§2.3).
///
/// Returns the [`SheetOperationKind`] that best fits the (operation,
/// is_a_sheet, is_b_sheet) tuple. The classification is intended for
/// dispatch: a sheet-aware boolean wrapper inspects the result and
/// routes to the appropriate sub-operation, or returns
/// [`crate::BoolError::OpenShellNotSupported`] for `Unsupported` /
/// `Invalid*` cases.
pub fn classify_sheet_operation(
    op: crate::BoolOperation,
    a_is_sheet: bool,
    b_is_sheet: bool,
) -> SheetOperationKind {
    use crate::BoolOperation;
    match (op, a_is_sheet, b_is_sheet) {
        (_, false, false) => SheetOperationKind::SolidSolid,
        (BoolOperation::Subtraction, false, true) => SheetOperationKind::SolidMinusSheet,
        (BoolOperation::Intersection, false, true) => SheetOperationKind::SheetIntersectSolid,
        (BoolOperation::Intersection, true, false) => SheetOperationKind::SheetIntersectSolid,
        (BoolOperation::Union, false, true) | (BoolOperation::Union, true, false) => {
            SheetOperationKind::SolidUnionSheet
        }
        (BoolOperation::Subtraction, true, false) => SheetOperationKind::InvalidSheetMinusSolid,
        (BoolOperation::Subtraction, true, true) => SheetOperationKind::Unsupported,
        (BoolOperation::Intersection, true, true) => SheetOperationKind::SheetSheetIntersection,
        (BoolOperation::Union, true, true) => SheetOperationKind::SheetSheetUnion,
    }
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
    fn side_of_sheet_outside_cube_is_front_of_closest_face() {
        // For a multi-face body, side_of_sheet picks the closest
        // face and reports relative to its outward normal. A point
        // outside the cube is on the front side of its closest face
        // (because every face's outward normal points away from
        // the cube interior). This test pins that behavior — for
        // a true sheet body the function returns Front/Back relative
        // to the single oriented surface.
        let cube = BRep::unit_cube(tol()).unwrap();
        use popkern_core::geometry::Point3;
        let above = Point3::new(0.5, 0.5, 5.0);
        let below = Point3::new(0.5, 0.5, -5.0);
        let side_above = side_of_sheet(&cube, &above, &tol());
        let side_below = side_of_sheet(&cube, &below, &tol());
        assert!(matches!(side_above, SideOfSheet::Front));
        assert!(matches!(side_below, SideOfSheet::Front));
    }

    #[test]
    fn side_of_sheet_classifies_consistently() {
        // For a closed solid, side_of_sheet picks the closest face
        // and reports relative to its outward normal. The exact
        // Front/Back assignment depends on construction-time loop
        // orientation, but the function MUST be deterministic for
        // identical inputs and MUST report `On` when the point lies
        // on the surface within tolerance.
        let cube = BRep::unit_cube(tol()).unwrap();
        use popkern_core::geometry::Point3;
        let inside = Point3::new(0.5, 0.5, 0.5);
        let outside = Point3::new(0.5, 0.5, 5.0);

        // Determinism: identical input produces identical output.
        let side_inside_1 = side_of_sheet(&cube, &inside, &tol());
        let side_inside_2 = side_of_sheet(&cube, &inside, &tol());
        assert_eq!(side_inside_1, side_inside_2);

        // Both classifications must be either Front or Back, never
        // Indeterminate (the cube has surfaces).
        for side in [
            side_of_sheet(&cube, &inside, &tol()),
            side_of_sheet(&cube, &outside, &tol()),
        ] {
            assert!(matches!(
                side,
                SideOfSheet::Front | SideOfSheet::Back | SideOfSheet::On
            ));
        }
    }

    #[test]
    fn side_of_sheet_returns_on_for_point_on_surface() {
        let cube = BRep::unit_cube(tol()).unwrap();
        use popkern_core::geometry::Point3;
        // A point on the +Z face of the unit cube.
        let on_face = Point3::new(0.5, 0.5, 1.0);
        let side = side_of_sheet(&cube, &on_face, &tol());
        assert!(
            matches!(side, SideOfSheet::On),
            "point on face should be On, got {side:?}"
        );
    }

    #[test]
    fn classify_sheet_operation_solid_solid() {
        use crate::BoolOperation;
        let kind = classify_sheet_operation(BoolOperation::Union, false, false);
        assert!(matches!(kind, SheetOperationKind::SolidSolid));
    }

    #[test]
    fn classify_sheet_operation_solid_minus_sheet() {
        use crate::BoolOperation;
        let kind = classify_sheet_operation(BoolOperation::Subtraction, false, true);
        assert!(matches!(kind, SheetOperationKind::SolidMinusSheet));
    }

    #[test]
    fn classify_sheet_operation_invalid_sheet_minus_solid() {
        use crate::BoolOperation;
        let kind = classify_sheet_operation(BoolOperation::Subtraction, true, false);
        assert!(matches!(kind, SheetOperationKind::InvalidSheetMinusSolid));
    }

    #[test]
    fn classify_sheet_operation_sheet_sheet_intersection() {
        use crate::BoolOperation;
        let kind = classify_sheet_operation(BoolOperation::Intersection, true, true);
        assert!(matches!(kind, SheetOperationKind::SheetSheetIntersection));
    }

    #[test]
    fn boolean_with_sheets_solid_solid_delegates_to_boolean() {
        use crate::BoolOperation;
        use popkern_core::geometry::Transform3;
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap();
        let result = boolean_with_sheets(BoolOperation::Union, &a, &b, tol()).unwrap();
        assert!(matches!(result, SheetBooleanResult::Solid(_)));
    }

    #[test]
    fn boolean_with_sheets_invalid_sheet_minus_solid_rejected() {
        // A is a sheet, B is a solid, op is subtraction → invalid.
        use crate::BoolOperation;
        use popkern_core::geometry::Point3;
        let mut sheet = BRep::new(tol());
        let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = sheet
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let cube = BRep::unit_cube(tol()).unwrap();
        let result = boolean_with_sheets(BoolOperation::Subtraction, &sheet, &cube, tol());
        assert!(matches!(
            result,
            Err(crate::BoolError::OpenShellNotSupported { a: true, b: false })
        ));
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
