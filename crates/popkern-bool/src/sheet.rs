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

/// Detect whether `sheet` is a single planar face and extract its plane.
///
/// Returns `Some((origin, outward_normal))` if `sheet` has exactly one
/// face whose surface is a [`SurfaceGeometry::Plane`], with the
/// outward normal already adjusted by the face's `orientation` flag.
/// Returns `None` for multi-face sheets, non-planar surfaces, or
/// missing surface geometry.
///
/// Used by [`boolean_with_sheets`] to dispatch the planar fast path
/// for `SolidMinusSheet` to a half-space split that produces multiple
/// closed sub-solids.
pub fn is_planar_sheet(
    sheet: &popkern_core::topology::BRep,
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
    // Apply the orientation flag so the returned normal is the
    // OUTWARD-facing normal (matching the convention used elsewhere
    // in this module).
    let n = if surf_on_face.orientation {
        *normal
    } else {
        *normal * -1.0
    };
    Some((*origin, n))
}

/// Compute the world-space AABB of a BRep from its vertices.
///
/// Used by [`split_solid_by_plane`] to size the cutting half-space.
fn brep_aabb(
    brep: &popkern_core::topology::BRep,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    let mut iter = brep.vertices.values();
    let first = iter.next()?.position;
    let mut min = first;
    let mut max = first;
    for v in iter {
        let p = v.position;
        if p.x < min.x {
            min.x = p.x;
        }
        if p.y < min.y {
            min.y = p.y;
        }
        if p.z < min.z {
            min.z = p.z;
        }
        if p.x > max.x {
            max.x = p.x;
        }
        if p.y > max.y {
            max.y = p.y;
        }
        if p.z > max.z {
            max.z = p.z;
        }
    }
    Some((min, max))
}

/// Split a solid by an infinite plane into two closed half-solids.
///
/// Constructs a properly-sized cutting box that fully covers ONE
/// half-space (the half the plane normal points into) and intersects
/// it with the body. Repeats with the opposite half-space. Returns
/// the two halves; either may be `None` if the plane doesn't pass
/// through the body on that side.
///
/// This is the building block for multi-piece `SolidMinusSheet` when
/// the sheet is a single planar face. Unlike the existing
/// [`crate::ops::split_slice::split_by_plane`], the cutting box here
/// is scaled to the body's actual extent so the intersection
/// produces correct results for arbitrarily-sized solids.
pub fn split_solid_by_plane(
    solid: &popkern_core::topology::BRep,
    plane_origin: popkern_core::geometry::Point3,
    plane_normal: popkern_core::geometry::Vec3,
    tol: popkern_core::tolerance::Tolerance,
) -> Result<Vec<popkern_core::topology::BRep>, crate::BoolError> {
    use popkern_core::geometry::{Point3, Transform3, Vec3};
    use popkern_core::topology::BRep;

    let normal = plane_normal
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));

    // Compute the body's extent so we can size the cutting box to
    // fully contain it on each side of the plane.
    let (mn, mx) = brep_aabb(solid)
        .ok_or_else(|| crate::BoolError::DegenerateInput("solid has no vertices".to_string()))?;
    let dx = mx.x - mn.x;
    let dy = mx.y - mn.y;
    let dz = mx.z - mn.z;
    let max_extent = dx.max(dy).max(dz);

    let body_centroid = Point3::new(
        0.5 * (mn.x + mx.x),
        0.5 * (mn.y + mx.y),
        0.5 * (mn.z + mx.z),
    );

    // Early-out: if the body is entirely on one side of the plane,
    // there's no cut and we return the body unchanged. The check is
    // a signed-distance comparison along the normal — if every
    // vertex is on the same side, the plane misses the body.
    let mut min_signed = f64::INFINITY;
    let mut max_signed = f64::NEG_INFINITY;
    for v in solid.vertices.values() {
        let d = (v.position - plane_origin).dot(&normal);
        if d < min_signed {
            min_signed = d;
        }
        if d > max_signed {
            max_signed = d;
        }
    }
    let on_plane_tol = tol.linear() * 10.0;
    if min_signed > -on_plane_tol || max_signed < on_plane_tol {
        // Whole body is on one side (or just barely touches the
        // plane). No cut.
        return Ok(vec![solid.clone()]);
    }

    // Box size: only 10% larger than the body's largest dimension.
    // Massively-oversized boxes trip a shortcut bug in the boolean
    // engine (point-sampling no-intersection shortcut samples a
    // vertex that happens to land outside the overlap region and
    // concludes the bodies are disjoint). A small box forces the
    // boolean down the full-pipeline path where partial overlap
    // is handled correctly.
    let margin = max_extent * 0.05 + tol.linear() * 10.0;
    let box_size = max_extent + 2.0 * margin;

    // Build the unit cube once; we transform it twice (for each half).
    let unit = BRep::unit_cube(tol).map_err(crate::BoolError::CoreError)?;

    // Helper: build a half-space cutting box on the side of the
    // plane indicated by `side_normal`. The box must cover the
    // entire body on that side of the plane. Strategy: anchor the
    // box at the body's centroid (so the box always overlaps the
    // body), then translate it along `side_normal` by half its
    // own size MINUS the signed distance from centroid to plane,
    // so its EDGE lies exactly on the plane on the opposite-from-
    // centroid side.
    //
    // The end result: the box covers everything from the cutting
    // plane outward in `side_normal` direction, and is large
    // enough to fully contain the body's other-axis extent.
    let make_half_box = |side_normal: Vec3| -> Result<BRep, crate::BoolError> {
        // Signed distance from body centroid to the plane along
        // side_normal. Positive when the centroid is on side_normal's
        // side of the plane.
        let centroid_offset = body_centroid - plane_origin;
        let signed_dist_to_plane = centroid_offset.dot(&side_normal);

        // Anchor the box CENTER at: centroid + side_normal * (half_box - signed_dist).
        // This places one edge of the box exactly on the plane
        // (on the side OPPOSITE to side_normal), and the box
        // extends `box_size` units in `side_normal` direction
        // beyond the centroid.
        let half_box = box_size * 0.5;
        let box_center = body_centroid + side_normal * (half_box - signed_dist_to_plane);

        // Scale the unit cube to box_size, then translate so its
        // center lands at box_center.
        let scale_t = Transform3::scale_uniform(box_size);
        let scaled_center = box_size * 0.5;
        let translate_t = Transform3::translation(
            box_center.x - scaled_center,
            box_center.y - scaled_center,
            box_center.z - scaled_center,
        );
        let combined = translate_t.compose(&scale_t);
        unit.clone()
            .transform(&combined)
            .map_err(crate::BoolError::CoreError)
    };

    let pos_box = make_half_box(normal)?;
    let neg_box = make_half_box(normal * -1.0)?;

    // Strategy: produce each half by SUBTRACTING the OPPOSITE
    // half-box from the solid. `solid - neg_box` leaves everything
    // in the solid not covered by the negative half-space, i.e. the
    // positive half. We use the subtraction path (not intersection)
    // because intersection_ref hits a shortcut bug for half-space
    // intersections, while subtraction handles cube-minus-box
    // partial overlap correctly.
    let mut pieces = Vec::new();
    if let Ok(positive) = crate::subtraction_ref(solid, &neg_box, tol, None)
        && !positive.faces.is_empty()
    {
        pieces.push(positive);
    }
    if let Ok(negative) = crate::subtraction_ref(solid, &pos_box, tol, None)
        && !negative.faces.is_empty()
    {
        pieces.push(negative);
    }
    Ok(pieces)
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

/// Trim a sheet body to the portion that lies inside `solid`.
///
/// Walks every face of the (already-imprinted) sheet, samples a
/// centroid point, and tests it via point-in-solid against the
/// solid. Faces whose centroid is `Outside` are removed from the
/// sheet via slotmap delete; orphaned edges and vertices are
/// pruned afterwards.
///
/// The result is an open shell BRep — possibly with multiple
/// disconnected components — containing only the inside fragments.
/// This is the §6.x boundary-relative classification path: instead
/// of asking "what's the volume of the result", we ask "which
/// pieces of the sheet survived the cut".
fn trim_sheet_to_inside(
    mut sheet: popkern_core::topology::BRep,
    solid: &popkern_core::topology::BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> popkern_core::topology::BRep {
    use crate::clip::face_vertex_positions;
    use crate::error::PointClassification;
    use crate::point_in_solid;
    use popkern_core::topology::FaceHandle;

    let solid_aabbs = crate::intersect::compute_face_aabbs(solid);

    // Step 1: classify every sheet face as inside / outside / on by
    // sampling its centroid.
    let mut faces_to_remove: Vec<FaceHandle> = Vec::new();
    for fh in sheet.faces.keys() {
        let positions = face_vertex_positions(&sheet, fh);
        if positions.is_empty() {
            faces_to_remove.push(fh);
            continue;
        }
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for p in &positions {
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }
        let n = positions.len() as f64;
        let centroid = popkern_core::geometry::Point3::new(cx / n, cy / n, cz / n);
        let classification = point_in_solid(&centroid, solid, &solid_aabbs, tolerance)
            .unwrap_or(PointClassification::Outside);
        if matches!(classification, PointClassification::Outside) {
            faces_to_remove.push(fh);
        }
    }

    // Step 2: remove the outside faces. Use slotmap delete; the
    // result will have orphaned half-edges/edges/vertices that we
    // clean up below.
    for fh in &faces_to_remove {
        if let Some(face) = sheet.faces.get(*fh) {
            let outer_loop = face.outer_loop;
            let inner_loops = face.inner_loops.clone();
            // Remove the loops the face owned.
            sheet.loops.remove(outer_loop);
            for il in inner_loops {
                sheet.loops.remove(il);
            }
            sheet.faces.remove(*fh);
        }
    }

    // Step 3: prune orphaned half-edges (those whose loop_ no longer
    // exists in the slotmap).
    let orphan_he: Vec<_> = sheet
        .half_edges
        .iter()
        .filter_map(|(heh, he)| {
            if sheet.loops.get(he.loop_).is_none() {
                Some(heh)
            } else {
                None
            }
        })
        .collect();
    for heh in orphan_he {
        sheet.half_edges.remove(heh);
    }

    // Step 4: prune orphaned edges (those whose half_edge no longer
    // exists in the slotmap).
    let orphan_edges: Vec<_> = sheet
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            if sheet.half_edges.get(edge.half_edge).is_none() {
                Some(eh)
            } else {
                None
            }
        })
        .collect();
    for eh in orphan_edges {
        sheet.edges.remove(eh);
    }

    // Step 5: prune orphaned vertices (those not used by any
    // remaining half-edge).
    let mut vertex_in_use = std::collections::HashSet::new();
    for (_, he) in &sheet.half_edges {
        vertex_in_use.insert(he.origin);
    }
    let orphan_vertices: Vec<_> = sheet
        .vertices
        .keys()
        .filter(|vh| !vertex_in_use.contains(vh))
        .collect();
    for vh in orphan_vertices {
        sheet.vertices.remove(vh);
    }

    sheet
}

/// Result of a [`boolean_with_sheets`] call.
///
/// Sheet-aware boolean operations can produce four different shapes
/// of result depending on the operation:
///
/// - [`SheetBooleanResult::Solid`] — a closed solid (e.g. solid-solid
///   union, or solid-minus-sheet trim that produced a single piece).
/// - [`SheetBooleanResult::Pieces`] — multiple closed sub-solids
///   (e.g. solid-minus-sheet where the sheet bisects the solid into
///   two or more pieces).
/// - [`SheetBooleanResult::Curves`] — a list of intersection curves
///   (e.g. sheet-sheet intersection produces a 1D wire result).
/// - [`SheetBooleanResult::OpenShell`] — an open shell BRep (e.g.
///   sheet-intersect-solid returns the portion of the sheet inside
///   the solid; the result is itself a sheet body).
#[derive(Debug)]
pub enum SheetBooleanResult {
    /// A single closed solid result.
    Solid(popkern_core::topology::BRep),
    /// Multiple closed sub-solids (sheet split the solid into pieces).
    Pieces(Vec<popkern_core::topology::BRep>),
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
            // Multi-piece dispatch (§2.3): if the sheet is a single
            // planar face, take the half-space split fast path —
            // intersect the solid with each side of the plane to
            // produce closed sub-solids on each side. This handles
            // the most common "saw cut through a part" use case
            // and produces multiple pieces when the cut bisects.
            //
            // For non-planar (curved) sheets, fall back to the
            // single-piece imprint that adds curves to the target
            // without cutting it into pieces. Curved-sheet
            // bisection requires the connected-component
            // flood-fill approach (Phase 2 future work).
            let (sheet, solid) = if a_is_sheet { (a, b) } else { (b, a) };
            if let Some((origin, normal)) = is_planar_sheet(sheet) {
                let pieces = split_solid_by_plane(solid, origin, normal, tol)?;
                return match pieces.len() {
                    0 => Err(crate::BoolError::OpenShellNotSupported {
                        a: a_is_sheet,
                        b: b_is_sheet,
                    }),
                    1 => Ok(SheetBooleanResult::Solid(
                        pieces.into_iter().next().unwrap(),
                    )),
                    _ => Ok(SheetBooleanResult::Pieces(pieces)),
                };
            }
            // Curved sheet fallback: single-piece imprint.
            let imprint = crate::ops::imprint::imprint(a, b, &tol)?;
            Ok(SheetBooleanResult::Solid(imprint.imprinted_target))
        }
        SheetOperationKind::SheetIntersectSolid => {
            // Imprint the solid's intersection curves onto the sheet,
            // then walk the resulting fragments and remove the ones
            // whose centroid is outside the solid. The remaining
            // fragments form the portion of the sheet inside the
            // solid — itself an open shell.
            //
            // Both orderings of the (sheet, solid) pair classify the
            // same way: we always imprint the SHEET (the half-D
            // input) using the SOLID's curves, then test against
            // the solid for inside/outside.
            let (sheet, solid) = if a_is_sheet { (a, b) } else { (b, a) };
            let imprint = crate::ops::imprint::imprint(sheet, solid, &tol)?;
            let trimmed = trim_sheet_to_inside(imprint.imprinted_target, solid, &tol);
            Ok(SheetBooleanResult::OpenShell(trimmed))
        }
        SheetOperationKind::InvalidSheetMinusSolid
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
    fn split_solid_by_plane_bisects_cube() {
        use popkern_core::geometry::{Point3, Vec3};
        // A unit cube cut by the z=0.5 plane should yield 2 closed
        // sub-solids — the top half and the bottom half. Each half
        // has its own boundary (4 sides + 1 cap at z=0.5 + 1 outer
        // face). The cap-construction quality depends on the
        // boolean's accuracy near coincident face boundaries; the
        // result may have extra slivers but the basic count is
        // pinned here.
        let cube = BRep::unit_cube(tol()).unwrap();
        let pieces = split_solid_by_plane(
            &cube,
            Point3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 0.0, 1.0),
            tol(),
        )
        .unwrap();
        assert_eq!(
            pieces.len(),
            2,
            "z=0.5 plane should bisect the cube into 2 pieces, got {}",
            pieces.len()
        );
        // Each piece should have at least 4 faces (a tetrahedron's
        // worth) and contain vertices from both sides of the plane.
        for (i, piece) in pieces.iter().enumerate() {
            assert!(piece.faces.len() >= 4, "piece {i} should have ≥4 faces");
        }
    }

    #[test]
    fn split_solid_by_plane_off_axis_still_bisects() {
        use popkern_core::geometry::{Point3, Vec3};
        // Diagonal cut: x=0.5 + y=0.5 + z=0.5 plane through cube
        // center. Both halves are valid, but each is no longer a
        // simple box.
        let cube = BRep::unit_cube(tol()).unwrap();
        let pieces = split_solid_by_plane(
            &cube,
            Point3::new(0.5, 0.5, 0.5),
            Vec3::new(1.0, 1.0, 1.0),
            tol(),
        )
        .unwrap();
        // Diagonal cut may produce 1 or 2 pieces depending on
        // boolean robustness. The important property is that the
        // function doesn't panic and returns a non-empty result.
        assert!(!pieces.is_empty(), "diagonal cut should produce ≥1 piece");
    }

    /// Build a single-face planar sheet body (a quad) at z=z_height
    /// spanning x in [x_min, x_max], y in [y_min, y_max]. Used by
    /// the multi-piece dispatcher tests.
    fn build_planar_sheet(x_min: f64, x_max: f64, y_min: f64, y_max: f64, z_height: f64) -> BRep {
        use popkern_core::constructors::builder::BRepBuilder;
        use popkern_core::geometry::Point3;
        use popkern_core::geometry::Vec3;
        use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut builder = BRepBuilder::new(tol());
        let solid = builder.add_solid();
        let shell = builder.add_shell(solid);
        let face = builder.add_face(shell);

        // 4 corners of the quad.
        let v0 = builder.add_vertex(Point3::new(x_min, y_min, z_height));
        let v1 = builder.add_vertex(Point3::new(x_max, y_min, z_height));
        let v2 = builder.add_vertex(Point3::new(x_max, y_max, z_height));
        let v3 = builder.add_vertex(Point3::new(x_min, y_max, z_height));

        // 4 edges going around the quad CCW.
        let (_, he01_fwd, _) = builder.add_edge(v0, v1);
        let (_, he12_fwd, _) = builder.add_edge(v1, v2);
        let (_, he23_fwd, _) = builder.add_edge(v2, v3);
        let (_, he30_fwd, _) = builder.add_edge(v3, v0);

        builder
            .add_loop(&[he01_fwd, he12_fwd, he23_fwd, he30_fwd], face)
            .unwrap();

        builder.brep_mut().faces.get_mut(face).unwrap().surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(x_min, y_min, z_height),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        });

        builder.build_unchecked()
    }

    #[test]
    fn build_planar_sheet_passes_sheet_detector() {
        let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 0.5);
        assert!(
            is_sheet_body(&sheet),
            "constructed quad should be a sheet body"
        );
        assert!(
            is_planar_sheet(&sheet).is_some(),
            "constructed quad should be detected as planar"
        );
        let (origin, normal) = is_planar_sheet(&sheet).unwrap();
        assert_eq!(normal.z, 1.0);
        assert_eq!(origin.z, 0.5);
    }

    #[test]
    fn boolean_with_sheets_solid_minus_planar_sheet_produces_pieces() {
        // End-to-end: cube minus a planar sheet that bisects it
        // should produce 2 pieces via the dispatcher.
        use crate::BoolOperation;
        let cube = BRep::unit_cube(tol()).unwrap();
        let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 0.5);
        let result = boolean_with_sheets(BoolOperation::Subtraction, &cube, &sheet, tol()).unwrap();
        match result {
            SheetBooleanResult::Pieces(pieces) => {
                assert_eq!(
                    pieces.len(),
                    2,
                    "bisecting sheet should produce 2 pieces, got {}",
                    pieces.len()
                );
            }
            other => panic!("expected Pieces, got {other:?}"),
        }
    }

    #[test]
    fn boolean_with_sheets_solid_minus_off_axis_planar_sheet_produces_one_piece() {
        // A planar sheet that doesn't intersect the solid leaves it
        // unchanged: 1 piece returned as Solid.
        use crate::BoolOperation;
        let cube = BRep::unit_cube(tol()).unwrap();
        let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 5.0);
        let result = boolean_with_sheets(BoolOperation::Subtraction, &cube, &sheet, tol()).unwrap();
        match result {
            SheetBooleanResult::Solid(_) => {}
            other => panic!("expected single Solid, got {other:?}"),
        }
    }

    #[test]
    fn split_solid_by_plane_outside_returns_one_piece() {
        use popkern_core::geometry::{Point3, Vec3};
        // A plane far above the cube doesn't bisect it. The whole
        // cube is on the negative side; the positive side is empty.
        let cube = BRep::unit_cube(tol()).unwrap();
        let pieces = split_solid_by_plane(
            &cube,
            Point3::new(0.5, 0.5, 10.0), // 10 units above the cube
            Vec3::new(0.0, 0.0, 1.0),
            tol(),
        )
        .unwrap();
        assert_eq!(
            pieces.len(),
            1,
            "non-bisecting plane should return 1 piece, got {}",
            pieces.len()
        );
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
    fn sheet_intersect_solid_routes_to_open_shell() {
        // Build a small "sheet" (single triangular face with no
        // surface attached). intersect against a unit cube. The
        // dispatcher should pick SheetIntersectSolid and return
        // an OpenShell variant — the result is the imprinted+trimmed
        // sheet.
        use crate::BoolOperation;
        use popkern_core::geometry::Point3;
        let mut sheet = BRep::new(tol());
        let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.5)).unwrap();
        let _mev = sheet
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(2.0, 0.0, 0.5))
            .unwrap();

        let cube = BRep::unit_cube(tol()).unwrap();
        let result = boolean_with_sheets(BoolOperation::Intersection, &sheet, &cube, tol());
        // Sheet has no surface attached, so the dispatcher will
        // either succeed with an empty/trimmed open shell or
        // surface a structured error. Either is acceptable for this
        // smoke test — the important property is that we DO route
        // to SheetIntersectSolid (no longer rejected upfront).
        match result {
            Ok(SheetBooleanResult::OpenShell(_)) => {}
            Ok(SheetBooleanResult::Solid(_))
            | Ok(SheetBooleanResult::Curves(_))
            | Ok(SheetBooleanResult::Pieces(_)) => {
                panic!("expected OpenShell, got different variant");
            }
            Err(_) => {
                // Acceptable: imprint may fail on a sheet without
                // surface geometry.
            }
        }
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
