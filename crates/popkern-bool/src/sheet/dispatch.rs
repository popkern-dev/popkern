//! Sheet-aware boolean dispatcher (§2.3).
//!
//! [`boolean_with_sheets`] is the public entry point for sheet-aware
//! boolean operations. It detects whether each input is a sheet body,
//! classifies the operation via [`crate::sheet::classify_sheet_operation`],
//! and routes to the appropriate sub-operation — producing a closed
//! solid, multiple closed sub-solids, an open shell, or a list of
//! intersection curves depending on the case.

use super::classify::{SheetOperationKind, classify_sheet_operation};
use super::detect::{is_approximately_planar_sheet, is_planar_sheet, is_sheet_body};
use super::split::split_solid_by_plane;
use super::trim::{KeepSide, trim_sheet};
use popkern_core::topology::BRep;

/// Result of a [`boolean_with_sheets`] call.
///
/// Sheet-aware boolean operations can produce several different shapes
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
/// - [`SheetBooleanResult::Mixed`] — a closed solid paired with an
///   open-shell "appendage" (e.g. solid-union-sheet returns the
///   original solid plus the portion of the sheet that lies outside
///   the solid).
#[derive(Debug)]
pub enum SheetBooleanResult {
    /// A single closed solid result.
    Solid(BRep),
    /// Multiple closed sub-solids (sheet split the solid into pieces).
    Pieces(Vec<BRep>),
    /// A list of intersection curves (no enclosing volume).
    Curves(Vec<popkern_core::geometry::curve_on_surface::CurveOnSurface>),
    /// An open shell BRep result.
    ///
    /// For `SheetSheetUnion` with a transversal crossing input, the
    /// resulting shell contains BOTH input shells concatenated into
    /// one BRep. The intersection curve is NOT welded into shared
    /// vertices/edges — both shells keep their own copy of the
    /// crossing topology, so the result may fail strict manifold
    /// validation at those points. Consumers that need welded
    /// intersections should post-process via the healing pipeline.
    OpenShell(BRep),
    /// A closed solid paired with an open-shell appendage. Used by
    /// `SolidUnionSheet` to return the original solid unchanged plus
    /// the trimmed portion of the sheet that lies outside the solid.
    /// `open_shell` is boxed so this variant does not inflate the
    /// stack size of the whole enum past the single-`BRep` variants.
    Mixed {
        /// The closed-solid portion of the result.
        solid: BRep,
        /// The open-shell portion — typically a trimmed sheet.
        open_shell: Box<BRep>,
    },
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
/// | `SolidMinusSheet` | Planar sheets: half-space split producing possibly multiple sub-solids. Curved sheets: single-piece imprint fallback (Phase 2). |
/// | `SheetIntersectSolid` | Imprints the solid's curves onto the sheet, then trims to the inside fragments. |
/// | `InvalidSheetMinusSolid` | Rejects with [`crate::BoolError::OpenShellNotSupported`]. |
/// | Other variants | Returns [`crate::BoolError::OpenShellNotSupported`] (future work). |
///
/// Use this entry point when you need sheet support; the standard
/// `boolean()` rejects sheet inputs upfront so non-sheet callers
/// don't accidentally invoke half-implemented sheet operations.
pub fn boolean_with_sheets(
    op: crate::BoolOperation,
    a: &BRep,
    b: &BRep,
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
            // Multi-piece dispatch (§2.3). Two paths:
            //
            // 1. Exact planar fast path: if every face of the sheet
            //    has a SurfaceGeometry::Plane surface on a single
            //    supporting plane, call split_solid_by_plane with
            //    that analytic plane.
            //
            // 2. Approximate planar fallback: if the sheet's vertex
            //    positions lie within `tol.linear() * 1e3` of a best-
            //    fit plane, use that best-fit plane for the split.
            //    This catches NURBS/spherical/cylindrical sheets that
            //    are actually near-flat (slight bow, triangulation
            //    noise). The 1e3 factor over tol.linear() is the
            //    "planarity budget" — larger deviations route to the
            //    single-piece imprint fallback so the user doesn't
            //    silently get a planar cut in place of a curved one.
            //
            // 3. Curved fallback: single-piece imprint that adds the
            //    cut curves to the solid's boundary without splitting
            //    it into sub-solids. Fully general curved splitting
            //    would need a sheet-aware merge pipeline (documented
            //    as future work in docs/reports/coverage/).
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
            let planarity_budget = tol.linear() * 1.0e3;
            if let Some((centroid, normal)) = is_approximately_planar_sheet(sheet, planarity_budget)
            {
                let pieces = split_solid_by_plane(solid, centroid, normal, tol)?;
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
            let (sheet, solid) = if a_is_sheet { (a, b) } else { (b, a) };
            let imprint = crate::ops::imprint::imprint(sheet, solid, &tol)?;
            let trimmed = trim_sheet(imprint.imprinted_target, solid, &tol, KeepSide::Inside);
            Ok(SheetBooleanResult::OpenShell(trimmed))
        }
        SheetOperationKind::SheetSheetUnion => {
            // §2.3 SheetSheetUnion: combine two open shells into one
            // open shell.
            //
            // Four cases:
            //
            // 1. Fully disjoint — no intersection curves, no
            //    coincident face overlaps: structural concatenation
            //    via `combine_breps`. Fast path.
            //
            // 2. Coincident overlap — face pairs that share a
            //    surface region. Drop the duplicate faces from `b`
            //    via `super::dedupe::drop_faces_and_prune` and then
            //    concatenate.
            //
            // 3. Transversal crossing — the two sheets properly
            //    intersect along a curve (no shared region).
            //    Concatenation produces a valid BRep with both
            //    shells intact. The intersection curve is NOT
            //    welded into shared vertices/edges — that would
            //    require a sheet-aware merge pipeline. The result
            //    is usable for visualization and mesh export but
            //    may contain duplicate topology on the crossing
            //    curve; this is documented on `SheetBooleanResult`.
            //
            // 4. Mixed (both curves AND coincident overlaps) —
            //    dedupe the coincident faces and then concatenate.
            let section = crate::ops::section::section(a, b, &tol)?;
            let aabbs_a = crate::intersect::compute_face_aabbs(a);
            let aabbs_b = crate::intersect::compute_face_aabbs(b);
            let bool_tol = crate::tolerance::BoolTolerance::from(tol);
            let coincidence = crate::coincidence::detect::detect_coincidences(
                a, b, &aabbs_a, &aabbs_b, &bool_tol,
            );

            // Collect any coincident `b` face handles so we can drop
            // their duplicates. Empty set when there are no overlaps.
            use std::collections::HashSet;
            let coincident_b_faces: HashSet<popkern_core::topology::FaceHandle> = coincidence
                .face_overlaps
                .iter()
                .map(|fo| fo.face_b)
                .collect();

            // Case 1: fully disjoint — direct concatenation.
            if section.curves.is_empty() && coincident_b_faces.is_empty() {
                let combined = crate::operation::helpers::combine_breps(a.clone(), b.clone())?;
                return Ok(SheetBooleanResult::OpenShell(combined));
            }

            // Cases 2, 3, 4 all funnel through the same code path:
            // dedupe any coincident faces (no-op for pure transversal)
            // and concatenate. Transversal crossings leave the cross
            // curve unwelded, which is documented as a known
            // limitation of this operation.
            let mut b_deduped = b.clone();
            if !coincident_b_faces.is_empty() {
                super::dedupe::drop_faces_and_prune(&mut b_deduped, &coincident_b_faces);
            }
            let combined = crate::operation::helpers::combine_breps(a.clone(), b_deduped)?;
            Ok(SheetBooleanResult::OpenShell(combined))
        }
        SheetOperationKind::SolidUnionSheet => {
            // §2.3 SolidUnionSheet: the solid stays unchanged, and the
            // portion of the sheet that lies OUTSIDE the solid is
            // attached as an open-shell appendage. Mirror of
            // SheetIntersectSolid — same imprint pipeline but the
            // opposite side of the trim result. The returned `Mixed`
            // variant pairs the original solid with the trimmed
            // outside-of-solid fragment.
            let (sheet, solid) = if a_is_sheet { (a, b) } else { (b, a) };
            let imprint = crate::ops::imprint::imprint(sheet, solid, &tol)?;
            let outside = trim_sheet(imprint.imprinted_target, solid, &tol, KeepSide::Outside);
            Ok(SheetBooleanResult::Mixed {
                solid: solid.clone(),
                open_shell: Box::new(outside),
            })
        }
        SheetOperationKind::InvalidSheetMinusSolid | SheetOperationKind::Unsupported => {
            Err(crate::BoolError::OpenShellNotSupported {
                a: a_is_sheet,
                b: b_is_sheet,
            })
        }
    }
}
