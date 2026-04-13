//! Sheet-aware boolean dispatcher (§2.3).
//!
//! [`boolean_with_sheets`] is the public entry point for sheet-aware
//! boolean operations. It detects whether each input is a sheet body,
//! classifies the operation via [`crate::sheet::classify_sheet_operation`],
//! and routes to the appropriate sub-operation — producing a closed
//! solid, multiple closed sub-solids, an open shell, or a list of
//! intersection curves depending on the case.

use super::classify::{SheetOperationKind, classify_sheet_operation};
use super::detect::{is_planar_sheet, is_sheet_body};
use super::split::split_solid_by_plane;
use super::trim::trim_sheet_to_inside;
use popkern_core::topology::BRep;

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
    Solid(BRep),
    /// Multiple closed sub-solids (sheet split the solid into pieces).
    Pieces(Vec<BRep>),
    /// A list of intersection curves (no enclosing volume).
    Curves(Vec<popkern_core::geometry::curve_on_surface::CurveOnSurface>),
    /// An open shell BRep result.
    OpenShell(BRep),
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
            // Multi-piece dispatch (§2.3): if the sheet is a single
            // planar face, take the half-space split fast path —
            // intersect the solid with each side of the plane to
            // produce closed sub-solids on each side.
            //
            // For non-planar (curved) sheets, fall back to the
            // single-piece imprint that adds curves to the target
            // without cutting it into pieces. Curved-sheet bisection
            // requires the Phase 2 approach in super::split_curved.
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
