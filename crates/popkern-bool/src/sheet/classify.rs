//! Sheet-aware operation matrix classification (§2.3).
//!
//! [`classify_sheet_operation`] maps a `(BoolOperation, a_is_sheet,
//! b_is_sheet)` tuple to a [`SheetOperationKind`] variant that the
//! dispatcher uses to route to the correct sub-operation. See
//! [`crate::sheet::boolean_with_sheets`] for the dispatcher itself.

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
