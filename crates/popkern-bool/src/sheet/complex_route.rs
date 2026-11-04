//! Sheet × solid booleans decided by the cell complex (PK-978).
//!
//! The sheet enters the partition as a lamina: sheet ∩ solid and
//! sheet − solid keep its faces inside or outside the solid
//! (`Partition::lamina_part`), solid − sheet keeps the solid's regions
//! as pieces (`split`). A decline is returned to the dispatcher, which
//! falls back to the planar split or the imprint-and-trim path where
//! one exists and records the fallback as a `PartitionFallback` shadow
//! line.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;

/// The faces of `sheet` inside (`inside`) or outside `solid`.
pub(super) fn sheet_part(
    sheet: &BRep,
    solid: &BRep,
    inside: bool,
    tolerance: Tolerance,
) -> Result<BRep, BoolError> {
    let p = crate::partition(&[solid, sheet], tolerance)?;
    p.lamina_part(1, |bodies| bodies.contains(0) == inside)
}

/// `solid` cut along `sheet`: its pieces, one when the sheet divides
/// nothing.
pub(super) fn solid_pieces(
    solid: &BRep,
    sheet: &BRep,
    tolerance: Tolerance,
) -> Result<Vec<BRep>, BoolError> {
    let cut = crate::split(solid, &[sheet], tolerance)?;
    Ok(cut.pieces.into_iter().map(|p| p.body).collect())
}

/// Record that `kind` fell back from the complex after `decline`.
pub(super) fn log_fallback(kind: super::SheetOperationKind, decline: &BoolError) {
    crate::complex_shadow::append_log_line(&format!("PartitionFallback sheet {kind:?}: {decline}"));
}

/// A one-face sheet for the constructors' two-face form (a front face
/// and its backside closing a zero-volume shell); `None` for any other
/// body.
pub(super) fn as_sheet(body: &BRep) -> Option<BRep> {
    super::lamina::single_face_sheet(body)
}
