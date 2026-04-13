//! Sheet body / open shell support (§2.3).
//!
//! A sheet body is a BRep whose boundary is not closed: at least one
//! edge is "free" — used by only one face. The boolean engine's
//! classification logic assumes a closed volume (point-set membership)
//! and produces undefined results on sheet bodies. This module
//! provides:
//!
//! - [`detect::is_sheet_body`] — the detector used upfront by the
//!   public boolean entry points to reject sheet inputs with a
//!   structured [`crate::BoolError::OpenShellNotSupported`] error
//!   rather than silently corrupting the result.
//! - [`side::side_of_sheet`] — boundary-relative classification for
//!   point-vs-sheet queries (§6.x).
//! - [`classify::classify_sheet_operation`] — operation-matrix
//!   classification (§2.3) of a (BoolOperation, a_is_sheet, b_is_sheet)
//!   tuple to a [`classify::SheetOperationKind`].
//! - [`dispatch::boolean_with_sheets`] — sheet-aware boolean
//!   operation that routes each classified kind to the right
//!   sub-operation: standard boolean for solid-solid, half-space
//!   split for planar solid-minus-sheet, imprint+trim for
//!   sheet-intersect-solid, and structured rejection for the
//!   not-yet-implemented cases.
//! - [`split::split_solid_by_plane`] — half-space split helper used
//!   by the multi-piece `SolidMinusSheet` fast path.

pub mod classify;
pub mod detect;
pub mod dispatch;
pub mod side;
pub mod split;
mod trim;

pub use classify::{SheetOperationKind, classify_sheet_operation};
pub use detect::{is_approximately_planar_sheet, is_closed_solid, is_planar_sheet, is_sheet_body};
pub use dispatch::{SheetBooleanResult, boolean_with_sheets};
pub use side::{SideOfSheet, side_of_sheet};
pub use split::split_solid_by_plane;

#[cfg(test)]
mod tests;
