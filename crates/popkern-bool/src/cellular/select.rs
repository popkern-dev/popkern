//! Cell selection predicates for boolean operations.
//!
//! Boolean operations reduce to selecting cells from a [`CellDecomposition`]
//! by a predicate on their membership bits. Union, intersection, subtraction,
//! and symmetric difference are built-in; custom predicates are also supported.

use crate::error::BoolOperation;

use super::cell::{CellDecomposition, CellId};

/// Select cells for a standard boolean operation.
pub fn select_for_operation(
    decomposition: &CellDecomposition,
    operation: BoolOperation,
) -> Vec<CellId> {
    match operation {
        BoolOperation::Union => decomposition.select_union(),
        BoolOperation::Intersection => decomposition.select_intersection(),
        BoolOperation::Subtraction => decomposition.select_subtraction(),
    }
}

/// Select cells using a custom predicate on membership bits.
///
/// The predicate receives a `u64` where bit 0 = body A, bit 1 = body B.
/// Returns the cell IDs where the predicate returns true.
///
/// # Example
///
/// ```ignore
/// // Symmetric difference: cells in A XOR B
/// let cells = select_custom(&decomposition, |m| (m & 1) ^ (m & 2) != 0);
/// ```
pub fn select_custom(
    decomposition: &CellDecomposition,
    predicate: impl Fn(u64) -> bool,
) -> Vec<CellId> {
    decomposition.select_cells(predicate)
}
