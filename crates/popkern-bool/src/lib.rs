#![warn(missing_docs)]

//! popkern-bool: Boolean set operations on BRep solids.
//!
//! Provides union, intersection, and subtraction operations on
//! [`BRep`](popkern_core::topology::BRep) solid models.
//!
//! # Pipeline
//!
//! Each boolean operation proceeds through five stages:
//!
//! 1. **Intersect** — compute face-face intersection curves.
//! 2. **Split** — trim faces along those curves.
//! 3. **Classify** — label each fragment as inside, outside, or on-boundary.
//! 4. **Select** — choose which fragments belong to the result.
//! 5. **Merge** — stitch fragments into a valid closed solid.
//!
//! # Example (placeholder — not yet implemented)
//!
//! ```rust,ignore
//! use popkern_bool::{BoolOperation, boolean};
//! use popkern_core::tolerance::Tolerance;
//!
//! let result = boolean(BoolOperation::Union, solid_a, solid_b, Tolerance::default())?;
//! ```

pub(crate) mod classify;
pub mod csg;
pub mod error;
pub(crate) mod intersect;
pub(crate) mod merge;
pub(crate) mod operation;
pub(crate) mod select;
pub(crate) mod split;

pub use classify::point_in_solid;
pub use csg::{CsgNode, csg_eval};

pub use error::{
    BoolDiagnostics, BoolError, BoolOperation, Classification, PointClassification, Stage,
};

#[cfg(test)]
mod integration_tests;

use popkern_core::{tolerance::Tolerance, topology::BRep};

/// Perform a boolean operation on two solids, consuming both inputs.
///
/// Returns the result solid on success, or a [`BoolError`] describing the failure.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
pub fn boolean(op: BoolOperation, a: BRep, b: BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    operation::execute_boolean(op, a, b, &tol)
}

/// Perform a boolean operation on two solids, borrowing both inputs.
///
/// Clones `a` and `b` internally before operating. Prefer [`boolean`] when
/// ownership can be transferred to avoid the clone cost.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
pub fn boolean_ref(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> Result<BRep, BoolError> {
    operation::execute_boolean(op, a.clone(), b.clone(), &tol)
}

/// Compute the union of two solids (A ∪ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Union, a, b, tol)`.
pub fn union(a: BRep, b: BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Union, a, b, tol)
}

/// Compute the union of two solids (A ∪ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Union, a, b, tol)`.
pub fn union_ref(a: &BRep, b: &BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Union, a, b, tol)
}

/// Compute the intersection of two solids (A ∩ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Intersection, a, b, tol)`.
pub fn intersection(a: BRep, b: BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Intersection, a, b, tol)
}

/// Compute the intersection of two solids (A ∩ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Intersection, a, b, tol)`.
pub fn intersection_ref(a: &BRep, b: &BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Intersection, a, b, tol)
}

/// Subtract `b` from `a` (A − B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Subtraction, a, b, tol)`.
pub fn subtraction(a: BRep, b: BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Subtraction, a, b, tol)
}

/// Subtract `b` from `a` (A − B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Subtraction, a, b, tol)`.
pub fn subtraction_ref(a: &BRep, b: &BRep, tol: Tolerance) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Subtraction, a, b, tol)
}
