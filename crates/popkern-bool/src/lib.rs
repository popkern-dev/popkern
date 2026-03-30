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
pub(crate) mod clip;
pub mod csg;
pub mod error;
pub(crate) mod intersect;
pub(crate) mod merge;
pub(crate) mod operation;
pub(crate) mod select;
pub(crate) mod split;
pub(crate) mod split_inner_loop;

pub use classify::point_in_solid;
pub use csg::{CsgNode, csg_eval};

pub use error::{
    BoolDiagnostics, BoolError, BoolOperation, Classification, PointClassification, Stage,
};

#[cfg(test)]
mod integration_tests;
#[cfg(test)]
mod integration_tests_curved;
#[cfg(test)]
mod integration_tests_curved_regression;
#[cfg(test)]
mod integration_tests_multicurve;
#[cfg(test)]
mod translate_helper_tests;

use popkern_core::history::TrackingContext;
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
pub fn boolean(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    let mut result = operation::execute_boolean(op, a, b, &tol)?;
    if let Some(ctx) = tracking {
        assign_boolean_names(&mut result, ctx);
    }
    Ok(result)
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
/// * `tracking` — optional tracking context for stable naming.
pub fn boolean_ref(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(op, a.clone(), b.clone(), tol, tracking)
}

/// Compute the union of two solids (A ∪ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the union of two solids (A ∪ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Assign stable names to all faces in a boolean result.
///
/// All surviving faces are assigned `Generated` names since the boolean
/// creates an entirely new solid topology.
fn assign_boolean_names(brep: &mut BRep, ctx: &mut TrackingContext) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let face_handles: Vec<_> = brep.faces.keys().collect();
    for (i, fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[*fh].name = Some(name_id);
        ctx.record_generated(name);
    }
}
