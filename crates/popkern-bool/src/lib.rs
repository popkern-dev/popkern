#![warn(missing_docs)]
#![allow(clippy::result_large_err)]

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

pub mod cellular;
pub(crate) mod classify;
pub mod classify_on_face;
pub mod classify_propagation;
pub mod classify_winding;
pub(crate) mod clip;
pub mod coincidence;
pub mod csg;
pub mod error;
pub mod imprint;
pub(crate) mod intersect;
pub(crate) mod intersect_hierarchy;
pub(crate) mod merge;
pub(crate) mod operation;
pub(crate) mod select;
pub(crate) mod split;
pub(crate) mod split_inner_loop;
pub mod tolerance;
pub mod validate;

pub use classify::point_in_solid;
pub use csg::{CsgNode, csg_eval};

pub use error::{
    BoolDiagnostics, BoolError, BoolOperation, Classification, PointClassification, Stage,
};

#[cfg(test)]
mod error_and_tracking_tests;
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

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
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
    // Collect input face data before the boolean consumes the inputs.
    let input_a_faces = if tracking.is_some() {
        collect_face_data(&a)
    } else {
        Vec::new()
    };
    let input_b_faces = if tracking.is_some() {
        collect_face_data(&b)
    } else {
        Vec::new()
    };

    let mut result = operation::execute_boolean(op, a, b, &tol)?;
    if let Some(ctx) = tracking {
        assign_boolean_names(&mut result, &input_a_faces, &input_b_faces, ctx, &tol);
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

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
///
/// For faces without an existing name, creates a placeholder name so that all
/// faces can participate in the matching step.
fn collect_face_data(brep: &BRep) -> Vec<(EntityName, Point3, Vec3)> {
    use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
    use popkern_core::naming_criteria::capture_face;

    let mut result = Vec::new();
    for (i, (fh, _face)) in brep.faces.iter().enumerate() {
        let Some(fc) = capture_face(brep, fh) else {
            continue;
        };
        let name = if let Some(name_id) = brep.faces[fh].name {
            if let Some(n) = brep.name_pool.get(name_id) {
                n.clone()
            } else {
                // Fallback placeholder for faces with an invalid name id.
                EntityName {
                    operation: OperationId(0),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: i as u32,
                    },
                    ancestors: Vec::new(),
                }
            }
        } else {
            // Placeholder for unnamed faces.
            EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            }
        };
        result.push((name, fc.centroid, fc.normal));
    }
    result
}

/// Assign stable names to faces in a boolean result by matching against inputs.
///
/// Output faces that geometrically match an input face (by centroid proximity
/// and normal alignment) receive `Modified` names with the input face as
/// ancestor. Genuinely new intersection faces receive `Generated` names.
fn assign_boolean_names(
    brep: &mut BRep,
    input_a_faces: &[(EntityName, Point3, Vec3)],
    input_b_faces: &[(EntityName, Point3, Vec3)],
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    use popkern_core::naming::{EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let all_inputs: Vec<&(EntityName, Point3, Vec3)> =
        input_a_faces.iter().chain(input_b_faces.iter()).collect();

    let face_handles: Vec<_> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &face_handles {
        let Some(fc) = capture_face(brep, *fh) else {
            continue;
        };
        let centroid = fc.centroid;
        let normal = fc.normal;

        // Find the closest input face by centroid + normal alignment.
        let generous_dist = tol.linear * 100.0;
        let best_match = all_inputs
            .iter()
            .filter(|(_, ic, in_)| {
                centroid.distance_to(ic) < generous_dist && normal.dot(in_).abs() > 0.9
            })
            .min_by(|a, b| {
                centroid
                    .distance_to(&a.1)
                    .partial_cmp(&centroid.distance_to(&b.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some((input_name, _, _)) = best_match {
            // Modified: this face descended from an input face.
            let ancestor_id = brep.name_pool.intern(input_name.clone(), None);
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors: vec![ancestor_id],
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_modified(input_name.clone(), new_name);
            mod_idx += 1;
        } else {
            // Generated: genuinely new intersection face.
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}
