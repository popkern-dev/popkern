//! Contact report and symmetric difference (PK-964 Stages 8 and 9).
//!
//! [`contact_report`] says where two bodies touch, and how: at a point,
//! along a curve, or over an area. It is the question behind every
//! non-manifold decline (an edge-touching union declines because the
//! bodies touch along an edge) asked directly, so a GUI can show the
//! contact before the user runs the boolean.
//!
//! [`symmetric_difference`] is A ⊕ B = (A − B) ∪ (B − A): the material
//! in exactly one operand. Its two pieces meet along every curve where
//! the operands' boundaries cross, so a ⊕ of overlapping bodies is
//! non-manifold by nature: it declines under the default policy and
//! serves as the pieces with recorded contacts under
//! [`NonManifoldPolicy::SplitToManifold`].

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolError, BoolOperation};
use crate::options::{BoolOptions, NonManifoldPolicy};
use crate::postprocess::shell_contact::{ContactExtent, shell_findings};

/// How far a touch between two bodies extends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TouchExtent {
    /// An isolated point.
    Point,
    /// A curve through the reported point.
    Curve,
    /// An area around the reported point (flush faces).
    Area,
}

/// One place two bodies touch: a face of each, a point of the touch,
/// and its extent.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyContact {
    /// The face of the first body.
    pub face_a: FaceHandle,
    /// The face of the second body.
    pub face_b: FaceHandle,
    /// A point where they touch.
    pub location: Point3,
    /// Point, curve or area.
    pub extent: TouchExtent,
}

/// Where bodies `a` and `b` touch (within the linear tolerance), one
/// entry per touching face pair, in face order. Bodies that overlap in
/// volume also report the faces that cross (their crossings are
/// curves); use [`crate::interference_check`] to tell overlap from
/// touch. Empty when they are apart.
///
/// # Errors
///
/// [`BoolError::MergeFailed`] when the bodies cannot be set side by
/// side for the query (a structurally broken operand).
pub fn contact_report(
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Result<Vec<BodyContact>, BoolError> {
    let (side_by_side, b_map) =
        crate::operation::helpers::copy_brep_into_mapped(a.clone(), b.clone(), false)?;
    // `a`'s handles survive the copy unchanged (versioned keys: a slot
    // `b`'s copy reuses carries a new version); `b`'s map back.
    let b_back: std::collections::HashMap<FaceHandle, FaceHandle> =
        b_map.into_iter().map(|(old, new)| (new, old)).collect();
    let back = |f: FaceHandle| -> Option<(bool, FaceHandle)> {
        if a.faces.contains_key(f) {
            Some((true, f))
        } else {
            b_back.get(&f).map(|&g| (false, g))
        }
    };
    let mut out = Vec::new();
    for finding in shell_findings(&side_by_side, tolerance) {
        let (Some(first), Some(second)) = (back(finding.faces.0), back(finding.faces.1)) else {
            continue;
        };
        // Only pairs across the two bodies; a body's own shells are its
        // own business.
        let (face_a, face_b) = match (first, second) {
            ((true, fa), (false, fb)) | ((false, fb), (true, fa)) => (fa, fb),
            _ => continue,
        };
        out.push(BodyContact {
            face_a,
            face_b,
            location: finding.contact.location,
            extent: match finding.extent {
                ContactExtent::Point => TouchExtent::Point,
                ContactExtent::Line => TouchExtent::Curve,
                ContactExtent::Area => TouchExtent::Area,
            },
        });
    }
    Ok(out)
}

/// A ⊕ B: the material in exactly one of `a` and `b` (see the module
/// docs for why its pieces meet with zero thickness).
///
/// # Errors
///
/// Either difference's failure; [`BoolError::NonManifoldResult`] when
/// the pieces touch and `options` declines non-manifold results (the
/// default), or when they cannot be served valid side by side.
pub fn symmetric_difference(
    a: BRep,
    b: BRep,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let a_only = crate::boolean_with_options(
        BoolOperation::Subtraction,
        a.clone(),
        b.clone(),
        tolerance,
        None,
        options,
    )?;
    let b_only =
        crate::boolean_with_options(BoolOperation::Subtraction, b, a, tolerance, None, options)?;
    if a_only.faces.is_empty() {
        return Ok(b_only);
    }
    if b_only.faces.is_empty() {
        return Ok(a_only);
    }
    let pieces = crate::operation::helpers::combine_breps(a_only, b_only)?;
    let contacts = crate::postprocess::shell_contact::shell_contacts(&pieces, &tolerance);
    if contacts.is_empty() {
        return Ok(pieces);
    }
    match options.non_manifold {
        NonManifoldPolicy::Decline => Err(BoolError::NonManifoldResult { contacts }),
        NonManifoldPolicy::SplitToManifold => {
            crate::postprocess::contact_record::accept_split(pieces, &tolerance)
        }
    }
}

#[cfg(test)]
#[path = "contact_tests.rs"]
mod tests;
