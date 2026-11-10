//! The body of [`crate::boolean`] and [`crate::boolean_with_options`]:
//! input checks, the pipeline under a non-manifold policy, and the
//! naming and lineage of the result.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;
use crate::{BoolOperation, TrackingContext, entity_naming_result, operation, options};

/// [`crate::boolean`] under `policy`, a non-manifold result as the
/// pipeline answers it: manifold pieces in one body under the split and
/// general policies (`crate::split_serve` serves them).
pub(crate) fn boolean_under(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
    policy: options::NonManifoldPolicy,
) -> Result<BRep, BoolError> {
    // A boolean of solids (empty operands go on to the pre-check's
    // "no faces" finding). A sheet is refused with its own error, whose
    // remediation names the sheet booleans (§2.3); any other kind (a
    // general body, a wire) with the kind gate's.
    // The general policy reads general operands as their manifold
    // pieces and runs split; the entry point glues the result.
    let (a, b) = if policy == options::NonManifoldPolicy::General {
        (
            crate::general::operand_pieces(a, 0)?,
            crate::general::operand_pieces(b, 1)?,
        )
    } else {
        (a, b)
    };
    let kinds = [a.body_kind(), b.body_kind()];
    let a_is_sheet = kinds[0] == BodyKind::Sheet;
    let b_is_sheet = kinds[1] == BodyKind::Sheet;
    if a_is_sheet || b_is_sheet {
        return Err(BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        });
    }
    for (operand, kind) in kinds.into_iter().enumerate() {
        if !matches!(kind, BodyKind::Solid | BodyKind::Empty) {
            return Err(BoolError::UnsupportedBodyKind {
                operation: "boolean",
                operand,
                kind,
                accepted: &[BodyKind::Solid],
            });
        }
    }

    // An operand thinner than the tolerance has no inside at its
    // resolution: the result is the other operand or empty (PK-1131),
    // reported as a sub-tolerance feature.
    for (operand, body) in [&a, &b].into_iter().enumerate() {
        if let Some(thin) = crate::sub_tolerance::thin_operand(body, &tol) {
            crate::sub_tolerance::record(crate::sub_tolerance::SubToleranceFeature {
                kind: crate::sub_tolerance::SubToleranceKind::OperandVanished { operand },
                location: thin.location,
                measured_width: thin.width,
                tolerance: tol.linear(),
            });
            return Ok(crate::sub_tolerance::resolve(op, &a, &b, operand));
        }
    }

    // No perturbation retry: plan Stage 7d retired it after it fired on
    // no boolean over three whole-suite runs (`boolean_with_perturbation`
    // remains an explicit opt-in).
    let mut result = operation::execute_boolean_with_policy(op, &a, &b, &tol, policy)?;
    if let Some(ctx) = tracking {
        entity_naming_result::name_result(&mut result, &a, &[&b], ctx, &tol);
    }
    Ok(result)
}

/// [`crate::boolean_with_options`] as a list of bodies: one body, except under
/// [`crate::NonManifoldPolicy::SplitToManifold`] a result whose manifold
/// pieces touch, each piece its own valid body (in slot order, named as
/// the one body would be).
///
/// # Errors
///
/// As [`crate::boolean_with_options`], never [`BoolError::SplitIntoBodies`].
pub fn boolean_bodies(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
    options: &crate::BoolOptions,
) -> Result<Vec<BRep>, BoolError> {
    let pieces = boolean_under(op, a, b, tol, tracking, options.non_manifold);
    crate::split_serve::serve_bodies(pieces, options.non_manifold, &tol)
}
