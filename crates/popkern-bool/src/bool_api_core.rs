//! The body of [`crate::boolean`] and [`crate::boolean_with_options`]:
//! input checks, the pipeline under a non-manifold policy, and the
//! naming and lineage of the result.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;
use crate::{BoolOperation, TrackingContext, entity_naming_result, operation, options};

/// [`crate::boolean`] under `policy`.
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

    // No perturbation retry: plan Stage 7d retired it after it fired on
    // no boolean over three whole-suite runs (`boolean_with_perturbation`
    // remains an explicit opt-in).
    // Operands served as pieces meeting at recorded contacts touch
    // themselves on purpose (PK-986): the result is served as pieces
    // too, its contacts recorded again, and declined only for a contact
    // no operand recorded.
    use crate::postprocess::contact_inherit::{inherited_or_declined, run_policy};
    let inherits = run_policy(policy, &[&a, &b]);
    let run = inherits.unwrap_or(policy);
    let mut result = operation::execute_boolean_with_policy(op, &a, &b, &tol, run)?;
    if inherits.is_some() {
        result = inherited_or_declined(result, &[&a, &b], &tol)?;
    }
    if let Some(ctx) = tracking {
        entity_naming_result::name_result(&mut result, &a, &[&b], ctx, &tol);
    }
    Ok(result)
}
