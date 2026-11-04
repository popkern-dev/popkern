//! The body of [`crate::boolean`] and [`crate::boolean_with_options`]:
//! input checks, the pipeline under a non-manifold policy, and the
//! naming and lineage of the result.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;
use crate::{BoolOperation, TrackingContext, entity_naming_result, operation, options, sheet};

/// [`crate::boolean`] under `policy`.
pub(crate) fn boolean_under(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
    policy: options::NonManifoldPolicy,
) -> Result<BRep, BoolError> {
    // §2.3: reject sheet body inputs with a structured error rather
    // than letting volume-based classification produce undefined
    // results. Future work will add boundary-relative semantics; for
    // now this is the only safe response.
    let a_is_sheet = sheet::is_sheet_body(&a);
    let b_is_sheet = sheet::is_sheet_body(&b);
    if a_is_sheet || b_is_sheet {
        return Err(BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        });
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
