//! The body-kind gate at a boolean operation's entry (PK-1101 S5): an
//! operand of a kind the operation does not accept is refused with
//! [`BoolError::UnsupportedBodyKind`], naming the operation, the
//! operand and the kinds it takes.

use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;

/// `body`'s kind when it is one of `accepted`, else the refusal naming
/// `operation` and the operand's index `operand`. Exported as
/// [`crate::require_operand_kind`] for operations built on booleans
/// (holes, hinges), which gate their target before building tools.
///
/// # Errors
///
/// [`BoolError::UnsupportedBodyKind`] when the kind is not accepted.
pub fn require_operand(
    body: &BRep,
    operation: &'static str,
    operand: usize,
    accepted: &'static [BodyKind],
) -> Result<BodyKind, BoolError> {
    let kind = body.body_kind();
    if accepted.contains(&kind) {
        Ok(kind)
    } else {
        Err(BoolError::UnsupportedBodyKind {
            operation,
            operand,
            kind,
            accepted,
        })
    }
}

#[cfg(test)]
#[path = "kind_gate_tests.rs"]
mod tests;
