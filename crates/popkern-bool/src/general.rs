//! `NonManifoldPolicy::General`: a boolean's non-manifold result served
//! as one general body (PK-1106, PK-1115).
//!
//! The split pipeline serves the result as manifold pieces that touch
//! without sharing (an intermediate); gluing them where they touch,
//! from their geometry alone (`BRep::glue_touching`), makes the general
//! body those pieces are. A general operand enters the same way in
//! reverse: taken apart into its manifold pieces
//! (`BRep::separate_at_contacts`), which the glue closes again.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;

/// `body` as an operand the split pipeline reads: a general body taken
/// apart into its manifold pieces, any other kind as it is.
pub(crate) fn operand_pieces(mut body: BRep, operand: usize) -> Result<BRep, BoolError> {
    if body.body_kind() != BodyKind::General {
        return Ok(body);
    }
    body.separate_at_contacts()
        .map_err(|_| BoolError::UnsupportedBodyKind {
            operation: "boolean",
            operand,
            kind: BodyKind::General,
            accepted: &[BodyKind::Solid, BodyKind::General],
        })?;
    Ok(body)
}

/// The pieces `pieces` glued where they touch into one general body,
/// served when it validates; pieces that touch nowhere as they are.
///
/// # Errors
///
/// [`BoolError::GeneralBodyDeclined`] with the glue's decline or the
/// glued body's first validation issue (a `GeneralDeclined` shadow
/// line too).
pub(crate) fn glue_result(pieces: BRep, tolerance: &Tolerance) -> Result<BRep, BoolError> {
    if crate::postprocess::split_pieces::touch_contacts(&pieces, tolerance).is_empty() {
        return Ok(pieces);
    }
    let mut body = pieces;
    let reason = match body.glue_touching() {
        Ok(_) => match body.validate_all().first() {
            None => return Ok(body),
            Some(issue) => format!("glued body invalid: {issue}"),
        },
        Err(e) => e.to_string(),
    };
    crate::complex_shadow::append_log_line(&format!("GeneralDeclined: {reason}"));
    Err(BoolError::GeneralBodyDeclined { reason })
}

#[cfg(test)]
#[path = "general_tests.rs"]
mod tests;
