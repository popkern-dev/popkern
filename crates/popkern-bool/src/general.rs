//! `NonManifoldPolicy::General`: a boolean's non-manifold result served
//! as one general body (PK-1106).
//!
//! The split policy's pieces meet at contacts recorded with the
//! coincident entities of every piece; gluing them with the core's
//! non-manifold operators (`BRep::glue_touching`, from their geometry) makes the general body
//! those pieces are, so the result is built from a path the split policy
//! already serves and validates. A general operand enters the same way
//! in reverse: split into manifold pieces meeting at recorded contacts
//! (`BRep::split_to_manifold`), which the boolean inherits and the glue
//! closes again.

use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;

/// `body` as an operand the split pipeline reads: a general body split
/// into its manifold pieces (contacts recorded), any other kind as it
/// is.
pub(crate) fn operand_pieces(mut body: BRep, operand: usize) -> Result<BRep, BoolError> {
    if body.body_kind() != BodyKind::General {
        return Ok(body);
    }
    body.split_to_manifold()
        .map_err(|_| BoolError::UnsupportedBodyKind {
            operation: "boolean",
            operand,
            kind: BodyKind::General,
            accepted: &[BodyKind::Solid, BodyKind::General],
        })?;
    Ok(body)
}

/// The split pieces `pieces` glued where they touch into one general
/// body, kept when it validates with nothing recorded; otherwise the
/// pieces, with a `GeneralFallback` shadow line naming why.
///
/// The glue reads the pieces' geometry alone (`BRep::glue_touching`,
/// PK-1117). Where it declines, the recorded contacts are glued as
/// before (`BRep::glue_contacts`, a `TouchGlueFallback` shadow line),
/// until the stored list retires (PK-1115).
pub(crate) fn glue_result(pieces: BRep) -> BRep {
    if pieces.topology_contacts.is_empty() {
        return pieces;
    }
    let mut body = pieces.clone();
    body.topology_contacts.clear();
    let touch = match body.glue_touching() {
        Ok(_) => {
            let issues = body.validate_all();
            if issues.is_empty() {
                return body;
            }
            format!("touch-glued body invalid: {}", issues[0])
        }
        Err(e) => e.to_string(),
    };
    crate::complex_shadow::append_log_line(&format!("TouchGlueFallback: {touch}"));
    let mut body = pieces.clone();
    let why = match body.glue_contacts() {
        Ok(_) => {
            let issues = body.validate_all();
            if issues.is_empty() {
                return body;
            }
            format!("glued body invalid: {}", issues[0])
        }
        Err(e) => e.to_string(),
    };
    crate::complex_shadow::append_log_line(&format!("GeneralFallback: {why}"));
    pieces
}

#[cfg(test)]
#[path = "general_tests.rs"]
mod tests;
