//! A non-manifold result as manifold pieces (PK-964, PK-1115,
//! `NonManifoldPolicy::SplitToManifold` and `General`).
//!
//! The pieces arrive already apart — their own copies of a contact
//! vertex or edge, or shells that touch — from whichever route served
//! them, in one body: an intermediate no single body may hold as its
//! answer. Acceptance requires every piece to validate on its own; a
//! touch with area (two shells coinciding over a neighbourhood) is no
//! contact (those pieces should have merged there), so it declines.
//! The public entry points then serve the pieces as separate bodies or
//! glue them into one general body (`crate::split_serve`).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::shell_contact::{ContactExtent, shell_findings};
use super::zero_thickness::zero_thickness_findings;
use crate::error::{BoolError, NonManifoldContact};
use crate::options::NonManifoldPolicy;

/// Where the pieces of `brep` touch, and whether a touch has area.
fn findings(brep: &BRep, tolerance: &Tolerance) -> (Vec<NonManifoldContact>, bool) {
    let mut all: Vec<NonManifoldContact> = zero_thickness_findings(brep, tolerance)
        .into_iter()
        .map(|(contact, _)| contact)
        .collect();
    let mut area = false;
    for finding in shell_findings(brep, tolerance) {
        area |= finding.extent == ContactExtent::Area;
        all.push(finding.contact);
    }
    (all, area)
}

/// Where the pieces of `brep` touch (none for a body whose lumps are
/// apart).
pub(crate) fn touch_contacts(brep: &BRep, tolerance: &Tolerance) -> Vec<NonManifoldContact> {
    findings(brep, tolerance).0
}

/// Accept `body` as manifold pieces: no touch with area, and under the
/// split policy each piece valid on its own (under the general policy a
/// piece may touch itself — a solid meeting itself along an edge — which
/// the glue closes and the glued body is validated whole), or decline
/// with the contacts.
pub(crate) fn accept_split(
    body: BRep,
    tolerance: &Tolerance,
    policy: NonManifoldPolicy,
) -> Result<BRep, BoolError> {
    let (contacts, area) = findings(&body, tolerance);
    if area {
        return Err(BoolError::NonManifoldResult { contacts });
    }
    if policy == NonManifoldPolicy::General {
        return Ok(body);
    }
    let pieces = body
        .piece_bodies()
        .map_err(|_| BoolError::NonManifoldResult {
            contacts: contacts.clone(),
        })?;
    for piece in &pieces {
        let issues = piece.validate_all();
        if let Some(first) = issues.first() {
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "[boolean] split piece declined: {} validation issue(s), first: {first}",
                    issues.len()
                );
            }
            return Err(BoolError::NonManifoldResult { contacts });
        }
    }
    Ok(body)
}
