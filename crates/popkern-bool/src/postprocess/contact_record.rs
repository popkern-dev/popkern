//! Serving a non-manifold result as manifold pieces with recorded
//! contacts (PK-964, `NonManifoldPolicy::SplitToManifold`).
//!
//! The pieces arrive already split — coincident copies of a contact
//! vertex or edge, or shells that touch — from whichever route served
//! them. Acceptance records every contact on the body (with the
//! entities that meet there) and then requires `validate_all` to pass:
//! recorded contacts are the only touches it excuses, so a split that
//! leaves any other defect still declines. A touch with area (two
//! shells coinciding over a neighbourhood) is not a contact to record:
//! those bodies should have merged there, so it declines.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, ContactEntity, TopologyContact, TopologyContactKind};

use super::shell_contact::{ContactExtent, shell_findings};
use super::zero_thickness::zero_thickness_findings;
use crate::error::{BoolError, ContactKind, NonManifoldContact};

/// Every contact of `brep` as a record, or `Err` with the contacts
/// when one has area.
fn contact_records(
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<Vec<TopologyContact>, Vec<NonManifoldContact>> {
    let mut records = Vec::new();
    let mut all = Vec::new();
    let mut area = false;
    let kind_of = |k: ContactKind| match k {
        ContactKind::Edge => TopologyContactKind::Edge,
        ContactKind::Vertex => TopologyContactKind::Vertex,
    };
    for (contact, entities) in zero_thickness_findings(brep, tolerance) {
        records.push(TopologyContact {
            kind: kind_of(contact.kind),
            location: contact.location,
            entities,
        });
        all.push(contact);
    }
    for finding in shell_findings(brep, tolerance) {
        area |= finding.extent == ContactExtent::Area;
        records.push(TopologyContact {
            kind: kind_of(finding.contact.kind),
            location: finding.contact.location,
            entities: vec![
                ContactEntity::Face(finding.faces.0),
                ContactEntity::Face(finding.faces.1),
            ],
        });
        all.push(finding.contact);
    }
    if area { Err(all) } else { Ok(records) }
}

/// Serve `body` as manifold pieces meeting at recorded contacts, or
/// decline with the contacts when it cannot be served valid.
pub(crate) fn accept_split(mut body: BRep, tolerance: &Tolerance) -> Result<BRep, BoolError> {
    let records = contact_records(&body, tolerance)
        .map_err(|contacts| BoolError::NonManifoldResult { contacts })?;
    let contacts: Vec<NonManifoldContact> = records
        .iter()
        .map(|r| NonManifoldContact {
            kind: match r.kind {
                TopologyContactKind::Edge => ContactKind::Edge,
                TopologyContactKind::Vertex => ContactKind::Vertex,
            },
            location: r.location,
        })
        .collect();
    body.topology_contacts = records;
    let issues = body.validate_all();
    if issues.is_empty() {
        Ok(body)
    } else {
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "[boolean] split result declined: {} validation issue(s), first: {}",
                issues.len(),
                issues[0]
            );
        }
        Err(BoolError::NonManifoldResult { contacts })
    }
}
