//! PK-201: reporting must not change what a route does.
//!
//! `accepted` in the facade requires the FAILURE-EVIDENCE vector to be
//! empty, so anything pushed there vetoes the route. That made
//! explaining yourself cost you the corner — measured: adding one
//! purely informational record flipped pyramid apexes from the Coons
//! fan to a fallback, with no other change.
//!
//! These gates pin the separation. They are structural rather than
//! end-to-end because the property is structural: informational
//! records join the report only after the accept decision is taken.

use super::pns_candidate::{PnsCandidateError, pns_substitution_note};
use crate::types::ApexConstructionDiagnostic;

#[test]
fn a_substitution_note_is_produced_for_both_outcomes() {
    let ok = pns_substitution_note(4, &Ok(4));
    let err = pns_substitution_note(
        4,
        &Err(PnsCandidateError::SectorNotJoined {
            sector: 2,
            gap: 1.5e-2,
        }),
    );
    for (note, expect) in [(ok, "substituted"), (err, "declined")] {
        let ApexConstructionDiagnostic::ExactBoundaryCap(d) = note else {
            panic!("substitution note must be an exact-boundary diagnostic");
        };
        assert!(
            d.stage.contains(expect),
            "the note must say which outcome it records, got stage {}",
            d.stage
        );
        assert!(
            d.detail_code.is_some(),
            "a note with no detail explains nothing"
        );
    }
}

/// The decline note carries the measured quantity, not just a label.
///
/// A note that says only "declined" reproduces the silence this
/// channel exists to remove — the reason has to be actionable.
#[test]
fn a_decline_note_carries_its_measurement() {
    let note = pns_substitution_note(
        3,
        &Err(PnsCandidateError::SectorNotJoined {
            sector: 1,
            gap: 2.361e-2,
        }),
    );
    let ApexConstructionDiagnostic::ExactBoundaryCap(d) = note else {
        panic!("expected an exact-boundary diagnostic");
    };
    let detail = d.detail_code.expect("decline notes carry detail");
    assert!(
        detail.contains("2.361e-2"),
        "the measured gap must survive into the record, got: {detail}"
    );
}
