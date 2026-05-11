use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationPurpose,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn empty_bridge_diagnostic_report_is_deterministic() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let request = BridgeDiagnosticRequest { operations: vec![] };

    let first = brep
        .bridge_diagnostic_report(request.clone())
        .expect("empty report builds");
    let second = brep
        .bridge_diagnostic_report(request)
        .expect("empty report builds");

    assert_eq!(first, second);
    assert_eq!(
        first.summary.diagnostic_label(),
        "bridge_diagnostics=operations=0 resources=0 conflicts=0 readiness=0 postconditions=0"
    );
}

#[test]
fn bridge_diagnostic_fingerprint_tracks_scoped_request_membership() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let loop_ = brep.loops.keys().next().expect("cube has a loop");
    let half_edge = brep.loops[loop_].half_edge;
    let next = brep.half_edges[half_edge].next;

    let first = brep
        .bridge_diagnostic_report(BridgeDiagnosticRequest {
            operations: vec![BridgeDiagnosticOperationRequest {
                purpose: DiagnosticOperationPurpose::BridgeCycle,
                first_loop: loop_,
                first_interval: vec![half_edge, next],
                second_loop: loop_,
                second_interval: vec![next, half_edge],
            }],
        })
        .expect("report builds");
    let second = brep
        .bridge_diagnostic_report(BridgeDiagnosticRequest {
            operations: vec![BridgeDiagnosticOperationRequest {
                purpose: DiagnosticOperationPurpose::BridgeCycle,
                first_loop: loop_,
                first_interval: vec![next, half_edge],
                second_loop: loop_,
                second_interval: vec![half_edge, next],
            }],
        })
        .expect("report builds");

    assert_ne!(
        first.graph.version.scoped_membership_checksum,
        second.graph.version.scoped_membership_checksum
    );
}
