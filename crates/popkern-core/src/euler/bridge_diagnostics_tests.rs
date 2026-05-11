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
        "bridge_diagnostics=operations=0 resources=0 conflicts=0 readiness=0 construction_candidates=0 postconditions=0 \
         conflicts_by_kind=duplicate_half_edge:0 duplicate_boundary_edge:0 shared_start:0 \
         shared_end:0 crossed_endpoint:0 aligned:0 nested:0 interior:0 repeated_edge:0 \
         already_boundary:0 paired_exposure:0 unsupported_owner:0 stale_handle:0 unknown:0 \
         readiness_by_status=ready:0 ready_after_normalization:0 missing_operator:0 \
         invalid_input:0 unknown:0 \
         candidate_families=disjoint_cross_loop:0 disjoint_same_loop:0 \
         paired_same_loop_twin_exposure:0 boundary_side_separator:0 unknown:0 \
         candidate_feasibility=executable_existing:0 executable_after_normalization:0 \
         requires_separator_subcomplex:0 requires_custom_aggregate_authoring:0 \
         requires_cycle_regrouping:0 invalid_non_manifold:0 unsupported_unknown_family:0 unknown:0 \
         candidate_requirements=live_revalidation:0 postcondition_verification:0 unknown:0 \
         postcondition_preview_candidates=0 \
         geometry_feasibility=facts:0 ready:0 unsupported_family:0 repair_required:0 \
         not_applicable:0 unknown:0 \
         validity_repair=facts:0 valid_as_input:0 valid_after_non_mutating_normalization:0 \
         requires_topology_repair:0 requires_geometry_repair:0 \
         rejected_design_intent_changing:0 not_applicable:0 unknown:0 \
         status_no_repair_required:0 status_non_mutating_normalization_available:0 \
         status_repair_required:0 status_rejected:0 status_unknown:0 \
         operation_dependencies=facts:0 executable:0 ambiguous:0 invalid:0 missing_primitive:0 \
         not_applicable:0 unknown:0 status_satisfied:0 status_not_required:0 \
         status_required_missing_primitive:0 status_blocked_invalid_input:0 \
         status_ambiguous_needs_resolver:0 status_unknown:0 \
         predicted_delta=dV=0 dE=0 dF=0 dL=0"
    );
    assert!(first.graph.candidates.is_empty());
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
