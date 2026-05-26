use crate::direct_apex_validate_helpers::{build_under_direct_apex, direct_apex_fallback_reasons};
use popkern_fillet::{
    ApexConstructionDiagnostic, ExactBoundaryCandidateOutcome,
    ExactBoundaryReportTransactionStatus, build_blend_repair_action_artifact,
};

/// Native cube-mixed direct-apex status. This fixture now completes through
/// the exact-boundary cap route, so the regression guard is the accepted
/// transaction and its structured route evidence rather than the old Gregory
/// cap-outer G1 blocker.
#[test]
fn direct_apex_5_cap_outer_cube_mixed_records_remaining_g1_blocker() {
    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "cube_mixed cap-outer gate must complete through native direct-apex; got:\n{}",
        reasons.join("\n"),
    );

    let routes = &result.repair_report.apex_routes;
    assert_eq!(routes.len(), 1, "cube_mixed should have one apex route");
    let route = &routes[0];
    assert_eq!(route.exact_boundary_sector_count, 3);
    assert_eq!(route.exact_boundary_pcurve_pair_count, 6);
    assert_eq!(route.exact_boundary_report.committed_faces.len(), 3);
    assert!(
        route.exact_boundary_report.candidate_status == ExactBoundaryCandidateOutcome::Accepted,
        "cube_mixed should complete through accepted exact-boundary construction"
    );
    assert!(
        route.exact_boundary_report.provenance.transaction_status
            == ExactBoundaryReportTransactionStatus::Committed,
        "cube_mixed exact-boundary topology transaction should commit"
    );
    let artifact = build_blend_repair_action_artifact(&result.repair_report);
    assert!(
        route.exact_boundary_report.diagnostics.is_empty(),
        "accepted cube_mixed exact-boundary route should not carry decline diagnostics"
    );
    assert!(
        artifact
            .actions
            .iter()
            .all(|action| action.stage != "exact_boundary_cap" || action.outcome != "declined"),
        "accepted cube_mixed exact-boundary route should not export decline actions"
    );
}

#[test]
fn cube_extreme_reaches_legacy_gregory_rewritten_boundary_gap_before_fallback() {
    let result = build_under_direct_apex(false, &[0.02, 0.10, 0.20]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons
            .iter()
            .any(|reason| reason.contains("Gregory FilletEnd[0] edge")
                && reason.contains("has stale cap-side end")),
        "cube_extreme must reach the legacy Gregory rewritten-boundary integration gap after exact-boundary sector ownership; got:\n{}",
        reasons.join("\n"),
    );

    let continuity = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => Some(diagnostic),
            _ => None,
        })
        .expect("cube_extreme fallback must promote cap/support continuity evidence");
    assert_eq!(
        continuity.continuity.decline_code.as_deref(),
        Some("cap_outer_g0_exceeds_gate")
    );
    assert!(continuity.continuity.worst_g0_gap.unwrap_or(0.0) > 5.0e-3);
}
