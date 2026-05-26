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
fn cube_extreme_reaches_native_gregory_cap_outer_g0_before_fallback() {
    let result = build_under_direct_apex(false, &[0.02, 0.10, 0.20]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.iter().any(|reason| reason
            .contains("native Gregory cap failed geometric acceptance: cap_outer_g0")),
        "cube_extreme must reach the native Gregory cap/support geometric gate after exact-boundary sector ownership; got:\n{}",
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

    let sphere_preflight = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic)
                if diagnostic.stage == "native_gregory_sphere_cap_preflight" =>
            {
                Some(diagnostic)
            }
            _ => None,
        })
        .expect("cube_extreme must classify the declined sphere-cap boundary preflight");
    assert_eq!(
        sphere_preflight.code,
        "sphere_plane_boundary_tangent_or_empty"
    );
    assert!(
        sphere_preflight.support_face.is_some(),
        "sphere-plane preflight diagnostic should identify the support face"
    );

    let boundary_owned_pcurve = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic)
                if diagnostic.stage == "native_gregory_boundary_owned_pcurve" =>
            {
                Some(diagnostic)
            }
            _ => None,
        })
        .expect("cube_extreme must classify the declined boundary-owned seam pcurve");
    assert!(
        matches!(
            boundary_owned_pcurve.code.as_str(),
            "boundary_owned_pcurve_endpoint_mismatch"
                | "boundary_owned_pcurve_not_on_canonical_boundary"
        ),
        "unexpected boundary-owned pcurve decline code: {}",
        boundary_owned_pcurve.code
    );
    assert!(
        boundary_owned_pcurve.edge.is_some(),
        "boundary-owned pcurve decline should identify the responsible edge"
    );
}
