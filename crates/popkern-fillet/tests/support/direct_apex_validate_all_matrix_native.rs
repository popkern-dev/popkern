use crate::direct_apex_validate_helpers::{
    build_under_direct_apex, direct_apex_fallback_reasons, native_cap_faces,
};
use crate::helpers::measure_cap_outer_boundary;
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
fn cube_mixed_exact_boundary_cap_outer_passes_commercial_gate() {
    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);

    assert!(
        metrics.edges_measured >= 1,
        "cube_mixed should expose measurable exact-boundary cap edges"
    );
    assert!(
        metrics.worst_g0_gap <= 5.0e-3,
        "cube_mixed cap_outer_g0 {:.3e} exceeds 5.0e-3",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity <= 5.0e-2,
        "cube_mixed cap_outer_g1 {:.3e} exceeds 5.0e-2",
        metrics.worst_g1_coplanarity
    );
}

#[test]
fn cube_extreme_completes_native_gregory_without_cap_outer_fallback() {
    let result = build_under_direct_apex(false, &[0.02, 0.10, 0.20]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "cube_extreme must complete through native Gregory after endpoint-preserving boundary construction; got:\n{}",
        reasons.join("\n"),
    );

    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(
        metrics.worst_g0_gap <= 5.0e-3 && metrics.worst_g1_coplanarity <= 5.0e-2,
        "cube_extreme exact-boundary cap_outer must pass after support-face subdivision; got g0={:.3e}, g1={:.3e}",
        metrics.worst_g0_gap,
        metrics.worst_g1_coplanarity,
    );

    let side_support = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::SideSupportContinuation(diagnostic) => Some(diagnostic),
            _ => None,
        })
        .expect("cube_extreme should report applied side-support topology continuation");
    assert_eq!(
        side_support.code,
        "side_support_continuation.outside_outside_cylinder_mediated.applied"
    );
    assert!(side_support.applied);
    assert_eq!(side_support.bridge_face_count, 3);
    assert_eq!(side_support.bridge_loop_count, 3);
}

#[test]
fn cube_extreme_exact_boundary_materializes_support_face_subdivision() {
    let result = build_under_direct_apex(false, &[0.02, 0.10, 0.20]);
    let route = result
        .repair_report
        .apex_routes
        .first()
        .expect("cube_extreme should record one direct-apex route");
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted,
        "support-face subdivision must let exact-boundary commit instead of declining"
    );
    assert!(
        !route.exact_boundary_report.committed_faces.is_empty(),
        "accepted exact-boundary route must record committed cap/support faces"
    );
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| diagnostic.code
                != "exact_boundary_cap.commit.support_face_subdivision_required"),
        "support-face subdivision-required is a pre-primitive decline and must not remain"
    );
    let support_scope_diagnostics = route
        .exact_boundary_report
        .diagnostics
        .iter()
        .filter(|diagnostic| diagnostic.stage == "exact_boundary_cap.support_scope")
        .collect::<Vec<_>>();
    assert!(
        support_scope_diagnostics
            .iter()
            .all(|diagnostic| diagnostic.support_half_edge.is_some()
                && diagnostic.support_face.is_some()),
        "support-scope diagnostics must preserve support ownership handles"
    );
}

#[test]
fn pyramid3_uniform_exact_boundary_declines_g1_instead_of_accepting_bad_cap() {
    let result = build_under_direct_apex(true, &[0.10, 0.10, 0.10]);
    let route = result
        .repair_report
        .apex_routes
        .first()
        .expect("pyramid3_uniform should record one direct-apex route");
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::Declined,
        "pyramid3 exact-boundary must not accept a support-normal-incompatible cap"
    );
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .any(|diagnostic| {
                diagnostic.stage == "exact_boundary_cap.validation"
                    && diagnostic.code
                        == "exact_boundary_cap.validation.support_normal_continuity_exceeded"
                    && diagnostic.edge.is_some()
                    && diagnostic.support_half_edge.is_some()
                    && diagnostic.support_face.is_some()
                    && diagnostic
                        .residual
                        .is_some_and(|residual| residual > 5.0e-2)
            }),
        "pyramid3 exact-boundary decline must preserve structured support-normal evidence"
    );
}
