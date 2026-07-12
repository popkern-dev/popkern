use crate::direct_apex_validate_helpers::{
    build_under_direct_apex, direct_apex_fallback_reasons, native_cap_faces,
};
use crate::helpers::measure_cap_outer_boundary;
use popkern_fillet::{
    ApexConstructionDiagnostic, ExactBoundaryCandidateOutcome,
    ExactBoundaryReportTransactionStatus, ExactBoundaryRouteDecisionOutcome,
    build_blend_repair_action_artifact,
};

/// Native cube-mixed direct-apex status. This fixture completes through the
/// exact-boundary cap route and resolves cap-outer residual diagnostics after
/// the final support-authority trim restore.
#[test]
fn direct_apex_5_cap_outer_cube_mixed_records_accepted_exact_boundary_route() {
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
    assert_eq!(route.exact_boundary_report.decision_chain.len(), 1);
    assert_eq!(
        route.exact_boundary_report.decision_chain[0].outcome,
        ExactBoundaryRouteDecisionOutcome::Accepted
    );
    let artifact = build_blend_repair_action_artifact(&result.repair_report);
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| {
                diagnostic.stage == "exact_boundary_cap.support_scope"
                    && diagnostic.code == "exact_boundary_cap.support_scope.residual_exceeded"
            }),
        "accepted cube_mixed exact-boundary route should not carry stale cap-outer residual diagnostics"
    );
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| diagnostic.stage != "exact_boundary_cap.cap_outer"),
        "accepted cube_mixed exact-boundary route should resolve cap-outer diagnostics"
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

fn assert_pyramid_exact_boundary_passes_commercial_cap_quality(
    fixture: &str,
    radii: &[f64],
    expect_multipatch_recovery: bool,
) {
    let result = build_under_direct_apex(true, radii);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "{fixture} exact-boundary route must stay native and not fall through to legacy; got:\n{}",
        reasons.join("\n"),
    );
    let route = result
        .repair_report
        .apex_routes
        .first()
        .expect("pyramid fixture should record one direct-apex route");
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted,
        "{fixture} exact-boundary route should satisfy live cap/support seams"
    );
    assert!(
        !route.exact_boundary_report.committed_faces.is_empty(),
        "accepted {fixture} exact-boundary route must report committed faces"
    );
    assert_eq!(
        route.exact_boundary_report.provenance.transaction_status,
        ExactBoundaryReportTransactionStatus::Committed,
        "accepted {fixture} exact-boundary transaction should commit"
    );
    let one_patch_step = route
        .exact_boundary_report
        .decision_chain
        .iter()
        .find(|step| step.route_label == "exact_boundary.one_patch")
        .expect("pyramid fixture route must preserve the one-patch decision");
    if expect_multipatch_recovery {
        assert_eq!(
            one_patch_step.outcome,
            ExactBoundaryRouteDecisionOutcome::Declined,
            "{fixture} should preserve the one-patch constrained-strip decline before multipatch recovery"
        );
        assert_eq!(
            one_patch_step.diagnostic_code.as_deref(),
            Some("exact_boundary_cap.validation.constrained_strip_band_unsatisfied"),
            "{fixture} should route to multipatch only for the measured constrained-strip gate"
        );
        let multipatch_step = route
            .exact_boundary_report
            .decision_chain
            .iter()
            .find(|step| step.route_label == "exact_boundary.multipatch")
            .expect("pyramid fixture constrained-strip decline must record the multipatch recovery decision");
        assert_eq!(
            multipatch_step.outcome,
            ExactBoundaryRouteDecisionOutcome::Accepted,
            "{fixture} multipatch route must accept only after live cap/support G1 passes"
        );
    } else {
        assert_eq!(
            one_patch_step.outcome,
            ExactBoundaryRouteDecisionOutcome::Accepted,
            "{fixture} should either accept directly or be promoted into the explicit multipatch recovery matrix"
        );
    }
    assert!(
        route
            .exact_boundary_report
            .continuity
            .decline_code
            .is_none(),
        "accepted {fixture} exact-boundary route must not carry stale cap/support decline evidence"
    );
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| {
                diagnostic.code != "cap_outer_g0_exceeds_gate"
                    && diagnostic.code != "cap_outer_g1_exceeds_gate"
                    && (diagnostic.code
                        != "exact_boundary_cap.multipatch.support_strip_band_unsatisfied"
                        || diagnostic
                            .residual
                            .map(|residual| residual.is_finite())
                            .unwrap_or(false))
            }),
        "{fixture} exact-boundary route must not regress to G0 or infinite support-strip diagnostics"
    );
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(
        metrics.worst_g0_gap <= 5.0e-3,
        "{fixture} exact-boundary cap_outer_g0 {:.3e} exceeds 5.0e-3",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity <= 5.0e-2,
        "{fixture} exact-boundary cap_outer_g1 {:.3e} exceeds 5.0e-2",
        metrics.worst_g1_coplanarity
    );
    let artifact = build_blend_repair_action_artifact(&result.repair_report);
    assert!(
        artifact
            .actions
            .iter()
            .all(|action| action.stage != "exact_boundary_cap" || action.outcome != "declined"),
        "accepted {fixture} exact-boundary route should not export exact-boundary decline actions"
    );
}

#[test]
fn pyramid3_exact_boundary_acceptance_matrix() {
    for (fixture, radii, expect_multipatch_recovery) in [
        ("pyramid3_uniform", [0.10, 0.10, 0.10], true),
        ("pyramid3_mixed", [0.08, 0.10, 0.12], false),
        ("pyramid3_extreme", [0.02, 0.10, 0.20], false),
    ] {
        assert_pyramid_exact_boundary_passes_commercial_cap_quality(
            fixture,
            &radii,
            expect_multipatch_recovery,
        );
    }
}

#[test]
fn pyramid4_exact_boundary_acceptance_matrix() {
    for (fixture, radii, expect_multipatch_recovery) in [
        ("pyramid4_uniform", [0.05, 0.05, 0.05, 0.05], false),
        ("pyramid4_mixed", [0.04, 0.06, 0.05, 0.07], false),
    ] {
        assert_pyramid_exact_boundary_passes_commercial_cap_quality(
            fixture,
            &radii,
            expect_multipatch_recovery,
        );
    }
}
