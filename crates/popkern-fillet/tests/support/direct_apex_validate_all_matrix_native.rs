use crate::direct_apex_validate_helpers::{build_under_direct_apex, direct_apex_fallback_reasons};
use crate::helpers::measure_cap_outer_boundary;
use popkern_fillet::{ApexConstructionDiagnostic, build_blend_repair_action_artifact};

/// Native cube-mixed direct-apex Gregory cap-outer status. This fixture must
/// complete natively and preserve cap/support G0 without relying on
/// the legacy fallback route. Cap/support G1 remains a recorded
/// construction-quality blocker; do not hide it with local support-frame
/// rewrites that break ring continuity or route ownership.
#[test]
fn direct_apex_5_cap_outer_cube_mixed_records_remaining_g1_blocker() {
    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "cube_mixed cap-outer gate must complete through native direct-apex; got:\n{}",
        reasons.join("\n"),
    );
    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<Vec<_>>();
    let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(m.edges_measured >= 1, "cube_mixed measured no cap seams");
    assert!(
        m.worst_g0_gap <= 5.0e-3,
        "cube_mixed cap_outer_g0 {:.3e} exceeds 5.0e-3",
        m.worst_g0_gap,
    );
    assert!(
        m.worst_g1_coplanarity > 5.0e-2,
        "cube_mixed should continue exposing the cap_outer_g1 blocker until construction-side support continuation is repaired; measured {:.3e}",
        m.worst_g1_coplanarity,
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                ApexConstructionDiagnostic::CapSupportContinuity(d)
                    if d.worst_g1_coplanarity > d.gate && d.edges_measured >= 1
            )),
        "cube_mixed cap_outer_g1 blocker must be promoted to a structured apex diagnostic"
    );
    let artifact = build_blend_repair_action_artifact(&result.repair_report);
    assert!(
        artifact.actions.iter().any(|action| {
            action.stage == "apex_cap_support_continuity"
                && action.outcome == "declined"
                && action.reason.as_deref() == Some("cap_outer_g1_exceeds_gate")
        }),
        "cube_mixed cap_outer_g1 blocker must be exported in the repair-action artifact"
    );
}
