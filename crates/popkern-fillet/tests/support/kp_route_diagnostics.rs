pub(crate) fn assert_route_declined(result: &popkern_fillet::BlendResult, stage: &str, code: &str) {
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| {
                matches!(
                    diagnostic,
                    popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                        if cap.stage == stage && cap.code == code
                )
            }),
        "expected route decline stage={stage} code={code}; diagnostics={:?}",
        exact_boundary_diagnostic_summary(result)
    );
}

pub(crate) fn assert_route_detail_contains(
    result: &popkern_fillet::BlendResult,
    stage: &str,
    code: &str,
    expected: &[&str],
) {
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| {
                matches!(
                    diagnostic,
                    popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                        if cap.stage == stage
                            && cap.code == code
                            && cap.detail_code.as_deref().is_some_and(|detail| {
                                expected.iter().all(|expected| detail.contains(expected))
                            })
                )
            }),
        "expected route detail stage={stage} code={code} detail containing {expected:?}; diagnostics={:?}",
        exact_boundary_diagnostic_summary(result)
    );
}

pub(crate) fn assert_support_cell_assembled_quality_recorded(result: &popkern_fillet::BlendResult) {
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                        if cap.stage == "kp20_polar_macro.dispatch.support_cell_assembled_quality"
                            && cap.code == "kp20_polar_macro.dispatch.fact"
                            && cap.detail_code.as_deref().is_some_and(|detail| {
                                detail.contains("selected=promoted")
                                    && detail.contains("selected_surface_count=33")
                                    && detail.contains("selected_cells_per_sector=4")
                                    && detail.contains("validation_status=valid")
                            })
            )),
        "expected KP20 support-cell assembled-quality fact; diagnostics={:?}",
        exact_boundary_diagnostic_summary(result)
    );
}

pub(crate) fn assert_selected_fallback_cap_support_continuity_recorded(
    result: &popkern_fillet::BlendResult,
) {
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::CapSupportContinuity(cap)
                    if cap.measurement_scope == "selected_fallback_live_topology"
                        && cap.edges_measured > 0
                        && cap.cap_uv.is_some()
                        && cap.support_uv.is_some()
                        && cap.cap_curve_gap.is_some()
                        && cap.support_curve_gap.is_some()
                        && cap.continuity.decline_code.is_some()
            )),
        "expected selected fallback cap/support continuity diagnostic; diagnostics={:?}",
        apex_diagnostic_summary(result)
    );
}

fn exact_boundary_diagnostic_summary(result: &popkern_fillet::BlendResult) -> Vec<String> {
    apex_diagnostic_summary(result)
        .into_iter()
        .filter(|summary| summary.starts_with("exact:"))
        .collect()
}

fn apex_diagnostic_summary(result: &popkern_fillet::BlendResult) -> Vec<String> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .map(|diagnostic| match diagnostic {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap) => format!(
                "exact:stage={};code={};detail={}",
                cap.stage,
                cap.code,
                cap.detail_code.as_deref().unwrap_or("<none>")
            ),
            popkern_fillet::ApexConstructionDiagnostic::CapSupportContinuity(cap) => format!(
                "cap_support:scope={};edges={};code={:?};cap_uv={:?};support_uv={:?};cap_curve_gap={:?};support_curve_gap={:?}",
                cap.measurement_scope,
                cap.edges_measured,
                cap.continuity.decline_code,
                cap.cap_uv,
                cap.support_uv,
                cap.cap_curve_gap,
                cap.support_curve_gap
            ),
            other => format!("{other:?}"),
        })
        .collect()
}
