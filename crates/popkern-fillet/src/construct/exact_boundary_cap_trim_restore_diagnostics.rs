//! Exact-boundary trim-restore diagnostic cleanup.

pub(super) fn prune_resolved_cap_outer_diagnostics(
    repair_report: &mut crate::types::BlendRepairReport,
    resolved: &[(String, String)],
) {
    if resolved.is_empty() {
        return;
    }
    for route in &mut repair_report.apex_routes {
        route
            .exact_boundary_report
            .diagnostics
            .retain(|diagnostic| !is_resolved_cap_outer_diagnostic(diagnostic, resolved));
    }
    repair_report.apex_diagnostics.retain(|diagnostic| {
        let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
            return true;
        };
        !resolved.iter().any(|(sector_name, outer_side)| {
            cap.stage == "exact_boundary_cap.cap_outer"
                && cap.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded"
                && &cap.sector_name == sector_name
                && cap.outer_side.as_ref() == Some(outer_side)
        })
    });
}

fn is_resolved_cap_outer_diagnostic(
    diagnostic: &crate::types::ExactBoundaryReportDiagnostic,
    resolved: &[(String, String)],
) -> bool {
    resolved.iter().any(|(sector_name, outer_side)| {
        diagnostic.stage == "exact_boundary_cap.cap_outer"
            && diagnostic.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded"
            && &diagnostic.sector_name == sector_name
            && diagnostic.outer_side.as_ref() == Some(outer_side)
    })
}
