//! Transactional acceptance for rational/KP three-quad dispatch results.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::three_quad_apply::ThreeQuadApplyResult;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};

use super::blend_insert_diagnostics::rational_gregory_dispatch_declined_diagnostic;
use super::blend_insert_support_cell_acceptance::{
    corner_patch_count_is_acceptable, partial_support_cell_coverage_detail,
    support_cell_patch_count_is_acceptable,
};

pub(super) enum ThreeQuadAcceptance {
    Accepted(ThreeQuadApplyResult),
    Declined(Vec<ApexConstructionDiagnostic>),
}

pub(super) fn accept_three_quad_result(
    brep: &mut BRep,
    mut result: ThreeQuadApplyResult,
    n_sides: usize,
) -> ThreeQuadAcceptance {
    if !corner_patch_count_is_acceptable(&result.patches, n_sides) {
        let mut diagnostics = std::mem::take(&mut result.diagnostics);
        let (reason, detail) =
            if let Some(detail) = partial_support_cell_coverage_detail(&result.patches, n_sides) {
                ("partial_support_cell_patch_coverage", Some(detail))
            } else {
                (
                    "patch_count_mismatch",
                    Some(format!(
                        "expected {n_sides} corner patches, got {}",
                        result.patches.len()
                    )),
                )
            };
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides, reason, detail,
        ));
        return ThreeQuadAcceptance::Declined(diagnostics);
    }
    if n_sides < 5 || support_cell_patch_count_is_acceptable(&result.patches, n_sides) {
        return ThreeQuadAcceptance::Accepted(result);
    }

    let cap_faces = result
        .patches
        .iter()
        .map(|patch| patch.face)
        .collect::<Vec<_>>();
    crate::construct::corner_patch_support_completion::complete_external_cap_support_pcurves(
        brep,
        &result.patches,
        &mut result.diagnostics,
    );
    let topology_error = brep.validate_all().into_iter().find(|error| {
        matches!(
            error,
            ValidationError::BrokenCycle { .. }
                | ValidationError::ShellGap { .. }
                | ValidationError::OpenLoop { .. }
        )
    });
    let continuity_error = crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep, &cap_faces, 5.0e-3, 5.0e-2,
    );
    if topology_error.is_none() && continuity_error.is_none() {
        return ThreeQuadAcceptance::Accepted(result);
    }

    let continuity_diagnostic =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep, &cap_faces, 5.0e-3, 5.0e-2,
        );
    let decline_context = continuity_diagnostic
        .as_ref()
        .and_then(|diagnostic| cap_support_decline_context(diagnostic, &result.patches));
    let mut diagnostics = std::mem::take(&mut result.diagnostics);
    diagnostics.push(
        crate::construct::three_quad_kp_diagnostics::support_cell_candidate_nonselection(
            "non_support_three_quad_acceptance_declined",
            n_sides,
            result.patches.len(),
        ),
    );
    if let Some(diagnostic) = continuity_diagnostic {
        diagnostics.push(diagnostic);
    }
    let mut detail = topology_error
        .map(|error| format!("topology_validation_failed={error}"))
        .or(continuity_error)
        .unwrap_or_else(|| "cap/support acceptance failed".to_string());
    if let Some(context) = decline_context {
        detail.push(';');
        detail.push_str(&context);
    }
    diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
        n_sides,
        "cap_support_acceptance",
        Some(detail),
    ));
    ThreeQuadAcceptance::Declined(diagnostics)
}

pub(super) fn record_fallback_cap_support_continuity(
    brep: &BRep,
    cap_face: FaceHandle,
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    if let Some(mut diagnostic) =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep,
            &[cap_face],
            5.0e-3,
            5.0e-2,
        )
    {
        if let ApexConstructionDiagnostic::CapSupportContinuity(cap) = &mut diagnostic {
            cap.measurement_scope = "selected_fallback_live_topology".to_string();
        }
        diagnostics.push(diagnostic);
    }
}

fn cap_support_decline_context(
    diagnostic: &ApexConstructionDiagnostic,
    patches: &[CornerPatchInfo],
) -> Option<String> {
    let ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) = diagnostic else {
        return None;
    };
    let cap_provenance = patches
        .iter()
        .find(|patch| patch.face == diagnostic.cap_face)
        .and_then(|patch| patch.sub_patch_provenance);
    let support_provenance = patches
        .iter()
        .find(|patch| patch.face == diagnostic.support_face)
        .and_then(|patch| patch.sub_patch_provenance);
    Some(format!(
        "cap_support_context=edge={:?};cap_face={:?};support_face={:?};sample={:.3};cap_uv={};support_uv={};cap_curve_gap={};support_curve_gap={};cap_provenance={:?};support_provenance={:?}",
        diagnostic.edge,
        diagnostic.cap_face,
        diagnostic.support_face,
        diagnostic.sample_parameter.unwrap_or_default(),
        format_optional_uv(diagnostic.cap_uv),
        format_optional_uv(diagnostic.support_uv),
        format_optional_metric(diagnostic.cap_curve_gap),
        format_optional_metric(diagnostic.support_curve_gap),
        cap_provenance,
        support_provenance,
    ))
}

fn format_optional_uv(uv: Option<(f64, f64)>) -> String {
    uv.map(|(u, v)| format!("({u:.3},{v:.3})"))
        .unwrap_or_else(|| "none".to_string())
}

fn format_optional_metric(metric: Option<f64>) -> String {
    metric
        .map(|value| format!("{value:.3e}"))
        .unwrap_or_else(|| "none".to_string())
}
