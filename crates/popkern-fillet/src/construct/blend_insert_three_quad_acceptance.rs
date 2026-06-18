//! Transactional acceptance for rational/KP three-quad dispatch results.

use popkern_core::error::ValidationError;
use popkern_core::topology::BRep;

use crate::construct::three_quad_apply::ThreeQuadApplyResult;
use crate::types::ApexConstructionDiagnostic;

use super::blend_insert_diagnostics::rational_gregory_dispatch_declined_diagnostic;
use super::blend_insert_support_cell_acceptance::{
    corner_patch_count_is_acceptable, support_cell_patch_count_is_acceptable,
};

pub(super) enum ThreeQuadAcceptance {
    Accepted(ThreeQuadApplyResult),
    Declined(Vec<ApexConstructionDiagnostic>),
}

pub(super) fn accept_three_quad_result(
    brep: &BRep,
    mut result: ThreeQuadApplyResult,
    n_sides: usize,
) -> ThreeQuadAcceptance {
    if !corner_patch_count_is_acceptable(&result.patches, n_sides) {
        let mut diagnostics = std::mem::take(&mut result.diagnostics);
        diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
            n_sides,
            "patch_count_mismatch",
            Some(format!(
                "expected {n_sides} corner patches, got {}",
                result.patches.len()
            )),
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

    let mut diagnostics = std::mem::take(&mut result.diagnostics);
    if let Some(diagnostic) =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep, &cap_faces, 5.0e-3, 5.0e-2,
        )
    {
        diagnostics.push(diagnostic);
    }
    let detail = topology_error
        .map(|error| format!("topology_validation_failed={error}"))
        .or(continuity_error)
        .unwrap_or_else(|| "cap/support acceptance failed".to_string());
    diagnostics.push(rational_gregory_dispatch_declined_diagnostic(
        n_sides,
        "cap_support_acceptance",
        Some(detail),
    ));
    ThreeQuadAcceptance::Declined(diagnostics)
}
