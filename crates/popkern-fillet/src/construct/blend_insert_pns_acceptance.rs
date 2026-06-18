//! Transactional PnS dispatch acceptance for vertex-blend insertion.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

use super::blend_insert_diagnostics::pns_dispatch_declined_diagnostic;

pub(super) struct PnsDispatchAcceptanceInput<'a> {
    pub brep: &'a mut BRep,
    pub blend_face: FaceHandle,
    pub boundary: &'a BlendBoundary,
    pub fillet_surfaces: &'a [FilletSurfaceResult],
    pub vertex: VertexHandle,
    pub classification: VertexBlendClassification,
    pub n_sides: usize,
    pub diagnostics: &'a mut Vec<ApexConstructionDiagnostic>,
}

pub(super) fn try_pns_dispatch_with_acceptance(
    input: PnsDispatchAcceptanceInput<'_>,
) -> Option<Vec<CornerPatchInfo>> {
    let PnsDispatchAcceptanceInput {
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        n_sides,
        diagnostics,
    } = input;
    let pns_snapshot = brep.clone();
    let Some(patches) = crate::blend::pns_promote::try_pns_dispatch(
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        boundary.bridge_fallbacks.clone(),
    ) else {
        return None;
    };

    let cap_faces = patches.iter().map(|patch| patch.face).collect::<Vec<_>>();
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
        return Some(patches);
    }

    if let Some(diagnostic) =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep, &cap_faces, 5.0e-3, 5.0e-2,
        )
    {
        diagnostics.push(diagnostic);
    }
    let detail_code = topology_error
        .map(|error| format!("topology_validation_failed={error}"))
        .or(continuity_error);
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "cap_support_acceptance",
        detail_code,
    ));
    *brep = pns_snapshot;
    None
}
