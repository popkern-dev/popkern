use popkern_core::geometry::gregory::GregoryPatch;

use crate::types::{ApexConstructionDiagnostic, CornerSubPatchProvenance};

pub(crate) struct Kp20PolarDispatchReport {
    pub(crate) patches: Option<Vec<GregoryPatch>>,
    pub(crate) support_cell_provenance: Option<Vec<CornerSubPatchProvenance>>,
    pub(crate) support_cell_floor_patches: Option<Vec<GregoryPatch>>,
    pub(crate) expanded_support_topology_accepted: bool,
    pub(crate) support_cell_topology_accepted: bool,
    pub(crate) diagnostics: Vec<ApexConstructionDiagnostic>,
}

impl Kp20PolarDispatchReport {
    pub(super) fn declined(diagnostics: Vec<ApexConstructionDiagnostic>) -> Self {
        Self {
            patches: None,
            support_cell_provenance: None,
            support_cell_floor_patches: None,
            expanded_support_topology_accepted: false,
            support_cell_topology_accepted: false,
            diagnostics,
        }
    }

    pub(super) fn accepted(
        patches: Vec<GregoryPatch>,
        diagnostics: Vec<ApexConstructionDiagnostic>,
    ) -> Self {
        Self {
            patches: Some(patches),
            support_cell_provenance: None,
            support_cell_floor_patches: None,
            expanded_support_topology_accepted: false,
            support_cell_topology_accepted: false,
            diagnostics,
        }
    }

    pub(super) fn accepted_expanded_support(
        patches: Vec<GregoryPatch>,
        diagnostics: Vec<ApexConstructionDiagnostic>,
    ) -> Self {
        Self {
            patches: Some(patches),
            support_cell_provenance: None,
            support_cell_floor_patches: None,
            expanded_support_topology_accepted: true,
            support_cell_topology_accepted: false,
            diagnostics,
        }
    }

    pub(super) fn accepted_support_cell(
        patches: Vec<GregoryPatch>,
        support_cell_provenance: Vec<CornerSubPatchProvenance>,
        support_cell_floor_patches: Option<Vec<GregoryPatch>>,
        diagnostics: Vec<ApexConstructionDiagnostic>,
    ) -> Self {
        Self {
            patches: Some(patches),
            support_cell_provenance: Some(support_cell_provenance),
            support_cell_floor_patches,
            expanded_support_topology_accepted: false,
            support_cell_topology_accepted: true,
            diagnostics,
        }
    }
}
