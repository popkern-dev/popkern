use popkern_core::geometry::gregory::GregoryPatch;

use crate::types::{ApexConstructionDiagnostic, CornerSubPatchProvenance};

pub(crate) struct Kp20PolarDispatchReport {
    pub(crate) patches: Option<Vec<GregoryPatch>>,
    /// `Some(relative g2 residual)` when the accepted paper candidate
    /// passed the COMMERCIAL tier rather than paper-exactness — the
    /// promotion ladder threads it into the corner's
    /// `relaxed_continuity_coplanarity` record so the containment is
    /// user-visible (CornerContinuityRelaxed).
    pub(crate) commercial_g2_relative: Option<f64>,
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
            commercial_g2_relative: None,
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
            commercial_g2_relative: None,
            support_cell_provenance: None,
            support_cell_floor_patches: None,
            expanded_support_topology_accepted: false,
            support_cell_topology_accepted: false,
            diagnostics,
        }
    }

    /// Paper candidate accepted at the COMMERCIAL tier (fold-free,
    /// corner-compatible, commercial derivative bounds) rather than
    /// paper-exactness; carries the measured relative g2 residual.
    pub(super) fn accepted_commercial(
        patches: Vec<GregoryPatch>,
        commercial_g2_relative: f64,
        diagnostics: Vec<ApexConstructionDiagnostic>,
    ) -> Self {
        Self {
            patches: Some(patches),
            commercial_g2_relative: Some(commercial_g2_relative),
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
            commercial_g2_relative: None,
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
            commercial_g2_relative: None,
            support_cell_provenance: Some(support_cell_provenance),
            support_cell_floor_patches,
            expanded_support_topology_accepted: false,
            support_cell_topology_accepted: true,
            diagnostics,
        }
    }
}
