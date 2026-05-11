//! Candidate coverage gate for cylinder-mediated continuation refusals.

use crate::error::CoreError;

use super::apply_cylinder_resource_blockers::CylinderResourceBlockerExposureResult;
use super::apply_cylinder_resources::BridgeResourceReconciliationStatus;
use super::apply_owner_spans::OwnerRepairExposureResult;
use super::apply_shared_start_bridge::SharedStartBridgePreflight;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CylinderMediatedCandidateCoverage {
    owner_exposure_candidates: usize,
    blocker_exposure_candidates: usize,
    resource_candidates: usize,
    shared_start_bridge_required: bool,
    shared_start_bridge_candidates: usize,
}

impl CylinderMediatedCandidateCoverage {
    pub(super) fn from_pipeline(
        owner_exposure: OwnerRepairExposureResult,
        blocker_exposure: CylinderResourceBlockerExposureResult,
        resource_status: BridgeResourceReconciliationStatus,
        shared_start_bridge: &SharedStartBridgePreflight,
    ) -> Self {
        Self {
            owner_exposure_candidates: owner_exposure.candidate_count(),
            blocker_exposure_candidates: blocker_exposure.candidate_count(),
            resource_candidates: resource_status.candidate_count(),
            shared_start_bridge_required: shared_start_bridge.requires_candidate_coverage(),
            shared_start_bridge_candidates: shared_start_bridge.candidate_count(),
        }
    }

    pub(super) fn require_complete(self) -> Result<Self, CoreError> {
        if self.owner_exposure_candidates == 0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated diagnostics missing owner exposure construction candidate"
                    .into(),
            ));
        }
        if self.blocker_exposure_candidates == 0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated diagnostics missing blocker exposure construction candidate"
                    .into(),
            ));
        }
        if self.resource_candidates == 0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated diagnostics missing resource construction candidate".into(),
            ));
        }
        if self.shared_start_bridge_required && self.shared_start_bridge_candidates == 0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated diagnostics missing shared-start bridge construction candidate"
                    .into(),
            ));
        }
        Ok(self)
    }

    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "cylinder_mediated_candidate_coverage=owner_exposure:{} blocker_exposure:{} resource:{} shared_start_required:{} shared_start_bridge:{} complete:{}",
            self.owner_exposure_candidates,
            self.blocker_exposure_candidates,
            self.resource_candidates,
            self.shared_start_bridge_required,
            self.shared_start_bridge_candidates,
            self.is_complete()
        )
    }

    fn is_complete(self) -> bool {
        self.owner_exposure_candidates > 0
            && self.blocker_exposure_candidates > 0
            && self.resource_candidates > 0
            && (!self.shared_start_bridge_required || self.shared_start_bridge_candidates > 0)
    }
}
