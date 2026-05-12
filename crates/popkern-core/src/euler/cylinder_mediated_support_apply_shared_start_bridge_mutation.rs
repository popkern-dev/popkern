//! Shared-start bridge mutation dispatch.

use crate::error::CoreError;
use crate::topology::BRep;

use super::CylinderMediatedSupportContinuationSegmentAuthoringPlan;
use super::apply_cylinder_resources::{BridgeResourcePairing, CylinderBridgeResourceSegment};
use super::apply_shared_start_bridge::{preflight_shared_start_bridge, shared_start_bridge_cycles};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SharedStartBridgeMutationStatus {
    NotApplicable,
    Committed,
}

pub(super) fn apply_shared_start_bridge_mutation(
    brep: &mut BRep,
    pairings: Option<&[BridgeResourcePairing]>,
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> Result<SharedStartBridgeMutationStatus, CoreError> {
    let Some(pairings) = pairings else {
        return Ok(SharedStartBridgeMutationStatus::NotApplicable);
    };
    if pairings.is_empty() {
        return Ok(SharedStartBridgeMutationStatus::NotApplicable);
    }

    let preflight = preflight_shared_start_bridge(brep, pairings, authoring, cylinder_segments);
    if !preflight.requires_candidate_coverage() || !preflight.candidate_gate_complete() {
        return Err(CoreError::PreconditionFailed(format!(
            "shared-start bridge mutation gate rejected {}",
            preflight.diagnostic_label()
        )));
    }

    let cycles = shared_start_bridge_cycles(brep, pairings, authoring, cylinder_segments)?;
    if let Ok(plan) = brep.preflight_face_family_bridge(cycles.clone()) {
        brep.face_family_bridge(&plan)?;
        return Ok(SharedStartBridgeMutationStatus::Committed);
    }
    let plan = brep.preflight_boundary_side_separator_bridge(
        crate::euler::BoundarySideSeparatorBridgeRequest { cycles },
    )?;
    brep.boundary_side_separator_bridge(&plan)?;
    Ok(SharedStartBridgeMutationStatus::Committed)
}
