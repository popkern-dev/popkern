//! Shared-start bridge mutation dispatch.

use crate::error::CoreError;
use crate::topology::BRep;

use super::CylinderMediatedSupportContinuationSegmentAuthoringPlan;
use crate::topology::{FaceHandle, LoopHandle};

use super::apply_cylinder_resources::{BridgeResourcePairing, CylinderBridgeResourceSegment};
use super::apply_shared_start_bridge::{preflight_shared_start_bridge, shared_start_bridge_cycles};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SharedStartBridgeMutationStatus {
    NotApplicable,
    Committed(SharedStartBridgeMutationResult),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SharedStartBridgeMutationResult {
    pub(super) bridge_faces: [FaceHandle; 3],
    pub(super) bridge_loops: [LoopHandle; 3],
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
        let result = brep.face_family_bridge(&plan)?;
        return Ok(SharedStartBridgeMutationStatus::Committed(bridge_result(
            result
                .cycles
                .iter()
                .map(|cycle| (cycle.bridge_face, cycle.bridge_loop)),
            result.cycles.len(),
        )?));
    }
    let plan = brep.preflight_boundary_side_separator_bridge(
        crate::euler::BoundarySideSeparatorBridgeRequest { cycles },
    )?;
    let result = brep.boundary_side_separator_bridge(&plan)?;
    Ok(SharedStartBridgeMutationStatus::Committed(bridge_result(
        result
            .bridge_faces
            .iter()
            .copied()
            .zip(result.bridge_loops.iter().copied()),
        result.bridge_faces.len(),
    )?))
}

fn bridge_result(
    items: impl Iterator<Item = (FaceHandle, LoopHandle)>,
    len: usize,
) -> Result<SharedStartBridgeMutationResult, CoreError> {
    let pairs = items.collect::<Vec<_>>();
    let [
        (first_face, first_loop),
        (second_face, second_loop),
        (third_face, third_loop),
    ] = pairs.as_slice()
    else {
        return Err(CoreError::PreconditionFailed(format!(
            "shared-start bridge mutation expected three bridge cycles, got {len}"
        )));
    };
    Ok(SharedStartBridgeMutationResult {
        bridge_faces: [*first_face, *second_face, *third_face],
        bridge_loops: [*first_loop, *second_loop, *third_loop],
    })
}
