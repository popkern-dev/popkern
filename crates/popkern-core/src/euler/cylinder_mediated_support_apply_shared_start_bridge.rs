//! Candidate bridge authoring for shared-start split cylinder resources.

use crate::error::CoreError;
use crate::euler::{FaceFamilyBridgeCycleRequest, FaceFamilyBridgePlan};
use crate::topology::{BRep, HalfEdgeHandle};

use super::CylinderMediatedSupportContinuationSegmentAuthoringPlan;
use super::apply_cylinder_resources::{BridgeResourcePairing, CylinderBridgeResourceSegment};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SharedStartBridgePreflight {
    status: SharedStartBridgePreflightStatus,
    pairings: usize,
    cycles: usize,
    expected_edges: i64,
    expected_faces: i64,
    expected_loops: i64,
    reason: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SharedStartBridgePreflightStatus {
    NotApplicable,
    Ready,
    Rejected,
}

impl SharedStartBridgePreflight {
    pub(super) fn not_applicable() -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::NotApplicable,
            pairings: 0,
            cycles: 0,
            expected_edges: 0,
            expected_faces: 0,
            expected_loops: 0,
            reason: "none",
        }
    }

    fn ready(pairings: usize, plan: &FaceFamilyBridgePlan) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Ready,
            pairings,
            cycles: plan.cycles.len(),
            expected_edges: plan.expected_delta.edges,
            expected_faces: plan.expected_delta.faces,
            expected_loops: plan.expected_delta.loops,
            reason: "none",
        }
    }

    fn rejected(pairings: usize, reason: &'static str) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Rejected,
            pairings,
            cycles: 0,
            expected_edges: 0,
            expected_faces: 0,
            expected_loops: 0,
            reason,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let status = match self.status {
            SharedStartBridgePreflightStatus::NotApplicable => "not-applicable",
            SharedStartBridgePreflightStatus::Ready => "ready",
            SharedStartBridgePreflightStatus::Rejected => "rejected",
        };
        format!(
            "shared_start_bridge_preflight=status={status} pairings={} cycles={} expected_delta=dV=0 dE={} dF={} dL={} reason={}",
            self.pairings,
            self.cycles,
            self.expected_edges,
            self.expected_faces,
            self.expected_loops,
            self.reason
        )
    }
}

pub(super) fn preflight_shared_start_bridge(
    brep: &BRep,
    pairings: &[BridgeResourcePairing],
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> SharedStartBridgePreflight {
    if pairings.is_empty() {
        return SharedStartBridgePreflight::not_applicable();
    }
    match shared_start_face_family_bridge_plan(brep, pairings, authoring, cylinder_segments) {
        Ok(plan) => SharedStartBridgePreflight::ready(pairings.len(), &plan),
        Err(err) => SharedStartBridgePreflight::rejected(pairings.len(), bridge_error_label(&err)),
    }
}

fn shared_start_face_family_bridge_plan(
    brep: &BRep,
    pairings: &[BridgeResourcePairing],
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> Result<FaceFamilyBridgePlan, CoreError> {
    let cycles = pairings
        .iter()
        .map(|pairing| bridge_cycle_request(brep, *pairing, authoring, cylinder_segments))
        .collect::<Result<Vec<_>, _>>()?;
    brep.preflight_face_family_bridge(cycles)
}

fn bridge_cycle_request(
    brep: &BRep,
    pairing: BridgeResourcePairing,
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> Result<FaceFamilyBridgeCycleRequest, CoreError> {
    let outside = authoring
        .get(pairing.outside_segment_index)
        .ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "shared-start bridge references missing outside segment {}",
                pairing.outside_segment_index
            ))
        })?;
    let cylinder = cylinder_segments
        .get(pairing.cylinder_segment_index)
        .ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "shared-start bridge references missing cylinder segment {}",
                pairing.cylinder_segment_index
            ))
        })?;
    let first_interval = interval_with_after_last(brep, &outside.half_edges, "outside")?;
    let second_interval = interval_with_after_last(brep, &cylinder.half_edges, "cylinder")?;
    let second_loop = brep
        .half_edges
        .get(second_interval[0])
        .map(|half_edge| half_edge.loop_)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "shared-start bridge cylinder interval has stale first half-edge".into(),
            )
        })?;
    Ok(FaceFamilyBridgeCycleRequest {
        first_loop: outside.outside_loop,
        first_interval,
        second_loop,
        second_interval,
    })
}

fn interval_with_after_last(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
    label: &'static str,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let first = half_edges.first().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!("shared-start bridge {label} interval is empty"))
    })?;
    let last = half_edges.last().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!("shared-start bridge {label} interval is empty"))
    })?;
    let after_last = brep
        .half_edges
        .get(last)
        .map(|half_edge| half_edge.next)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "shared-start bridge {label} interval has stale last half-edge"
            ))
        })?;
    let mut interval = Vec::with_capacity(half_edges.len() + 1);
    interval.push(first);
    for pair in half_edges.windows(2) {
        if brep.half_edges[pair[0]].next != pair[1] {
            return Err(CoreError::PreconditionFailed(format!(
                "shared-start bridge {label} interval is not directed"
            )));
        }
        interval.push(pair[1]);
    }
    interval.push(after_last);
    Ok(interval)
}

fn bridge_error_label(err: &CoreError) -> &'static str {
    let detail = err.to_string();
    if detail.contains("requires distinct source loops") {
        "face-family-bridge-cycle-requires-distinct-source-loops"
    } else if detail.contains("interval is not directed") {
        "face-family-bridge-cycle-interval-not-directed"
    } else if detail.contains("references missing outside segment") {
        "shared-start-bridge-references-missing-outside-segment"
    } else if detail.contains("references missing cylinder segment") {
        "shared-start-bridge-references-missing-cylinder-segment"
    } else {
        "face-family-bridge-preflight-rejected"
    }
}
