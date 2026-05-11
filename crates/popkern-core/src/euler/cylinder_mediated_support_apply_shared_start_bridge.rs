//! Candidate bridge authoring for shared-start split cylinder resources.

use crate::error::CoreError;
use crate::euler::FaceFamilyBridgeCycleRequest;
use crate::topology::{BRep, HalfEdgeHandle};

use super::CylinderMediatedSupportContinuationSegmentAuthoringPlan;
use super::apply_cylinder_resources::{BridgeResourcePairing, CylinderBridgeResourceSegment};
use super::apply_shared_start_cycle_mix::{
    SharedStartCycleMixSummary, SharedStartCycleShape, preflight_shared_start_cycle_mix,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SharedStartBridgePreflight {
    status: SharedStartBridgePreflightStatus,
    pairings: usize,
    cycles: usize,
    cross_loop_cycles: usize,
    same_loop_cycles: usize,
    same_loop_overlap_cycles: usize,
    same_loop_overlapping_edges: usize,
    same_loop_shared_start_overlaps: usize,
    same_loop_shared_end_overlaps: usize,
    same_loop_crossed_endpoint_overlaps: usize,
    same_loop_interior_overlaps: usize,
    one_edge_interior_overlap_candidates: usize,
    repeated_one_edge_overlap_occurrences: usize,
    one_edge_interior_overlap_plans: usize,
    one_edge_interior_overlap_boundary_edges: usize,
    one_edge_interior_overlap_normalizations: usize,
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
            cross_loop_cycles: 0,
            same_loop_cycles: 0,
            same_loop_overlap_cycles: 0,
            same_loop_overlapping_edges: 0,
            same_loop_shared_start_overlaps: 0,
            same_loop_shared_end_overlaps: 0,
            same_loop_crossed_endpoint_overlaps: 0,
            same_loop_interior_overlaps: 0,
            one_edge_interior_overlap_candidates: 0,
            repeated_one_edge_overlap_occurrences: 0,
            one_edge_interior_overlap_plans: 0,
            one_edge_interior_overlap_boundary_edges: 0,
            one_edge_interior_overlap_normalizations: 0,
            expected_edges: 0,
            expected_faces: 0,
            expected_loops: 0,
            reason: "none",
        }
    }

    fn ready(pairings: usize, summary: SharedStartCycleMixSummary) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Ready,
            pairings,
            cycles: summary.cycles,
            cross_loop_cycles: summary.cross_loop_cycles,
            same_loop_cycles: summary.same_loop_cycles,
            same_loop_overlap_cycles: summary.same_loop_overlap_cycles,
            same_loop_overlapping_edges: summary.same_loop_overlapping_edges,
            same_loop_shared_start_overlaps: summary.same_loop_shared_start_overlaps,
            same_loop_shared_end_overlaps: summary.same_loop_shared_end_overlaps,
            same_loop_crossed_endpoint_overlaps: summary.same_loop_crossed_endpoint_overlaps,
            same_loop_interior_overlaps: summary.same_loop_interior_overlaps,
            one_edge_interior_overlap_candidates: summary.one_edge_interior_overlap_candidates,
            repeated_one_edge_overlap_occurrences: summary.repeated_one_edge_overlap_occurrences,
            one_edge_interior_overlap_plans: summary.one_edge_interior_overlap_plans,
            one_edge_interior_overlap_boundary_edges: summary
                .one_edge_interior_overlap_boundary_edges,
            one_edge_interior_overlap_normalizations: summary
                .one_edge_interior_overlap_normalizations,
            expected_edges: summary.expected_edges,
            expected_faces: summary.expected_faces,
            expected_loops: summary.expected_loops,
            reason: "none",
        }
    }

    fn rejected(pairings: usize, shape: SharedStartCycleShape, reason: &'static str) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Rejected,
            pairings,
            cycles: shape.cycles,
            cross_loop_cycles: shape.cross_loop_cycles,
            same_loop_cycles: shape.same_loop_cycles,
            same_loop_overlap_cycles: shape.same_loop_overlap_cycles,
            same_loop_overlapping_edges: shape.same_loop_overlapping_edges,
            same_loop_shared_start_overlaps: shape.same_loop_shared_start_overlaps,
            same_loop_shared_end_overlaps: shape.same_loop_shared_end_overlaps,
            same_loop_crossed_endpoint_overlaps: shape.same_loop_crossed_endpoint_overlaps,
            same_loop_interior_overlaps: shape.same_loop_interior_overlaps,
            one_edge_interior_overlap_candidates: shape.one_edge_interior_overlap_candidates,
            repeated_one_edge_overlap_occurrences: shape.repeated_one_edge_overlap_occurrences,
            one_edge_interior_overlap_plans: shape.one_edge_interior_overlap_plans,
            one_edge_interior_overlap_boundary_edges: shape
                .one_edge_interior_overlap_boundary_edges,
            one_edge_interior_overlap_normalizations: 0,
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
            "shared_start_bridge_preflight=status={status} pairings={} cycles={} cross_loop_cycles={} same_loop_cycles={} same_loop_overlap_cycles={} same_loop_overlapping_edges={} same_loop_shared_start_overlaps={} same_loop_shared_end_overlaps={} same_loop_crossed_endpoint_overlaps={} same_loop_interior_overlaps={} one_edge_interior_overlap_candidates={} repeated_one_edge_overlap_occurrences={} one_edge_interior_overlap_plans={} one_edge_interior_overlap_boundary_edges={} one_edge_interior_overlap_normalizations={} expected_delta=dV=0 dE={} dF={} dL={} reason={}",
            self.pairings,
            self.cycles,
            self.cross_loop_cycles,
            self.same_loop_cycles,
            self.same_loop_overlap_cycles,
            self.same_loop_overlapping_edges,
            self.same_loop_shared_start_overlaps,
            self.same_loop_shared_end_overlaps,
            self.same_loop_crossed_endpoint_overlaps,
            self.same_loop_interior_overlaps,
            self.one_edge_interior_overlap_candidates,
            self.repeated_one_edge_overlap_occurrences,
            self.one_edge_interior_overlap_plans,
            self.one_edge_interior_overlap_boundary_edges,
            self.one_edge_interior_overlap_normalizations,
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
    let cycles = match shared_start_bridge_cycles(brep, pairings, authoring, cylinder_segments) {
        Ok(cycles) => cycles,
        Err(err) => {
            return SharedStartBridgePreflight::rejected(
                pairings.len(),
                SharedStartCycleShape::empty(),
                bridge_error_label(&err),
            );
        }
    };
    let cycle_shape = SharedStartCycleShape::from_cycles(&cycles, brep);
    match preflight_shared_start_cycle_mix(brep, cycles) {
        Ok(summary) => SharedStartBridgePreflight::ready(pairings.len(), summary),
        Err(err) => SharedStartBridgePreflight::rejected(
            pairings.len(),
            cycle_shape,
            bridge_error_label(&err),
        ),
    }
}

fn shared_start_bridge_cycles(
    brep: &BRep,
    pairings: &[BridgeResourcePairing],
    authoring: &[CylinderMediatedSupportContinuationSegmentAuthoringPlan],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> Result<Vec<FaceFamilyBridgeCycleRequest>, CoreError> {
    let cycles = pairings
        .iter()
        .map(|pairing| bridge_cycle_request(brep, *pairing, authoring, cylinder_segments))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(cycles)
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
    } else if detail.contains("same-loop bridge connector would be degenerate") {
        "same-loop-bridge-connector-would-be-degenerate"
    } else if detail.contains("same-loop bridge intervals overlap") {
        "same-loop-bridge-intervals-overlap"
    } else if detail.contains("same-loop edge is already boundary") {
        "same-loop-overlap-edge-already-boundary"
    } else if detail.contains("face-family bridge cycles overlap") {
        "face-family-bridge-cycles-overlap"
    } else if detail.contains("is not an unowned boundary side") {
        "bridge-cycle-side-not-unowned-boundary"
    } else if detail.contains("source faces are in different shells") {
        "face-family-bridge-source-faces-different-shells"
    } else {
        "face-family-bridge-preflight-rejected"
    }
}
