//! Candidate bridge authoring for shared-start split cylinder resources.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::euler::{FaceFamilyBridgeCycleRequest, SameLoopFaceBridgeRequest};
use crate::topology::{BRep, HalfEdgeHandle};

use super::CylinderMediatedSupportContinuationSegmentAuthoringPlan;
use super::apply_cylinder_resources::{BridgeResourcePairing, CylinderBridgeResourceSegment};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SharedStartBridgePreflight {
    status: SharedStartBridgePreflightStatus,
    pairings: usize,
    cycles: usize,
    cross_loop_cycles: usize,
    same_loop_cycles: usize,
    same_loop_overlap_cycles: usize,
    same_loop_overlapping_edges: usize,
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
            "shared_start_bridge_preflight=status={status} pairings={} cycles={} cross_loop_cycles={} same_loop_cycles={} same_loop_overlap_cycles={} same_loop_overlapping_edges={} expected_delta=dV=0 dE={} dF={} dL={} reason={}",
            self.pairings,
            self.cycles,
            self.cross_loop_cycles,
            self.same_loop_cycles,
            self.same_loop_overlap_cycles,
            self.same_loop_overlapping_edges,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SharedStartCycleShape {
    cycles: usize,
    cross_loop_cycles: usize,
    same_loop_cycles: usize,
    same_loop_overlap_cycles: usize,
    same_loop_overlapping_edges: usize,
}

impl SharedStartCycleShape {
    fn empty() -> Self {
        Self {
            cycles: 0,
            cross_loop_cycles: 0,
            same_loop_cycles: 0,
            same_loop_overlap_cycles: 0,
            same_loop_overlapping_edges: 0,
        }
    }

    fn from_cycles(cycles: &[FaceFamilyBridgeCycleRequest], brep: &BRep) -> Self {
        let mut same_loop_cycles = 0;
        let mut same_loop_overlap_cycles = 0;
        let mut same_loop_overlapping_edges = 0;
        for cycle in cycles
            .iter()
            .filter(|cycle| cycle.first_loop == cycle.second_loop)
        {
            same_loop_cycles += 1;
            let overlap_edges = same_loop_overlap_edge_count(brep, cycle);
            if overlap_edges > 0 {
                same_loop_overlap_cycles += 1;
                same_loop_overlapping_edges += overlap_edges;
            }
        }
        Self {
            cycles: cycles.len(),
            cross_loop_cycles: cycles.len() - same_loop_cycles,
            same_loop_cycles,
            same_loop_overlap_cycles,
            same_loop_overlapping_edges,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SharedStartCycleMixSummary {
    cycles: usize,
    cross_loop_cycles: usize,
    same_loop_cycles: usize,
    same_loop_overlap_cycles: usize,
    same_loop_overlapping_edges: usize,
    expected_edges: i64,
    expected_faces: i64,
    expected_loops: i64,
}

fn preflight_shared_start_cycle_mix(
    brep: &BRep,
    cycles: Vec<FaceFamilyBridgeCycleRequest>,
) -> Result<SharedStartCycleMixSummary, CoreError> {
    let cycle_shape = SharedStartCycleShape::from_cycles(&cycles, brep);
    let mut cross_loop_cycles = Vec::new();
    let mut same_loop_cycles = Vec::new();
    for cycle in cycles {
        if cycle.first_loop == cycle.second_loop {
            same_loop_cycles.push(cycle);
        } else {
            cross_loop_cycles.push(cycle);
        }
    }
    let mut expected_edges = 0;
    let mut expected_faces = 0;
    let mut expected_loops = 0;
    let cross_loop_count = cross_loop_cycles.len();
    let same_loop_count = same_loop_cycles.len();
    if !cross_loop_cycles.is_empty() {
        let plan = brep.preflight_face_family_bridge(cross_loop_cycles)?;
        expected_edges += plan.expected_delta.edges;
        expected_faces += plan.expected_delta.faces;
        expected_loops += plan.expected_delta.loops;
    }
    for cycle in &same_loop_cycles {
        let plan = brep.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_: cycle.first_loop,
            first_interval: cycle.first_interval.clone(),
            second_interval: cycle.second_interval.clone(),
        })?;
        expected_edges += plan.expected_delta.edges;
        expected_faces += plan.expected_delta.faces;
        expected_loops += plan.expected_delta.loops;
    }
    Ok(SharedStartCycleMixSummary {
        cycles: cross_loop_count + same_loop_count,
        cross_loop_cycles: cross_loop_count,
        same_loop_cycles: same_loop_count,
        same_loop_overlap_cycles: cycle_shape.same_loop_overlap_cycles,
        same_loop_overlapping_edges: cycle_shape.same_loop_overlapping_edges,
        expected_edges,
        expected_faces,
        expected_loops,
    })
}

fn same_loop_overlap_edge_count(brep: &BRep, cycle: &FaceFamilyBridgeCycleRequest) -> usize {
    let first_edges = match bridge_path_edges(brep, &cycle.first_interval) {
        Ok(edges) => edges,
        Err(_) => return 0,
    };
    let second_edges = match bridge_path_edges(brep, &cycle.second_interval) {
        Ok(edges) => edges,
        Err(_) => return 0,
    };
    first_edges.intersection(&second_edges).count()
}

fn bridge_path_edges(
    brep: &BRep,
    interval: &[HalfEdgeHandle],
) -> Result<BTreeSet<crate::topology::EdgeHandle>, CoreError> {
    let mut edges = BTreeSet::new();
    for &half_edge in half_open_bridge_path(interval)? {
        let edge = brep
            .half_edges
            .get(half_edge)
            .map(|data| data.edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle("shared-start bridge interval has stale half-edge".into())
            })?;
        edges.insert(edge);
    }
    Ok(edges)
}

fn half_open_bridge_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "shared-start bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
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
    } else {
        "face-family-bridge-preflight-rejected"
    }
}
