//! Candidate bridge authoring for shared-start split cylinder resources.

use crate::error::CoreError;
use crate::euler::{
    BoundarySideSeparatorBridgeRequest, BridgeDiagnosticOperationRequest, BridgeDiagnosticReport,
    BridgeDiagnosticRequest, BridgeDiagnosticSummary, CandidateMissingFacet,
    DiagnosticOperationPurpose, FaceFamilyBridgeCycleRequest,
};
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
    one_edge_interior_overlap_candidates: usize,
    repeated_one_edge_overlap_occurrences: usize,
    one_edge_interior_overlap_plans: usize,
    one_edge_interior_overlap_boundary_edges: usize,
    one_edge_interior_overlap_normalizations: usize,
    graph_summary: BridgeDiagnosticSummary,
    expected_edges: i64,
    expected_faces: i64,
    expected_loops: i64,
    candidate_gate_complete: bool,
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
            one_edge_interior_overlap_candidates: 0,
            repeated_one_edge_overlap_occurrences: 0,
            one_edge_interior_overlap_plans: 0,
            one_edge_interior_overlap_boundary_edges: 0,
            one_edge_interior_overlap_normalizations: 0,
            graph_summary: empty_bridge_diagnostic_summary(),
            expected_edges: 0,
            expected_faces: 0,
            expected_loops: 0,
            candidate_gate_complete: true,
            reason: "none",
        }
    }

    fn ready(
        pairings: usize,
        summary: SharedStartCycleMixSummary,
        graph_summary: BridgeDiagnosticSummary,
    ) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Ready,
            pairings,
            cycles: summary.cycles,
            cross_loop_cycles: summary.cross_loop_cycles,
            same_loop_cycles: summary.same_loop_cycles,
            one_edge_interior_overlap_candidates: summary.one_edge_interior_overlap_candidates,
            repeated_one_edge_overlap_occurrences: summary.repeated_one_edge_overlap_occurrences,
            one_edge_interior_overlap_plans: summary.one_edge_interior_overlap_plans,
            one_edge_interior_overlap_boundary_edges: summary
                .one_edge_interior_overlap_boundary_edges,
            one_edge_interior_overlap_normalizations: summary
                .one_edge_interior_overlap_normalizations,
            graph_summary,
            expected_edges: summary.expected_edges,
            expected_faces: summary.expected_faces,
            expected_loops: summary.expected_loops,
            candidate_gate_complete: true,
            reason: "none",
        }
    }

    fn rejected(
        pairings: usize,
        shape: SharedStartCycleShape,
        graph_summary: BridgeDiagnosticSummary,
        candidate_gate_complete: bool,
        reason: &'static str,
    ) -> Self {
        Self {
            status: SharedStartBridgePreflightStatus::Rejected,
            pairings,
            cycles: shape.cycles,
            cross_loop_cycles: shape.cross_loop_cycles,
            same_loop_cycles: shape.same_loop_cycles,
            one_edge_interior_overlap_candidates: shape.one_edge_interior_overlap_candidates,
            repeated_one_edge_overlap_occurrences: shape.repeated_one_edge_overlap_occurrences,
            one_edge_interior_overlap_plans: shape.one_edge_interior_overlap_plans,
            one_edge_interior_overlap_boundary_edges: shape
                .one_edge_interior_overlap_boundary_edges,
            one_edge_interior_overlap_normalizations: 0,
            graph_summary,
            expected_edges: 0,
            expected_faces: 0,
            expected_loops: 0,
            candidate_gate_complete,
            reason,
        }
    }

    fn rejected_with_delta(
        pairings: usize,
        shape: SharedStartCycleShape,
        graph_summary: BridgeDiagnosticSummary,
        expected_delta: SharedStartExpectedDelta,
        candidate_gate_complete: bool,
        reason: &'static str,
    ) -> Self {
        Self {
            expected_edges: expected_delta.edges,
            expected_faces: expected_delta.faces,
            expected_loops: expected_delta.loops,
            ..Self::rejected(
                pairings,
                shape,
                graph_summary,
                candidate_gate_complete,
                reason,
            )
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        let status = match self.status {
            SharedStartBridgePreflightStatus::NotApplicable => "not-applicable",
            SharedStartBridgePreflightStatus::Ready => "ready",
            SharedStartBridgePreflightStatus::Rejected => "rejected",
        };
        format!(
            "shared_start_bridge_preflight=status={status} pairings={} cycles={} cross_loop_cycles={} same_loop_cycles={} one_edge_interior_overlap_candidates={} repeated_one_edge_overlap_occurrences={} one_edge_interior_overlap_plans={} one_edge_interior_overlap_boundary_edges={} one_edge_interior_overlap_normalizations={} {} expected_delta=dV=0 dE={} dF={} dL={} candidate_gate_complete={} reason={}",
            self.pairings,
            self.cycles,
            self.cross_loop_cycles,
            self.same_loop_cycles,
            self.one_edge_interior_overlap_candidates,
            self.repeated_one_edge_overlap_occurrences,
            self.one_edge_interior_overlap_plans,
            self.one_edge_interior_overlap_boundary_edges,
            self.one_edge_interior_overlap_normalizations,
            self.graph_summary.diagnostic_label(),
            self.expected_edges,
            self.expected_faces,
            self.expected_loops,
            self.candidate_gate_complete,
            self.reason
        )
    }

    pub(super) fn requires_candidate_coverage(&self) -> bool {
        self.status != SharedStartBridgePreflightStatus::NotApplicable
    }

    pub(super) fn candidate_count(&self) -> usize {
        if self.requires_candidate_coverage() {
            self.graph_summary.construction_candidates
        } else {
            0
        }
    }

    pub(super) fn candidate_gate_complete(&self) -> bool {
        self.candidate_gate_complete
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SharedStartExpectedDelta {
    edges: i64,
    faces: i64,
    loops: i64,
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
                empty_bridge_diagnostic_summary(),
                false,
                bridge_error_label(&err),
            );
        }
    };
    let cycle_shape = SharedStartCycleShape::from_cycles(&cycles, brep);
    let report = match bridge_diagnostic_report(brep, &cycles) {
        Ok(report) => report,
        Err(err) => {
            return SharedStartBridgePreflight::rejected(
                pairings.len(),
                cycle_shape,
                empty_bridge_diagnostic_summary(),
                false,
                bridge_error_label(&err),
            );
        }
    };
    let graph_summary = report.summary;
    if let Err(reason) = validate_shared_start_candidate_gate(&report) {
        return SharedStartBridgePreflight::rejected(
            pairings.len(),
            cycle_shape,
            graph_summary,
            false,
            reason,
        );
    }
    match preflight_shared_start_cycle_mix(brep, cycles.clone()) {
        Ok(summary) => SharedStartBridgePreflight::ready(pairings.len(), summary, graph_summary),
        Err(err) => match brep
            .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest { cycles })
        {
            Ok(plan) => SharedStartBridgePreflight::rejected_with_delta(
                pairings.len(),
                cycle_shape,
                graph_summary,
                SharedStartExpectedDelta {
                    edges: plan.expected_delta.edges,
                    faces: plan.expected_delta.faces,
                    loops: plan.expected_delta.loops,
                },
                true,
                "boundary-side-separator-bridge-mutation-required",
            ),
            Err(_) => SharedStartBridgePreflight::rejected(
                pairings.len(),
                cycle_shape,
                graph_summary,
                true,
                bridge_error_label(&err),
            ),
        },
    }
}

fn bridge_diagnostic_report(
    brep: &BRep,
    cycles: &[FaceFamilyBridgeCycleRequest],
) -> Result<BridgeDiagnosticReport, CoreError> {
    brep.bridge_diagnostic_report(BridgeDiagnosticRequest {
        operations: cycles
            .iter()
            .map(|cycle| BridgeDiagnosticOperationRequest {
                purpose: DiagnosticOperationPurpose::BridgeCycle,
                first_loop: cycle.first_loop,
                first_interval: cycle.first_interval.clone(),
                second_loop: cycle.second_loop,
                second_interval: cycle.second_interval.clone(),
            })
            .collect(),
    })
}

fn validate_shared_start_candidate_gate(
    report: &BridgeDiagnosticReport,
) -> Result<(), &'static str> {
    if report.graph.candidates.is_empty() {
        return Err("shared-start-bridge-candidate-gate-missing-candidate");
    }
    for operation in &report.graph.operations {
        if !report
            .graph
            .candidates
            .iter()
            .any(|candidate| candidate.operations.contains(&operation.id))
        {
            return Err("shared-start-bridge-candidate-gate-missing-operation-candidate");
        }
    }
    for candidate in &report.graph.candidates {
        if let Err(err) = candidate.validate_mutation_gate() {
            return Err(candidate_gate_error_label(err.missing_facet));
        }
    }
    Ok(())
}

fn candidate_gate_error_label(missing_facet: CandidateMissingFacet) -> &'static str {
    match missing_facet {
        CandidateMissingFacet::Family => "shared-start-bridge-candidate-gate-missing-family",
        CandidateMissingFacet::Feasibility => {
            "shared-start-bridge-candidate-gate-missing-feasibility"
        }
        CandidateMissingFacet::TopologyResource => {
            "shared-start-bridge-candidate-gate-missing-topology-resource"
        }
        CandidateMissingFacet::PostconditionPreview => {
            "shared-start-bridge-candidate-gate-missing-postcondition-preview"
        }
        CandidateMissingFacet::ExpectedDelta => {
            "shared-start-bridge-candidate-gate-missing-expected-delta"
        }
        CandidateMissingFacet::GeometryFeasibility => {
            "shared-start-bridge-candidate-gate-missing-geometry-feasibility"
        }
        CandidateMissingFacet::ValidityRepair => {
            "shared-start-bridge-candidate-gate-missing-validity-repair"
        }
        CandidateMissingFacet::OperationDependency => {
            "shared-start-bridge-candidate-gate-missing-operation-dependency"
        }
        CandidateMissingFacet::SemanticHistoryDependency => {
            "shared-start-bridge-candidate-gate-missing-semantic-history"
        }
        CandidateMissingFacet::DownstreamReadiness => {
            "shared-start-bridge-candidate-gate-missing-downstream-readiness"
        }
        CandidateMissingFacet::LiveTopologyRevalidation => {
            "shared-start-bridge-candidate-gate-missing-live-revalidation"
        }
        CandidateMissingFacet::PostconditionVerification => {
            "shared-start-bridge-candidate-gate-missing-postcondition-verification"
        }
    }
}

fn empty_bridge_diagnostic_summary() -> BridgeDiagnosticSummary {
    BridgeDiagnosticSummary::empty()
}

pub(super) fn shared_start_bridge_cycles(
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

#[cfg(test)]
#[path = "cylinder_mediated_support_apply_shared_start_bridge_tests.rs"]
mod tests;
