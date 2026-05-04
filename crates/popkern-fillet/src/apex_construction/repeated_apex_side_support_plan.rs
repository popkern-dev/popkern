//! First-class selected plans for repeated-apex side-support region rewrites.

use std::collections::BTreeMap;

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_region_candidates::{
    CandidateRejection, CandidateSetSnapshot, CandidateStatus, RegionOperationCandidate,
    RegionOperationKind,
};

#[derive(Clone, Debug)]
pub(super) enum SideSupportRegionPlanSelection {
    Selected(SelectedSideSupportRegionPlan),
    Blocked(SelectedSideSupportRegionPlan, SideSupportRegionRefusal),
    Refused(SideSupportRegionRefusal),
}

#[derive(Clone, Debug)]
pub(super) struct SelectedSideSupportRegionPlan {
    operations: Vec<PlannedTopologyOperation>,
    expected_delta: ExpectedTopologyDelta,
}

#[derive(Clone, Debug)]
pub(super) enum PlannedTopologyOperation {
    SplitFaceWithEdge { fillet_index: usize },
    BridgeSupportRegion(BridgeSupportRegionPlan),
}

#[derive(Clone, Debug)]
pub(super) struct BridgeSupportRegionPlan {
    support_face: FaceHandle,
    trims: Vec<BridgeSupportTrimPlan>,
    expected_delta: ExpectedTopologyDelta,
    preconditions: BridgeSupportRegionPreconditions,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeSupportTrimPlan {
    fillet_index: usize,
    support_loop: LoopHandle,
    support_half_edge: HalfEdgeHandle,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    outside_half_edge: HalfEdgeHandle,
    start: VertexHandle,
    end: VertexHandle,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct ExpectedTopologyDelta {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeSupportRegionPreconditions {
    exactly_one_endpoint_on_support_loop: bool,
    other_endpoint_has_known_membership: bool,
    support_owner_is_unique: bool,
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionRefusal {
    pub(super) stage: PlannerStage,
    pub(super) reason: SideSupportRegionRefusalReason,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum PlannerStage {
    SelectPlan,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum SideSupportRegionRefusalReason {
    UnknownEndpointMembership,
    AmbiguousEndpointMembership,
    UnsupportedBridgeTrimCount,
    MissingValidatedBridgeExecutor,
    NoAcceptedOperationCandidate,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BridgeSupportRegionPlanRejection {
    EndpointMembershipUnknown,
    UnsupportedTrimCount,
}

impl SideSupportRegionPlanSelection {
    pub(super) fn select(candidates: &CandidateSetSnapshot) -> Self {
        let mut operations = Vec::new();
        let mut bridge_candidates_by_support = BTreeMap::<FaceHandle, Vec<_>>::new();
        let mut saw_unknown_endpoint = false;
        let mut saw_ambiguous_endpoint = false;
        let mut saw_unsupported_bridge_region = false;
        let mut selected_unexecutable_bridge = false;
        for candidate in &candidates.candidates {
            match (&candidate.kind, &candidate.status) {
                (RegionOperationKind::SplitFaceWithEdge, CandidateStatus::Accepted) => {
                    operations.push(PlannedTopologyOperation::SplitFaceWithEdge {
                        fillet_index: candidate.fillet_index,
                    });
                }
                (RegionOperationKind::BridgeSupportRegion, CandidateStatus::Accepted) => {
                    if candidate.endpoints.bridge_loop_membership().is_some() {
                        bridge_candidates_by_support
                            .entry(candidate.endpoints.support_face)
                            .or_default()
                            .push(candidate);
                        selected_unexecutable_bridge = true;
                    } else {
                        saw_unknown_endpoint = true;
                    }
                }
                (
                    _,
                    CandidateStatus::Rejected {
                        reason: CandidateRejection::EndpointMembershipUnknown,
                    },
                ) => saw_unknown_endpoint = true,
                (
                    _,
                    CandidateStatus::Rejected {
                        reason: CandidateRejection::EndpointMembershipAmbiguous,
                    },
                ) => saw_ambiguous_endpoint = true,
                _ => {}
            }
        }
        for (support_face, bridge_candidates) in bridge_candidates_by_support {
            match BridgeSupportRegionPlan::from_candidates(support_face, &bridge_candidates) {
                Ok(plan) => operations.push(PlannedTopologyOperation::BridgeSupportRegion(plan)),
                Err(BridgeSupportRegionPlanRejection::EndpointMembershipUnknown) => {
                    saw_unknown_endpoint = true;
                }
                Err(BridgeSupportRegionPlanRejection::UnsupportedTrimCount) => {
                    saw_unsupported_bridge_region = true;
                }
            }
        }
        if !operations.is_empty()
            && !saw_unknown_endpoint
            && !saw_ambiguous_endpoint
            && !saw_unsupported_bridge_region
            && !selected_unexecutable_bridge
        {
            return SideSupportRegionPlanSelection::Selected(SelectedSideSupportRegionPlan::new(
                operations,
            ));
        }
        if !operations.is_empty()
            && !saw_unknown_endpoint
            && !saw_ambiguous_endpoint
            && !saw_unsupported_bridge_region
            && selected_unexecutable_bridge
        {
            return SideSupportRegionPlanSelection::Blocked(
                SelectedSideSupportRegionPlan::new(operations),
                SideSupportRegionRefusal {
                    stage: PlannerStage::SelectPlan,
                    reason: SideSupportRegionRefusalReason::MissingValidatedBridgeExecutor,
                },
            );
        }
        let reason = if saw_unknown_endpoint {
            SideSupportRegionRefusalReason::UnknownEndpointMembership
        } else if saw_ambiguous_endpoint {
            SideSupportRegionRefusalReason::AmbiguousEndpointMembership
        } else if saw_unsupported_bridge_region {
            SideSupportRegionRefusalReason::UnsupportedBridgeTrimCount
        } else {
            SideSupportRegionRefusalReason::NoAcceptedOperationCandidate
        };
        SideSupportRegionPlanSelection::Refused(SideSupportRegionRefusal {
            stage: PlannerStage::SelectPlan,
            reason,
        })
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionPlanSelection::Selected(plan) => plan.diagnostic_label(),
            SideSupportRegionPlanSelection::Blocked(plan, refusal) => {
                format!(
                    "{}; {}",
                    plan.diagnostic_label(),
                    refusal.diagnostic_label()
                )
            }
            SideSupportRegionPlanSelection::Refused(refusal) => refusal.diagnostic_label(),
        }
    }

    pub(super) fn execution_readiness(&self) -> SideSupportRegionExecutionReadiness {
        match self {
            SideSupportRegionPlanSelection::Selected(_) => {
                SideSupportRegionExecutionReadiness::Ready
            }
            SideSupportRegionPlanSelection::Blocked(_, refusal)
            | SideSupportRegionPlanSelection::Refused(refusal) => {
                SideSupportRegionExecutionReadiness::Blocked {
                    planner_stage: refusal.stage,
                    planner_reason: refusal.reason,
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum SideSupportRegionExecutionReadiness {
    Ready,
    Blocked {
        planner_stage: PlannerStage,
        planner_reason: SideSupportRegionRefusalReason,
    },
}

#[derive(Clone, Debug)]
pub(super) enum SideSupportRegionLivePreconditionFailure {
    EndpointHalfEdgeMismatch {
        support_half_edge: HalfEdgeHandle,
        outside_half_edge: HalfEdgeHandle,
        start: VertexHandle,
        end: VertexHandle,
    },
    MissingFace(FaceHandle),
    MissingHalfEdge(HalfEdgeHandle),
    MissingLoop(LoopHandle),
    MissingVertex(VertexHandle),
    HalfEdgeLoopMismatch {
        half_edge: HalfEdgeHandle,
        expected_loop: LoopHandle,
        actual_loop: LoopHandle,
    },
}

impl SideSupportRegionLivePreconditionFailure {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionLivePreconditionFailure::EndpointHalfEdgeMismatch {
                support_half_edge,
                outside_half_edge,
                start,
                end,
            } => {
                format!(
                    "endpoint-half-edge-mismatch support={support_half_edge:?} outside={outside_half_edge:?} endpoints={start:?}->{end:?}"
                )
            }
            SideSupportRegionLivePreconditionFailure::MissingFace(face) => {
                format!("missing-face {face:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingHalfEdge(half_edge) => {
                format!("missing-half-edge {half_edge:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingLoop(loop_handle) => {
                format!("missing-loop {loop_handle:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingVertex(vertex) => {
                format!("missing-vertex {vertex:?}")
            }
            SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                half_edge,
                expected_loop,
                actual_loop,
            } => {
                format!(
                    "half-edge-loop-mismatch {half_edge:?} expected={expected_loop:?} actual={actual_loop:?}"
                )
            }
        }
    }
}

impl SelectedSideSupportRegionPlan {
    fn new(operations: Vec<PlannedTopologyOperation>) -> Self {
        let expected_delta = operations
            .iter()
            .fold(ExpectedTopologyDelta::zero(), |acc, operation| {
                acc.plus(operation.expected_delta())
            });
        Self {
            operations,
            expected_delta,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "selected_plan total_delta={} operations=[{}]",
            self.expected_delta.diagnostic_label(),
            self.operations
                .iter()
                .map(PlannedTopologyOperation::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl SideSupportRegionPlanSelection {
    pub(super) fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        match self {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => plan.validate_live_operands(brep),
            SideSupportRegionPlanSelection::Refused(_) => Ok(()),
        }
    }
}

impl SelectedSideSupportRegionPlan {
    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        for operation in &self.operations {
            operation.validate_live_operands(brep)?;
        }
        Ok(())
    }
}

impl PlannedTopologyOperation {
    fn expected_delta(&self) -> ExpectedTopologyDelta {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { .. } => ExpectedTopologyDelta {
                vertices: 0,
                edges: 1,
                faces: 1,
                loops: 1,
            },
            PlannedTopologyOperation::BridgeSupportRegion(plan) => plan.expected_delta,
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { fillet_index } => {
                format!("FilletEnd[{fillet_index}] split-face-with-edge")
            }
            PlannedTopologyOperation::BridgeSupportRegion(plan) => plan.diagnostic_label(),
        }
    }

    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { .. } => Ok(()),
            PlannedTopologyOperation::BridgeSupportRegion(plan) => {
                plan.validate_live_operands(brep)
            }
        }
    }
}

impl BridgeSupportRegionPlan {
    fn from_candidates(
        support_face: FaceHandle,
        candidates: &[&RegionOperationCandidate],
    ) -> Result<Self, BridgeSupportRegionPlanRejection> {
        if candidates.len() != 2 {
            return Err(BridgeSupportRegionPlanRejection::UnsupportedTrimCount);
        }
        let mut trims = Vec::new();
        for candidate in candidates {
            let bridge_membership = candidate
                .endpoints
                .bridge_loop_membership()
                .ok_or(BridgeSupportRegionPlanRejection::EndpointMembershipUnknown)?;
            trims.push(BridgeSupportTrimPlan {
                fillet_index: candidate.fillet_index,
                support_loop: bridge_membership.support_loop,
                support_half_edge: bridge_membership.support_half_edge,
                outside_face: bridge_membership.outside_face,
                outside_loop: bridge_membership.outside_loop,
                outside_half_edge: bridge_membership.outside_half_edge,
                start: candidate.endpoints.start,
                end: candidate.endpoints.end,
            });
        }
        let expected_delta = ExpectedTopologyDelta {
            vertices: 0,
            edges: trims.len() as i64,
            faces: trims.len() as i64,
            loops: trims.len() as i64,
        };
        Ok(Self {
            support_face,
            trims,
            expected_delta,
            preconditions: BridgeSupportRegionPreconditions {
                exactly_one_endpoint_on_support_loop: true,
                other_endpoint_has_known_membership: true,
                support_owner_is_unique: true,
            },
        })
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "bridge-support-region support={:?} trims=[{}] delta={} preconditions={}",
            self.support_face,
            self.trims
                .iter()
                .map(BridgeSupportTrimPlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.expected_delta.diagnostic_label(),
            self.preconditions.diagnostic_label()
        )
    }

    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        if !brep.faces.contains_key(self.support_face) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.support_face,
            ));
        }
        for trim in &self.trims {
            trim.validate_live_operands(brep)?;
        }
        Ok(())
    }
}

impl BridgeSupportTrimPlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] support_loop={:?}@{:?} outside={:?}/{:?}@{:?} {:?}->{:?}",
            self.fillet_index,
            self.support_loop,
            self.support_half_edge,
            self.outside_face,
            self.outside_loop,
            self.outside_half_edge,
            self.start,
            self.end
        )
    }

    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        if !brep.loops.contains_key(self.support_loop) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.support_loop,
            ));
        }
        let Some(support_half_edge) = brep.half_edges.get(self.support_half_edge) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.support_half_edge,
            ));
        };
        if support_half_edge.loop_ != self.support_loop {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.support_half_edge,
                    expected_loop: self.support_loop,
                    actual_loop: support_half_edge.loop_,
                },
            );
        }
        if !brep.faces.contains_key(self.outside_face) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.outside_face,
            ));
        }
        if !brep.loops.contains_key(self.outside_loop) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.outside_loop,
            ));
        }
        let Some(outside_half_edge) = brep.half_edges.get(self.outside_half_edge) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.outside_half_edge,
            ));
        };
        if outside_half_edge.loop_ != self.outside_loop {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.outside_half_edge,
                    expected_loop: self.outside_loop,
                    actual_loop: outside_half_edge.loop_,
                },
            );
        }
        if !brep.vertices.contains_key(self.start) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.start,
            ));
        }
        if !brep.vertices.contains_key(self.end) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.end,
            ));
        }
        let support_origin = support_half_edge.origin;
        let outside_origin = outside_half_edge.origin;
        if !((support_origin == self.start && outside_origin == self.end)
            || (support_origin == self.end && outside_origin == self.start))
        {
            return Err(
                SideSupportRegionLivePreconditionFailure::EndpointHalfEdgeMismatch {
                    support_half_edge: self.support_half_edge,
                    outside_half_edge: self.outside_half_edge,
                    start: self.start,
                    end: self.end,
                },
            );
        }
        Ok(())
    }
}

impl ExpectedTopologyDelta {
    fn zero() -> Self {
        Self {
            vertices: 0,
            edges: 0,
            faces: 0,
            loops: 0,
        }
    }

    fn plus(self, other: Self) -> Self {
        Self {
            vertices: self.vertices + other.vertices,
            edges: self.edges + other.edges,
            faces: self.faces + other.faces,
            loops: self.loops + other.loops,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "dV={} dE={} dF={} dL={}",
            self.vertices, self.edges, self.faces, self.loops
        )
    }
}

impl BridgeSupportRegionPreconditions {
    fn diagnostic_label(&self) -> String {
        format!(
            "one_endpoint_on_support_loop={} other_endpoint_known={} unique_support_owner={}",
            self.exactly_one_endpoint_on_support_loop,
            self.other_endpoint_has_known_membership,
            self.support_owner_is_unique
        )
    }
}

impl SideSupportRegionRefusal {
    fn diagnostic_label(&self) -> String {
        format!("refusal stage={:?} reason={:?}", self.stage, self.reason)
    }
}
