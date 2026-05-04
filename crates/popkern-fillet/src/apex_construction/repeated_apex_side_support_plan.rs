//! First-class selected plans for repeated-apex side-support region rewrites.

use popkern_core::topology::{FaceHandle, LoopHandle, VertexHandle};

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
    fillet_index: usize,
    support_face: FaceHandle,
    support_loop: LoopHandle,
    start: VertexHandle,
    end: VertexHandle,
    expected_delta: ExpectedTopologyDelta,
    preconditions: BridgeSupportRegionPreconditions,
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
    MissingValidatedBridgeExecutor,
    NoAcceptedOperationCandidate,
}

impl SideSupportRegionPlanSelection {
    pub(super) fn select(candidates: &CandidateSetSnapshot) -> Self {
        let mut operations = Vec::new();
        let mut saw_unknown_endpoint = false;
        let mut selected_unexecutable_bridge = false;
        for candidate in &candidates.candidates {
            match (&candidate.kind, &candidate.status) {
                (RegionOperationKind::SplitFaceWithEdge, CandidateStatus::Accepted) => {
                    operations.push(PlannedTopologyOperation::SplitFaceWithEdge {
                        fillet_index: candidate.fillet_index,
                    });
                }
                (RegionOperationKind::BridgeSupportRegion, CandidateStatus::Accepted) => {
                    if let Some(plan) = BridgeSupportRegionPlan::from_candidate(candidate) {
                        operations.push(PlannedTopologyOperation::BridgeSupportRegion(plan));
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
                _ => {}
            }
        }
        if !operations.is_empty() && !saw_unknown_endpoint && !selected_unexecutable_bridge {
            return SideSupportRegionPlanSelection::Selected(SelectedSideSupportRegionPlan::new(
                operations,
            ));
        }
        if !operations.is_empty() && !saw_unknown_endpoint && selected_unexecutable_bridge {
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
}

impl BridgeSupportRegionPlan {
    fn from_candidate(candidate: &RegionOperationCandidate) -> Option<Self> {
        let (support_face, support_loop) = candidate.endpoints.support_loop_membership()?;
        Some(Self {
            fillet_index: candidate.fillet_index,
            support_face,
            support_loop,
            start: candidate.endpoints.start,
            end: candidate.endpoints.end,
            expected_delta: ExpectedTopologyDelta {
                vertices: 0,
                edges: 1,
                faces: 1,
                loops: 1,
            },
            preconditions: BridgeSupportRegionPreconditions {
                exactly_one_endpoint_on_support_loop: true,
                other_endpoint_has_known_membership: true,
                support_owner_is_unique: support_face == candidate.endpoints.support_face,
            },
        })
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] bridge-support-region support={:?}/{:?} {:?}->{:?} delta={} preconditions={}",
            self.fillet_index,
            self.support_face,
            self.support_loop,
            self.start,
            self.end,
            self.expected_delta.diagnostic_label(),
            self.preconditions.diagnostic_label()
        )
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
