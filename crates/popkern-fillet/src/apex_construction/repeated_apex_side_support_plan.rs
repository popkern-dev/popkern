//! First-class selected plans for repeated-apex side-support region rewrites.

use std::collections::BTreeMap;

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

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
    OutsideOutsideSupportContinuation(OutsideOutsideSupportContinuationPlan),
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
    cylinder_face: FaceHandle,
    support_loop: LoopHandle,
    support_half_edge: HalfEdgeHandle,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    outside_half_edge: HalfEdgeHandle,
    start: VertexHandle,
    end: VertexHandle,
}

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportContinuationPlan {
    support_face: FaceHandle,
    trims: Vec<OutsideOutsideSupportTrimPlan>,
    expected_delta: ExpectedTopologyDelta,
    preconditions: OutsideOutsideSupportPreconditions,
}

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportTrimPlan {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start_face: FaceHandle,
    start_loop: LoopHandle,
    start_half_edge: HalfEdgeHandle,
    end_face: FaceHandle,
    end_loop: LoopHandle,
    end_half_edge: HalfEdgeHandle,
    start: VertexHandle,
    end: VertexHandle,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum ExpectedTopologyDelta {
    Known {
        vertices: i64,
        edges: i64,
        faces: i64,
        loops: i64,
    },
    Unknown {
        reason: &'static str,
    },
}

#[derive(Clone, Debug)]
pub(super) struct BridgeSupportRegionPreconditions {
    exactly_one_endpoint_on_support_loop: bool,
    other_endpoint_has_known_membership: bool,
    support_owner_is_unique: bool,
    face_family_bridge_operands_known: bool,
}

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportPreconditions {
    both_endpoints_outside_support_family: bool,
    endpoint_membership_known: bool,
    support_owner_is_unique: bool,
    topology_operator_required: bool,
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionRefusal {
    stage: PlannerStage,
    reason: SideSupportRegionRefusalReason,
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
    MissingValidatedFaceFamilyBridgePrimitive,
    MissingValidatedOutsideOutsideSupportContinuationPrimitive,
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
        let mut outside_outside_candidates_by_support = BTreeMap::<FaceHandle, Vec<_>>::new();
        let mut saw_unknown_endpoint = false;
        let mut saw_ambiguous_endpoint = false;
        let mut saw_unsupported_bridge_region = false;
        let mut selected_unexecutable_bridge = false;
        let mut selected_unexecutable_outside_outside = false;
        for candidate in candidates.iter() {
            match (candidate.kind(), candidate.status()) {
                (RegionOperationKind::SplitFaceWithEdge, CandidateStatus::Accepted) => {
                    operations.push(PlannedTopologyOperation::SplitFaceWithEdge {
                        fillet_index: candidate.fillet_index(),
                    });
                }
                (RegionOperationKind::BridgeSupportRegion, CandidateStatus::Accepted) => {
                    if candidate.endpoints().bridge_loop_membership().is_some() {
                        bridge_candidates_by_support
                            .entry(candidate.endpoints().support_face())
                            .or_default()
                            .push(candidate);
                        selected_unexecutable_bridge = true;
                    } else {
                        saw_unknown_endpoint = true;
                    }
                }
                (
                    RegionOperationKind::OutsideOutsideSupportContinuation,
                    CandidateStatus::Accepted,
                ) => {
                    if candidate.endpoints().outside_outside_membership().is_some() {
                        outside_outside_candidates_by_support
                            .entry(candidate.endpoints().support_face())
                            .or_default()
                            .push(candidate);
                        selected_unexecutable_outside_outside = true;
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
        for (support_face, outside_outside_candidates) in outside_outside_candidates_by_support {
            operations.push(PlannedTopologyOperation::OutsideOutsideSupportContinuation(
                OutsideOutsideSupportContinuationPlan::from_candidates(
                    support_face,
                    &outside_outside_candidates,
                ),
            ));
        }
        if !operations.is_empty()
            && !saw_unknown_endpoint
            && !saw_ambiguous_endpoint
            && !saw_unsupported_bridge_region
            && !selected_unexecutable_bridge
            && !selected_unexecutable_outside_outside
        {
            return SideSupportRegionPlanSelection::Selected(SelectedSideSupportRegionPlan::new(
                operations,
            ));
        }
        if !operations.is_empty()
            && !saw_unknown_endpoint
            && !saw_ambiguous_endpoint
            && !saw_unsupported_bridge_region
            && selected_unexecutable_outside_outside
        {
            return SideSupportRegionPlanSelection::Blocked(
                SelectedSideSupportRegionPlan::new(operations),
                SideSupportRegionRefusal {
                    stage: PlannerStage::SelectPlan,
                    reason:
                        SideSupportRegionRefusalReason::MissingValidatedOutsideOutsideSupportContinuationPrimitive,
                },
            );
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
                    reason:
                        SideSupportRegionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive,
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
                    planner_stage: refusal.stage(),
                    planner_reason: refusal.reason(),
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

    pub(super) fn operations(&self) -> impl Iterator<Item = &PlannedTopologyOperation> {
        self.operations.iter()
    }
}

impl PlannedTopologyOperation {
    fn expected_delta(&self) -> ExpectedTopologyDelta {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { .. } => ExpectedTopologyDelta::Known {
                vertices: 0,
                edges: 1,
                faces: 1,
                loops: 1,
            },
            PlannedTopologyOperation::BridgeSupportRegion(plan) => plan.expected_delta.clone(),
            PlannedTopologyOperation::OutsideOutsideSupportContinuation(plan) => {
                plan.expected_delta.clone()
            }
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { fillet_index } => {
                format!("FilletEnd[{fillet_index}] split-face-with-edge")
            }
            PlannedTopologyOperation::BridgeSupportRegion(plan) => plan.diagnostic_label(),
            PlannedTopologyOperation::OutsideOutsideSupportContinuation(plan) => {
                plan.diagnostic_label()
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
                .endpoints()
                .bridge_loop_membership()
                .ok_or(BridgeSupportRegionPlanRejection::EndpointMembershipUnknown)?;
            trims.push(BridgeSupportTrimPlan {
                fillet_index: candidate.fillet_index(),
                cylinder_face: bridge_membership.cylinder_face(),
                support_loop: bridge_membership.support_loop(),
                support_half_edge: bridge_membership.support_half_edge(),
                outside_face: bridge_membership.outside_face(),
                outside_loop: bridge_membership.outside_loop(),
                outside_half_edge: bridge_membership.outside_half_edge(),
                start: candidate.endpoints().start(),
                end: candidate.endpoints().end(),
            });
        }
        let expected_delta = ExpectedTopologyDelta::Known {
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
                face_family_bridge_operands_known: true,
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

    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn trims(&self) -> &[BridgeSupportTrimPlan] {
        &self.trims
    }
}

impl OutsideOutsideSupportContinuationPlan {
    fn from_candidates(support_face: FaceHandle, candidates: &[&RegionOperationCandidate]) -> Self {
        let trims = candidates
            .iter()
            .filter_map(|candidate| {
                let membership = candidate.endpoints().outside_outside_membership()?;
                Some(OutsideOutsideSupportTrimPlan {
                    fillet_index: candidate.fillet_index(),
                    cylinder_face: membership.cylinder_face(),
                    support_face: membership.support_face(),
                    start_face: membership.start_face(),
                    start_loop: membership.start_loop(),
                    start_half_edge: membership.start_half_edge(),
                    end_face: membership.end_face(),
                    end_loop: membership.end_loop(),
                    end_half_edge: membership.end_half_edge(),
                    start: candidate.endpoints().start(),
                    end: candidate.endpoints().end(),
                })
            })
            .collect::<Vec<_>>();
        Self {
            support_face,
            trims,
            expected_delta: ExpectedTopologyDelta::Unknown {
                reason: "outside-outside-support-continuation-operator-required",
            },
            preconditions: OutsideOutsideSupportPreconditions {
                both_endpoints_outside_support_family: true,
                endpoint_membership_known: true,
                support_owner_is_unique: true,
                topology_operator_required: true,
            },
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "outside-outside-support-continuation support={:?} trims=[{}] delta={} preconditions={}",
            self.support_face,
            self.trims
                .iter()
                .map(OutsideOutsideSupportTrimPlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.expected_delta.diagnostic_label(),
            self.preconditions.diagnostic_label()
        )
    }

    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn trims(&self) -> &[OutsideOutsideSupportTrimPlan] {
        &self.trims
    }
}

impl OutsideOutsideSupportTrimPlan {
    pub(super) fn cylinder_face(&self) -> FaceHandle {
        self.cylinder_face
    }

    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn start_face(&self) -> FaceHandle {
        self.start_face
    }

    pub(super) fn start_loop(&self) -> LoopHandle {
        self.start_loop
    }

    pub(super) fn start_half_edge(&self) -> HalfEdgeHandle {
        self.start_half_edge
    }

    pub(super) fn end_face(&self) -> FaceHandle {
        self.end_face
    }

    pub(super) fn end_loop(&self) -> LoopHandle {
        self.end_loop
    }

    pub(super) fn end_half_edge(&self) -> HalfEdgeHandle {
        self.end_half_edge
    }

    pub(super) fn start(&self) -> VertexHandle {
        self.start
    }

    pub(super) fn end(&self) -> VertexHandle {
        self.end
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} support={:?} start={:?}/{:?}@{:?} end={:?}/{:?}@{:?} {:?}->{:?}",
            self.fillet_index,
            self.cylinder_face,
            self.support_face,
            self.start_face,
            self.start_loop,
            self.start_half_edge,
            self.end_face,
            self.end_loop,
            self.end_half_edge,
            self.start,
            self.end
        )
    }
}

impl BridgeSupportTrimPlan {
    pub(super) fn fillet_index(&self) -> usize {
        self.fillet_index
    }

    pub(super) fn cylinder_face(&self) -> FaceHandle {
        self.cylinder_face
    }

    pub(super) fn support_loop(&self) -> LoopHandle {
        self.support_loop
    }

    pub(super) fn support_half_edge(&self) -> HalfEdgeHandle {
        self.support_half_edge
    }

    pub(super) fn outside_face(&self) -> FaceHandle {
        self.outside_face
    }

    pub(super) fn outside_loop(&self) -> LoopHandle {
        self.outside_loop
    }

    pub(super) fn outside_half_edge(&self) -> HalfEdgeHandle {
        self.outside_half_edge
    }

    pub(super) fn start(&self) -> VertexHandle {
        self.start
    }

    pub(super) fn end(&self) -> VertexHandle {
        self.end
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} support_loop={:?}@{:?} outside={:?}/{:?}@{:?} {:?}->{:?}",
            self.fillet_index,
            self.cylinder_face,
            self.support_loop,
            self.support_half_edge,
            self.outside_face,
            self.outside_loop,
            self.outside_half_edge,
            self.start,
            self.end
        )
    }
}

impl ExpectedTopologyDelta {
    fn zero() -> Self {
        Self::Known {
            vertices: 0,
            edges: 0,
            faces: 0,
            loops: 0,
        }
    }

    fn plus(self, other: Self) -> Self {
        match (self, other) {
            (
                Self::Known {
                    vertices,
                    edges,
                    faces,
                    loops,
                },
                Self::Known {
                    vertices: other_vertices,
                    edges: other_edges,
                    faces: other_faces,
                    loops: other_loops,
                },
            ) => Self::Known {
                vertices: vertices + other_vertices,
                edges: edges + other_edges,
                faces: faces + other_faces,
                loops: loops + other_loops,
            },
            (Self::Unknown { reason }, _) | (_, Self::Unknown { reason }) => {
                Self::Unknown { reason }
            }
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
            Self::Known {
                vertices,
                edges,
                faces,
                loops,
            } => {
                format!("dV={vertices} dE={edges} dF={faces} dL={loops}")
            }
            Self::Unknown { reason } => format!("unknown reason={reason}"),
        }
    }
}

impl BridgeSupportRegionPreconditions {
    fn diagnostic_label(&self) -> String {
        format!(
            "one_endpoint_on_support_loop={} other_endpoint_known={} unique_support_owner={} face_family_bridge_operands_known={}",
            self.exactly_one_endpoint_on_support_loop,
            self.other_endpoint_has_known_membership,
            self.support_owner_is_unique,
            self.face_family_bridge_operands_known
        )
    }
}

impl OutsideOutsideSupportPreconditions {
    fn diagnostic_label(&self) -> String {
        format!(
            "both_endpoints_outside_support_family={} endpoint_membership_known={} unique_support_owner={} topology_operator_required={}",
            self.both_endpoints_outside_support_family,
            self.endpoint_membership_known,
            self.support_owner_is_unique,
            self.topology_operator_required
        )
    }
}

impl SideSupportRegionRefusal {
    pub(super) fn stage(&self) -> PlannerStage {
        self.stage
    }

    pub(super) fn reason(&self) -> SideSupportRegionRefusalReason {
        self.reason
    }

    fn diagnostic_label(&self) -> String {
        format!("refusal stage={:?} reason={:?}", self.stage, self.reason)
    }
}
