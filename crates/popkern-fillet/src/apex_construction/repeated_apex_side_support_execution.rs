//! Execution boundary for selected repeated-apex side-support plans.

use popkern_core::topology::BRep;

use super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitiveDerivation;
use super::repeated_apex_side_support_plan::{
    SideSupportRegionExecutionReadiness, SideSupportRegionPlanSelection,
    SideSupportRegionRefusalReason,
};

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionExecutionReport {
    status: SideSupportRegionExecutionStatus,
}

#[derive(Clone, Debug)]
enum SideSupportRegionExecutionStatus {
    ReadyForExecutor,
    Refused(SideSupportRegionExecutionRefusal),
}

#[derive(Clone, Debug)]
struct SideSupportRegionExecutionRefusal {
    stage: ExecutionStage,
    reason: SideSupportRegionExecutionRefusalReason,
    detail: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ExecutionStage {
    Preflight,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SideSupportRegionExecutionRefusalReason {
    StaleKnowledgeGraph,
    PlannerRefused,
    MissingLiveTopologyOperand,
    MissingValidatedFaceFamilyBridgePrimitive,
    MissingFaceFamilyBridgeTopologyOperator,
}

impl SideSupportRegionExecutionReport {
    pub(super) fn preflight(
        selection: &SideSupportRegionPlanSelection,
        brep: &BRep,
        knowledge_is_fresh: bool,
    ) -> Self {
        if !knowledge_is_fresh {
            return Self::refused(SideSupportRegionExecutionRefusalReason::StaleKnowledgeGraph);
        }
        if let Err(failure) = selection.validate_live_operands(brep) {
            return Self::refused_with_detail(
                SideSupportRegionExecutionRefusalReason::MissingLiveTopologyOperand,
                failure.diagnostic_label(),
            );
        }
        match selection.execution_readiness() {
            SideSupportRegionExecutionReadiness::Ready => Self {
                status: SideSupportRegionExecutionStatus::ReadyForExecutor,
            },
            SideSupportRegionExecutionReadiness::Blocked { planner_reason, .. } => {
                let reason = match planner_reason {
                    SideSupportRegionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive => {
                        let bridge_primitive =
                            FaceFamilyBridgePrimitiveDerivation::build(selection, brep);
                        let reason = if bridge_primitive.is_complete() {
                            SideSupportRegionExecutionRefusalReason::MissingFaceFamilyBridgeTopologyOperator
                        } else {
                            SideSupportRegionExecutionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive
                        };
                        return Self::refused_with_detail(
                            reason,
                            format!(
                                "{}; {}",
                                bridge_primitive.diagnostic_label(),
                                bridge_primitive
                                    .topology_operator_report(brep)
                                    .diagnostic_label()
                            ),
                        );
                    }
                    SideSupportRegionRefusalReason::UnknownEndpointMembership
                    | SideSupportRegionRefusalReason::AmbiguousEndpointMembership
                    | SideSupportRegionRefusalReason::UnsupportedBridgeTrimCount
                    | SideSupportRegionRefusalReason::NoAcceptedOperationCandidate => {
                        SideSupportRegionExecutionRefusalReason::PlannerRefused
                    }
                };
                Self::refused(reason)
            }
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match &self.status {
            SideSupportRegionExecutionStatus::ReadyForExecutor => {
                "execution ready_for_executor".to_string()
            }
            SideSupportRegionExecutionStatus::Refused(refusal) => {
                let mut label = format!(
                    "execution refused stage={:?} reason={:?}",
                    refusal.stage, refusal.reason
                );
                if let Some(detail) = &refusal.detail {
                    label.push_str(&format!(" detail={detail}"));
                }
                label
            }
        }
    }

    fn refused(reason: SideSupportRegionExecutionRefusalReason) -> Self {
        Self {
            status: SideSupportRegionExecutionStatus::Refused(SideSupportRegionExecutionRefusal {
                stage: ExecutionStage::Preflight,
                reason,
                detail: None,
            }),
        }
    }

    fn refused_with_detail(
        reason: SideSupportRegionExecutionRefusalReason,
        detail: String,
    ) -> Self {
        Self {
            status: SideSupportRegionExecutionStatus::Refused(SideSupportRegionExecutionRefusal {
                stage: ExecutionStage::Preflight,
                reason,
                detail: Some(detail),
            }),
        }
    }
}
