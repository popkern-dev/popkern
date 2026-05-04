//! Execution boundary for selected repeated-apex side-support plans.

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
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ExecutionStage {
    Preflight,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SideSupportRegionExecutionRefusalReason {
    StaleKnowledgeGraph,
    PlannerRefused,
    MissingValidatedBridgeExecutor,
}

impl SideSupportRegionExecutionReport {
    pub(super) fn preflight(
        selection: &SideSupportRegionPlanSelection,
        knowledge_is_fresh: bool,
    ) -> Self {
        if !knowledge_is_fresh {
            return Self::refused(SideSupportRegionExecutionRefusalReason::StaleKnowledgeGraph);
        }
        match selection.execution_readiness() {
            SideSupportRegionExecutionReadiness::Ready => Self {
                status: SideSupportRegionExecutionStatus::ReadyForExecutor,
            },
            SideSupportRegionExecutionReadiness::Blocked { planner_reason, .. } => {
                let reason = match planner_reason {
                    SideSupportRegionRefusalReason::MissingValidatedBridgeExecutor => {
                        SideSupportRegionExecutionRefusalReason::MissingValidatedBridgeExecutor
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
                format!(
                    "execution refused stage={:?} reason={:?}",
                    refusal.stage, refusal.reason
                )
            }
        }
    }

    fn refused(reason: SideSupportRegionExecutionRefusalReason) -> Self {
        Self {
            status: SideSupportRegionExecutionStatus::Refused(SideSupportRegionExecutionRefusal {
                stage: ExecutionStage::Preflight,
                reason,
            }),
        }
    }
}
