//! Execution boundary for selected repeated-apex side-support plans.

use popkern_core::topology::BRep;

use super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitiveDerivation;
use super::repeated_apex_side_support_merge::SupportRegionMergePrimitiveDerivation;
use super::repeated_apex_side_support_outside::OutsideOutsideSupportPrimitiveDerivation;
use super::repeated_apex_side_support_plan::{
    SideSupportRegionExecutionReadiness, SideSupportRegionPlanSelection,
    SideSupportRegionRefusalReason,
};
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionExecutionReport {
    status: SideSupportRegionExecutionStatus,
}

#[derive(Clone, Debug)]
enum SideSupportRegionExecutionStatus {
    ReadyForExecutor { detail: Option<String> },
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
    MissingValidatedOutsideOutsideSupportContinuationPrimitive,
    MissingFaceFamilyBridgeTopologyOperator,
    MissingCylinderMediatedSupportContinuationTopologyOperator,
    MissingValidatedSupportRegionMergePrimitive,
    MissingExactSideGeneratorCurveAuthoring,
}

impl SideSupportRegionExecutionReport {
    pub(super) fn preflight(
        selection: &SideSupportRegionPlanSelection,
        knowledge: &TopologyKnowledgeGraph,
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
                status: SideSupportRegionExecutionStatus::ReadyForExecutor { detail: None },
            },
            SideSupportRegionExecutionReadiness::Blocked { planner_reason, .. } => {
                let reason = match planner_reason {
                    SideSupportRegionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive => {
                        let bridge_primitive =
                            FaceFamilyBridgePrimitiveDerivation::build(selection, brep);
                        let topology_report = bridge_primitive.topology_operator_report(brep);
                        if bridge_primitive.is_complete() && topology_report.is_ready() {
                            let mut candidate = brep.clone();
                            if let Err(err) =
                                bridge_primitive.apply_topology_operator(&mut candidate)
                            {
                                return Self::refused_with_detail(
                                    SideSupportRegionExecutionRefusalReason::MissingFaceFamilyBridgeTopologyOperator,
                                    err,
                                );
                            }
                            return Self {
                                status: SideSupportRegionExecutionStatus::ReadyForExecutor {
                                    detail: Some(format!(
                                        "{}; {}",
                                        bridge_primitive.diagnostic_label(),
                                        topology_report.diagnostic_label()
                                    )),
                                },
                            };
                        }
                        let merge_primitive =
                            SupportRegionMergePrimitiveDerivation::build(selection, brep);
                        if merge_primitive.exact_geometry_ready() {
                            let mut candidate = brep.clone();
                            if let Err(err) =
                                merge_primitive.apply_topology_and_geometry(&mut candidate)
                            {
                                return Self::refused_with_detail(
                                    SideSupportRegionExecutionRefusalReason::MissingExactSideGeneratorCurveAuthoring,
                                    format!(
                                        "{}; {}; {}; apply_refusal={}",
                                        bridge_primitive.diagnostic_label(),
                                        topology_report.diagnostic_label(),
                                        merge_primitive.diagnostic_label(),
                                        err
                                    ),
                                );
                            }
                            return Self {
                                status: SideSupportRegionExecutionStatus::ReadyForExecutor {
                                    detail: Some(format!(
                                        "{}; {}; {}",
                                        bridge_primitive.diagnostic_label(),
                                        topology_report.diagnostic_label(),
                                        merge_primitive.diagnostic_label()
                                    )),
                                },
                            };
                        }
                        if merge_primitive.is_ready() {
                            return Self::refused_with_detail(
                                SideSupportRegionExecutionRefusalReason::MissingExactSideGeneratorCurveAuthoring,
                                format!(
                                    "{}; {}; {}",
                                    bridge_primitive.diagnostic_label(),
                                    topology_report.diagnostic_label(),
                                    merge_primitive.diagnostic_label()
                                ),
                            );
                        }
                        let reason = if bridge_primitive.is_complete() {
                            SideSupportRegionExecutionRefusalReason::MissingFaceFamilyBridgeTopologyOperator
                        } else if !merge_primitive.has_regions() {
                            SideSupportRegionExecutionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive
                        } else {
                            SideSupportRegionExecutionRefusalReason::MissingValidatedSupportRegionMergePrimitive
                        };
                        return Self::refused_with_detail(
                            reason,
                            format!(
                                "{}; {}; {}",
                                bridge_primitive.diagnostic_label(),
                                topology_report.diagnostic_label(),
                                merge_primitive.diagnostic_label()
                            ),
                        );
                    }
                    SideSupportRegionRefusalReason::MissingValidatedOutsideOutsideSupportContinuationPrimitive => {
                        let outside_primitive =
                            OutsideOutsideSupportPrimitiveDerivation::build(
                                selection, knowledge, brep,
                            );
                        let topology_report = outside_primitive.topology_operator_report(brep);
                        let reason = if outside_primitive.is_complete()
                            && !topology_report.is_ready()
                        {
                            SideSupportRegionExecutionRefusalReason::MissingCylinderMediatedSupportContinuationTopologyOperator
                        } else {
                            SideSupportRegionExecutionRefusalReason::MissingValidatedOutsideOutsideSupportContinuationPrimitive
                        };
                        return Self::refused_with_detail(
                            reason,
                            format!(
                                "{}; {}; {}",
                                selection.diagnostic_label(),
                                outside_primitive.diagnostic_label(),
                                topology_report.diagnostic_label()
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
            SideSupportRegionExecutionStatus::ReadyForExecutor { detail } => {
                let mut label = "execution ready_for_executor".to_string();
                if let Some(detail) = detail {
                    label.push_str(&format!(" detail={detail}"));
                }
                label
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
