//! Execution preflight for cylinder-mediated support-continuation topology.

use super::repeated_apex_side_support_outside::OutsideOutsideSupportPrimitivePlan;

#[derive(Clone, Debug)]
pub(super) struct CylinderMediatedSupportContinuationTopologyOperatorReport {
    requests: usize,
    status: CylinderMediatedSupportContinuationTopologyOperatorStatus,
}

#[derive(Clone, Debug)]
enum CylinderMediatedSupportContinuationTopologyOperatorStatus {
    Refused(CylinderMediatedSupportContinuationTopologyOperatorGap),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CylinderMediatedSupportContinuationTopologyOperatorGap {
    MissingPrimitivePlan,
    EmptyPrimitivePlan,
    MissingMutatingEulerOperator,
}

pub(super) fn topology_operator_report(
    primitive: Option<&OutsideOutsideSupportPrimitivePlan>,
) -> CylinderMediatedSupportContinuationTopologyOperatorReport {
    let Some(primitive) = primitive else {
        return CylinderMediatedSupportContinuationTopologyOperatorReport::refused(
            0,
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingPrimitivePlan,
        );
    };
    let requests = primitive.requests().len();
    if requests == 0 {
        return CylinderMediatedSupportContinuationTopologyOperatorReport::refused(
            requests,
            CylinderMediatedSupportContinuationTopologyOperatorGap::EmptyPrimitivePlan,
        );
    }
    CylinderMediatedSupportContinuationTopologyOperatorReport::refused(
        requests,
        CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator,
    )
}

impl CylinderMediatedSupportContinuationTopologyOperatorReport {
    fn refused(
        requests: usize,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
    ) -> Self {
        Self {
            requests,
            status: CylinderMediatedSupportContinuationTopologyOperatorStatus::Refused(reason),
        }
    }

    pub(super) fn is_ready(&self) -> bool {
        false
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "topology_operator_ready=false requests={} expected_delta=unknown reason={}",
            self.requests,
            self.status.diagnostic_label()
        )
    }
}

impl CylinderMediatedSupportContinuationTopologyOperatorStatus {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            CylinderMediatedSupportContinuationTopologyOperatorStatus::Refused(reason) => {
                reason.diagnostic_label()
            }
        }
    }
}

impl CylinderMediatedSupportContinuationTopologyOperatorGap {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingPrimitivePlan => {
                "missing-primitive-plan"
            }
            CylinderMediatedSupportContinuationTopologyOperatorGap::EmptyPrimitivePlan => {
                "empty-primitive-plan"
            }
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator => {
                "missing-mutating-euler-operator"
            }
        }
    }
}
