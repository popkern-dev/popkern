//! Execution preflight for cylinder-mediated support-continuation topology.

use popkern_core::euler::CylinderMediatedOutsideIntervalInteraction;
use popkern_core::topology::BRep;

use super::repeated_apex_side_support_outside::OutsideOutsideSupportPrimitivePlan;

#[derive(Clone, Debug)]
pub(super) struct CylinderMediatedSupportContinuationTopologyOperatorReport {
    requests: usize,
    status: CylinderMediatedSupportContinuationTopologyOperatorStatus,
}

#[derive(Clone, Debug)]
enum CylinderMediatedSupportContinuationTopologyOperatorStatus {
    Refused(CylinderMediatedSupportContinuationTopologyOperatorGap),
    RefusedWithDetail(
        CylinderMediatedSupportContinuationTopologyOperatorGap,
        String,
    ),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CylinderMediatedSupportContinuationTopologyOperatorGap {
    MissingPrimitivePlan,
    EmptyPrimitivePlan,
    CoreBatchPreflightRejected,
    OverlappingOutsideIntervals,
    MissingMutatingEulerOperator,
}

pub(super) fn topology_operator_report(
    primitive: Option<&OutsideOutsideSupportPrimitivePlan>,
    brep: &BRep,
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
    let batch = match brep
        .preflight_cylinder_mediated_support_continuation_batch(primitive.requests().to_vec())
    {
        Ok(batch) => batch,
        Err(err) => {
            return CylinderMediatedSupportContinuationTopologyOperatorReport::refused_with_detail(
                requests,
                CylinderMediatedSupportContinuationTopologyOperatorGap::CoreBatchPreflightRejected,
                err.to_string(),
            );
        }
    };
    if let CylinderMediatedOutsideIntervalInteraction::Overlapping { shared_half_edges } =
        batch.outside_interval_interaction
    {
        return CylinderMediatedSupportContinuationTopologyOperatorReport::refused_with_detail(
            requests,
            CylinderMediatedSupportContinuationTopologyOperatorGap::OverlappingOutsideIntervals,
            format!("shared_outside_half_edges={}", shared_half_edges.len()),
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

    fn refused_with_detail(
        requests: usize,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
        detail: String,
    ) -> Self {
        Self {
            requests,
            status: CylinderMediatedSupportContinuationTopologyOperatorStatus::RefusedWithDetail(
                reason, detail,
            ),
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
    fn diagnostic_label(&self) -> String {
        match self {
            CylinderMediatedSupportContinuationTopologyOperatorStatus::Refused(reason) => {
                reason.diagnostic_label().to_string()
            }
            CylinderMediatedSupportContinuationTopologyOperatorStatus::RefusedWithDetail(
                reason,
                detail,
            ) => {
                format!("{} detail={detail}", reason.diagnostic_label())
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
            CylinderMediatedSupportContinuationTopologyOperatorGap::CoreBatchPreflightRejected => {
                "core-batch-preflight-rejected"
            }
            CylinderMediatedSupportContinuationTopologyOperatorGap::OverlappingOutsideIntervals => {
                "overlapping-outside-intervals"
            }
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator => {
                "missing-mutating-euler-operator"
            }
        }
    }
}
