//! Execution preflight for cylinder-mediated support-continuation topology.

use popkern_core::euler::CylinderMediatedOutsideIntervalInteraction;
use popkern_core::euler::CylinderMediatedOutsideIntervalPartition;
use popkern_core::euler::CylinderMediatedSupportContinuationBatchPlan;
use popkern_core::topology::BRep;

use super::repeated_apex_side_support_outside::OutsideOutsideSupportPrimitivePlan;

#[derive(Clone, Debug)]
pub(super) struct CylinderMediatedSupportContinuationTopologyOperatorReport {
    requests: usize,
    expected_delta_label: String,
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
    let expected_delta_label = batch.expected_delta.diagnostic_label();
    let mut candidate = brep.clone();
    let Err(err) = candidate.cylinder_mediated_support_continuation_batch(&batch) else {
        return CylinderMediatedSupportContinuationTopologyOperatorReport::refused_with_delta(
            requests,
            expected_delta_label,
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator,
        );
    };
    CylinderMediatedSupportContinuationTopologyOperatorReport::refused_with_delta_detail(
        requests,
        expected_delta_label,
        CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator,
        format!("{} core_refusal={err}", batch_interaction_detail(&batch)),
    )
}

fn batch_interaction_detail(batch: &CylinderMediatedSupportContinuationBatchPlan) -> String {
    match &batch.outside_interval_interaction {
        CylinderMediatedOutsideIntervalInteraction::Disjoint => {
            "outside_interval_interaction=disjoint".to_string()
        }
        CylinderMediatedOutsideIntervalInteraction::Overlapping { shared_half_edges } => format!(
            "outside_interval_interaction=overlapping shared_outside_half_edges={} {}",
            shared_half_edges.len(),
            partition_overlap_detail(&batch.outside_interval_partition)
        ),
    }
}

fn partition_overlap_detail(partition: &CylinderMediatedOutsideIntervalPartition) -> String {
    let mut shared_segments = 0;
    let mut max_shared_segment_len = 0;
    for segment in partition
        .loops
        .iter()
        .flat_map(|loop_partition| loop_partition.segments.iter())
    {
        if segment.request_indices.len() < 2 {
            continue;
        }
        shared_segments += 1;
        max_shared_segment_len = max_shared_segment_len.max(segment.half_edges.len());
    }
    format!("shared_segments={shared_segments} max_shared_segment_len={max_shared_segment_len}")
}

impl CylinderMediatedSupportContinuationTopologyOperatorReport {
    fn refused(
        requests: usize,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
    ) -> Self {
        Self::refused_with_delta(requests, "unknown".to_string(), reason)
    }

    fn refused_with_delta(
        requests: usize,
        expected_delta_label: String,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
    ) -> Self {
        Self {
            requests,
            expected_delta_label,
            status: CylinderMediatedSupportContinuationTopologyOperatorStatus::Refused(reason),
        }
    }

    fn refused_with_detail(
        requests: usize,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
        detail: String,
    ) -> Self {
        Self::refused_with_delta_detail(requests, "unknown".to_string(), reason, detail)
    }

    fn refused_with_delta_detail(
        requests: usize,
        expected_delta_label: String,
        reason: CylinderMediatedSupportContinuationTopologyOperatorGap,
        detail: String,
    ) -> Self {
        Self {
            requests,
            expected_delta_label,
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
            "topology_operator_ready=false requests={} expected_delta={} reason={}",
            self.requests,
            self.expected_delta_label,
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
            CylinderMediatedSupportContinuationTopologyOperatorGap::MissingMutatingEulerOperator => {
                "missing-mutating-euler-operator"
            }
        }
    }
}
