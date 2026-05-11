//! Preconditions for cylinder-mediated support continuation application.

use crate::error::CoreError;

use super::super::cylinder_mediated_support::CylinderMediatedSupportContinuationAggregateDelta;
use super::super::cylinder_mediated_support_partition::CylinderMediatedOutsideIntervalPartition;

pub(super) fn require_known_aggregate_delta(
    delta: &CylinderMediatedSupportContinuationAggregateDelta,
) -> Result<(), CoreError> {
    match delta {
        CylinderMediatedSupportContinuationAggregateDelta::Known { .. } => Ok(()),
        CylinderMediatedSupportContinuationAggregateDelta::Unknown { reason } => {
            Err(CoreError::PreconditionFailed(format!(
                "cylinder-mediated support continuation aggregate delta is unknown: {reason}"
            )))
        }
    }
}

pub(super) fn require_non_empty_partition(
    partition: &CylinderMediatedOutsideIntervalPartition,
) -> Result<(), CoreError> {
    if partition.loops.is_empty()
        || partition
            .loops
            .iter()
            .all(|loop_partition| loop_partition.segments.is_empty())
    {
        return Err(CoreError::PreconditionFailed(
            "cylinder-mediated support continuation requires a non-empty outside interval partition"
                .into(),
        ));
    }
    Ok(())
}
