//! Candidate selection for repeated-apex partition execution.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::partition_plan::RepeatedApexPartitionPlan;
use super::{ApexBuildError, plan_repeated_apex_loop_normalization};

pub(super) fn try_plan_repeated_apex_partition(
    brep: &BRep,
    cap_face: FaceHandle,
    apex_vertex: VertexHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[Option<EdgeHandle>],
    n_fillets: usize,
) -> Result<Option<RepeatedApexPartitionPlan>, ApexBuildError> {
    let Some(loop_plan) = plan_repeated_apex_loop_normalization(
        brep,
        cap_face,
        apex_vertex,
        apex_vertex_handles,
        n_fillets,
    ) else {
        return Ok(None);
    };
    if !loop_plan.is_commercial_grade_candidate(n_fillets) {
        return Ok(None);
    }

    RepeatedApexPartitionPlan::build(
        brep,
        cap_face,
        &loop_plan,
        apex_vertex_handles,
        face_edge_handles,
        n_fillets,
    )
    .map(Some)
}
