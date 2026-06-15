//! Constrained normal solve composed on exact/partial boundary ownership.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    assigned_support_cell_slots, boundary_composed, full_boundary_trims, normal_solve_status,
    reauthor_support_cell_shared_edge_curves, top_g1_normal_solve,
};

pub(crate) fn reconcile_assigned_exact_partial_top_g1_normal_solve_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let base_status =
        boundary_composed::reconcile_assigned_ranked_endpoint_exact_then_partial_status(
            brep, faces,
        )?;
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    let solve_report =
        top_g1_normal_solve::reconcile_top_g1_normal_samples_if_improves(brep, faces, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell exact/partial constrained normal solve produced invalid BRep: {error}"
        ))
    })?;
    Ok(format!(
        "exact_partial_top_g1_normal_solve;base_status={base_status};{}",
        normal_solve_status::worst_g1_normal_solve_status(&solve_report)
    ))
}
