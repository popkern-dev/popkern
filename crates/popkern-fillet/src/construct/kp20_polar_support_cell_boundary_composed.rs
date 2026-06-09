//! Composed KP20 support-cell boundary repair candidates.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    assigned_support_cell_slots, boundary_cluster, boundary_prefix, full_boundary_trims,
    reauthor_support_cell_shared_edge_curves,
};

pub(super) fn reconcile_ranked_prefix_then_worst_endpoint_cluster(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    boundary_prefix::reconcile_ranked_boundary_prefix(brep, &trims)?;
    let endpoint_status = boundary_cluster::reconcile_worst_endpoint_cluster(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell ranked-prefix endpoint-cluster reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_endpoint_status(endpoint_status))
}

fn composed_endpoint_status(endpoint_status: &'static str) -> &'static str {
    match endpoint_status {
        "endpoint_cluster_endpoint_band_fit_kept" => "ranked_prefix_then_endpoint_band_fit_kept",
        "endpoint_cluster_adjacent_controls_kept" => "ranked_prefix_then_adjacent_controls_kept",
        "endpoint_cluster_corners_only_kept" => "ranked_prefix_then_endpoint_corners_kept",
        _ => "ranked_prefix_then_endpoint_cluster_declined",
    }
}
