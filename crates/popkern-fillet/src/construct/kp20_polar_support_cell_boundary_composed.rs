//! Composed KP20 support-cell boundary repair candidates.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    assigned_support_cell_slots, boundary_cluster, boundary_prefix, boundary_row_group,
    full_boundary_trims, group_reconcile, reauthor_support_cell_shared_edge_curves, trim_gap,
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

pub(super) fn reconcile_ranked_prefix_endpoint_cluster_then_row_group(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    boundary_prefix::reconcile_ranked_boundary_prefix(brep, &trims)?;
    let endpoint_status = boundary_cluster::reconcile_worst_endpoint_cluster(brep, &trims)?;
    let row_status = reconcile_row_group_if_improving(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell ranked-prefix endpoint row-group reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_row_group_status(endpoint_status, row_status))
}

fn composed_endpoint_status(endpoint_status: &'static str) -> &'static str {
    match endpoint_status {
        "endpoint_cluster_endpoint_band_fit_kept" => "ranked_prefix_then_endpoint_band_fit_kept",
        "endpoint_cluster_adjacent_controls_kept" => "ranked_prefix_then_adjacent_controls_kept",
        "endpoint_cluster_corners_only_kept" => "ranked_prefix_then_endpoint_corners_kept",
        _ => "ranked_prefix_then_endpoint_cluster_declined",
    }
}

fn reconcile_row_group_if_improving(
    brep: &mut BRep,
    trims: &[super::FullBoundaryTrim],
) -> Result<bool, CoreError> {
    let row_group = boundary_row_group::worst_boundary_row_group_trims(brep, trims)?;
    let before_gap = worst_gap(brep, trims)?;
    let before = brep.clone();
    group_reconcile::reconcile_grouped_boundary_rows(brep, &row_group)?;
    if worst_gap(brep, trims)? + 1.0e-12 < before_gap {
        Ok(true)
    } else {
        *brep = before;
        Ok(false)
    }
}

fn worst_gap(brep: &BRep, trims: &[super::FullBoundaryTrim]) -> Result<f64, CoreError> {
    let mut worst = 0.0_f64;
    for trim in trims {
        let patch_a = super::rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = super::rectangular_binonic_on_face(brep, trim.face_b)?;
        worst = worst.max(trim_gap(&patch_a, &patch_b, *trim));
    }
    Ok(worst)
}

fn composed_row_group_status(endpoint_status: &'static str, row_kept: bool) -> &'static str {
    match (endpoint_status, row_kept) {
        ("endpoint_cluster_endpoint_band_fit_kept", true) => {
            "ranked_prefix_endpoint_band_then_row_group_kept"
        }
        ("endpoint_cluster_endpoint_band_fit_kept", false) => {
            "ranked_prefix_endpoint_band_then_row_group_declined"
        }
        (_, true) => "ranked_prefix_endpoint_cluster_then_row_group_kept",
        _ => "ranked_prefix_endpoint_cluster_then_row_group_declined",
    }
}
