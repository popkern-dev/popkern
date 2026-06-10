//! Composed KP20 support-cell boundary repair candidates.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{
    assigned_support_cell_slots, boundary_cluster, boundary_prefix, boundary_row_group,
    exact_boundary_row, full_boundary_trims, group_reconcile,
    reauthor_support_cell_shared_edge_curves, tangent_row, trim_gap,
};

pub(crate) fn reconcile_assigned_ranked_endpoint_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_then_worst_endpoint_cluster(brep, faces)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_row_group_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_two_row_groups_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 2)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_three_row_groups_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 3)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_row_group_tangent_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_ranked_prefix_endpoint_row_group_then_tangent(brep, faces)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_exact_boundary_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_row_group_then_exact_boundary(brep, faces)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_partial_boundary_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_row_group_then_partial_boundary(brep, faces)
}

pub(crate) fn reconcile_assigned_ranked_endpoint_exact_then_partial_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    reconcile_ranked_prefix_endpoint_row_group_then_exact_and_partial_boundary(brep, faces)
}

fn reconcile_ranked_prefix_then_worst_endpoint_cluster(
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

fn reconcile_ranked_prefix_endpoint_cluster_then_row_groups(
    brep: &mut BRep,
    faces: &[FaceHandle],
    max_row_group_passes: usize,
) -> Result<&'static str, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    boundary_prefix::reconcile_ranked_boundary_prefix(brep, &trims)?;
    let endpoint_status = boundary_cluster::reconcile_worst_endpoint_cluster(brep, &trims)?;
    let kept_pass_count = reconcile_row_group_passes(brep, &trims, max_row_group_passes)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell ranked-prefix endpoint row-group reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_row_group_status(endpoint_status, kept_pass_count))
}

fn reconcile_ranked_prefix_endpoint_row_group_then_tangent(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let status = reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)?;
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    let tangent_status =
        tangent_row::reconcile_worst_row_group_if_normal_improves(brep, faces, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell tangent-row reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_tangent_status(status, tangent_status))
}

fn reconcile_ranked_prefix_endpoint_row_group_then_exact_boundary(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let base_status = reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)?;
    let slots = assigned_support_cell_slots(brep, faces);
    let exact_status = exact_boundary_row::exact_full_span_shared_boundary_rows(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell exact-boundary composed reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_exact_boundary_status(base_status, exact_status))
}

fn reconcile_ranked_prefix_endpoint_row_group_then_partial_boundary(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let base_status = reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)?;
    let slots = assigned_support_cell_slots(brep, faces);
    let partial_status = exact_boundary_row::fit_partial_shared_boundary_segments(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell partial-boundary composed reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_partial_boundary_status(
        base_status,
        partial_status,
    ))
}

fn reconcile_ranked_prefix_endpoint_row_group_then_exact_and_partial_boundary(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let base_status = reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)?;
    let slots = assigned_support_cell_slots(brep, faces);
    let exact_status = exact_boundary_row::exact_full_span_shared_boundary_rows(brep, &slots)?;
    let partial_status = exact_boundary_row::fit_partial_shared_boundary_segments(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell exact-plus-partial boundary reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_exact_and_partial_boundary_status(
        base_status,
        exact_status,
        partial_status,
    ))
}

pub(crate) fn reconcile_assigned_ranked_endpoint_exact_partial_tangent_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let base_status = reconcile_ranked_prefix_endpoint_cluster_then_row_groups(brep, faces, 1)?;
    let slots = assigned_support_cell_slots(brep, faces);
    let exact_status = exact_boundary_row::exact_full_span_shared_boundary_rows(brep, &slots)?;
    let partial_status = exact_boundary_row::fit_partial_shared_boundary_segments(brep, &slots)?;
    let tangent_trims = full_boundary_trims(brep, &slots)?;
    let tangent_status =
        tangent_row::reconcile_segment_owned_if_normal_improves(brep, faces, &tangent_trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell exact-partial tangent reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(composed_exact_partial_tangent_status(
        base_status,
        exact_status,
        partial_status,
        tangent_status,
    ))
}

fn composed_endpoint_status(endpoint_status: &'static str) -> &'static str {
    match endpoint_status {
        "endpoint_cluster_endpoint_band_fit_kept" => "ranked_prefix_then_endpoint_band_fit_kept",
        "endpoint_cluster_adjacent_controls_kept" => "ranked_prefix_then_adjacent_controls_kept",
        "endpoint_cluster_corners_only_kept" => "ranked_prefix_then_endpoint_corners_kept",
        _ => "ranked_prefix_then_endpoint_cluster_declined",
    }
}

fn reconcile_row_group_passes(
    brep: &mut BRep,
    trims: &[super::FullBoundaryTrim],
    max_row_group_passes: usize,
) -> Result<usize, CoreError> {
    let mut kept_pass_count = 0;
    for _ in 0..max_row_group_passes {
        if reconcile_row_group_if_improving(brep, trims)? {
            kept_pass_count += 1;
        } else {
            break;
        }
    }
    Ok(kept_pass_count)
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

fn composed_row_group_status(
    endpoint_status: &'static str,
    kept_pass_count: usize,
) -> &'static str {
    match (endpoint_status, kept_pass_count) {
        ("endpoint_cluster_endpoint_band_fit_kept", 3..) => {
            "ranked_prefix_endpoint_band_then_three_row_groups_kept"
        }
        ("endpoint_cluster_endpoint_band_fit_kept", 2..) => {
            "ranked_prefix_endpoint_band_then_two_row_groups_kept"
        }
        ("endpoint_cluster_endpoint_band_fit_kept", 1) => {
            "ranked_prefix_endpoint_band_then_row_group_kept"
        }
        ("endpoint_cluster_endpoint_band_fit_kept", 0) => {
            "ranked_prefix_endpoint_band_then_row_group_declined"
        }
        (_, 3..) => "ranked_prefix_endpoint_cluster_then_three_row_groups_kept",
        (_, 2..) => "ranked_prefix_endpoint_cluster_then_two_row_groups_kept",
        (_, 1) => "ranked_prefix_endpoint_cluster_then_row_group_kept",
        _ => "ranked_prefix_endpoint_cluster_then_row_group_declined",
    }
}

fn composed_tangent_status(
    base_status: &'static str,
    tangent_status: Option<tangent_row::TangentRowUpdate>,
) -> String {
    match tangent_status {
        Some(update) if base_status.contains("row_group_kept") => format!(
            "ranked_prefix_endpoint_row_group_{}_tangent_kept;tangent_alignment={};tangent_damping={:.15e}",
            update.alignment_label(),
            update.alignment_label(),
            update.damping
        ),
        Some(update) => format!(
            "ranked_prefix_endpoint_{}_tangent_kept;tangent_alignment={};tangent_damping={:.15e}",
            update.alignment_label(),
            update.alignment_label(),
            update.damping
        ),
        None if base_status.contains("row_group_kept") => {
            "ranked_prefix_endpoint_row_group_tangent_declined;tangent_alignment=none;tangent_damping=none".to_string()
        }
        None => {
            "ranked_prefix_endpoint_tangent_declined;tangent_alignment=none;tangent_damping=none"
                .to_string()
        }
    }
}

fn composed_exact_boundary_status(
    base_status: &'static str,
    exact_status: &'static str,
) -> &'static str {
    match exact_status {
        "exact_boundary_full_span_rows_authored" if base_status.contains("row_group_kept") => {
            "ranked_prefix_endpoint_row_group_exact_boundary_rows_authored"
        }
        "exact_boundary_full_span_rows_authored" => {
            "ranked_prefix_endpoint_exact_boundary_rows_authored"
        }
        _ if base_status.contains("row_group_kept") => {
            "ranked_prefix_endpoint_row_group_exact_boundary_declined"
        }
        _ => "ranked_prefix_endpoint_exact_boundary_declined",
    }
}

fn composed_partial_boundary_status(
    base_status: &'static str,
    partial_status: &'static str,
) -> &'static str {
    if segment_topology_already_segmented(partial_status) && base_status.contains("row_group_kept")
    {
        "ranked_prefix_endpoint_row_group_segment_topology_already_segmented_bounded_partial_boundary_segments_fit"
    } else if segment_topology_already_segmented(partial_status) {
        "ranked_prefix_endpoint_segment_topology_already_segmented_bounded_partial_boundary_segments_fit"
    } else if partial_boundary_segments_authored(partial_status)
        && base_status.contains("row_group_kept")
    {
        "ranked_prefix_endpoint_row_group_partial_boundary_segments_fit"
    } else if partial_boundary_segments_authored(partial_status) {
        "ranked_prefix_endpoint_partial_boundary_segments_fit"
    } else if base_status.contains("row_group_kept") {
        "ranked_prefix_endpoint_row_group_partial_boundary_declined"
    } else {
        "ranked_prefix_endpoint_partial_boundary_declined"
    }
}

fn composed_exact_and_partial_boundary_status(
    base_status: &'static str,
    exact_status: &'static str,
    partial_status: &'static str,
) -> &'static str {
    match (
        base_status.contains("row_group_kept"),
        exact_status,
        partial_status,
    ) {
        (true, "exact_boundary_full_span_rows_authored", partial_status)
            if segment_topology_already_segmented(partial_status) =>
        {
            "ranked_prefix_endpoint_row_group_segment_topology_already_segmented_exact_and_partial_boundary_authored"
        }
        (true, "exact_boundary_full_span_rows_authored", partial_status)
            if partial_boundary_segments_authored(partial_status) =>
        {
            "ranked_prefix_endpoint_row_group_exact_and_partial_boundary_authored"
        }
        (true, "exact_boundary_full_span_rows_authored", _) => {
            "ranked_prefix_endpoint_row_group_exact_boundary_rows_authored"
        }
        (true, _, partial_status) if partial_boundary_segments_authored(partial_status) => {
            "ranked_prefix_endpoint_row_group_partial_boundary_segments_fit"
        }
        _ if base_status.contains("row_group_kept") => {
            "ranked_prefix_endpoint_row_group_exact_and_partial_boundary_declined"
        }
        _ => "ranked_prefix_endpoint_exact_and_partial_boundary_declined",
    }
}

fn composed_exact_partial_tangent_status(
    base_status: &'static str,
    exact_status: &'static str,
    partial_status: &'static str,
    tangent_status: Option<tangent_row::TangentRowUpdate>,
) -> String {
    match (
        base_status.contains("row_group_kept"),
        exact_status,
        segment_topology_already_segmented(partial_status),
        partial_boundary_segments_authored(partial_status),
        tangent_status,
    ) {
        (
            true,
            "exact_boundary_full_span_rows_authored",
            true,
            _,
            Some(update),
        ) => {
            format!(
                "ranked_prefix_endpoint_row_group_segment_topology_already_segmented_exact_partial_{}_tangent_kept;tangent_alignment={};tangent_damping={:.15e}",
                update.alignment_label(),
                update.alignment_label(),
                update.damping
            )
        }
        (true, "exact_boundary_full_span_rows_authored", true, _, None) => {
            "ranked_prefix_endpoint_row_group_segment_topology_already_segmented_exact_and_partial_boundary_authored_tangent_declined;tangent_alignment=none;tangent_damping=none".to_string()
        }
        (true, "exact_boundary_full_span_rows_authored", _, true, Some(update)) => {
            format!(
                "ranked_prefix_endpoint_row_group_exact_partial_{}_tangent_kept;tangent_alignment={};tangent_damping={:.15e}",
                update.alignment_label(),
                update.alignment_label(),
                update.damping
            )
        }
        (true, "exact_boundary_full_span_rows_authored", _, true, None) => {
            "ranked_prefix_endpoint_row_group_exact_and_partial_boundary_authored_tangent_declined;tangent_alignment=none;tangent_damping=none".to_string()
        }
        _ => composed_exact_and_partial_boundary_status(base_status, exact_status, partial_status)
            .to_string(),
    }
}

fn segment_topology_already_segmented(status: &str) -> bool {
    status.contains("segment_topology_already_segmented")
}

fn partial_boundary_segments_authored(status: &str) -> bool {
    matches!(
        status,
        "partial_boundary_segments_fit" | "bounded_partial_boundary_segments_fit"
    ) || status.contains("then_bounded_partial_boundary_segments_fit")
}
