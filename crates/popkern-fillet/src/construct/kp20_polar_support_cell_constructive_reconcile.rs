//! Two-pass constructive support-cell reconciliation wrappers.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

use super::{
    assigned_support_cell_slots, chain_fit, full_boundary_trims,
    reauthor_support_cell_shared_edge_curves, tangent_row, worst_g1_normal_solve,
};

type ChainConstructiveAuthor =
    fn(&mut BRep, &[Kp20ExpandedSupportSurfaceSlot]) -> Result<String, CoreError>;

pub(crate) fn reconcile_assigned_chain_constructive_cells_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_assigned_chain_two_pass_status(
        brep,
        faces,
        chain_fit::author_chain_network_with_constructive_cells,
        "kp20 support-cell constructive cell re-authoring produced invalid BRep",
        "max_twist_conflict_distance=",
        "first_pass_max_twist_conflict_distance",
    )
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_cells_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_assigned_chain_two_pass_status(
        brep,
        faces,
        chain_fit::author_chain_network_with_constructive_dual_twist_cells,
        "kp20 support-cell dual-twist constructive re-authoring produced invalid BRep",
        "max_dual_twist_pair_separation=",
        "first_pass_max_dual_twist_pair_separation",
    )
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_beta_envelope_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_assigned_chain_two_pass_status(
        brep,
        faces,
        chain_fit::author_chain_network_with_constructive_dual_twist_beta_envelope,
        "kp20 support-cell dual-twist beta-envelope solve produced invalid BRep",
        "max_beta_envelope_control_delta=",
        "first_pass_max_beta_envelope_control_delta",
    )
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_pair_normal_solve_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_assigned_chain_two_pass_status(
        brep,
        faces,
        chain_fit::author_chain_network_with_constructive_dual_twist_pair_normal_solve,
        "kp20 support-cell dual-twist pair-normal solve produced invalid BRep",
        "max_pair_normal_control_delta=",
        "first_pass_max_pair_normal_control_delta",
    )
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_then_worst_g1_tangent_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let first = chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    let status = chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    let trims = full_boundary_trims(brep, &slots)?;
    let tangent_status =
        tangent_row::reconcile_worst_g1_trim_if_normal_improves(brep, faces, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell dual-twist worst-G1 tangent solve produced invalid BRep: {error}"
        ))
    })?;
    let first_value = first
        .split("max_dual_twist_pair_separation=")
        .nth(1)
        .and_then(|tail| tail.split(';').next())
        .unwrap_or("unknown");
    Ok(format!(
        "{status};constructive_pass_count=2;first_pass_max_dual_twist_pair_separation={first_value};{}",
        worst_g1_tangent_status(tangent_status)
    ))
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_then_worst_g1_normal_solve_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let first = chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    let status = chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    let trims = full_boundary_trims(brep, &slots)?;
    let solve_report =
        worst_g1_normal_solve::reconcile_worst_g1_normal_sample_if_improves(brep, faces, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell dual-twist worst-G1 normal solve produced invalid BRep: {error}"
        ))
    })?;
    let first_value = first
        .split("max_dual_twist_pair_separation=")
        .nth(1)
        .and_then(|tail| tail.split(';').next())
        .unwrap_or("unknown");
    Ok(format!(
        "{status};constructive_pass_count=2;first_pass_max_dual_twist_pair_separation={first_value};{}",
        worst_g1_normal_solve_status(&solve_report)
    ))
}

pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_joint_leg_solve_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    reconcile_assigned_chain_two_pass_status(
        brep,
        faces,
        chain_fit::author_chain_network_with_constructive_dual_twist_joint_leg_solve,
        "kp20 support-cell dual-twist joint leg solve produced invalid BRep",
        "max_joint_endpoint_residual_after=",
        "first_pass_max_joint_endpoint_residual_after",
    )
}

fn reconcile_assigned_chain_two_pass_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
    author: ChainConstructiveAuthor,
    invalid_message: &'static str,
    first_metric_key: &'static str,
    first_metric_status_key: &'static str,
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let first = author(brep, &slots)?;
    let status = author(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate()
        .map_err(|error| CoreError::PreconditionFailed(format!("{invalid_message}: {error}")))?;
    let first_value = first
        .split(first_metric_key)
        .nth(1)
        .and_then(|tail| tail.split(';').next())
        .unwrap_or("unknown");
    Ok(format!(
        "{status};constructive_pass_count=2;{first_metric_status_key}={first_value}"
    ))
}

fn worst_g1_tangent_status(tangent_status: Option<tangent_row::TangentRowUpdate>) -> String {
    match tangent_status {
        Some(update) => format!(
            "worst_g1_tangent_status=kept;tangent_alignment={};tangent_damping={:.15e}",
            update.alignment_label(),
            update.damping
        ),
        None => "worst_g1_tangent_status=declined;tangent_alignment=none;tangent_damping=none"
            .to_string(),
    }
}

fn worst_g1_normal_solve_status(
    report: &worst_g1_normal_solve::WorstG1NormalSolveReport,
) -> String {
    format!(
        "worst_g1_normal_solve_status={};worst_g1_normal_sample_count={};worst_g1_normal_variable_count={};worst_g1_normal_limited_control_count={};max_worst_g1_normal_control_delta={:.6e};worst_g1_normal_damping={};worst_g1_before={:.15e};worst_g1_after={:.15e};local_worst_g1_before={:.15e};local_worst_g1_after={:.15e}",
        report.status,
        report.sample_count,
        report.variable_count,
        report.limited_control_count,
        report.max_control_delta,
        report
            .damping
            .map(|damping| format!("{damping:.15e}"))
            .unwrap_or_else(|| "none".to_string()),
        report.before_worst_g1,
        report.after_worst_g1,
        report.local_disagreement_before,
        report.local_disagreement_after
    )
}
