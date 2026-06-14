//! Two-pass constructive support-cell reconciliation wrappers.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

use super::{assigned_support_cell_slots, chain_fit, reauthor_support_cell_shared_edge_curves};

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
