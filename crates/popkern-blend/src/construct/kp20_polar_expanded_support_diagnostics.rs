//! Diagnostics for expanded KP20 support topology.

use popkern_core::geometry::gregory::GregoryPatch;

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_expanded_support_topology::Kp20ExpandedSupportTopologyReport;
use crate::types::{CornerSubPatchProvenance, CornerSubPatchRole};

#[derive(Debug, Clone, Default)]
pub(super) struct Kp20ExpandedSupportFoldLocalization {
    pub(super) macro_fragment_folds: usize,
    pub(super) support_boundary_lane_folds: usize,
    pub(super) support_interior_lane_folds: usize,
    pub(super) worst_fold_role: &'static str,
    pub(super) worst_fold_sector: Option<usize>,
    pub(super) worst_fold_patch: Option<usize>,
    pub(super) worst_fold_normal_cos: Option<f64>,
}

impl Kp20ExpandedSupportFoldLocalization {
    pub(super) fn total_folds(&self) -> usize {
        self.macro_fragment_folds
            + self.support_boundary_lane_folds
            + self.support_interior_lane_folds
    }

    pub(super) fn support_lane_folds(&self) -> usize {
        self.support_boundary_lane_folds + self.support_interior_lane_folds
    }
}

pub(super) fn localize_expanded_support_folds(
    patches: &[GregoryPatch],
    provenance: &[CornerSubPatchProvenance],
) -> Kp20ExpandedSupportFoldLocalization {
    let mut localization = Kp20ExpandedSupportFoldLocalization::default();
    for (patch_index, (patch, provenance)) in patches.iter().zip(provenance.iter()).enumerate() {
        let report = patch.self_intersection_report(SELF_INT_TOL);
        if !report.intersects() {
            continue;
        }
        match provenance.role {
            CornerSubPatchRole::SupportBoundaryLane => {
                localization.support_boundary_lane_folds += 1;
            }
            CornerSubPatchRole::SupportInteriorLane => {
                localization.support_interior_lane_folds += 1;
            }
            role if role.is_support_cell() => {
                localization.support_interior_lane_folds += 1;
            }
            _ => localization.macro_fragment_folds += 1,
        }
        if localization
            .worst_fold_normal_cos
            .is_none_or(|worst| report.min_jacobian_normal_cos < worst)
        {
            localization.worst_fold_role = role_label(provenance.role);
            localization.worst_fold_sector = Some(provenance.sector_index);
            localization.worst_fold_patch = Some(patch_index);
            localization.worst_fold_normal_cos = Some(report.min_jacobian_normal_cos);
        }
    }
    localization
}

pub(super) fn expanded_support_topology_detail(
    label: &str,
    report: &Kp20ExpandedSupportTopologyReport,
) -> String {
    format!(
        ";{label}_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_measured_folds={};{label}_retained_macro_folds={};{label}_expanded_patch_folds={};{label}_target_piece={:?};{label}_target_sector={};{label}_emitted_patch_count={};{label}_measured_patch_count={};{label}_topology_face_count={};{label}_topology_pentagon_count={};{label}_topology_alignment={};{label}_split_parameter={:.3};{label}_fragment_rebuild={};{label}_support_lane_reauthor={};{label}_support_lane_reauthor_scale={};{label}_support_cell_primitive_evaluated={};{label}_support_cell_primitive_target_count={};{label}_support_cell_primitive_patch_count={};{label}_support_cell_primitive_folded_patch_count={};{label}_support_cell_primitive_cells_u={};{label}_support_cell_primitive_cells_v={};{label}_support_cell_primitive_topology_path_g0={:.6e};{label}_support_cell_primitive_selection_policy={};{label}_support_cell_primitive_mode={};{label}_support_cell_primitive_target_piece={:?};{label}_support_cell_primitive_decline_reason={};{label}_macro_fragment_folds={};{label}_support_boundary_lane_folds={};{label}_support_interior_lane_folds={};{label}_worst_fold_role={};{label}_worst_fold_sector={};{label}_worst_fold_patch={};{label}_worst_fold_normal_cos={};{label}_decline_reason={}",
        report.evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.measured_folds,
        report.retained_macro_folds,
        report.expanded_patch_folds,
        report.target_piece,
        option_usize_label(report.target_sector),
        report.emitted_patch_count,
        report.measured_patch_count,
        report.topology_face_count,
        report.topology_pentagon_count,
        report.topology_alignment,
        report.split_parameter,
        report.fragment_rebuild,
        report.support_lane_reauthor,
        option_f64_label(report.support_lane_reauthor_scale),
        report.support_cell_primitive.evaluated,
        report.support_cell_primitive.target_count,
        report.support_cell_primitive.patch_count,
        report.support_cell_primitive.folded_patch_count,
        report.support_cell_primitive.cells_u,
        report.support_cell_primitive.cells_v,
        report.support_cell_primitive.topology_path_g0,
        report.support_cell_primitive.selection_policy,
        report.support_cell_primitive.construction_mode,
        report.support_cell_primitive.target_piece,
        report.support_cell_primitive.decline_reason,
        report.fold_localization.macro_fragment_folds,
        report.fold_localization.support_boundary_lane_folds,
        report.fold_localization.support_interior_lane_folds,
        report.fold_localization.worst_fold_role,
        option_usize_label(report.fold_localization.worst_fold_sector),
        option_usize_label(report.fold_localization.worst_fold_patch),
        report
            .fold_localization
            .worst_fold_normal_cos
            .map(|value| format!("{value:.3e}"))
            .unwrap_or_else(|| "none".to_string()),
        report.decline_reason,
    )
}

fn option_usize_label(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "none".to_string())
}

fn option_f64_label(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.3}"))
        .unwrap_or_else(|| "none".to_string())
}

fn role_label(role: CornerSubPatchRole) -> &'static str {
    match role {
        CornerSubPatchRole::NearOriginLeft => "NearOriginLeft",
        CornerSubPatchRole::NearOriginRight => "NearOriginRight",
        CornerSubPatchRole::NearMidpointRight => "NearMidpointRight",
        CornerSubPatchRole::NearMidpointLeft => "NearMidpointLeft",
        CornerSubPatchRole::SupportBoundaryLane => "SupportBoundaryLane",
        CornerSubPatchRole::SupportInteriorLane => "SupportInteriorLane",
        CornerSubPatchRole::SupportCell => "SupportCell",
        CornerSubPatchRole::SupportCellNearMidpointRight => "SupportCellNearMidpointRight",
        CornerSubPatchRole::SupportCellNearMidpointLeft => "SupportCellNearMidpointLeft",
        CornerSubPatchRole::Sector => "Sector",
    }
}
