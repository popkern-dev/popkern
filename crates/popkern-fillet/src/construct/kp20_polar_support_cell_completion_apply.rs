//! Transactional support-cell pcurve/boundary completion application.

use std::borrow::Cow;

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{SupportCellCompletionApplication, support_cell_pcurve_author};
use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

pub(super) fn boundary_after_pcurve_if_quality_safe(
    brep: &mut BRep,
    faces: &[FaceHandle],
    reconcile_boundary: impl FnOnce(&mut BRep, &[FaceHandle]) -> Result<(), CoreError>,
) -> SupportCellCompletionApplication {
    boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
        reconcile_boundary(brep, faces).map(|_| Cow::Borrowed("boundary_kept"))
    })
}

pub(super) fn boundary_after_pcurve_with_status_if_quality_safe<S>(
    brep: &mut BRep,
    faces: &[FaceHandle],
    reconcile_boundary: impl FnOnce(&mut BRep, &[FaceHandle]) -> Result<S, CoreError>,
) -> SupportCellCompletionApplication
where
    S: Into<Cow<'static, str>>,
{
    if support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_err() {
        return SupportCellCompletionApplication {
            operator_succeeded: false,
            operation_status: Cow::Borrowed("pcurve_operator_failed"),
        };
    }
    let after_pcurve =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let before_boundary = brep.clone();
    let operation_status = match reconcile_boundary(brep, faces) {
        Ok(status) => status.into(),
        Err(_) => {
            *brep = before_boundary;
            return SupportCellCompletionApplication {
                operator_succeeded: true,
                operation_status: Cow::Borrowed("boundary_operator_failed_rolled_back"),
            };
        }
    };
    let after_boundary =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    if brep.validate().is_ok() && after_boundary.is_no_worse_than(after_pcurve) {
        SupportCellCompletionApplication {
            operator_succeeded: true,
            operation_status,
        }
    } else if brep.validate().is_err() {
        *brep = before_boundary;
        SupportCellCompletionApplication {
            operator_succeeded: true,
            operation_status: Cow::Borrowed("boundary_invalid_rolled_back"),
        }
    } else {
        *brep = before_boundary;
        // Carry the regression evidence in the status so a
        // rolled-back candidate is diagnosable from the fact
        // alone (which metric moved, by how much).
        SupportCellCompletionApplication {
            operator_succeeded: true,
            operation_status: Cow::Owned(format!(
                "boundary_quality_regressed_rolled_back;regressed_g0_after_pcurve={:.15e};regressed_g0_after_boundary={:.15e};regressed_out_of_unit_after_pcurve={};regressed_out_of_unit_after_boundary={};regressed_degenerate_after_pcurve={};regressed_degenerate_after_boundary={}",
                after_pcurve.worst_g0_gap,
                after_boundary.worst_g0_gap,
                after_pcurve.pcurve_out_of_unit_count,
                after_boundary.pcurve_out_of_unit_count,
                after_pcurve.degenerate_pcurve_count,
                after_boundary.degenerate_pcurve_count
            )),
        }
    }
}

pub(super) fn boundary_after_pcurve_with_status_if_g1_safe_against_baseline<S, B>(
    brep: &mut BRep,
    faces: &[FaceHandle],
    reconcile_baseline: impl FnOnce(&mut BRep, &[FaceHandle]) -> Result<B, CoreError>,
    reconcile_boundary: impl FnOnce(&mut BRep, &[FaceHandle]) -> Result<S, CoreError>,
) -> SupportCellCompletionApplication
where
    S: Into<Cow<'static, str>>,
    B: Into<Cow<'static, str>>,
{
    const G1_COMPARISON_EPS: f64 = 1.0e-12;
    if support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_err() {
        return SupportCellCompletionApplication {
            operator_succeeded: false,
            operation_status: Cow::Borrowed("pcurve_operator_failed"),
        };
    }
    let after_pcurve =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let mut baseline_brep = brep.clone();
    if reconcile_baseline(&mut baseline_brep, faces).is_err() || baseline_brep.validate().is_err() {
        return SupportCellCompletionApplication {
            operator_succeeded: false,
            operation_status: Cow::Borrowed("g1_safe_baseline_failed"),
        };
    }
    let baseline_quality =
        shared_trim_quality::measure_support_cell_shared_trim_quality(&baseline_brep, faces, 9);
    let before_boundary = brep.clone();
    let operation_status = match reconcile_boundary(brep, faces) {
        Ok(status) => status.into(),
        Err(_) => {
            *brep = before_boundary;
            return SupportCellCompletionApplication {
                operator_succeeded: false,
                operation_status: Cow::Borrowed("g1_safe_operator_failed_rolled_back"),
            };
        }
    };
    let after_boundary =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let g1_safe = after_boundary.worst_normal_disagreement
        <= baseline_quality.worst_normal_disagreement + G1_COMPARISON_EPS;
    if brep.validate().is_ok()
        && after_boundary.is_no_worse_than(after_pcurve)
        && after_boundary.is_no_worse_than(baseline_quality)
        && g1_safe
    {
        return SupportCellCompletionApplication {
            operator_succeeded: true,
            operation_status: Cow::Owned(format!(
                "{};joint_baseline_worst_g0={:.15e};joint_baseline_worst_g1={:.15e}",
                operation_status,
                baseline_quality.worst_g0_gap,
                baseline_quality.worst_normal_disagreement
            )),
        };
    }
    *brep = before_boundary;
    SupportCellCompletionApplication {
        operator_succeeded: false,
        operation_status: Cow::Owned(format!(
            "g1_safe_quality_declined_rolled_back;declined_operation_status=[{}];baseline_g0={:.15e};candidate_g0={:.15e};baseline_g1={:.15e};candidate_g1={:.15e};baseline_out_of_unit={};candidate_out_of_unit={};baseline_degenerate={};candidate_degenerate={}",
            operation_status,
            baseline_quality.worst_g0_gap,
            after_boundary.worst_g0_gap,
            baseline_quality.worst_normal_disagreement,
            after_boundary.worst_normal_disagreement,
            baseline_quality.pcurve_out_of_unit_count,
            after_boundary.pcurve_out_of_unit_count,
            baseline_quality.degenerate_pcurve_count,
            after_boundary.degenerate_pcurve_count
        )),
    }
}
