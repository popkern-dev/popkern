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
