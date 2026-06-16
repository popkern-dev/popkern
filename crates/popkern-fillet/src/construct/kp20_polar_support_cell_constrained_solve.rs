//! Constrained normal solve composed on exact/partial boundary ownership.

use popkern_core::error::CoreError;
use popkern_core::geometry::Vec3;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

use super::{
    assigned_support_cell_slots, boundary_composed, full_boundary_trims, normal_solve_status,
    reauthor_support_cell_shared_edge_curves, top_g1_normal_solve,
};

const PARTIAL_ENDPOINT_AUTHORITY_BAND: f64 = 0.25;

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
    let solve_status = reconcile_classified_top_g1_normal_authority(brep, faces, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell exact/partial constrained normal solve produced invalid BRep: {error}"
        ))
    })?;
    Ok(format!(
        "exact_partial_top_g1_normal_solve;base_status={base_status};{}",
        solve_status
    ))
}

fn reconcile_classified_top_g1_normal_authority(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[super::FullBoundaryTrim],
) -> Result<String, CoreError> {
    let primary_report =
        top_g1_normal_solve::reconcile_top_g1_normal_samples_if_improves(brep, faces, trims)?;
    let quality = shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    if !partial_endpoint_g1_authority_needed(quality) {
        return Ok(format!(
            "g1_authority_class=interior_adjacent_rows;{}",
            normal_solve_status::worst_g1_normal_solve_status(&primary_report)
        ));
    }
    let corner_report =
        top_g1_normal_solve::reconcile_top_g1_corner_adjacent_normal_samples_if_improves(
            brep, faces, trims,
        )?;
    if normal_solve_status::report_was_kept(&corner_report) {
        return Ok(format!(
            "g1_authority_class=partial_endpoint_adjacent_corner;primary_normal_solve_status={};{}",
            primary_report.status,
            normal_solve_status::worst_g1_normal_solve_status(&corner_report)
        ));
    }
    Ok(format!(
        "g1_authority_class=partial_endpoint_adjacent_rows;{};{}",
        normal_solve_status::worst_g1_normal_solve_status(&primary_report),
        normal_solve_status::declined_authority_summary("corner_authority", &corner_report)
    ))
}

fn partial_endpoint_g1_authority_needed(
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> bool {
    quality.top_g1_samples.iter().flatten().any(|sample| {
        !sample.trim.full_span
            && (uv_is_near_boundary_endpoint(sample.uv_a)
                || uv_is_near_boundary_endpoint(sample.uv_b))
    })
}

fn uv_is_near_boundary_endpoint(uv: Vec3) -> bool {
    let parameter = if near_unit_endpoint(uv.x) {
        uv.y
    } else if near_unit_endpoint(uv.y) {
        uv.x
    } else {
        return false;
    };
    parameter <= PARTIAL_ENDPOINT_AUTHORITY_BAND
        || parameter >= 1.0 - PARTIAL_ENDPOINT_AUTHORITY_BAND
}

fn near_unit_endpoint(value: f64) -> bool {
    value.abs() <= 1.0e-9 || (1.0 - value).abs() <= 1.0e-9
}
