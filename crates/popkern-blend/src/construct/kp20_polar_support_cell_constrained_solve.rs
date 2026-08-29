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
const FINAL_BOUNDARY_RESTORE_COMMERCIAL_G0_TARGET: f64 = 1.0e-10;
const FINAL_BOUNDARY_RESTORE_EPS: f64 = 1.0e-12;
const FINAL_BOUNDARY_RESTORE_G2_REGRESSION_BAND: f64 = 5.0e-2;

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
        let restore_report = restore_exact_partial_boundary_after_normal_if_quality_safe(
            brep,
            faces,
            "interior_adjacent_rows",
        )?;
        return Ok(format!(
            "g1_authority_class=interior_adjacent_rows;{};{}",
            normal_solve_status::worst_g1_normal_solve_status(&primary_report),
            restore_report
        ));
    }
    let corner_report =
        top_g1_normal_solve::reconcile_top_g1_corner_adjacent_normal_samples_if_improves(
            brep, faces, trims,
        )?;
    if normal_solve_status::report_was_kept(&corner_report) {
        let restore_report = restore_exact_partial_boundary_after_normal_if_quality_safe(
            brep,
            faces,
            "partial_endpoint_adjacent_corner",
        )?;
        return Ok(format!(
            "g1_authority_class=partial_endpoint_adjacent_corner;primary_normal_solve_status={};{};{}",
            primary_report.status,
            normal_solve_status::worst_g1_normal_solve_status(&corner_report),
            restore_report
        ));
    }
    let restore_report = restore_exact_partial_boundary_after_normal_if_quality_safe(
        brep,
        faces,
        "partial_endpoint_adjacent_rows",
    )?;
    Ok(format!(
        "g1_authority_class=partial_endpoint_adjacent_rows;{};{};{}",
        normal_solve_status::worst_g1_normal_solve_status(&primary_report),
        normal_solve_status::declined_authority_summary("corner_authority", &corner_report),
        restore_report
    ))
}

fn restore_exact_partial_boundary_after_normal_if_quality_safe(
    brep: &mut BRep,
    faces: &[FaceHandle],
    authority_class: &'static str,
) -> Result<String, CoreError> {
    let before_quality =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let before = brep.clone();
    let boundary_status =
        boundary_composed::reconcile_assigned_ranked_endpoint_exact_then_partial_status(
            brep, faces,
        )?;
    let after_quality =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let g0_reaches_commercial =
        after_quality.worst_g0_gap <= FINAL_BOUNDARY_RESTORE_COMMERCIAL_G0_TARGET;
    let aligned_g2_non_regressing = after_quality.worst_normal_aligned_g2_residual
        <= before_quality.worst_normal_aligned_g2_residual + FINAL_BOUNDARY_RESTORE_EPS;
    let aligned_g2_bounded_for_g0_frontier = g0_reaches_commercial
        && after_quality.worst_normal_aligned_g2_residual
            <= before_quality.worst_normal_aligned_g2_residual
                * (1.0 + FINAL_BOUNDARY_RESTORE_G2_REGRESSION_BAND);
    let keep = after_quality.pcurve_out_of_unit_count <= before_quality.pcurve_out_of_unit_count
        && after_quality.degenerate_pcurve_count <= before_quality.degenerate_pcurve_count
        && after_quality.worst_g0_gap + FINAL_BOUNDARY_RESTORE_EPS < before_quality.worst_g0_gap
        && after_quality.worst_normal_disagreement
            <= before_quality.worst_normal_disagreement + FINAL_BOUNDARY_RESTORE_EPS
        && (aligned_g2_non_regressing || aligned_g2_bounded_for_g0_frontier)
        && brep.validate().is_ok();
    if keep {
        return Ok(final_boundary_restore_status(
            "kept",
            authority_class,
            boundary_status,
            before_quality,
            after_quality,
        ));
    }
    *brep = before;
    Ok(final_boundary_restore_status(
        "declined",
        authority_class,
        boundary_status,
        before_quality,
        after_quality,
    ))
}

fn final_boundary_restore_status(
    status: &'static str,
    authority_class: &'static str,
    boundary_status: &'static str,
    before: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
    after: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> String {
    format!(
        "final_boundary_restore_status={status};final_boundary_restore_authority_class={authority_class};final_boundary_restore_boundary_status={boundary_status};final_boundary_restore_g0_before={:.15e};final_boundary_restore_g0_after={:.15e};final_boundary_restore_g1_before={:.15e};final_boundary_restore_g1_after={:.15e};final_boundary_restore_aligned_g2_before={:.15e};final_boundary_restore_aligned_g2_after={:.15e};final_boundary_restore_g2_regression_band={:.15e};final_boundary_restore_pcurve_out_before={};final_boundary_restore_pcurve_out_after={};final_boundary_restore_degenerate_before={};final_boundary_restore_degenerate_after={}",
        before.worst_g0_gap,
        after.worst_g0_gap,
        before.worst_normal_disagreement,
        after.worst_normal_disagreement,
        before.worst_normal_aligned_g2_residual,
        after.worst_normal_aligned_g2_residual,
        FINAL_BOUNDARY_RESTORE_G2_REGRESSION_BAND,
        before.pcurve_out_of_unit_count,
        after.pcurve_out_of_unit_count,
        before.degenerate_pcurve_count,
        after.degenerate_pcurve_count
    )
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
