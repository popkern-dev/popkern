//! Structured comparison diagnostics for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20PieceId, TensorBorderBi3,
    build_kp20_guide_ring, sample_kp20_macro_patch_from_orientation_valid_maps,
    sample_kp20_macro_patch_from_tau_maps,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_fold_summary::{Kp20SamplerFoldSummary, candidate_summary_from_sampled};
use super::kp20_polar_support_image::support_image_composed_candidate;
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, Kp20SupportNormalStep, normal_lift_candidate_summary,
};
use super::kp20_polar_support_surface::{
    Kp20SupportSurfaceReport, support_surface_candidate_summary,
};
use crate::types::ApexConstructionDiagnostic;
use std::time::Instant;

pub(super) fn sampler_comparison_diagnostic(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    bounded_affine: Kp20SamplerFoldSummary,
) -> ApexConstructionDiagnostic {
    let total_start = Instant::now();
    let (orientation_valid, orientation_valid_ms) =
        timed(|| orientation_valid_candidate_summary(star_point, tensor_border));
    let (normal_lift, normal_lift_ms) =
        timed(|| normal_lift_candidate_summary(star_point, tensor_border));
    let (support_image_composed, support_image_composed_ms) =
        timed(|| support_image_composed_candidate(star_point, tensor_border));
    let support_image = support_image_result(&support_image_composed);
    let support_image_lift = support_image_lift_result(&support_image_composed);
    let support_surface =
        support_surface_result(star_point, tensor_border, &support_image_composed);
    let (sector_local_tau, sector_local_tau_ms) =
        timed(|| tau_candidate_summary(star_point, tensor_border));
    let total_ms = elapsed_ms(total_start);
    let mut detail = format!(
        "sampler_total_ms={total_ms};orientation_valid_ms={orientation_valid_ms};support_normal_lift_ms={normal_lift_ms};support_image_composed_ms={support_image_composed_ms};sector_local_tau_ms={sector_local_tau_ms};bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={}",
        bounded_affine.fold_count,
        bounded_affine.min_jacobian_normal_cos,
        piece_counts_label(bounded_affine.piece_fold_counts)
    );
    append_candidate_detail(&mut detail, "orientation_valid", &orientation_valid);
    append_support_lift_detail(&mut detail, "support_normal_lift", &normal_lift);
    append_candidate_detail(&mut detail, "support_image", &support_image);
    append_support_lift_detail(&mut detail, "support_image_lift", &support_image_lift);
    append_support_surface_detail(&mut detail, "support_surface", &support_surface);
    append_candidate_detail(&mut detail, "sector_local_tau", &sector_local_tau);
    comparison_diagnostic(
        tensor_border.sectors.len(),
        detail,
        sector_local_tau
            .as_ref()
            .ok()
            .map(|(summary, _)| summary.min_jacobian_normal_cos),
    )
}

fn support_image_result(
    report: &Result<super::kp20_polar_support_image::Kp20SupportImageSearchCandidate, String>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    report
        .as_ref()
        .map(|candidate| (candidate.image.summary.clone(), candidate.image.validation))
        .map_err(Clone::clone)
}

fn support_image_lift_result(
    report: &Result<super::kp20_polar_support_image::Kp20SupportImageSearchCandidate, String>,
) -> Result<Kp20SupportNormalLiftReport, String> {
    report
        .as_ref()
        .map(|candidate| candidate.lift.clone())
        .map_err(Clone::clone)
}

fn support_surface_result(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    report: &Result<super::kp20_polar_support_image::Kp20SupportImageSearchCandidate, String>,
) -> Result<Kp20SupportSurfaceReport, String> {
    report.as_ref().map_err(Clone::clone).and_then(|candidate| {
        support_surface_candidate_summary(star_point, tensor_border, &candidate.image)
    })
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, u128) {
    let start = Instant::now();
    let value = f();
    (value, elapsed_ms(start))
}

fn elapsed_ms(start: Instant) -> u128 {
    start.elapsed().as_millis()
}

fn tau_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .enumerate()
        .map(|(sector, guide_sector)| {
            sample_kp20_macro_patch_from_tau_maps(guide_sector, guide.n, sector)
        })
        .collect::<Result<Vec<_>, _>>()?;
    candidate_summary_from_sampled(tensor_border, sampled_macro_patches)
}

fn orientation_valid_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(sample_kp20_macro_patch_from_orientation_valid_maps)
        .collect::<Vec<_>>();
    candidate_summary_from_sampled(tensor_border, sampled_macro_patches)
}

fn append_candidate_detail(
    detail: &mut String,
    label: &str,
    summary: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
) {
    match summary {
        Ok((folds, validation)) => {
            detail.push_str(&format!(
                ";{label}_status=sampled;{label}_folds={};{label}_min_jacobian_cos={:.3e};{label}_piece_folds={};{label}_fold_details={};{label}_tensor_border_residual={:.3e};{label}_internal_c2_position={:.3e};{label}_internal_c2_first={:.3e};{label}_origin_g2_strip={:.3e}",
                folds.fold_count,
                folds.min_jacobian_normal_cos,
                piece_counts_label(folds.piece_fold_counts),
                fold_details_label(&folds.fold_details),
                validation.max_tensor_border_residual,
                validation.max_internal_c2_position_residual,
                validation.max_internal_c2_first_residual,
                validation.max_origin_g2_strip_residual,
            ));
        }
        Err(err) => {
            detail.push_str(&format!(";{label}_status=declined;{label}_error={err}"));
        }
    }
}

fn append_support_lift_detail(
    detail: &mut String,
    label: &str,
    report: &Result<Kp20SupportNormalLiftReport, String>,
) {
    match report {
        Ok(report) => {
            detail.push_str(&format!(
                ";{label}_status=sampled;{label}_folds={};{label}_min_jacobian_cos={:.3e};{label}_piece_folds={};{label}_fold_details={};{label}_tensor_border_residual={:.3e};{label}_internal_c2_position={:.3e};{label}_internal_c2_first={:.3e};{label}_origin_g2_strip={:.3e};{label}_optimizer_active_fold_sites={};{label}_optimizer_row_variables={};{label}_optimizer_evaluated={};{label}_optimizer_accepted={};{label}_optimizer_rejected={};{label}_accepted_steps={};{label}_step_trace={}",
                report.summary.fold_count,
                report.summary.min_jacobian_normal_cos,
                piece_counts_label(report.summary.piece_fold_counts),
                fold_details_label(&report.summary.fold_details),
                report.validation.max_tensor_border_residual,
                report.validation.max_internal_c2_position_residual,
                report.validation.max_internal_c2_first_residual,
                report.validation.max_origin_g2_strip_residual,
                report.optimizer_stats.active_fold_sites,
                report.optimizer_stats.row_variable_count,
                report.optimizer_stats.evaluated_candidates,
                report.optimizer_stats.accepted_candidates,
                report.optimizer_stats.rejected_candidates,
                report.accepted_steps.len(),
                step_trace_label(&report.accepted_steps),
            ));
        }
        Err(err) => {
            detail.push_str(&format!(";{label}_status=declined;{label}_error={err}"));
        }
    }
}

fn append_support_surface_detail(
    detail: &mut String,
    label: &str,
    report: &Result<Kp20SupportSurfaceReport, String>,
) {
    match report {
        Ok(report) => {
            detail.push_str(&format!(
                ";{label}_status=sampled;{label}_folds={};{label}_min_jacobian_cos={:.3e};{label}_piece_folds={};{label}_fold_details={};{label}_tensor_border_residual={:.3e};{label}_internal_c2_position={:.3e};{label}_internal_c2_first={:.3e};{label}_origin_g2_strip={:.3e};{label}_active_fold_pieces={};{label}_orientation_rows={};{label}_evaluated={};{label}_accepted={};{label}_rejected={};{label}_decline_reason={}",
                report.summary.fold_count,
                report.summary.min_jacobian_normal_cos,
                piece_counts_label(report.summary.piece_fold_counts),
                fold_details_label(&report.summary.fold_details),
                report.validation.max_tensor_border_residual,
                report.validation.max_internal_c2_position_residual,
                report.validation.max_internal_c2_first_residual,
                report.validation.max_origin_g2_strip_residual,
                report.stats.active_fold_pieces,
                report.stats.orientation_rows,
                report.stats.evaluated_candidates,
                report.stats.accepted_candidates,
                report.stats.rejected_candidates,
                report.decline_reason,
            ));
        }
        Err(err) => {
            detail.push_str(&format!(";{label}_status=declined;{label}_error={err}"));
        }
    }
}

fn comparison_diagnostic(
    n_sides: usize,
    detail: String,
    residual: Option<f64>,
) -> ApexConstructionDiagnostic {
    super::kp20_polar_diagnostics::dispatch_fact(
        "sampler_comparison",
        n_sides,
        Some(detail),
        residual,
        Some(SELF_INT_TOL),
    )
}

fn piece_counts_label(counts: [usize; 4]) -> String {
    let pieces = Kp20PieceId::all();
    format!(
        "{:?}:{}|{:?}:{}|{:?}:{}|{:?}:{}",
        pieces[0], counts[0], pieces[1], counts[1], pieces[2], counts[2], pieces[3], counts[3]
    )
}

fn fold_details_label(details: &[String]) -> String {
    if details.is_empty() {
        return "none".to_string();
    }
    details.join("|")
}

fn step_trace_label(steps: &[Kp20SupportNormalStep]) -> String {
    if steps.is_empty() {
        return "none".to_string();
    }
    steps
        .iter()
        .enumerate()
        .map(|(idx, step)| {
            format!(
                "{}:{}@{:.3}->{}:{}:{:?}",
                idx,
                piece_mask_label(step.piece_mask),
                step.scale,
                step.fold_count,
                piece_counts_label(step.piece_fold_counts),
                step.mode
            )
        })
        .collect::<Vec<_>>()
        .join("/")
}

fn piece_mask_label(mask: [bool; 4]) -> String {
    let pieces = Kp20PieceId::all();
    pieces
        .into_iter()
        .zip(mask)
        .filter_map(|(piece, enabled)| enabled.then_some(format!("{piece:?}")))
        .collect::<Vec<_>>()
        .join("+")
}
