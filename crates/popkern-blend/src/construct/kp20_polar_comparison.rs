//! Structured comparison diagnostics for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20PieceId, Kp20TauMap, TensorBorderBi3,
    kp20_c2_stitched_sector_piece_maps, kp20_coupled_stitched_sector_piece_maps,
    kp20_default_piece_map_domain_report, kp20_paper_oriented_sector_piece_maps,
    kp20_piece_map_domain_report, kp20_piece_map_stitch_report,
    kp20_positive_nonoverlap_c2_sector_piece_maps,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_multipatch_labels::{
    multipatch_evaluation_label, multipatch_inner_label, multipatch_interfaces_label,
    multipatch_inventory_label, multipatch_strips_label, multipatch_validation_label,
};
use super::kp20_polar_support_image::support_image_composed_candidate;
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, Kp20SupportNormalStep, normal_lift_candidate_summary,
};
use super::kp20_polar_support_surface::{
    Kp20SupportSurfaceReport, support_surface_candidate_summary,
};
use crate::types::ApexConstructionDiagnostic;
use candidates::{
    lifted_mapped_candidate_summary, mapped_candidate_summary, tau_candidate_summary,
};

#[path = "kp20_polar_comparison_candidates.rs"]
mod candidates;

pub(super) fn sampler_comparison_diagnostic(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    bounded_affine: Kp20SamplerFoldSummary,
) -> ApexConstructionDiagnostic {
    let total_start = popkern_core::time::Stopwatch::start();
    let (orientation_valid, orientation_valid_ms) = timed(|| {
        mapped_candidate_summary(
            star_point,
            tensor_border,
            popkern_core::geometry::gregory::kp20_orientation_valid_sector_piece_maps(),
        )
    });
    let (paper_oriented, paper_oriented_ms) = timed(|| {
        mapped_candidate_summary(
            star_point,
            tensor_border,
            kp20_paper_oriented_sector_piece_maps(),
        )
    });
    let (c2_stitched, c2_stitched_ms) = timed(|| {
        mapped_candidate_summary(
            star_point,
            tensor_border,
            kp20_c2_stitched_sector_piece_maps(),
        )
    });
    let (coupled_stitched, coupled_stitched_ms) = timed(|| {
        mapped_candidate_summary(
            star_point,
            tensor_border,
            kp20_coupled_stitched_sector_piece_maps(),
        )
    });
    let (nonoverlap_coupled, nonoverlap_coupled_ms) = timed(|| {
        mapped_candidate_summary(
            star_point,
            tensor_border,
            kp20_positive_nonoverlap_c2_sector_piece_maps(),
        )
    });
    let (nonoverlap_coupled_lift, nonoverlap_coupled_lift_ms) = timed(|| {
        lifted_mapped_candidate_summary(
            star_point,
            tensor_border,
            kp20_positive_nonoverlap_c2_sector_piece_maps(),
        )
    });
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
        "sampler_total_ms={total_ms};orientation_valid_ms={orientation_valid_ms};paper_oriented_ms={paper_oriented_ms};c2_stitched_ms={c2_stitched_ms};coupled_stitched_ms={coupled_stitched_ms};nonoverlap_coupled_ms={nonoverlap_coupled_ms};nonoverlap_coupled_lift_ms={nonoverlap_coupled_lift_ms};support_normal_lift_ms={normal_lift_ms};support_image_composed_ms={support_image_composed_ms};sector_local_tau_ms={sector_local_tau_ms};bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={}",
        bounded_affine.fold_count,
        bounded_affine.min_jacobian_normal_cos,
        piece_counts_label(bounded_affine.piece_fold_counts)
    );
    append_candidate_detail(&mut detail, "orientation_valid", &orientation_valid);
    append_candidate_detail(&mut detail, "paper_oriented", &paper_oriented);
    append_candidate_detail(&mut detail, "c2_stitched", &c2_stitched);
    append_candidate_detail(&mut detail, "coupled_stitched", &coupled_stitched);
    append_candidate_detail(&mut detail, "nonoverlap_coupled", &nonoverlap_coupled);
    append_support_lift_detail(
        &mut detail,
        "nonoverlap_coupled_lift",
        &nonoverlap_coupled_lift,
    );
    append_support_lift_detail(&mut detail, "support_normal_lift", &normal_lift);
    append_candidate_detail(&mut detail, "support_image", &support_image);
    append_support_lift_detail(&mut detail, "support_image_lift", &support_image_lift);
    append_support_surface_detail(&mut detail, "support_surface", &support_surface);
    append_candidate_detail(&mut detail, "sector_local_tau", &sector_local_tau);
    append_paper_alignment_detail(
        &mut detail,
        tensor_border,
        &bounded_affine,
        &paper_oriented,
        &c2_stitched,
        &coupled_stitched,
        &nonoverlap_coupled,
        &sector_local_tau,
    );
    comparison_diagnostic(
        tensor_border.sectors.len(),
        detail,
        sector_local_tau
            .as_ref()
            .ok()
            .map(|(summary, _)| summary.min_jacobian_normal_cos),
    )
}

#[allow(clippy::too_many_arguments)]
fn append_paper_alignment_detail(
    detail: &mut String,
    tensor_border: &TensorBorderBi3,
    bounded_affine: &Kp20SamplerFoldSummary,
    paper_oriented: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
    c2_stitched: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
    coupled_stitched: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
    nonoverlap_coupled: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
    sector_local_tau: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
) {
    let domain = kp20_default_piece_map_domain_report();
    let paper_domain = kp20_piece_map_domain_report(kp20_paper_oriented_sector_piece_maps());
    let default_stitch = kp20_piece_map_stitch_report(
        popkern_core::geometry::gregory::kp20_default_sector_piece_maps(),
    );
    let paper_stitch = kp20_piece_map_stitch_report(kp20_paper_oriented_sector_piece_maps());
    let c2_domain = kp20_piece_map_domain_report(kp20_c2_stitched_sector_piece_maps());
    let c2_stitch = kp20_piece_map_stitch_report(kp20_c2_stitched_sector_piece_maps());
    let coupled_domain = kp20_piece_map_domain_report(kp20_coupled_stitched_sector_piece_maps());
    let coupled_stitch = kp20_piece_map_stitch_report(kp20_coupled_stitched_sector_piece_maps());
    let nonoverlap_domain =
        kp20_piece_map_domain_report(kp20_positive_nonoverlap_c2_sector_piece_maps());
    let nonoverlap_stitch =
        kp20_piece_map_stitch_report(kp20_positive_nonoverlap_c2_sector_piece_maps());
    let tau_available = (0..tensor_border.sectors.len())
        .all(|sector| Kp20TauMap::for_sector_local(tensor_border.sectors.len(), sector).is_ok());
    let tau_label = candidate_label(sector_local_tau);
    let paper_label = candidate_label(paper_oriented);
    let c2_label = candidate_label(c2_stitched);
    let coupled_label = candidate_label(coupled_stitched);
    let nonoverlap_label = candidate_label(nonoverlap_coupled);
    detail.push_str(&format!(
        ";kp20_paper_alignment=valence:{}|pieces_per_sector:4|piece_degree:bi6|diagnostic_degree:bi9|paper_domain_inside:{}|outside_corners:{}|max_simplex_excess:{:.3e}|positive_orientation_maps:{}|non_positive_orientation_maps:{}|paper_oriented_inside:{}|paper_oriented_non_positive_maps:{}|c2_stitched_inside:{}|c2_stitched_non_positive_maps:{}|coupled_stitched_inside:{}|coupled_stitched_non_positive_maps:{}|nonoverlap_coupled_inside:{}|nonoverlap_coupled_non_positive_maps:{}|default_c2_edge:{:.3e}|paper_c2_edge:{:.3e}|c2_stitched_c2_edge:{:.3e}|c2_stitched_origin_strip_edge:{:.3e}|coupled_stitched_c2_edge:{:.3e}|coupled_stitched_origin_strip_edge:{:.3e}|coupled_stitched_min_center_sep:{:.3e}|nonoverlap_coupled_c2_edge:{:.3e}|nonoverlap_coupled_origin_strip_edge:{:.3e}|nonoverlap_coupled_origin_reversed:{}|nonoverlap_coupled_min_center_sep:{:.3e}|tau_available:{}|bounded_affine_folds:{}|bounded_affine_piece_folds:{}|paper_oriented:{}|c2_stitched:{}|coupled_stitched:{}|nonoverlap_coupled:{}|sector_local_tau:{}",
        tensor_border.sectors.len(),
        domain.stays_inside_triangle(),
        domain.outside_triangle_corners,
        domain.max_simplex_excess,
        4_usize.saturating_sub(domain.non_positive_orientation_maps),
        domain.non_positive_orientation_maps,
        paper_domain.stays_inside_triangle(),
        paper_domain.non_positive_orientation_maps,
        c2_domain.stays_inside_triangle(),
        c2_domain.non_positive_orientation_maps,
        coupled_domain.stays_inside_triangle(),
        coupled_domain.non_positive_orientation_maps,
        nonoverlap_domain.stays_inside_triangle(),
        nonoverlap_domain.non_positive_orientation_maps,
        default_stitch.max_internal_c2_edge_mismatch,
        paper_stitch.max_internal_c2_edge_mismatch,
        c2_stitch.max_internal_c2_edge_mismatch,
        c2_stitch.max_origin_strip_edge_mismatch,
        coupled_stitch.max_internal_c2_edge_mismatch,
        coupled_stitch.max_origin_strip_edge_mismatch,
        coupled_stitch.min_piece_center_separation,
        nonoverlap_stitch.max_internal_c2_edge_mismatch,
        nonoverlap_stitch.max_origin_strip_edge_mismatch,
        nonoverlap_stitch.origin_strip_reversed,
        nonoverlap_stitch.min_piece_center_separation,
        tau_available,
        bounded_affine.fold_count,
        piece_counts_label(bounded_affine.piece_fold_counts),
        paper_label,
        c2_label,
        coupled_label,
        nonoverlap_label,
        tau_label,
    ));
}

fn candidate_label(
    candidate: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
) -> String {
    candidate.as_ref().map_or_else(
        |err| format!("declined:{err}"),
        |(summary, validation)| {
            format!(
                "sampled|folds:{}|piece_folds:{}|validation:{}",
                summary.fold_count,
                piece_counts_label(summary.piece_fold_counts),
                validation_label(validation)
            )
        },
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
        support_surface_candidate_summary(star_point, tensor_border, &candidate.lift)
    })
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, u128) {
    let start = popkern_core::time::Stopwatch::start();
    let value = f();
    (value, elapsed_ms(start))
}

fn elapsed_ms(start: popkern_core::time::Stopwatch) -> u128 {
    start.elapsed_or_zero_ms()
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
                ";{label}_status=sampled;{label}_folds={};{label}_min_jacobian_cos={:.3e};{label}_piece_folds={};{label}_fold_details={};{label}_tensor_border_residual={:.3e};{label}_internal_c2_position={:.3e};{label}_internal_c2_first={:.3e};{label}_origin_g2_strip={:.3e};{label}_active_fold_pieces={};{label}_orientation_rows={};{label}_evaluated={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_stage_trace={};{label}_multipatch_inventory={};{label}_multipatch_strips={};{label}_multipatch_interfaces={};{label}_multipatch_inner={};{label}_multipatch_validation={};{label}_multipatch_evaluation={};{label}_decline_reason={}",
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
                report.stats.exact_evaluated_candidates,
                report.stats.accepted_candidates,
                report.stats.rejected_candidates,
                support_surface_stage_trace_label(&report.stage_trace),
                multipatch_inventory_label(&report.multipatch_inventory),
                multipatch_strips_label(&report.multipatch_strips),
                multipatch_interfaces_label(&report.multipatch_interfaces),
                multipatch_inner_label(&report.multipatch_inner),
                multipatch_validation_label(&report.multipatch_validation),
                multipatch_evaluation_label(&report.multipatch_evaluation),
                report.decline_reason,
            ));
        }
        Err(err) => {
            detail.push_str(&format!(";{label}_status=declined;{label}_error={err}"));
        }
    }
}

fn support_surface_stage_trace_label(
    stages: &[super::kp20_polar_support_surface::Kp20SupportSurfaceStageReport],
) -> String {
    if stages.is_empty() {
        return "none".to_string();
    }
    stages
        .iter()
        .map(|stage| {
            format!(
                "{}:eval={},fast={},exact={},exact_folds={},accepted={}",
                stage.mode,
                stage.evaluated_candidates,
                optional_usize_label(stage.fast_finalist_folds),
                stage.exact_evaluated,
                optional_usize_label(stage.exact_folds),
                stage.accepted,
            )
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn optional_usize_label(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "none".to_string())
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

fn validation_label(validation: &Kp20Bi6ValidationReport) -> String {
    format!(
        "tensor:{:.3e}|c2_pos:{:.3e}|c2_first:{:.3e}|origin_g2:{:.3e}",
        validation.max_tensor_border_residual,
        validation.max_internal_c2_position_residual,
        validation.max_internal_c2_first_residual,
        validation.max_origin_g2_strip_residual,
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
