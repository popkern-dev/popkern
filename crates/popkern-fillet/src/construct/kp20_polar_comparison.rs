//! Structured comparison diagnostics for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput,
    Kp20MacroPatch, Kp20PieceId, RectangularGregoryBiNonic, TensorBorderBi3,
    assemble_kp20_bi6_macro_cap, build_kp20_guide_ring, elevate_bi6_grid_to_bi9,
    measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_orientation_valid_maps,
    sample_kp20_macro_patch_from_tau_maps,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_support_image::support_image_candidate_summary;
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, Kp20SupportNormalStep, normal_lift_candidate_summary,
};
use crate::types::ApexConstructionDiagnostic;

#[derive(Debug, Clone)]
pub(super) struct Kp20SamplerFoldSummary {
    pub fold_count: usize,
    pub min_jacobian_normal_cos: f64,
    pub piece_fold_counts: [usize; 4],
    pub fold_details: Vec<String>,
    pub(super) fold_sites: Vec<Kp20FoldSite>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20FoldSite {
    pub(super) sector: usize,
    pub(super) piece: Kp20PieceId,
    pub(super) u: f64,
    pub(super) v: f64,
    pub(super) min_jacobian_ratio: f64,
    pub(super) min_jacobian_normal_cos: f64,
}

pub(super) fn sampler_comparison_diagnostic(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    bounded_affine: Kp20SamplerFoldSummary,
) -> ApexConstructionDiagnostic {
    let orientation_valid = orientation_valid_candidate_summary(star_point, tensor_border);
    let normal_lift = normal_lift_candidate_summary(star_point, tensor_border);
    let support_image = support_image_candidate_summary(star_point, tensor_border);
    let sector_local_tau = tau_candidate_summary(star_point, tensor_border);
    let mut detail = format!(
        "bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={}",
        bounded_affine.fold_count,
        bounded_affine.min_jacobian_normal_cos,
        piece_counts_label(bounded_affine.piece_fold_counts)
    );
    append_candidate_detail(&mut detail, "orientation_valid", &orientation_valid);
    append_support_lift_detail(&mut detail, &normal_lift);
    append_candidate_detail(&mut detail, "support_image", &support_image);
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

pub(super) fn summarize_intersections(
    intersections: &[(
        usize,
        popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    )],
) -> Kp20SamplerFoldSummary {
    let mut piece_fold_counts = [0_usize; 4];
    let mut fold_details = Vec::new();
    for (idx, report) in intersections {
        piece_fold_counts[idx % 4] += 1;
        let sector = idx / 4;
        let piece = Kp20PieceId::all()[idx % 4];
        fold_details.push(format!(
            "sector={sector},piece={piece:?},ratio={:.3e},cos={:.3e}",
            report.min_jacobian_ratio, report.min_jacobian_normal_cos
        ));
    }
    let min_jacobian_normal_cos = intersections
        .iter()
        .map(|(_, report)| report.min_jacobian_normal_cos)
        .fold(f64::INFINITY, f64::min);
    Kp20SamplerFoldSummary {
        fold_count: intersections.len(),
        min_jacobian_normal_cos,
        piece_fold_counts,
        fold_details,
        fold_sites: Vec::new(),
    }
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

pub(super) fn candidate_summary_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<popkern_core::geometry::gregory::Kp20MacroPatch>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    candidate_summary_from_cap(tensor_border, cap.macro_patches)
}

pub(super) fn is_better_summary(
    candidate: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    best: Option<&(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport)>,
) -> bool {
    let Some(best) = best else {
        return true;
    };
    candidate.0.fold_count < best.0.fold_count
        || (candidate.0.fold_count == best.0.fold_count
            && candidate.0.min_jacobian_normal_cos > best.0.min_jacobian_normal_cos)
}

pub(super) fn validation_is_clean(validation: Kp20Bi6ValidationReport) -> bool {
    validation.max_tensor_border_residual <= 1.0e-8
        && validation.max_internal_c2_position_residual <= 1.0e-8
        && validation.max_internal_c2_first_residual <= 1.0e-8
        && validation.max_origin_g2_strip_residual <= 1.0e-8
}

pub(super) fn candidate_summary_from_cap(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = popkern_core::geometry::gregory::Kp20Bi6MacroCapOutput {
        n: tensor_border.sectors.len(),
        patches_bi6: flatten_macro_patches(&macro_patches),
        macro_patches,
        strip_projections: Vec::new(),
    };
    let validation = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    let candidate = cap
        .patches_bi6
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect::<Vec<_>>();
    let intersections = candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect::<Vec<_>>();
    Ok((
        summarize_candidate_intersections(&candidate, &intersections),
        validation,
    ))
}

pub(super) fn candidate_summary_from_changed_patches(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    unchanged_summary: &Kp20SamplerFoldSummary,
    changed_patches: &[[bool; 4]],
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = popkern_core::geometry::gregory::Kp20Bi6MacroCapOutput {
        n: tensor_border.sectors.len(),
        patches_bi6: flatten_macro_patches(&macro_patches),
        macro_patches,
        strip_projections: Vec::new(),
    };
    let validation = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    let mut fold_sites = unchanged_summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| !changed_patches[site.sector][site.piece.index()])
        .collect::<Vec<_>>();
    for (idx, grid) in cap.patches_bi6.into_iter().enumerate() {
        let sector = idx / 4;
        let piece = Kp20PieceId::all()[idx % 4];
        if !changed_patches[sector][piece.index()] {
            continue;
        }
        let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
            grid: elevate_bi6_grid_to_bi9(grid),
        }));
        let report = patch.self_intersection_report(SELF_INT_TOL);
        if report.intersects() {
            let sample = worst_jacobian_sample(&patch);
            fold_sites.push(Kp20FoldSite {
                sector,
                piece,
                u: sample.u,
                v: sample.v,
                min_jacobian_ratio: report.min_jacobian_ratio,
                min_jacobian_normal_cos: report.min_jacobian_normal_cos,
            });
        }
    }
    Ok((summary_from_fold_sites(fold_sites), validation))
}

fn summarize_candidate_intersections(
    candidate: &[GregoryPatch],
    intersections: &[(
        usize,
        popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    )],
) -> Kp20SamplerFoldSummary {
    let mut summary = summarize_intersections(intersections);
    summary.fold_sites = intersections
        .iter()
        .map(|(idx, report)| {
            let sector = idx / 4;
            let piece = Kp20PieceId::all()[idx % 4];
            let sample = worst_jacobian_sample(&candidate[*idx]);
            Kp20FoldSite {
                sector,
                piece,
                u: sample.u,
                v: sample.v,
                min_jacobian_ratio: report.min_jacobian_ratio,
                min_jacobian_normal_cos: report.min_jacobian_normal_cos,
            }
        })
        .collect();
    summary.fold_details = summary
        .fold_sites
        .iter()
        .map(|site| {
            format!(
                "sector={},piece={:?},u={:.3},v={:.3},ratio={:.3e},cos={:.3e}",
                site.sector,
                site.piece,
                site.u,
                site.v,
                site.min_jacobian_ratio,
                site.min_jacobian_normal_cos
            )
        })
        .collect();
    summary
}

fn summary_from_fold_sites(fold_sites: Vec<Kp20FoldSite>) -> Kp20SamplerFoldSummary {
    let mut piece_fold_counts = [0_usize; 4];
    for site in &fold_sites {
        piece_fold_counts[site.piece.index()] += 1;
    }
    let min_jacobian_normal_cos = fold_sites
        .iter()
        .map(|site| site.min_jacobian_normal_cos)
        .fold(f64::INFINITY, f64::min);
    let fold_details = fold_sites
        .iter()
        .map(|site| {
            format!(
                "sector={},piece={:?},u={:.3},v={:.3},ratio={:.3e},cos={:.3e}",
                site.sector,
                site.piece,
                site.u,
                site.v,
                site.min_jacobian_ratio,
                site.min_jacobian_normal_cos
            )
        })
        .collect();
    Kp20SamplerFoldSummary {
        fold_count: fold_sites.len(),
        min_jacobian_normal_cos,
        piece_fold_counts,
        fold_details,
        fold_sites,
    }
}

#[derive(Debug, Clone, Copy)]
struct WorstJacobianSample {
    u: f64,
    v: f64,
    min_jacobian_normal_cos: f64,
}

fn worst_jacobian_sample(patch: &GregoryPatch) -> WorstJacobianSample {
    let n = 8_usize;
    let mut samples = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let u = (i as f64 + 0.5) / n as f64;
            let v = (j as f64 + 0.5) / n as f64;
            let jacobian = patch.partial_u(u, v).cross(&patch.partial_v(u, v));
            let mag = jacobian.length();
            if mag >= 1.0e-14 {
                samples.push((u, v, jacobian, mag));
            }
        }
    }
    let Some((_, _, reference, max_mag)) = samples
        .iter()
        .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
        .copied()
    else {
        return WorstJacobianSample {
            u: 0.5,
            v: 0.5,
            min_jacobian_normal_cos: 1.0,
        };
    };
    let mag_threshold = 0.01 * max_mag;
    let mut worst = WorstJacobianSample {
        u: 0.5,
        v: 0.5,
        min_jacobian_normal_cos: 1.0,
    };
    for (u, v, jacobian, mag) in samples {
        if mag < mag_threshold {
            continue;
        }
        let cos = jacobian.dot(&reference) / (mag * max_mag);
        if cos < worst.min_jacobian_normal_cos {
            worst.u = u;
            worst.v = v;
            worst.min_jacobian_normal_cos = cos;
        }
    }
    worst
}

fn flatten_macro_patches(macro_patches: &[Kp20MacroPatch]) -> Vec<[[Point3; 7]; 7]> {
    let mut out = Vec::with_capacity(4 * macro_patches.len());
    for macro_patch in macro_patches {
        for id in Kp20PieceId::all() {
            out.push(macro_patch.piece(id).cps);
        }
    }
    out
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
    report: &Result<Kp20SupportNormalLiftReport, String>,
) {
    match report {
        Ok(report) => {
            detail.push_str(&format!(
                ";support_normal_lift_status=sampled;support_normal_lift_folds={};support_normal_lift_min_jacobian_cos={:.3e};support_normal_lift_piece_folds={};support_normal_lift_fold_details={};support_normal_lift_tensor_border_residual={:.3e};support_normal_lift_internal_c2_position={:.3e};support_normal_lift_internal_c2_first={:.3e};support_normal_lift_origin_g2_strip={:.3e};support_normal_lift_optimizer_active_fold_sites={};support_normal_lift_optimizer_row_variables={};support_normal_lift_optimizer_evaluated={};support_normal_lift_optimizer_accepted={};support_normal_lift_optimizer_rejected={};support_normal_lift_accepted_steps={};support_normal_lift_step_trace={}",
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
            detail.push_str(&format!(
                ";support_normal_lift_status=declined;support_normal_lift_error={err}"
            ));
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
