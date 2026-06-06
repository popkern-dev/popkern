//! Support-image reparameterization candidates for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    Kp20PieceId, Kp20SectorPieceMaps, TensorBorderBi3, assemble_kp20_bi6_macro_cap,
    build_kp20_guide_ring, kp20_default_sector_piece_maps,
    kp20_positive_nonoverlap_c2_sector_piece_maps, sample_kp20_macro_patch_from_piece_maps,
};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_cap,
    fast_candidate_summary_from_changed_piece_sites, is_better_summary, validation_is_clean,
};
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_support_image_macro_patches,
};
use super::kp20_polar_support_reparameterization::Kp20SupportReparameterizationReport;
#[path = "kp20_polar_support_image_params.rs"]
mod params;
use params::Kp20SupportImageParams;
use std::time::Instant;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportImageCandidate {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) macro_patches: Vec<Kp20MacroPatch>,
    pub(super) params: Kp20SupportImageParams,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportImageSearchCandidate {
    pub(super) image: Kp20SupportImageCandidate,
    pub(super) lift: Kp20SupportNormalLiftReport,
}

struct Kp20AssembledSupportImageCandidate {
    macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SupportImageParams,
}

type LiftReport = Kp20SupportNormalLiftReport;
type ReparamReport = Kp20SupportReparameterizationReport;
type SupportReparamResult = Result<(LiftReport, ReparamReport), String>;

pub(super) fn refine_support_side_reparameterization(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> SupportReparamResult {
    let target_piece = Kp20PieceId::NearMidpointLeft;
    let mut report = Kp20SupportReparameterizationReport {
        evaluated: 0,
        fast_accepted: 0,
        fast_rejected: 0,
        lift_evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.summary.fold_count,
        end_folds: baseline.summary.fold_count,
        target_piece,
        params: "none".to_string(),
        decline_reason: "support_reparameterization_no_candidate",
    };
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let baseline_summary = (baseline.summary.clone(), baseline.validation);
    let mut best_fast: Option<(
        Kp20AssembledSupportImageCandidate,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    )> = None;
    for params in Kp20SupportImageParams::support_side_reparameterization_family() {
        report.evaluated += 1;
        let mut changed_piece_sites = Vec::new();
        let sampled_macro_patches = guide
            .sectors
            .iter()
            .enumerate()
            .map(|(sector, guide_sector)| {
                let (maps, changed_sites) = targeted_support_image_piece_maps(
                    kp20_positive_nonoverlap_c2_sector_piece_maps(),
                    sector,
                    params,
                    &baseline.summary.fold_sites,
                    target_piece,
                );
                changed_piece_sites.extend(changed_sites);
                sample_kp20_macro_patch_from_piece_maps(guide_sector, maps)
            })
            .collect::<Vec<_>>();
        if changed_piece_sites.is_empty() {
            report.rejected += 1;
            continue;
        }
        let candidate =
            assembled_candidate_from_sampled(tensor_border, sampled_macro_patches, params)?;
        let fast = fast_candidate_summary_from_changed_piece_sites(
            tensor_border,
            candidate.macro_patches.clone(),
            &baseline.summary,
            &changed_piece_sites,
            baseline.summary.fold_count,
        )?;
        if !validation_is_clean(fast.1) || fast.0.fold_count > baseline.summary.fold_count {
            report.fast_rejected += 1;
            report.rejected += 1;
            continue;
        }
        report.fast_accepted += 1;
        if best_fast
            .as_ref()
            .is_none_or(|(_, best)| is_better_summary(&fast, Some(best)))
        {
            best_fast = Some((candidate, fast));
        } else {
            report.rejected += 1;
        }
    }
    let Some((candidate, _)) = best_fast else {
        report.rejected = report.evaluated;
        return Ok((baseline, report));
    };
    report.lift_evaluated = 1;
    let params = candidate.params;
    let lift = lift_candidate_from_support_image_macro_patches(
        star_point,
        tensor_border,
        candidate.macro_patches,
    )?;
    if !validation_is_clean(lift.validation) {
        report.rejected = report.evaluated;
        report.decline_reason = "support_reparameterization_lift_validation_not_clean";
        return Ok((baseline, report));
    }
    let lift_summary = (lift.summary.clone(), lift.validation);
    if !is_better_summary(&lift_summary, Some(&baseline_summary)) {
        report.rejected = report.evaluated;
        report.decline_reason = "support_reparameterization_no_lift_improvement";
        return Ok((baseline, report));
    }
    report.end_folds = lift.summary.fold_count;
    report.params = params.label();
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.decline_reason = "none";
    Ok((lift, report))
}

pub(super) fn support_image_composed_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20SupportImageSearchCandidate, String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let default_sampled = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(sector, kp20_default_sector_piece_maps())
        })
        .collect::<Vec<_>>();
    let baseline = candidate_from_sampled(
        tensor_border,
        default_sampled,
        Kp20SupportImageParams::identity(),
    )?;
    let image_search_start = Instant::now();
    let mut image_candidates: Vec<Kp20SupportImageCandidate> = Vec::new();
    let mut evaluated = 0_usize;
    let mut fast_evaluated = 0_usize;
    let mut exact_evaluated = 0_usize;
    let mut rejected = 0_usize;
    let mut fast_rejected = 0_usize;
    let mut exact_rejected = 0_usize;
    let baseline_summary = (baseline.summary.clone(), baseline.validation);
    for params in Kp20SupportImageParams::coordinate_search_family() {
        evaluated += 1;
        let mut changed_piece_sites = Vec::new();
        let sampled_macro_patches = guide
            .sectors
            .iter()
            .enumerate()
            .map(|(sector, guide_sector)| {
                let (maps, changed_sites) = support_image_piece_maps(
                    kp20_default_sector_piece_maps(),
                    sector,
                    params,
                    &baseline.summary.fold_sites,
                );
                changed_piece_sites.extend(changed_sites);
                sample_kp20_macro_patch_from_piece_maps(guide_sector, maps)
            })
            .collect::<Vec<_>>();
        let candidate =
            assembled_candidate_from_sampled(tensor_border, sampled_macro_patches, params)?;
        fast_evaluated += 1;
        let fast_summary = fast_candidate_summary_from_changed_piece_sites(
            tensor_border,
            candidate.macro_patches.clone(),
            &baseline.summary,
            &changed_piece_sites,
            baseline.summary.fold_count,
        )?;
        if !validation_is_clean(fast_summary.1) {
            rejected += 1;
            fast_rejected += 1;
            continue;
        }
        if !is_better_summary(&fast_summary, Some(&baseline_summary)) {
            rejected += 1;
            fast_rejected += 1;
            continue;
        }
        if image_candidates.first().is_some_and(|best| {
            !is_better_summary(
                &fast_summary,
                Some(&(best.summary.clone(), best.validation)),
            )
        }) {
            rejected += 1;
            fast_rejected += 1;
            continue;
        }
        exact_evaluated += 1;
        let candidate = exact_summary_candidate(tensor_border, candidate)?;
        if !validation_is_clean(candidate.validation) {
            rejected += 1;
            exact_rejected += 1;
            continue;
        }
        let candidate_summary = (candidate.summary.clone(), candidate.validation);
        if is_better_summary(&candidate_summary, Some(&baseline_summary)) {
            image_candidates.push(candidate);
            image_candidates.sort_by(|a, b| compare_summaries(&a.summary, &b.summary));
        } else {
            rejected += 1;
            exact_rejected += 1;
        }
    }
    let image_search_ms = elapsed_ms(image_search_start);
    image_candidates.sort_by(|a, b| compare_summaries(&a.summary, &b.summary));
    let mut best: Option<Kp20SupportImageSearchCandidate> = None;
    let mut composed_evaluated = 0_usize;
    let mut lift_budget_skipped = 0_usize;
    let mut targeted_evaluated = 0_usize;
    let mut targeted_accepted = 0_usize;
    let lift_start = Instant::now();
    for candidate in image_candidates {
        if composed_evaluated >= COUPLED_LIFT_FINALIST_BUDGET {
            lift_budget_skipped += 1;
            continue;
        }
        composed_evaluated += 1;
        let lift = lift_candidate_from_support_image_macro_patches(
            star_point,
            tensor_border,
            candidate.macro_patches.clone(),
        )?;
        if !validation_is_clean(lift.validation) {
            rejected += 1;
            continue;
        }
        let candidate_lift_summary = (lift.summary.clone(), lift.validation);
        let best_lift_summary = best
            .as_ref()
            .map(|candidate| (candidate.lift.summary.clone(), candidate.lift.validation));
        if is_better_summary(&candidate_lift_summary, best_lift_summary.as_ref()) {
            best = Some(Kp20SupportImageSearchCandidate {
                image: candidate,
                lift,
            });
        } else {
            rejected += 1;
        }
    }
    if let Some(current) = best.clone() {
        let fold_sites = current.lift.summary.fold_sites.clone();
        for params in Kp20SupportImageParams::lift_informed_family(current.image.params) {
            if fold_sites.is_empty() {
                break;
            }
            evaluated += 1;
            targeted_evaluated += 1;
            let mut changed_piece_sites = Vec::new();
            let sampled_macro_patches = guide
                .sectors
                .iter()
                .enumerate()
                .map(|(sector, guide_sector)| {
                    let (maps, changed_sites) = support_image_piece_maps(
                        kp20_default_sector_piece_maps(),
                        sector,
                        params,
                        &fold_sites,
                    );
                    changed_piece_sites.extend(changed_sites);
                    sample_kp20_macro_patch_from_piece_maps(guide_sector, maps)
                })
                .collect::<Vec<_>>();
            let candidate =
                assembled_candidate_from_sampled(tensor_border, sampled_macro_patches, params)?;
            fast_evaluated += 1;
            let fast_summary = fast_candidate_summary_from_changed_piece_sites(
                tensor_border,
                candidate.macro_patches.clone(),
                &baseline.summary,
                &changed_piece_sites,
                baseline.summary.fold_count,
            )?;
            if !validation_is_clean(fast_summary.1)
                || fast_summary.0.fold_count > baseline.summary.fold_count
            {
                rejected += 1;
                fast_rejected += 1;
                continue;
            }
            exact_evaluated += 1;
            let candidate = exact_summary_candidate(tensor_border, candidate)?;
            if !validation_is_clean(candidate.validation)
                || candidate.summary.fold_count > baseline.summary.fold_count
            {
                rejected += 1;
                exact_rejected += 1;
                continue;
            }
            let lift = lift_candidate_from_support_image_macro_patches(
                star_point,
                tensor_border,
                candidate.macro_patches.clone(),
            )?;
            if !validation_is_clean(lift.validation) {
                rejected += 1;
                continue;
            }
            let candidate_lift_summary = (lift.summary.clone(), lift.validation);
            let best_lift_summary = best
                .as_ref()
                .map(|candidate| (candidate.lift.summary.clone(), candidate.lift.validation));
            if is_better_summary(&candidate_lift_summary, best_lift_summary.as_ref()) {
                targeted_accepted += 1;
                best = Some(Kp20SupportImageSearchCandidate {
                    image: candidate,
                    lift,
                });
            } else {
                rejected += 1;
            }
        }
    }
    let composed_lift_ms = elapsed_ms(lift_start);
    best.map(|mut candidate| {
        candidate
            .image
            .summary
            .fold_details
            .push(format!("map_params={}", candidate.image.params.label()));
        candidate.image.summary.fold_details.push(format!(
            "map_search=evaluated:{evaluated},fast_evaluated:{fast_evaluated},exact_evaluated:{exact_evaluated},composed:{composed_evaluated},lift_budget:{COUPLED_LIFT_FINALIST_BUDGET},lift_budget_skipped:{lift_budget_skipped},targeted_evaluated:{targeted_evaluated},targeted_accepted:{targeted_accepted},accepted:1,rejected:{rejected},fast_rejected:{fast_rejected},exact_rejected:{exact_rejected},image_ms:{image_search_ms},lift_ms:{composed_lift_ms}"
        ));
        candidate.lift.summary.fold_details.push(format!(
            "map_params={};map_search=evaluated:{evaluated},fast_evaluated:{fast_evaluated},exact_evaluated:{exact_evaluated},composed:{composed_evaluated},lift_budget:{COUPLED_LIFT_FINALIST_BUDGET},lift_budget_skipped:{lift_budget_skipped},targeted_evaluated:{targeted_evaluated},targeted_accepted:{targeted_accepted},accepted:1,rejected:{rejected},fast_rejected:{fast_rejected},exact_rejected:{exact_rejected},image_ms:{image_search_ms},lift_ms:{composed_lift_ms}",
            candidate.image.params.label()
        ));
        candidate
    })
    .ok_or_else(|| {
        "support-image reparameterization found no residual-clean improvement".to_string()
    })
}

const COUPLED_LIFT_FINALIST_BUDGET: usize = 1;

fn elapsed_ms(start: Instant) -> u128 {
    start.elapsed().as_millis()
}

fn compare_summaries(a: &Kp20SamplerFoldSummary, b: &Kp20SamplerFoldSummary) -> std::cmp::Ordering {
    a.fold_count.cmp(&b.fold_count).then_with(|| {
        b.min_jacobian_normal_cos
            .partial_cmp(&a.min_jacobian_normal_cos)
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}

fn candidate_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SupportImageParams,
) -> Result<Kp20SupportImageCandidate, String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let (summary, validation) =
        candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
    Ok(Kp20SupportImageCandidate {
        summary,
        validation,
        macro_patches: cap.macro_patches,
        params,
    })
}

fn assembled_candidate_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SupportImageParams,
) -> Result<Kp20AssembledSupportImageCandidate, String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    Ok(Kp20AssembledSupportImageCandidate {
        macro_patches: cap.macro_patches,
        params,
    })
}

fn exact_summary_candidate(
    tensor_border: &TensorBorderBi3,
    candidate: Kp20AssembledSupportImageCandidate,
) -> Result<Kp20SupportImageCandidate, String> {
    let (summary, validation) =
        candidate_summary_from_cap(tensor_border, candidate.macro_patches.clone())?;
    Ok(Kp20SupportImageCandidate {
        summary,
        validation,
        macro_patches: candidate.macro_patches,
        params: candidate.params,
    })
}

fn support_image_piece_maps(
    mut maps: Kp20SectorPieceMaps,
    sector: usize,
    params: Kp20SupportImageParams,
    fold_sites: &[Kp20FoldSite],
) -> (Kp20SectorPieceMaps, Vec<(usize, Kp20PieceId)>) {
    let mut changed_sites = Vec::new();
    for (piece, corner_maps) in &mut maps {
        if let Some(site) = worst_site_for_piece(fold_sites, sector, *piece) {
            *corner_maps = params.reparameterized_maps(*corner_maps, site);
            changed_sites.push((sector, *piece));
        }
    }
    (maps, changed_sites)
}

fn targeted_support_image_piece_maps(
    mut maps: Kp20SectorPieceMaps,
    sector: usize,
    params: Kp20SupportImageParams,
    fold_sites: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> (Kp20SectorPieceMaps, Vec<(usize, Kp20PieceId)>) {
    let mut changed_sites = Vec::new();
    for (piece, corner_maps) in &mut maps {
        if *piece != target_piece {
            continue;
        }
        if let Some(site) = worst_site_for_piece(fold_sites, sector, *piece) {
            *corner_maps = params.reparameterized_maps(*corner_maps, site);
            changed_sites.push((sector, *piece));
        }
    }
    (maps, changed_sites)
}

fn worst_site_for_piece(
    fold_sites: &[Kp20FoldSite],
    sector: usize,
    piece: Kp20PieceId,
) -> Option<Kp20FoldSite> {
    fold_sites
        .iter()
        .copied()
        .filter(|site| site.sector == sector && site.piece == piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}
