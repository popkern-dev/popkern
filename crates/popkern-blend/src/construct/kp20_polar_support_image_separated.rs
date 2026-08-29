//! Fold-site separated KP20 support-image reparameterization candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    Kp20PieceId, Kp20SectorPieceMaps, TensorBorderBi3, assemble_kp20_bi6_macro_cap,
    build_kp20_guide_ring, kp20_positive_nonoverlap_c2_sector_piece_maps,
    sample_kp20_macro_patch_from_piece_maps,
};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, fast_candidate_summary_from_changed_piece_sites,
    is_better_summary, validation_is_clean,
};
use super::kp20_polar_support_image_params::{
    Kp20SeparatedSupportImageParams, Kp20SupportImageParams,
};
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_support_image_macro_patches,
};
use super::kp20_polar_support_reparameterization::Kp20SupportReparameterizationReport;
use super::kp20_polar_support_surface_split::support_owner_piece;

type FoldSummary = (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport);
type SeparatedReparamResult = Result<
    (
        Kp20SupportNormalLiftReport,
        Kp20SupportReparameterizationReport,
    ),
    String,
>;

struct Kp20SeparatedImageCandidate {
    macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SeparatedSupportImageParams,
}

pub(super) fn refine_fold_site_separated_reparameterization(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> SeparatedReparamResult {
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
        params: "fold_site_separated:none".to_string(),
        decline_reason: "fold_site_separated_reparameterization_no_candidate",
    };
    if !has_support_owned_fold_site(&baseline.summary.fold_sites, target_piece) {
        report.decline_reason = "fold_site_separated_reparameterization_no_owned_fold_site";
        return Ok((baseline, report));
    }
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let baseline_summary = (baseline.summary.clone(), baseline.validation);
    let mut best_fast: Option<(Kp20SeparatedImageCandidate, FoldSummary)> = None;
    for params in Kp20SupportImageParams::fold_site_separated_family() {
        report.evaluated += 1;
        let mut changed_piece_sites = Vec::new();
        let sampled_macro_patches = guide
            .sectors
            .iter()
            .enumerate()
            .map(|(sector, guide_sector)| {
                let (maps, changed_sites) = separated_image_piece_maps(
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
            separated_candidate_from_sampled(tensor_border, sampled_macro_patches, params)?;
        let fast = fast_candidate_summary_from_changed_piece_sites(
            tensor_border,
            candidate.macro_patches.clone(),
            &baseline.summary,
            &changed_piece_sites,
            baseline.summary.fold_count,
        )?;
        if !validation_is_clean(fast.1) || fast.0.fold_count >= baseline.summary.fold_count {
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
        report.decline_reason = "fold_site_separated_reparameterization_lift_validation_not_clean";
        return Ok((baseline, report));
    }
    let lift_summary = (lift.summary.clone(), lift.validation);
    if !is_better_summary(&lift_summary, Some(&baseline_summary)) {
        report.rejected = report.evaluated;
        report.decline_reason = "fold_site_separated_reparameterization_no_lift_improvement";
        return Ok((baseline, report));
    }
    report.end_folds = lift.summary.fold_count;
    report.params = params.label();
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.decline_reason = "none";
    Ok((lift, report))
}

fn separated_candidate_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SeparatedSupportImageParams,
) -> Result<Kp20SeparatedImageCandidate, String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    Ok(Kp20SeparatedImageCandidate {
        macro_patches: cap.macro_patches,
        params,
    })
}

fn separated_image_piece_maps(
    mut maps: Kp20SectorPieceMaps,
    sector: usize,
    params: Kp20SeparatedSupportImageParams,
    fold_sites: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> (Kp20SectorPieceMaps, Vec<(usize, Kp20PieceId)>) {
    let mut changed_sites = Vec::new();
    for (piece, corner_maps) in &mut maps {
        if *piece != target_piece {
            continue;
        }
        if let Some(site) = worst_site_for_support_owner(fold_sites, sector, target_piece) {
            *corner_maps = params.reparameterized_maps_for_piece(*corner_maps, target_piece, site);
            changed_sites.push((sector, *piece));
        }
    }
    (maps, changed_sites)
}

fn has_support_owned_fold_site(fold_sites: &[Kp20FoldSite], target_piece: Kp20PieceId) -> bool {
    fold_sites
        .iter()
        .any(|site| support_owner_piece(site.piece) == target_piece)
}

fn worst_site_for_support_owner(
    fold_sites: &[Kp20FoldSite],
    sector: usize,
    target_piece: Kp20PieceId,
) -> Option<Kp20FoldSite> {
    fold_sites
        .iter()
        .copied()
        .filter(|site| site.sector == sector && support_owner_piece(site.piece) == target_piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}
