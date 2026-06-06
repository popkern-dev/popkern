//! Initial support-domain basis candidate for bounded high-valence KP20 lift.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    TensorBorderBi3, assemble_kp20_bi6_macro_cap, build_kp20_guide_ring,
    kp20_positive_nonoverlap_c2_sector_piece_maps, sample_kp20_macro_patch_from_piece_maps,
};

use super::super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, fast_candidate_summary_from_cap,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::super::super::kp20_polar_support_reparameterization::Kp20SupportDomainReport;
use super::TARGET_PIECE;
use super::guide_repair::transactional_kp20_guide_repair;
use super::maps::{Kp20SupportDomainParams, changed_target_pieces, support_domain_piece_maps};

pub(crate) fn support_domain_initial_basis_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Vec<Kp20MacroPatch>, Kp20SupportDomainReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let baseline_sampled = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_positive_nonoverlap_c2_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let baseline_cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches: baseline_sampled,
    })
    .map_err(|e| e.to_string())?;
    let baseline =
        fast_candidate_summary_from_cap(tensor_border, baseline_cap.macro_patches.clone())?;
    let (guide, _, _) =
        transactional_kp20_guide_repair(&guide, tensor_border, star_point, &baseline)?;
    let mut report = report_from_baseline(&baseline);
    let mut best_fast = None;
    let changed_patches = changed_target_pieces(tensor_border.sectors.len());
    for params in Kp20SupportDomainParams::family() {
        report.evaluated += 1;
        let sampled = guide
            .sectors
            .iter()
            .map(|sector| {
                sample_kp20_macro_patch_from_piece_maps(sector, support_domain_piece_maps(params))
            })
            .collect::<Vec<_>>();
        let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border: tensor_border.clone(),
            sampled_macro_patches: sampled,
        })
        .map_err(|e| e.to_string())?;
        let fast = fast_candidate_summary_from_changed_patches(
            tensor_border,
            cap.macro_patches.clone(),
            &baseline.0,
            &changed_patches,
            baseline.0.fold_count,
        )?;
        if !validation_is_clean(fast.1) || !is_better_summary(&fast, Some(&baseline)) {
            report.rejected += 1;
            report.decline_reason = "support_domain_initial_fast_declined";
            continue;
        }
        if best_fast
            .as_ref()
            .is_none_or(|(_, best_summary, _)| is_better_summary(&fast, Some(best_summary)))
        {
            best_fast = Some((cap.macro_patches, fast, params));
        } else {
            report.rejected += 1;
        }
    }
    let Some((macro_patches, _, params)) = best_fast else {
        report.rejected = report.evaluated;
        return Ok((baseline_cap.macro_patches, report));
    };
    let exact = candidate_summary_from_cap(tensor_border, macro_patches.clone())?;
    if !validation_is_clean(exact.1) || !is_better_summary(&exact, Some(&baseline)) {
        report.rejected = report.evaluated;
        report.decline_reason = "support_domain_initial_exact_declined";
        return Ok((baseline_cap.macro_patches, report));
    }
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = exact.0.fold_count;
    report.origin_u = params.origin[0];
    report.origin_v = params.origin[1];
    report.span_u = params.span_u;
    report.span_v = params.span_v;
    report.skew = params.skew;
    report.decline_reason = "none";
    Ok((macro_patches, report))
}

fn report_from_baseline(
    baseline: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
) -> Kp20SupportDomainReport {
    Kp20SupportDomainReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.0.fold_count,
        end_folds: baseline.0.fold_count,
        target_piece: TARGET_PIECE,
        origin_u: 0.0,
        origin_v: 0.0,
        span_u: 0.0,
        span_v: 0.0,
        skew: 0.0,
        decline_reason: "support_domain_no_candidate",
    }
}
