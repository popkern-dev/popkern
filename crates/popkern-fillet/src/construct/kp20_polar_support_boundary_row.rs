//! Support boundary-adjacent row refinement for KP20 lifted candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};

use super::super::super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::{Kp20InnerRowReport, Kp20OriginInnerMode, origin_inner_direction, sector_span};

pub(super) fn refine_support_boundary_adjacent_row(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20InnerRowReport), String> {
    let start_folds = lift.summary.fold_count;
    let mut best: Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20OriginInnerMode,
        f64,
    )> = None;
    let mut stats = Kp20InnerRowReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds,
        end_folds: start_folds,
        mode: "support_boundary_row_none",
        scale: 0.0,
        pre_evaluated: 0,
        pre_accepted: 0,
        pre_rejected: 0,
        pre_mode: "none",
    };
    for mode in Kp20OriginInnerMode::all() {
        for &scale in &[-0.32, -0.24, -0.16, -0.08, 0.08, 0.16, 0.24, 0.32] {
            let mut candidate = lift.macro_patches.clone();
            let changed = apply_support_boundary_row_update(
                &mut candidate,
                star_point,
                tensor_border,
                mode,
                scale,
                &lift.summary.fold_sites,
            );
            stats.evaluated += 1;
            let fast = fast_candidate_summary_from_changed_patches(
                tensor_border,
                candidate.clone(),
                &lift.summary,
                &changed,
                lift.summary.fold_count,
            )?;
            if !validation_is_clean(fast.1)
                || !is_better_summary(&fast, Some(&(lift.summary.clone(), lift.validation)))
            {
                stats.rejected += 1;
                continue;
            }
            let exact = candidate_summary_from_changed_patches(
                tensor_border,
                candidate.clone(),
                &lift.summary,
                &changed,
                lift.summary.fold_count,
            )?;
            if !validation_is_clean(exact.1)
                || !is_better_summary(&exact, Some(&(lift.summary.clone(), lift.validation)))
            {
                stats.rejected += 1;
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, summary, _, _)| is_better_summary(&exact, Some(summary)))
            {
                best = Some((candidate, exact, mode, scale));
            } else {
                stats.rejected += 1;
            }
        }
    }
    let Some((macro_patches, summary, mode, scale)) = best else {
        return Ok((lift, stats));
    };
    stats.accepted = 1;
    stats.end_folds = summary.0.fold_count;
    stats.mode = support_boundary_row_label(mode);
    stats.scale = scale;
    let mut optimizer_stats = lift.optimizer_stats.clone();
    optimizer_stats.evaluated_candidates += stats.evaluated;
    optimizer_stats.accepted_candidates += stats.accepted;
    optimizer_stats.rejected_candidates += stats.rejected;
    Ok((
        Kp20SupportNormalLiftReport {
            summary: summary.0,
            validation: summary.1,
            macro_patches,
            accepted_steps: lift.accepted_steps,
            optimizer_stats,
        },
        stats,
    ))
}

fn apply_support_boundary_row_update(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    mode: Kp20OriginInnerMode,
    scale: f64,
    fold_sites: &[Kp20FoldSite],
) -> Vec<[bool; 4]> {
    let mut changed = vec![[false; 4]; macro_patches.len()];
    for (sector, macro_patch) in macro_patches.iter_mut().enumerate() {
        let direction = origin_inner_direction(star_point, tensor_border, sector, mode);
        let span = sector_span(star_point, tensor_border, sector);
        let sites = fold_sites
            .iter()
            .filter(|site| {
                site.sector == sector
                    && site.piece == Kp20PieceId::NearMidpointLeft
                    && site.v <= 0.25
            })
            .copied()
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let piece = macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft);
        for i in 1..=5 {
            for j in 1..=2 {
                let weight = support_boundary_row_weight(i, j, &sites);
                if weight == 0.0 {
                    continue;
                }
                piece.cps[i][j] = piece.cps[i][j] + direction * (scale * span * weight);
                changed[sector][Kp20PieceId::NearMidpointLeft.index()] = true;
            }
        }
    }
    changed
}

fn support_boundary_row_label(mode: Kp20OriginInnerMode) -> &'static str {
    match mode {
        Kp20OriginInnerMode::Normal => "support_boundary_row_normal",
        Kp20OriginInnerMode::Radial => "support_boundary_row_radial",
        Kp20OriginInnerMode::Chord => "support_boundary_row_chord",
    }
}

fn support_boundary_row_weight(i: usize, j: usize, sites: &[Kp20FoldSite]) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let row_fade = if j == 1 { 1.0 } else { 0.35 };
    let boundary_fade = (core::f64::consts::PI * u).sin().max(0.0);
    sites
        .iter()
        .map(|site| {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.32_f64 * 0.32_f64;
            let locality = (1.0 - r2 / radius2).max(0.0).powi(2);
            let severity = (1.0_f64 - site.min_jacobian_normal_cos).max(0.0);
            locality * severity * row_fade * boundary_fade
        })
        .fold(0.0, f64::max)
        .clamp(0.0, 1.0)
}
