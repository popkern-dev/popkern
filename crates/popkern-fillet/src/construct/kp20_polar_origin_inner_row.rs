//! Bounded origin-owned interior-row refinement for KP20 lifted candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};
use popkern_core::geometry::{Point3, Vec3};

use super::super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20OriginInnerRowReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) mode: &'static str,
    pub(super) scale: f64,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20PairedInnerRowReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) mode: &'static str,
    pub(super) scale: f64,
}

pub(super) fn refine_origin_inner_interior_row(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20OriginInnerRowReport), String> {
    let start_folds = lift.summary.fold_count;
    let mut best: Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20OriginInnerMode,
        f64,
    )> = None;
    let mut stats = Kp20OriginInnerRowReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds,
        end_folds: start_folds,
        mode: "none",
        scale: 0.0,
    };
    for mode in Kp20OriginInnerMode::all() {
        for scale in [-0.20, -0.12, 0.12, 0.20] {
            let mut candidate = lift.macro_patches.clone();
            let changed = apply_origin_inner_row_update(
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
    stats.mode = mode.label();
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

pub(super) fn refine_paired_inner_interior_row(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20PairedInnerRowReport), String> {
    let start_folds = lift.summary.fold_count;
    let mut best: Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20OriginInnerMode,
        f64,
    )> = None;
    let mut stats = Kp20PairedInnerRowReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds,
        end_folds: start_folds,
        mode: "none",
        scale: 0.0,
    };
    for mode in Kp20OriginInnerMode::all() {
        for scale in [-0.16, -0.08, 0.08, 0.16] {
            let mut candidate = lift.macro_patches.clone();
            let changed = apply_inner_row_update(
                &mut candidate,
                star_point,
                tensor_border,
                mode,
                scale,
                &lift.summary.fold_sites,
                &Kp20PieceId::all(),
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
    stats.mode = mode.label();
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

#[derive(Debug, Clone, Copy)]
enum Kp20OriginInnerMode {
    Normal,
    Radial,
    Chord,
}

impl Kp20OriginInnerMode {
    const fn all() -> [Self; 3] {
        [Self::Normal, Self::Radial, Self::Chord]
    }

    const fn label(self) -> &'static str {
        match self {
            Self::Normal => "normal",
            Self::Radial => "radial",
            Self::Chord => "chord",
        }
    }
}

fn apply_origin_inner_row_update(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    mode: Kp20OriginInnerMode,
    scale: f64,
    fold_sites: &[Kp20FoldSite],
) -> Vec<[bool; 4]> {
    apply_inner_row_update(
        macro_patches,
        star_point,
        tensor_border,
        mode,
        scale,
        fold_sites,
        &[Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight],
    )
}

fn apply_inner_row_update(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    mode: Kp20OriginInnerMode,
    scale: f64,
    fold_sites: &[Kp20FoldSite],
    target_pieces: &[Kp20PieceId],
) -> Vec<[bool; 4]> {
    let mut changed = vec![[false; 4]; macro_patches.len()];
    for (sector, macro_patch) in macro_patches.iter_mut().enumerate() {
        let direction = origin_inner_direction(star_point, tensor_border, sector, mode);
        let span = sector_span(star_point, tensor_border, sector);
        for &piece_id in target_pieces {
            let sites = fold_sites
                .iter()
                .filter(|site| {
                    site.sector == sector && site.piece == piece_id && interior_site(site)
                })
                .copied()
                .collect::<Vec<_>>();
            if sites.is_empty() {
                continue;
            }
            let piece = macro_patch.piece_mut(piece_id);
            for i in 2..=4 {
                for j in 2..=4 {
                    let weight = origin_inner_weight(i, j, &sites);
                    if weight == 0.0 {
                        continue;
                    }
                    piece.cps[i][j] = piece.cps[i][j] + direction * (scale * span * weight);
                    changed[sector][piece_id.index()] = true;
                }
            }
        }
    }
    changed
}

fn origin_inner_direction(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    sector: usize,
    mode: Kp20OriginInnerMode,
) -> Vec3 {
    let row = tensor_border.sectors[sector][0];
    let spoke_a = row[0] - star_point;
    let spoke_b = row[3] - star_point;
    let normal = spoke_a
        .cross(&spoke_b)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    match mode {
        Kp20OriginInnerMode::Normal => normal,
        Kp20OriginInnerMode::Radial => (spoke_a + spoke_b).normalized().unwrap_or(normal),
        Kp20OriginInnerMode::Chord => (spoke_b - spoke_a).normalized().unwrap_or(normal),
    }
}

fn sector_span(star_point: Point3, tensor_border: &TensorBorderBi3, sector: usize) -> f64 {
    let row = tensor_border.sectors[sector][0];
    0.5 * ((row[0] - star_point).length() + (row[3] - star_point).length())
}

fn interior_site(site: &Kp20FoldSite) -> bool {
    site.u < 0.80 && site.v < 0.80
}

fn origin_inner_weight(i: usize, j: usize, sites: &[Kp20FoldSite]) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    sites
        .iter()
        .map(|site| {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.35_f64 * 0.35_f64;
            (1.0 - r2 / radius2).max(0.0).powi(2)
        })
        .fold(0.0, f64::max)
}
