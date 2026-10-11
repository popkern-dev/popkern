//! Fold-site-local interior-field refinements for KP20 lifted candidates.
#![allow(clippy::needless_range_loop)] // coefficient-grid index algebra (PK-795)
use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};
use popkern_math::PortableMath;

use super::super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::origin_inner_row::{
    Kp20InnerRowReport, Kp20OriginInnerMode, Kp20OriginInnerRowReport, Kp20SupportInnerRowReport,
    origin_inner_direction, sector_span,
};

type Kp20InnerFieldSelection = (
    Vec<Kp20MacroPatch>,
    (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    Kp20OriginInnerMode,
    f64,
);

pub(super) fn refine_origin_corner_field(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20OriginInnerRowReport), String> {
    refine_field(
        star_point,
        tensor_border,
        lift,
        "origin_corner_none",
        &[-0.30, -0.20, -0.12, -0.06, 0.06, 0.12, 0.20, 0.30],
        Kp20InnerFieldKind::OriginCorner,
    )
}

pub(super) fn refine_constrained_interior_solve(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20OriginInnerRowReport), String> {
    refine_field(
        star_point,
        tensor_border,
        lift,
        "constrained_interior_solve_none",
        &[-0.20, -0.12, -0.06, 0.06, 0.12, 0.20],
        Kp20InnerFieldKind::ConstrainedInteriorSolve,
    )
}

pub(super) fn refine_coupled_origin_support_field(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20OriginInnerRowReport), String> {
    refine_field(
        star_point,
        tensor_border,
        lift,
        "coupled_origin_support_none",
        &[-0.24, -0.16, -0.08, 0.08, 0.16, 0.24],
        Kp20InnerFieldKind::CoupledOriginSupport,
    )
}

pub(super) fn refine_support_interior_field(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20SupportInnerRowReport), String> {
    refine_field(
        star_point,
        tensor_border,
        lift,
        "none",
        &[-0.24, -0.16, -0.08, -0.04, 0.04, 0.08, 0.16, 0.24],
        Kp20InnerFieldKind::SupportInterior,
    )
}

fn refine_field(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    lift: Kp20SupportNormalLiftReport,
    decline_mode: &'static str,
    scales: &[f64],
    kind: Kp20InnerFieldKind,
) -> Result<(Kp20SupportNormalLiftReport, Kp20InnerRowReport), String> {
    let start_folds = lift.summary.fold_count;
    let mut stats = Kp20InnerRowReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds,
        end_folds: start_folds,
        mode: decline_mode,
        scale: 0.0,
        pre_evaluated: 0,
        pre_accepted: 0,
        pre_rejected: 0,
        pre_mode: "none",
    };
    let mut best: Option<Kp20InnerFieldSelection> = None;
    for mode in Kp20OriginInnerMode::all() {
        for &scale in scales {
            let mut candidate = lift.macro_patches.clone();
            let changed = apply_field_update(
                &mut candidate,
                star_point,
                tensor_border,
                mode,
                scale,
                &lift.summary.fold_sites,
                kind,
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
    stats.mode = kind.label(mode);
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
enum Kp20InnerFieldKind {
    ConstrainedInteriorSolve,
    CoupledOriginSupport,
    OriginCorner,
    SupportInterior,
}

impl Kp20InnerFieldKind {
    fn label(self, mode: Kp20OriginInnerMode) -> &'static str {
        match self {
            Self::ConstrainedInteriorSolve => mode.constrained_solve_label(),
            Self::CoupledOriginSupport => mode.coupled_field_label(),
            Self::OriginCorner => mode.origin_corner_label(),
            Self::SupportInterior => mode.field_label(),
        }
    }
}

fn apply_field_update(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    mode: Kp20OriginInnerMode,
    scale: f64,
    fold_sites: &[Kp20FoldSite],
    kind: Kp20InnerFieldKind,
) -> Vec<[bool; 4]> {
    match kind {
        Kp20InnerFieldKind::ConstrainedInteriorSolve => apply_constrained_interior_solve(
            macro_patches,
            star_point,
            tensor_border,
            mode,
            scale,
            fold_sites,
        ),
        Kp20InnerFieldKind::CoupledOriginSupport => apply_coupled_origin_support_update(
            macro_patches,
            star_point,
            tensor_border,
            mode,
            scale,
            fold_sites,
        ),
        Kp20InnerFieldKind::OriginCorner => apply_origin_corner_field_update(
            macro_patches,
            star_point,
            tensor_border,
            mode,
            scale,
            fold_sites,
        ),
        Kp20InnerFieldKind::SupportInterior => apply_support_field_update(
            macro_patches,
            star_point,
            tensor_border,
            mode,
            scale,
            fold_sites,
        ),
    }
}

fn apply_constrained_interior_solve(
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
        for piece_id in Kp20PieceId::all() {
            let sites = fold_sites
                .iter()
                .filter(|site| site.sector == sector && site.piece == piece_id)
                .copied()
                .collect::<Vec<_>>();
            if sites.is_empty() {
                continue;
            }
            let piece = macro_patch.piece_mut(piece_id);
            for i in 1..=5 {
                for j in 1..=5 {
                    let Some(displacement) =
                        constrained_control_delta(i, j, direction, scale * span, &sites)
                    else {
                        continue;
                    };
                    piece.cps[i][j] = piece.cps[i][j] + displacement;
                    changed[sector][piece_id.index()] = true;
                }
            }
        }
    }
    changed
}

fn constrained_control_delta(
    i: usize,
    j: usize,
    direction: Vec3,
    scale: f64,
    sites: &[Kp20FoldSite],
) -> Option<Vec3> {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let mut numerator = 0.0;
    let mut denominator = 2.0e-2;
    for site in sites {
        let du = u - site.u;
        let dv = v - site.v;
        let r2 = du * du + dv * dv;
        let radius2 = 0.42_f64 * 0.42_f64;
        let basis = (1.0 - r2 / radius2).max(0.0).powi(2);
        if basis == 0.0 {
            continue;
        }
        let severity = (1.0 - site.min_jacobian_normal_cos).max(0.0);
        numerator += basis * severity;
        denominator += basis * basis;
    }
    if numerator == 0.0 {
        None
    } else {
        let boundary_fade = (core::f64::consts::PI * u).psin() * (core::f64::consts::PI * v).psin();
        Some(direction * (scale * boundary_fade.max(0.0) * numerator / denominator))
    }
}

fn apply_coupled_origin_support_update(
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
        let origin_sites = fold_sites
            .iter()
            .filter(|site| {
                site.sector == sector
                    && matches!(
                        site.piece,
                        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight
                    )
            })
            .copied()
            .collect::<Vec<_>>();
        let support_sites = fold_sites
            .iter()
            .filter(|site| site.sector == sector && site.piece == Kp20PieceId::NearMidpointLeft)
            .copied()
            .collect::<Vec<_>>();
        if origin_sites.is_empty() && support_sites.is_empty() {
            continue;
        }
        for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            let sites = origin_sites
                .iter()
                .filter(|site| site.piece == piece_id)
                .copied()
                .collect::<Vec<_>>();
            if sites.is_empty() {
                continue;
            }
            let piece = macro_patch.piece_mut(piece_id);
            for i in 1..=3 {
                for j in 1..=3 {
                    let weight = origin_corner_weight(i, j, &sites);
                    if weight > 0.0 {
                        piece.cps[i][j] = piece.cps[i][j] + direction * (scale * span * weight);
                        changed[sector][piece_id.index()] = true;
                    }
                }
            }
        }
        if !support_sites.is_empty() {
            let piece = macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft);
            for i in 1..=5 {
                for j in 1..=5 {
                    let weight = 0.65 * support_field_weight(i, j, &support_sites);
                    if weight > 0.0 {
                        piece.cps[i][j] = piece.cps[i][j] + direction * (-scale * span * weight);
                        changed[sector][Kp20PieceId::NearMidpointLeft.index()] = true;
                    }
                }
            }
        }
    }
    changed
}

fn apply_origin_corner_field_update(
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
        for piece_id in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            let sites = fold_sites
                .iter()
                .filter(|site| site.sector == sector && site.piece == piece_id)
                .copied()
                .collect::<Vec<_>>();
            if sites.is_empty() {
                continue;
            }
            let piece = macro_patch.piece_mut(piece_id);
            for i in 1..=3 {
                for j in 1..=3 {
                    let weight = origin_corner_weight(i, j, &sites);
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

fn apply_support_field_update(
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
            .filter(|site| site.sector == sector && site.piece == Kp20PieceId::NearMidpointLeft)
            .copied()
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let piece = macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft);
        for i in 1..=5 {
            for j in 1..=5 {
                let weight = support_field_weight(i, j, &sites);
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

fn origin_corner_weight(i: usize, j: usize, sites: &[Kp20FoldSite]) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let locality = sites
        .iter()
        .map(|site| {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.30_f64 * 0.30_f64;
            let locality = (1.0 - r2 / radius2).max(0.0).powi(2);
            let severity = (1.0 - site.min_jacobian_normal_cos).max(0.0);
            locality * severity
        })
        .fold(0.0, f64::max);
    if locality == 0.0 {
        0.0
    } else {
        let corner_fade = ((1.0 - u).max(0.0) * (1.0 - v).max(0.0)).sqrt();
        (locality * corner_fade).clamp(0.0, 1.0)
    }
}

fn support_field_weight(i: usize, j: usize, sites: &[Kp20FoldSite]) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let fold_locality = sites
        .iter()
        .map(|site| {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.48_f64 * 0.48_f64;
            let locality = (1.0 - r2 / radius2).max(0.0).powi(2);
            let severity = (1.0 - site.min_jacobian_normal_cos).max(0.0);
            locality * severity
        })
        .fold(0.0, f64::max);
    if fold_locality == 0.0 {
        0.0
    } else {
        let boundary_fade = (core::f64::consts::PI * u).psin() * (core::f64::consts::PI * v).psin();
        (fold_locality * boundary_fade.max(0.0)).clamp(0.0, 1.0)
    }
}
