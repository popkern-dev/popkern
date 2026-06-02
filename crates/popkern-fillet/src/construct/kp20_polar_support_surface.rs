//! Support-surface-owned KP20 polar candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6Patch, Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};
use popkern_core::geometry::{Point3, Vec3};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::kp20_polar_multipatch_evaluate::{
    Kp20MultiPatchEvaluationReport, evaluate_multipatch_candidate,
};
use super::kp20_polar_multipatch_inner::{
    Kp20MultiPatchInnerBuildReport, build_multipatch_inner_patches,
};
use super::kp20_polar_multipatch_interface::{
    Kp20MultiPatchInterfaceBuildReport, build_multipatch_interfaces,
};
use super::kp20_polar_multipatch_inventory::{
    Kp20MultiPatchInventoryReport, multipatch_inventory_from_fold_summary,
};
use super::kp20_polar_multipatch_strip::{
    Kp20MultiPatchStripBuildReport, build_multipatch_strip_surfaces,
};
use super::kp20_polar_multipatch_validate::{
    Kp20MultiPatchValidationReport, validate_multipatch_candidate,
};
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::kp20_polar_support_surface_split::{
    evaluate_bi6_normal, evaluate_bi6_patch, propagate_support_pole_interfaces,
    support_owner_piece, support_side_index,
};

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceReport {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) stats: Kp20SupportSurfaceStats,
    pub(super) stage_trace: Vec<Kp20SupportSurfaceStageReport>,
    pub(super) multipatch_inventory: Kp20MultiPatchInventoryReport,
    pub(super) multipatch_strips: Kp20MultiPatchStripBuildReport,
    pub(super) multipatch_interfaces: Kp20MultiPatchInterfaceBuildReport,
    pub(super) multipatch_inner: Kp20MultiPatchInnerBuildReport,
    pub(super) multipatch_validation: Kp20MultiPatchValidationReport,
    pub(super) multipatch_evaluation: Kp20MultiPatchEvaluationReport,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceStats {
    pub(super) active_fold_pieces: usize,
    pub(super) orientation_rows: usize,
    pub(super) evaluated_candidates: usize,
    pub(super) exact_evaluated_candidates: usize,
    pub(super) accepted_candidates: usize,
    pub(super) rejected_candidates: usize,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceStageReport {
    pub(super) mode: &'static str,
    pub(super) evaluated_candidates: usize,
    pub(super) fast_finalist_folds: Option<usize>,
    pub(super) exact_evaluated: bool,
    pub(super) exact_folds: Option<usize>,
    pub(super) accepted: bool,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SupportSurfaceTrial {
    mode: Kp20SupportSurfaceMode,
    scale: f64,
    piece: Kp20PieceId,
}

#[derive(Debug, Clone, Copy)]
enum Kp20SupportSurfaceMode {
    SupportNormal,
    RadialBisector,
    SectorChord,
    CoupledSupportRing,
    FoldAwareSurfaceFrame,
    FoldAwareNormalFrame,
    SupportPoleInterfaceSplit,
}

pub(super) fn support_surface_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    support_lift: &Kp20SupportNormalLiftReport,
) -> Result<Kp20SupportSurfaceReport, String> {
    let active_fold_pieces = active_fold_piece_count(&support_lift.summary);
    let mut accepted_patches = support_lift.macro_patches.clone();
    let mut accepted_summary = (support_lift.summary.clone(), support_lift.validation);
    let mut stats = Kp20SupportSurfaceStats {
        active_fold_pieces,
        orientation_rows: active_fold_pieces,
        evaluated_candidates: 0,
        exact_evaluated_candidates: 0,
        accepted_candidates: 0,
        rejected_candidates: 0,
    };
    let mut stage_trace = Vec::new();
    for mode in Kp20SupportSurfaceMode::stages() {
        let mut stage = Kp20SupportSurfaceStageReport {
            mode: mode.label(),
            evaluated_candidates: 0,
            fast_finalist_folds: None,
            exact_evaluated: false,
            exact_folds: None,
            accepted: false,
        };
        let mut best_fast: Option<(
            Vec<Kp20MacroPatch>,
            Vec<[bool; 4]>,
            (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        )> = None;
        for trial in support_surface_trials(mode, &accepted_summary.0) {
            stats.evaluated_candidates += 1;
            stage.evaluated_candidates += 1;
            let mut candidate = accepted_patches.clone();
            let changed_patches = apply_orientation_row_trial(
                &mut candidate,
                star_point,
                tensor_border,
                &accepted_summary.0.fold_sites,
                trial,
            );
            let fast_summary = fast_candidate_summary_from_changed_patches(
                tensor_border,
                candidate.clone(),
                &accepted_summary.0,
                &changed_patches,
                accepted_summary.0.fold_count,
            )?;
            if !validation_is_clean(fast_summary.1)
                || !is_fold_reducing_summary(&fast_summary, &accepted_summary)
            {
                stats.rejected_candidates += 1;
                continue;
            }
            if best_fast
                .as_ref()
                .is_some_and(|(_, _, best)| !is_better_summary(&fast_summary, Some(best)))
            {
                stats.rejected_candidates += 1;
                continue;
            }
            best_fast = Some((candidate, changed_patches, fast_summary));
        }
        let Some((candidate, changed_patches, fast_summary)) = best_fast else {
            stage_trace.push(stage);
            continue;
        };
        stage.fast_finalist_folds = Some(fast_summary.0.fold_count);
        stage.exact_evaluated = true;
        stats.exact_evaluated_candidates += 1;
        let summary = candidate_summary_from_changed_patches(
            tensor_border,
            candidate.clone(),
            &accepted_summary.0,
            &changed_patches,
            accepted_summary.0.fold_count,
        )?;
        stage.exact_folds = Some(summary.0.fold_count);
        if validation_is_clean(summary.1) && is_fold_reducing_summary(&summary, &accepted_summary) {
            accepted_patches = candidate;
            accepted_summary = summary;
            stats.accepted_candidates += 1;
            stage.accepted = true;
        } else {
            stats.rejected_candidates += 1;
        }
        stage_trace.push(stage);
    }
    let decline_reason = if stats.accepted_candidates > 0 {
        "support_surface_orientation_rows_accepted"
    } else {
        "support_surface_no_fold_reducing_orientation_row"
    };
    let multipatch_inventory = multipatch_inventory_from_fold_summary(&accepted_summary.0);
    let multipatch_strips =
        build_multipatch_strip_surfaces(&multipatch_inventory, &accepted_patches);
    let multipatch_interfaces = build_multipatch_interfaces(&multipatch_strips);
    let multipatch_inner =
        build_multipatch_inner_patches(&multipatch_inventory, &multipatch_interfaces);
    let multipatch_validation = validate_multipatch_candidate(
        &multipatch_inventory,
        &multipatch_strips,
        &multipatch_interfaces,
        &multipatch_inner,
    );
    let multipatch_evaluation = evaluate_multipatch_candidate(
        tensor_border,
        &accepted_patches,
        &accepted_summary.0,
        &multipatch_validation,
        &multipatch_inner,
    )?;
    Ok(Kp20SupportSurfaceReport {
        multipatch_inventory,
        multipatch_strips,
        multipatch_interfaces,
        multipatch_inner,
        multipatch_validation,
        multipatch_evaluation,
        summary: accepted_summary.0,
        validation: accepted_summary.1,
        stats,
        stage_trace,
        decline_reason,
    })
}

fn is_fold_reducing_summary(
    candidate: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    accepted: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
) -> bool {
    candidate.0.fold_count < accepted.0.fold_count
}

fn active_fold_piece_count(summary: &Kp20SamplerFoldSummary) -> usize {
    active_fold_pieces(summary).len()
}

fn active_fold_pieces(summary: &Kp20SamplerFoldSummary) -> Vec<Kp20PieceId> {
    let mut active = [false; 4];
    for site in &summary.fold_sites {
        active[site.piece.index()] = true;
    }
    Kp20PieceId::all()
        .into_iter()
        .filter(|piece| active[piece.index()])
        .collect()
}

impl Kp20SupportSurfaceMode {
    fn stages() -> [Self; 7] {
        [
            Self::SupportNormal,
            Self::RadialBisector,
            Self::SectorChord,
            Self::CoupledSupportRing,
            Self::FoldAwareSurfaceFrame,
            Self::FoldAwareNormalFrame,
            Self::SupportPoleInterfaceSplit,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::SupportNormal => "support_normal",
            Self::RadialBisector => "radial_bisector",
            Self::SectorChord => "sector_chord",
            Self::CoupledSupportRing => "coupled_support_ring",
            Self::FoldAwareSurfaceFrame => "fold_aware_surface_frame",
            Self::FoldAwareNormalFrame => "fold_aware_normal_frame",
            Self::SupportPoleInterfaceSplit => "support_pole_interface_split",
        }
    }
}

fn support_surface_trials(
    mode: Kp20SupportSurfaceMode,
    summary: &Kp20SamplerFoldSummary,
) -> Vec<Kp20SupportSurfaceTrial> {
    let mut trials = Vec::new();
    for scale in [0.06, -0.06] {
        for piece in active_trial_pieces(mode, summary) {
            trials.push(Kp20SupportSurfaceTrial { mode, scale, piece });
        }
    }
    trials
}

fn active_trial_pieces(
    mode: Kp20SupportSurfaceMode,
    summary: &Kp20SamplerFoldSummary,
) -> Vec<Kp20PieceId> {
    let mut pieces = active_fold_pieces(summary);
    if matches!(mode, Kp20SupportSurfaceMode::SupportPoleInterfaceSplit) {
        pieces = pieces.into_iter().map(support_owner_piece).collect();
        pieces.sort_by_key(|piece| piece.index());
        pieces.dedup();
    }
    pieces
}

fn apply_orientation_row_trial(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    fold_sites: &[Kp20FoldSite],
    trial: Kp20SupportSurfaceTrial,
) -> Vec<[bool; 4]> {
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    let mut split_support_changed = vec![[false; 2]; macro_patches.len()];
    for site in fold_sites
        .iter()
        .filter(|site| trial_affects_site(trial, **site))
        .copied()
    {
        let row = tensor_border.sectors[site.sector][0];
        let spoke_a = row[0] - star_point;
        let spoke_b = row[3] - star_point;
        let normal = spoke_a
            .cross(&spoke_b)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let direction = support_surface_direction(
            trial.mode,
            spoke_a,
            spoke_b,
            normal,
            star_point,
            site,
            macro_patches[site.sector].piece(site.piece),
        );
        let span = 0.5 * (spoke_a.length() + spoke_b.length());
        let target_piece = target_piece_for_trial(trial, site.piece);
        let piece = macro_patches[site.sector].piece_mut(target_piece);
        for i in 0..7 {
            for j in 0..7 {
                if !can_adjust_orientation_row(trial.mode, target_piece, i, j) {
                    continue;
                }
                let weight = orientation_row_weight(trial.mode, site, target_piece, i, j);
                if weight == 0.0 {
                    continue;
                }
                piece.cps[i][j] = piece.cps[i][j] + direction * (trial.scale * span * weight);
                changed_patches[site.sector][target_piece.index()] = true;
                if matches!(
                    trial.mode,
                    Kp20SupportSurfaceMode::SupportPoleInterfaceSplit
                ) {
                    split_support_changed[site.sector][support_side_index(target_piece)] = true;
                }
            }
        }
    }
    if matches!(
        trial.mode,
        Kp20SupportSurfaceMode::SupportPoleInterfaceSplit
    ) {
        propagate_support_pole_interfaces(
            macro_patches,
            &mut changed_patches,
            split_support_changed,
        );
    }
    changed_patches
}

fn trial_affects_site(trial: Kp20SupportSurfaceTrial, site: Kp20FoldSite) -> bool {
    if matches!(
        trial.mode,
        Kp20SupportSurfaceMode::SupportPoleInterfaceSplit
    ) {
        support_owner_piece(site.piece) == trial.piece
    } else {
        site.piece == trial.piece
    }
}

fn target_piece_for_trial(trial: Kp20SupportSurfaceTrial, site_piece: Kp20PieceId) -> Kp20PieceId {
    if matches!(
        trial.mode,
        Kp20SupportSurfaceMode::SupportPoleInterfaceSplit
    ) {
        support_owner_piece(site_piece)
    } else {
        site_piece
    }
}

fn support_surface_direction(
    mode: Kp20SupportSurfaceMode,
    spoke_a: Vec3,
    spoke_b: Vec3,
    normal: Vec3,
    star_point: Point3,
    site: Kp20FoldSite,
    piece: &Kp20Bi6Patch,
) -> Vec3 {
    match mode {
        Kp20SupportSurfaceMode::SupportNormal | Kp20SupportSurfaceMode::CoupledSupportRing => {
            normal
        }
        Kp20SupportSurfaceMode::RadialBisector => {
            (spoke_a + spoke_b).normalized().unwrap_or(normal)
        }
        Kp20SupportSurfaceMode::SectorChord => (spoke_b - spoke_a).normalized().unwrap_or(normal),
        Kp20SupportSurfaceMode::FoldAwareSurfaceFrame => {
            let point = evaluate_bi6_patch(piece, site.u, site.v);
            (point - star_point).normalized().unwrap_or(normal)
        }
        Kp20SupportSurfaceMode::FoldAwareNormalFrame
        | Kp20SupportSurfaceMode::SupportPoleInterfaceSplit => {
            evaluate_bi6_normal(piece, site.u, site.v)
                .map(|local| {
                    if local.dot(&normal) < 0.0 {
                        -local
                    } else {
                        local
                    }
                })
                .unwrap_or(normal)
        }
    }
}

fn can_adjust_orientation_row(
    mode: Kp20SupportSurfaceMode,
    piece: Kp20PieceId,
    i: usize,
    j: usize,
) -> bool {
    match piece {
        Kp20PieceId::NearMidpointLeft | Kp20PieceId::NearMidpointRight => {
            let row_ok = match mode {
                Kp20SupportSurfaceMode::CoupledSupportRing
                | Kp20SupportSurfaceMode::FoldAwareSurfaceFrame
                | Kp20SupportSurfaceMode::FoldAwareNormalFrame
                | Kp20SupportSurfaceMode::SupportPoleInterfaceSplit => (2..=4).contains(&i),
                _ => i == 3,
            };
            row_ok && (2..=4).contains(&j)
        }
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => {
            let row_ok = match mode {
                Kp20SupportSurfaceMode::CoupledSupportRing
                | Kp20SupportSurfaceMode::FoldAwareSurfaceFrame
                | Kp20SupportSurfaceMode::FoldAwareNormalFrame
                | Kp20SupportSurfaceMode::SupportPoleInterfaceSplit => (2..=4).contains(&j),
                _ => j == 3,
            };
            row_ok && (2..=4).contains(&i)
        }
    }
}

fn orientation_row_weight(
    mode: Kp20SupportSurfaceMode,
    site: Kp20FoldSite,
    target_piece: Kp20PieceId,
    i: usize,
    j: usize,
) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let du = u - site.u;
    let dv = v - site.v;
    let r2 = du * du + dv * dv;
    let locality = (1.0 - r2 / (0.42_f64 * 0.42_f64)).max(0.0);
    let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
    let row_weight = match mode {
        Kp20SupportSurfaceMode::CoupledSupportRing
        | Kp20SupportSurfaceMode::FoldAwareSurfaceFrame
        | Kp20SupportSurfaceMode::FoldAwareNormalFrame
        | Kp20SupportSurfaceMode::SupportPoleInterfaceSplit => match target_piece {
            Kp20PieceId::NearMidpointLeft | Kp20PieceId::NearMidpointRight => {
                coupled_row_weight(site, i)
            }
            Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => {
                coupled_row_weight(site, j)
            }
        },
        _ => 1.0,
    };
    locality * locality * severity.max(0.25) * row_weight
}

fn coupled_row_weight(site: Kp20FoldSite, row: usize) -> f64 {
    let center = match site.piece {
        Kp20PieceId::NearMidpointLeft | Kp20PieceId::NearMidpointRight => 6.0 * site.u,
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => 6.0 * site.v,
    };
    match (row as f64 - center).abs() {
        distance if distance <= 0.75 => 1.0,
        distance if distance <= 1.75 => 0.45,
        _ => 0.0,
    }
}
