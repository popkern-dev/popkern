//! Bounded support-strip reauthoring for KP20 limited-lift candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};

use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
use super::kp20_polar_multipatch_strip::{
    Kp20MultiPatchStripBuildReport, Kp20MultiPatchStripSurface, build_multipatch_strip_surfaces,
};
use super::kp20_polar_multipatch_strip_optimize::{
    Kp20StripInteriorOptimizerMode, optimized_strip_surface,
};
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::kp20_polar_support_surface_split::{
    propagate_support_pole_interfaces, support_owner_piece, support_side_index,
};

type FoldSummary = (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport);

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportStripReauthorReport {
    pub(super) evaluated: usize,
    pub(super) fast_accepted: usize,
    pub(super) fast_rejected: usize,
    pub(super) exact_evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) target_sector: Option<usize>,
    pub(super) mode: &'static str,
    pub(super) coupling: &'static str,
    pub(super) coupling_family: &'static str,
    pub(super) scale: f64,
    pub(super) worst_before: String,
    pub(super) worst_after: String,
    pub(super) decline_reason: &'static str,
}

pub(super) fn refine_support_strip_reauthoring(
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> Result<(Kp20SupportNormalLiftReport, Kp20SupportStripReauthorReport), String> {
    let target_piece = Kp20PieceId::NearMidpointLeft;
    let mut report = Kp20SupportStripReauthorReport {
        evaluated: 0,
        fast_accepted: 0,
        fast_rejected: 0,
        exact_evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.summary.fold_count,
        end_folds: baseline.summary.fold_count,
        target_piece,
        target_sector: worst_target_sector(&baseline.summary, target_piece),
        mode: "none",
        coupling: "none",
        coupling_family: Kp20SupportStripCouplingMode::family_label(),
        scale: 0.0,
        worst_before: worst_site_label(&baseline.summary),
        worst_after: worst_site_label(&baseline.summary),
        decline_reason: "support_strip_reauthor_no_candidate",
    };
    if baseline.summary.fold_count == 0 {
        report.decline_reason = "support_strip_reauthor_no_folds";
        return Ok((baseline, report));
    }
    let inventory = multipatch_inventory_from_fold_summary(&baseline.summary);
    let strips =
        build_multipatch_strip_surfaces(&inventory, &baseline.macro_patches, &baseline.summary);
    if strips.decline_reason.is_some() {
        report.decline_reason = "support_strip_reauthor_strip_build_declined";
        return Ok((baseline, report));
    }
    let mut best_exact: Option<Kp20SelectedStripTrial> = None;
    for strip in selected_target_strips(&strips, &baseline.summary, target_piece) {
        for mode in Kp20StripInteriorOptimizerMode::exact_guided_all() {
            for coupling in Kp20SupportStripCouplingMode::all() {
                for scale in [1.0, 0.75, 0.5] {
                    report.evaluated += 1;
                    let applied = selected_strip_candidate(
                        &baseline.macro_patches,
                        &baseline.summary,
                        strip,
                        mode,
                        coupling,
                        scale,
                    );
                    let fast = fast_candidate_summary_from_changed_patches(
                        tensor_border,
                        applied.macro_patches.clone(),
                        &baseline.summary,
                        &applied.changed_patches,
                        baseline.summary.fold_count,
                    )?;
                    if !validation_is_clean(fast.1)
                        || fast.0.fold_count >= baseline.summary.fold_count
                    {
                        report.fast_rejected += 1;
                        report.rejected += 1;
                        continue;
                    }
                    report.fast_accepted += 1;
                    report.exact_evaluated += 1;
                    let exact = candidate_summary_from_changed_patches(
                        tensor_border,
                        applied.macro_patches.clone(),
                        &baseline.summary,
                        &applied.changed_patches,
                        baseline.summary.fold_count,
                    )?;
                    if !validation_is_clean(exact.1)
                        || exact.0.fold_count >= baseline.summary.fold_count
                    {
                        report.rejected += 1;
                        continue;
                    }
                    let trial = Kp20SelectedStripTrial {
                        macro_patches: applied.macro_patches,
                        exact,
                        mode,
                        coupling,
                        scale,
                    };
                    if best_exact
                        .as_ref()
                        .is_none_or(|best| trial_is_better(&trial.exact, &best.exact))
                    {
                        best_exact = Some(trial);
                    }
                }
            }
        }
    }
    let Some(best) = best_exact else {
        if report.evaluated == 0 {
            report.decline_reason = "support_strip_reauthor_no_target_strip";
        }
        return Ok((baseline, report));
    };
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = best.exact.0.fold_count;
    report.mode = best.mode.label();
    report.coupling = best.coupling.label();
    report.scale = best.scale;
    report.worst_after = worst_site_label(&best.exact.0);
    report.decline_reason = "none";
    let mut refined = baseline;
    refined.summary = best.exact.0;
    refined.validation = best.exact.1;
    refined.macro_patches = best.macro_patches;
    Ok((refined, report))
}

#[derive(Debug, Clone)]
struct Kp20SelectedStripApplied {
    macro_patches: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
}

#[derive(Debug, Clone)]
struct Kp20SelectedStripTrial {
    macro_patches: Vec<Kp20MacroPatch>,
    exact: FoldSummary,
    mode: Kp20StripInteriorOptimizerMode,
    coupling: Kp20SupportStripCouplingMode,
    scale: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kp20SupportStripCouplingMode {
    PropagateSupportPole,
    SharedFrameExpand,
    SharedFrameContract,
}

impl Kp20SupportStripCouplingMode {
    fn all() -> [Self; 3] {
        [
            Self::PropagateSupportPole,
            Self::SharedFrameExpand,
            Self::SharedFrameContract,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::PropagateSupportPole => "propagate_support_pole",
            Self::SharedFrameExpand => "shared_frame_expand",
            Self::SharedFrameContract => "shared_frame_contract",
        }
    }

    fn family_label() -> &'static str {
        "propagate_support_pole|shared_frame_expand|shared_frame_contract"
    }
}

fn selected_target_strips<'a>(
    strips: &'a Kp20MultiPatchStripBuildReport,
    summary: &Kp20SamplerFoldSummary,
    target_piece: Kp20PieceId,
) -> Vec<&'a Kp20MultiPatchStripSurface> {
    let Some(worst) = summary
        .fold_sites
        .iter()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    else {
        return Vec::new();
    };
    strips
        .strip_surfaces
        .iter()
        .filter(|strip| strip.strip_piece == target_piece && strip.sector_index == worst.sector)
        .collect()
}

fn selected_strip_candidate(
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    coupling: Kp20SupportStripCouplingMode,
    scale: f64,
) -> Kp20SelectedStripApplied {
    let mut macro_patches = accepted_patches.to_vec();
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    if strip.sector_index < macro_patches.len() {
        *macro_patches[strip.sector_index].piece_mut(strip.strip_piece) =
            optimized_strip_surface(accepted_patches, accepted_summary, strip, mode, scale);
        changed_patches[strip.sector_index][strip.strip_piece.index()] = true;
        apply_support_strip_coupling(
            &mut macro_patches,
            &mut changed_patches,
            accepted_summary,
            strip,
            coupling,
            scale,
        );
    }
    Kp20SelectedStripApplied {
        macro_patches,
        changed_patches,
    }
}

fn apply_support_strip_coupling(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    coupling: Kp20SupportStripCouplingMode,
    scale: f64,
) {
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    support_changed[strip.sector_index][support_side_index(strip.strip_piece)] = true;
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    if matches!(coupling, Kp20SupportStripCouplingMode::PropagateSupportPole) {
        return;
    }
    let Some(site) = worst_strip_site(accepted_summary, strip) else {
        return;
    };
    reauthor_shared_support_inner_frame(
        macro_patches,
        changed_patches,
        strip.sector_index,
        strip.strip_piece,
        site.v,
        matches!(coupling, Kp20SupportStripCouplingMode::SharedFrameExpand),
        scale,
    );
}

fn worst_strip_site(
    summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
) -> Option<super::kp20_polar_fold_summary::Kp20FoldSite> {
    summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| {
            site.sector == strip.sector_index
                && support_owner_piece(site.piece) == strip.strip_piece
        })
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

fn reauthor_shared_support_inner_frame(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    sector: usize,
    support_piece: Kp20PieceId,
    site_v: f64,
    expand: bool,
    scale: f64,
) {
    let origin_piece = origin_piece_for_strip(support_piece);
    let target_column = if site_v <= 0.5 { 1_usize } else { 5_usize };
    for column in 0..7 {
        let weight = coupled_shoulder_column_weight(column, target_column);
        if weight <= 0.0 {
            continue;
        }
        let support_current = macro_patches[sector].piece(support_piece).cps[5][column];
        let origin_current = macro_patches[sector].piece(origin_piece).cps[1][column];
        let interface = macro_patches[sector].piece(support_piece).cps[6][column];
        let half_step = (origin_current - support_current) * 0.5;
        let frame_scale = if expand { 1.15 } else { 0.75 };
        let support_target = interface - half_step * frame_scale;
        let origin_target = interface + half_step * frame_scale;
        let factor = (weight * scale).clamp(0.0, 1.0);
        macro_patches[sector].piece_mut(support_piece).cps[5][column] =
            support_current + (support_target - support_current) * factor;
        macro_patches[sector].piece_mut(origin_piece).cps[1][column] =
            origin_current + (origin_target - origin_current) * factor;
    }
    changed_patches[sector][support_piece.index()] = true;
    changed_patches[sector][origin_piece.index()] = true;
}

fn origin_piece_for_strip(strip_piece: Kp20PieceId) -> Kp20PieceId {
    match strip_piece {
        Kp20PieceId::NearMidpointLeft => Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearMidpointRight => Kp20PieceId::NearOriginRight,
        _ => unreachable!("support-strip reauthoring targets support-owned pieces"),
    }
}

fn coupled_shoulder_column_weight(column: usize, target_column: usize) -> f64 {
    match column.abs_diff(target_column) {
        0 => 1.0,
        1 => 0.45,
        _ => 0.0,
    }
}

fn trial_is_better(candidate: &FoldSummary, best: &FoldSummary) -> bool {
    candidate.0.fold_count < best.0.fold_count
        || (candidate.0.fold_count == best.0.fold_count
            && candidate.0.min_jacobian_normal_cos > best.0.min_jacobian_normal_cos)
}

fn worst_target_sector(
    summary: &Kp20SamplerFoldSummary,
    target_piece: Kp20PieceId,
) -> Option<usize> {
    summary
        .fold_sites
        .iter()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|site| site.sector)
}

fn worst_site_label(summary: &Kp20SamplerFoldSummary) -> String {
    summary.worst_fold_site().map_or_else(
        || "none".to_string(),
        |site| {
            format!(
                "sector={},piece={:?},u={:.3},v={:.3},cos={:.3e}",
                site.sector, site.piece, site.u, site.v, site.min_jacobian_normal_cos
            )
        },
    )
}
