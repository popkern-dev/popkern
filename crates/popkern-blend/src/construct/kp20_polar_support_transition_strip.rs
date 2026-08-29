//! Support-owned transition-strip candidates for KP20 limited-lift patches.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::kp20_polar_support_surface_split::{
    evaluate_bi6_normal, propagate_support_pole_interfaces, support_owner_piece, support_side_index,
};

type FoldSummary = (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport);
type TransitionResult = Result<
    (
        Kp20SupportNormalLiftReport,
        Kp20SupportTransitionStripReport,
    ),
    String,
>;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportTransitionStripReport {
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
    pub(super) band_width: usize,
    pub(super) mode: &'static str,
    pub(super) scale: f64,
    pub(super) fold_u: f64,
    pub(super) fold_v: f64,
    pub(super) decline_reason: &'static str,
}

pub(super) fn refine_support_transition_strip(
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> TransitionResult {
    let target_piece = Kp20PieceId::NearMidpointLeft;
    let target = worst_support_site(&baseline.summary, target_piece);
    let mut report = Kp20SupportTransitionStripReport {
        evaluated: 0,
        fast_accepted: 0,
        fast_rejected: 0,
        exact_evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.summary.fold_count,
        end_folds: baseline.summary.fold_count,
        target_piece,
        target_sector: target.map(|site| site.sector),
        band_width: 0,
        mode: "none",
        scale: 0.0,
        fold_u: target.map_or(0.0, |site| site.u),
        fold_v: target.map_or(0.0, |site| site.v),
        decline_reason: "support_transition_strip_no_candidate",
    };
    let Some(target) = target else {
        report.decline_reason = "support_transition_strip_no_target_site";
        return Ok((baseline, report));
    };
    let mut best_exact: Option<Kp20TransitionTrial> = None;
    for band_width in [2_usize, 3] {
        for mode in Kp20TransitionStripMode::all() {
            for scale in [0.45, 0.70, 1.0] {
                report.evaluated += 1;
                let applied = transition_candidate(
                    &baseline.macro_patches,
                    &baseline.summary,
                    target,
                    band_width,
                    mode,
                    scale,
                );
                let fast = fast_candidate_summary_from_changed_patches(
                    tensor_border,
                    applied.macro_patches.clone(),
                    &baseline.summary,
                    &applied.changed_patches,
                    baseline.summary.fold_count,
                )?;
                if !validation_is_clean(fast.1) || fast.0.fold_count >= baseline.summary.fold_count
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
                let trial = Kp20TransitionTrial {
                    macro_patches: applied.macro_patches,
                    exact,
                    band_width,
                    mode,
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
    let Some(best) = best_exact else {
        return Ok((baseline, report));
    };
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = best.exact.0.fold_count;
    report.band_width = best.band_width;
    report.mode = best.mode.label();
    report.scale = best.scale;
    report.decline_reason = "none";
    let mut refined = baseline;
    refined.summary = best.exact.0;
    refined.validation = best.exact.1;
    refined.macro_patches = best.macro_patches;
    Ok((refined, report))
}

#[derive(Debug, Clone)]
struct Kp20TransitionApplied {
    macro_patches: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
}

#[derive(Debug, Clone)]
struct Kp20TransitionTrial {
    macro_patches: Vec<Kp20MacroPatch>,
    exact: FoldSummary,
    band_width: usize,
    mode: Kp20TransitionStripMode,
    scale: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kp20TransitionStripMode {
    FoldNormalRelief,
    FoldNormalOppose,
    ParameterExpand,
    ParameterContract,
    SupportImageBoundaryExpand,
    SupportImageBoundaryContract,
}

impl Kp20TransitionStripMode {
    fn all() -> [Self; 6] {
        [
            Self::FoldNormalRelief,
            Self::FoldNormalOppose,
            Self::ParameterExpand,
            Self::ParameterContract,
            Self::SupportImageBoundaryExpand,
            Self::SupportImageBoundaryContract,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::FoldNormalRelief => "fold_normal_relief",
            Self::FoldNormalOppose => "fold_normal_oppose",
            Self::ParameterExpand => "parameter_expand",
            Self::ParameterContract => "parameter_contract",
            Self::SupportImageBoundaryExpand => "support_image_boundary_expand",
            Self::SupportImageBoundaryContract => "support_image_boundary_contract",
        }
    }
}

fn transition_candidate(
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    target: Kp20FoldSite,
    band_width: usize,
    mode: Kp20TransitionStripMode,
    scale: f64,
) -> Kp20TransitionApplied {
    let mut macro_patches = accepted_patches.to_vec();
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    if target.sector >= macro_patches.len() {
        return Kp20TransitionApplied {
            macro_patches,
            changed_patches,
        };
    }
    let support_piece = support_owner_piece(target.piece);
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    let base = *macro_patches[target.sector].piece(support_piece);
    *macro_patches[target.sector].piece_mut(support_piece) = transition_surface(
        &base,
        accepted_summary,
        target,
        support_piece,
        band_width,
        mode,
        scale,
    );
    changed_patches[target.sector][support_piece.index()] = true;
    support_changed[target.sector][support_side_index(support_piece)] = true;
    propagate_support_pole_interfaces(&mut macro_patches, &mut changed_patches, support_changed);
    Kp20TransitionApplied {
        macro_patches,
        changed_patches,
    }
}

fn transition_surface(
    base: &popkern_core::geometry::gregory::Kp20Bi6Patch,
    summary: &Kp20SamplerFoldSummary,
    target: Kp20FoldSite,
    support_piece: Kp20PieceId,
    band_width: usize,
    mode: Kp20TransitionStripMode,
    scale: f64,
) -> popkern_core::geometry::gregory::Kp20Bi6Patch {
    if matches!(
        mode,
        Kp20TransitionStripMode::SupportImageBoundaryExpand
            | Kp20TransitionStripMode::SupportImageBoundaryContract
    ) {
        return support_image_boundary_strip_surface(
            base,
            summary,
            target,
            support_piece,
            mode,
            scale,
        );
    }
    let mut surface = *base;
    let span = transition_span(base);
    let normal = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let row_range = match band_width {
        2 => 2..=3,
        _ => 2..=4,
    };
    for column in 0..7 {
        let support_anchor = base.cps[1][column];
        let interior_anchor = base.cps[5][column];
        let column_v = column as f64 / 6.0;
        for row in row_range.clone() {
            let row_u = row as f64 / 6.0;
            let t = transition_parameter(row_u, target.u, mode, scale);
            let fair = support_anchor + (interior_anchor - support_anchor) * t;
            let relief = transition_relief(
                summary,
                target.sector,
                support_piece,
                row_u,
                column_v,
                normal,
                mode,
            ) * (span * scale);
            surface.cps[row][column] = fair + relief;
        }
        if matches!(
            mode,
            Kp20TransitionStripMode::FoldNormalRelief | Kp20TransitionStripMode::FoldNormalOppose
        ) {
            let row_u = 5.0 / 6.0;
            let relief = transition_relief(
                summary,
                target.sector,
                support_piece,
                row_u,
                column_v,
                normal,
                mode,
            ) * (0.35 * span * scale);
            surface.cps[5][column] = base.cps[5][column] + relief;
        }
    }
    surface
}

fn support_image_boundary_strip_surface(
    base: &popkern_core::geometry::gregory::Kp20Bi6Patch,
    summary: &Kp20SamplerFoldSummary,
    target: Kp20FoldSite,
    support_piece: Kp20PieceId,
    mode: Kp20TransitionStripMode,
    scale: f64,
) -> popkern_core::geometry::gregory::Kp20Bi6Patch {
    let mut surface = *base;
    let span = transition_span(base);
    let normal = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let direction = match mode {
        Kp20TransitionStripMode::SupportImageBoundaryContract => -normal,
        _ => normal,
    };
    for row in 1..=5 {
        let boundary_anchor = base.cps[row][0];
        let interior_anchor = base.cps[row][5];
        let row_u = row as f64 / 6.0;
        for column in 1..=2 {
            let column_v = column as f64 / 6.0;
            let t = support_image_boundary_parameter(column_v, target.v, mode, scale);
            let fair = boundary_anchor + (interior_anchor - boundary_anchor) * t;
            let relief = transition_relief(
                summary,
                target.sector,
                support_piece,
                row_u,
                column_v,
                direction,
                mode,
            ) * (0.55 * span * scale);
            surface.cps[row][column] = fair + relief;
        }
    }
    surface
}

fn support_image_boundary_parameter(
    column_v: f64,
    site_v: f64,
    mode: Kp20TransitionStripMode,
    scale: f64,
) -> f64 {
    let base = (column_v / (5.0 / 6.0)).clamp(0.0, 1.0);
    let locality = (1.0 - ((column_v - site_v).abs() / 0.24)).max(0.0);
    let sign = match mode {
        Kp20TransitionStripMode::SupportImageBoundaryExpand => 1.0,
        Kp20TransitionStripMode::SupportImageBoundaryContract => -1.0,
        _ => 0.0,
    };
    (base + sign * 0.12 * locality * scale).clamp(0.04, 0.92)
}

fn transition_parameter(row_u: f64, site_u: f64, mode: Kp20TransitionStripMode, scale: f64) -> f64 {
    let base = ((row_u - (1.0 / 6.0)) / (4.0 / 6.0)).clamp(0.0, 1.0);
    let local = (1.0 - ((row_u - site_u).abs() / 0.38)).max(0.0);
    let sign = match mode {
        Kp20TransitionStripMode::ParameterExpand => {
            if base <= 0.5 {
                -1.0
            } else {
                1.0
            }
        }
        Kp20TransitionStripMode::ParameterContract => {
            if base <= 0.5 {
                1.0
            } else {
                -1.0
            }
        }
        _ => 0.0,
    };
    (base + sign * 0.14 * local * scale).clamp(0.06, 0.94)
}

fn transition_relief(
    summary: &Kp20SamplerFoldSummary,
    sector: usize,
    support_piece: Kp20PieceId,
    row_u: f64,
    column_v: f64,
    normal: Vec3,
    mode: Kp20TransitionStripMode,
) -> Vec3 {
    let mut displacement = Vec3::zero();
    for site in summary
        .fold_sites
        .iter()
        .filter(|site| site.sector == sector && support_owner_piece(site.piece) == support_piece)
    {
        let du = row_u - site.u;
        let dv = column_v - site.v;
        let locality = (1.0 - (du * du + dv * dv) / (0.34_f64 * 0.34_f64)).max(0.0);
        let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
        displacement = displacement + normal * (locality * locality * severity);
    }
    let displacement = capped_displacement(displacement, 0.10);
    match mode {
        Kp20TransitionStripMode::FoldNormalOppose => -displacement,
        _ => displacement,
    }
}

fn capped_displacement(displacement: Vec3, max_length: f64) -> Vec3 {
    let length = displacement.length();
    if length > max_length
        && let Some(direction) = displacement.normalized()
    {
        return direction * max_length;
    }
    displacement
}

fn transition_span(surface: &popkern_core::geometry::gregory::Kp20Bi6Patch) -> f64 {
    let support_width = surface.cps[1][0].distance_to(&surface.cps[1][6]);
    let interface_width = surface.cps[5][0].distance_to(&surface.cps[5][6]);
    let cross_depth = surface.cps[1][3].distance_to(&surface.cps[5][3]);
    ((support_width + interface_width + cross_depth) / 3.0).max(1.0e-9)
}

fn worst_support_site(
    summary: &Kp20SamplerFoldSummary,
    target_piece: Kp20PieceId,
) -> Option<Kp20FoldSite> {
    summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

fn trial_is_better(candidate: &FoldSummary, best: &FoldSummary) -> bool {
    candidate.0.fold_count < best.0.fold_count
        || (candidate.0.fold_count == best.0.fold_count
            && candidate.0.min_jacobian_normal_cos > best.0.min_jacobian_normal_cos)
}
