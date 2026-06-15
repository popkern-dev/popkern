//! Alternate support-owned macro-patch construction for KP20 limited lift.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    Kp20Bi6Patch, Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
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
type MacroConstructionResult = Result<
    (
        Kp20SupportNormalLiftReport,
        Kp20SupportMacroConstructionReport,
    ),
    String,
>;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportMacroConstructionReport {
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
    pub(super) scale: f64,
    pub(super) fold_u: f64,
    pub(super) fold_v: f64,
    pub(super) decline_reason: &'static str,
}

pub(super) fn refine_support_macro_construction(
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> MacroConstructionResult {
    let target_piece = Kp20PieceId::NearMidpointLeft;
    let target = worst_support_site(&baseline.summary, target_piece);
    let mut report = Kp20SupportMacroConstructionReport {
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
        mode: "none",
        scale: 0.0,
        fold_u: target.map_or(0.0, |site| site.u),
        fold_v: target.map_or(0.0, |site| site.v),
        decline_reason: "support_macro_construction_no_candidate",
    };
    let Some(target) = target else {
        report.decline_reason = "support_macro_construction_no_target_site";
        return Ok((baseline, report));
    };
    let mut best_exact: Option<Kp20MacroConstructionTrial> = None;
    for mode in Kp20SupportMacroConstructionMode::all() {
        for scale in [0.35, 0.60, 0.85] {
            report.evaluated += 1;
            let applied =
                macro_construction_candidate(&baseline.macro_patches, target, mode, scale);
            let fast = fast_candidate_summary_from_changed_patches(
                tensor_border,
                applied.macro_patches.clone(),
                &baseline.summary,
                &applied.changed_patches,
                baseline.summary.fold_count,
            )?;
            if !validation_is_clean(fast.1) || fast.0.fold_count >= baseline.summary.fold_count {
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
            if !validation_is_clean(exact.1) || exact.0.fold_count >= baseline.summary.fold_count {
                report.rejected += 1;
                continue;
            }
            let trial = Kp20MacroConstructionTrial {
                macro_patches: applied.macro_patches,
                exact,
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
    let Some(best) = best_exact else {
        return Ok((baseline, report));
    };
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = best.exact.0.fold_count;
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
struct Kp20MacroConstructionApplied {
    macro_patches: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
}

#[derive(Debug, Clone)]
struct Kp20MacroConstructionTrial {
    macro_patches: Vec<Kp20MacroPatch>,
    exact: FoldSummary,
    mode: Kp20SupportMacroConstructionMode,
    scale: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
enum Kp20SupportMacroConstructionMode {
    SupportNormalFrame,
    RadialChordFrame,
    FoldReliefFrame,
    PairedSupportFrame,
}

impl Kp20SupportMacroConstructionMode {
    fn all() -> [Self; 4] {
        [
            Self::SupportNormalFrame,
            Self::RadialChordFrame,
            Self::FoldReliefFrame,
            Self::PairedSupportFrame,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::SupportNormalFrame => "support_normal_frame",
            Self::RadialChordFrame => "radial_chord_frame",
            Self::FoldReliefFrame => "fold_relief_frame",
            Self::PairedSupportFrame => "paired_support_frame",
        }
    }
}

fn macro_construction_candidate(
    accepted_patches: &[Kp20MacroPatch],
    target: Kp20FoldSite,
    mode: Kp20SupportMacroConstructionMode,
    scale: f64,
) -> Kp20MacroConstructionApplied {
    let mut macro_patches = accepted_patches.to_vec();
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    if target.sector >= macro_patches.len() {
        return Kp20MacroConstructionApplied {
            macro_patches,
            changed_patches,
        };
    }
    let support_piece = support_owner_piece(target.piece);
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    if mode == Kp20SupportMacroConstructionMode::PairedSupportFrame {
        for piece in [
            Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointRight,
        ] {
            let base = *macro_patches[target.sector].piece(piece);
            let local_target = support_piece_target(target, piece);
            *macro_patches[target.sector].piece_mut(piece) =
                constructed_support_macro_patch(&base, local_target, mode, scale);
            changed_patches[target.sector][piece.index()] = true;
            support_changed[target.sector][support_side_index(piece)] = true;
        }
    } else {
        let base = *macro_patches[target.sector].piece(support_piece);
        *macro_patches[target.sector].piece_mut(support_piece) =
            constructed_support_macro_patch(&base, target, mode, scale);
        changed_patches[target.sector][support_piece.index()] = true;
        support_changed[target.sector][support_side_index(support_piece)] = true;
    }
    propagate_support_pole_interfaces(&mut macro_patches, &mut changed_patches, support_changed);
    Kp20MacroConstructionApplied {
        macro_patches,
        changed_patches,
    }
}

fn support_piece_target(target: Kp20FoldSite, piece: Kp20PieceId) -> Kp20FoldSite {
    Kp20FoldSite { piece, ..target }
}

fn constructed_support_macro_patch(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    mode: Kp20SupportMacroConstructionMode,
    scale: f64,
) -> Kp20Bi6Patch {
    let mut surface = *base;
    let span = support_span(base);
    let normal = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    for row in 2..=5 {
        let row_u = row as f64 / 6.0;
        for column in 0..7 {
            let column_v = column as f64 / 6.0;
            let support = base.cps[1][column];
            let interface = base.cps[6][column];
            let base_t = ((row_u - (1.0 / 6.0)) / (5.0 / 6.0)).clamp(0.0, 1.0);
            let fair = support + (interface - support) * base_t;
            let frame = macro_frame_displacement(base, target, row_u, column_v, normal, mode);
            surface.cps[row][column] =
                fair + frame * (span * scale * macro_frame_weight(target, row_u, column_v));
        }
    }
    surface
}

fn macro_frame_displacement(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    row_u: f64,
    column_v: f64,
    normal: Vec3,
    mode: Kp20SupportMacroConstructionMode,
) -> Vec3 {
    match mode {
        Kp20SupportMacroConstructionMode::SupportNormalFrame
        | Kp20SupportMacroConstructionMode::PairedSupportFrame => normal,
        Kp20SupportMacroConstructionMode::RadialChordFrame => {
            let chord = base.cps[1][6] - base.cps[1][0];
            let depth = base.cps[6][column_index(column_v)] - base.cps[1][column_index(column_v)];
            (chord + depth).normalized().unwrap_or(normal)
        }
        Kp20SupportMacroConstructionMode::FoldReliefFrame => {
            let local = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(normal);
            if row_u < target.u { local } else { -local }
        }
    }
}

fn macro_frame_weight(target: Kp20FoldSite, row_u: f64, column_v: f64) -> f64 {
    let du = row_u - target.u;
    let dv = column_v - target.v;
    let locality = (1.0 - (du * du + dv * dv) / (0.42_f64 * 0.42_f64)).max(0.0);
    let severity = (1.0 - target.min_jacobian_normal_cos).clamp(0.0, 1.0);
    locality * locality * severity.max(0.25)
}

fn column_index(column_v: f64) -> usize {
    ((column_v * 6.0).round() as usize).clamp(0, 6)
}

fn support_span(surface: &Kp20Bi6Patch) -> f64 {
    let support_width = surface.cps[1][0].distance_to(&surface.cps[1][6]);
    let interface_width = surface.cps[6][0].distance_to(&surface.cps[6][6]);
    let cross_depth = surface.cps[1][3].distance_to(&surface.cps[6][3]);
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
