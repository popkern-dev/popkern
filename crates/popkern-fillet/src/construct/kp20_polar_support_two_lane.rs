//! Candidate-only two-lane support patch construction for KP20 limited lift.

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
type TwoLaneResult = Result<(Kp20SupportNormalLiftReport, Kp20SupportTwoLaneReport), String>;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportTwoLaneReport {
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
    pub(super) lane_count: usize,
    pub(super) split_parameter: f64,
    pub(super) max_interface_residual: f64,
    pub(super) mode: &'static str,
    pub(super) scale: f64,
    pub(super) decline_reason: &'static str,
}

pub(super) fn refine_support_two_lane(
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> TwoLaneResult {
    let target_piece = Kp20PieceId::NearMidpointLeft;
    let target = worst_support_site(&baseline.summary, target_piece);
    let mut report = Kp20SupportTwoLaneReport {
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
        lane_count: 2,
        split_parameter: 0.0,
        max_interface_residual: 0.0,
        mode: "none",
        scale: 0.0,
        decline_reason: "support_two_lane_no_candidate",
    };
    let Some(target) = target else {
        report.decline_reason = "support_two_lane_no_target_site";
        return Ok((baseline, report));
    };
    let mut best_exact: Option<Kp20TwoLaneTrial> = None;
    for split_parameter in [0.42, 0.50, 0.58] {
        for mode in Kp20TwoLaneMode::all() {
            for scale in [0.45, 0.75] {
                report.evaluated += 1;
                let applied = two_lane_candidate(
                    &baseline.macro_patches,
                    target,
                    split_parameter,
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
                let trial = Kp20TwoLaneTrial {
                    macro_patches: applied.macro_patches,
                    exact,
                    split_parameter,
                    max_interface_residual: applied.max_interface_residual,
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
    report.split_parameter = best.split_parameter;
    report.max_interface_residual = best.max_interface_residual;
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
struct Kp20TwoLaneApplied {
    macro_patches: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
    max_interface_residual: f64,
}

#[derive(Debug, Clone)]
struct Kp20TwoLaneTrial {
    macro_patches: Vec<Kp20MacroPatch>,
    exact: FoldSummary,
    split_parameter: f64,
    max_interface_residual: f64,
    mode: Kp20TwoLaneMode,
    scale: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kp20TwoLaneMode {
    SupportNormalRelief,
    ChordRelief,
    InterfaceContract,
}

impl Kp20TwoLaneMode {
    fn all() -> [Self; 3] {
        [
            Self::SupportNormalRelief,
            Self::ChordRelief,
            Self::InterfaceContract,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::SupportNormalRelief => "support_normal_relief",
            Self::ChordRelief => "chord_relief",
            Self::InterfaceContract => "interface_contract",
        }
    }
}

fn two_lane_candidate(
    accepted_patches: &[Kp20MacroPatch],
    target: Kp20FoldSite,
    split_parameter: f64,
    mode: Kp20TwoLaneMode,
    scale: f64,
) -> Kp20TwoLaneApplied {
    let mut macro_patches = accepted_patches.to_vec();
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    if target.sector >= macro_patches.len() {
        return Kp20TwoLaneApplied {
            macro_patches,
            changed_patches,
            max_interface_residual: 0.0,
        };
    }
    let support_piece = support_owner_piece(target.piece);
    let base = *macro_patches[target.sector].piece(support_piece);
    let (surface, max_interface_residual) =
        two_lane_surface(&base, target, split_parameter, mode, scale);
    *macro_patches[target.sector].piece_mut(support_piece) = surface;
    changed_patches[target.sector][support_piece.index()] = true;
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    support_changed[target.sector][support_side_index(support_piece)] = true;
    propagate_support_pole_interfaces(&mut macro_patches, &mut changed_patches, support_changed);
    Kp20TwoLaneApplied {
        macro_patches,
        changed_patches,
        max_interface_residual,
    }
}

fn two_lane_surface(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    split_parameter: f64,
    mode: Kp20TwoLaneMode,
    scale: f64,
) -> (Kp20Bi6Patch, f64) {
    let mut surface = *base;
    let span = support_span(base);
    let normal = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let interface_row = support_lane_interface(base, target, split_parameter, normal, mode, scale);
    surface.cps[3] = interface_row;
    for column in 0..7 {
        surface.cps[2][column] = midpoint(base.cps[1][column], interface_row[column], 0.5);
        surface.cps[4][column] = midpoint(interface_row[column], base.cps[6][column], 0.34);
        surface.cps[5][column] = midpoint(interface_row[column], base.cps[6][column], 0.67);
        if matches!(
            mode,
            Kp20TwoLaneMode::SupportNormalRelief | Kp20TwoLaneMode::ChordRelief
        ) {
            let v = column as f64 / 6.0;
            let local = (1.0 - ((v - target.v).abs() / 0.34)).max(0.0);
            let relief = two_lane_direction(base, normal, mode, v) * (span * scale * 0.08 * local);
            surface.cps[4][column] = surface.cps[4][column] + relief;
            surface.cps[5][column] = surface.cps[5][column] + relief * 0.5;
        }
    }
    (
        surface,
        max_interface_residual(&surface.cps[3], &interface_row),
    )
}

fn support_lane_interface(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    split_parameter: f64,
    normal: Vec3,
    mode: Kp20TwoLaneMode,
    scale: f64,
) -> [popkern_core::geometry::Point3; 7] {
    let mut interface = base.cps[3];
    let span = support_span(base);
    for (column, out) in interface.iter_mut().enumerate() {
        let v = column as f64 / 6.0;
        let support = base.cps[1][column];
        let interior = base.cps[6][column];
        let mut point = support + (interior - support) * split_parameter;
        let local = (1.0 - ((v - target.v).abs() / 0.38)).max(0.0);
        let direction = two_lane_direction(base, normal, mode, v);
        let sign = if matches!(mode, Kp20TwoLaneMode::InterfaceContract) {
            -1.0
        } else {
            1.0
        };
        point = point + direction * (sign * span * scale * 0.10 * local * local);
        *out = point;
    }
    interface
}

fn two_lane_direction(base: &Kp20Bi6Patch, normal: Vec3, mode: Kp20TwoLaneMode, v: f64) -> Vec3 {
    match mode {
        Kp20TwoLaneMode::SupportNormalRelief => normal,
        Kp20TwoLaneMode::ChordRelief | Kp20TwoLaneMode::InterfaceContract => {
            let column = ((v * 6.0).round() as usize).clamp(0, 6);
            let chord = base.cps[1][6] - base.cps[1][0];
            let depth = base.cps[6][column] - base.cps[1][column];
            (chord + depth).normalized().unwrap_or(normal)
        }
    }
}

fn midpoint(
    a: popkern_core::geometry::Point3,
    b: popkern_core::geometry::Point3,
    t: f64,
) -> popkern_core::geometry::Point3 {
    a + (b - a) * t
}

fn max_interface_residual(
    row: &[popkern_core::geometry::Point3; 7],
    interface: &[popkern_core::geometry::Point3; 7],
) -> f64 {
    row.iter()
        .zip(interface)
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0, f64::max)
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
