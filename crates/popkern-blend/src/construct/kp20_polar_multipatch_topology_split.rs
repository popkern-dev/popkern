//! Bounded multipatch topology-split candidates for KP20 limited-lift patches.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6Patch, Kp20Bi6ValidationReport, Kp20MacroPatch, RectangularGregoryBiNonic,
    TensorBorderBi3, elevate_bi6_grid_to_bi9,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_apply::{Kp20InnerConstructionMode, apply_multipatch_candidate};
use super::kp20_polar_multipatch_inner::build_multipatch_inner_patches;
use super::kp20_polar_multipatch_interface::build_multipatch_interfaces;
use super::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
use super::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;
use super::kp20_polar_multipatch_strip_optimize::Kp20StripInteriorOptimizerMode;
use super::kp20_polar_multipatch_validate::validate_multipatch_candidate;
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;

type FoldSummary = (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport);
type TopologySplitResult = Result<(Kp20SupportNormalLiftReport, Kp20TopologySplitReport), String>;

#[derive(Debug, Clone)]
pub(super) struct Kp20TopologySplitReport {
    pub(super) evaluated: usize,
    pub(super) fast_accepted: usize,
    pub(super) fast_rejected: usize,
    pub(super) exact_evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) strip_count: usize,
    pub(super) inner_count: usize,
    pub(super) direct_patch_count: usize,
    pub(super) direct_patch_folds: usize,
    pub(super) changed_piece_count: usize,
    pub(super) strip_mode: &'static str,
    pub(super) inner_mode: &'static str,
    pub(super) scale: f64,
    pub(super) validation_decline_reason: &'static str,
    pub(super) decline_reason: &'static str,
}

pub(super) fn refine_multipatch_topology_split(
    tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> TopologySplitResult {
    let mut report = Kp20TopologySplitReport {
        evaluated: 0,
        fast_accepted: 0,
        fast_rejected: 0,
        exact_evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.summary.fold_count,
        end_folds: baseline.summary.fold_count,
        strip_count: 0,
        inner_count: 0,
        direct_patch_count: 0,
        direct_patch_folds: 0,
        changed_piece_count: 0,
        strip_mode: "none",
        inner_mode: "none",
        scale: 0.0,
        validation_decline_reason: "none",
        decline_reason: "multipatch_topology_split_no_candidate",
    };
    if baseline.summary.fold_count == 0 {
        report.decline_reason = "multipatch_topology_split_no_folds";
        return Ok((baseline, report));
    }
    let inventory = multipatch_inventory_from_fold_summary(&baseline.summary);
    let strips =
        build_multipatch_strip_surfaces(&inventory, &baseline.macro_patches, &baseline.summary);
    let interfaces = build_multipatch_interfaces(&strips);
    let inner = build_multipatch_inner_patches(&inventory, &interfaces);
    let validation = validate_multipatch_candidate(&inventory, &strips, &interfaces, &inner);
    report.strip_count = strips.built_strip_patch_count;
    report.inner_count = inner.built_inner_patch_count;
    let direct = measure_direct_multipatch_surfaces(&strips, &inner);
    report.direct_patch_count = direct.patch_count;
    report.direct_patch_folds = direct.fold_count;
    if let Some(reason) = validation.decline_reason {
        report.validation_decline_reason = reason;
        report.decline_reason = "multipatch_topology_split_validation_declined";
        return Ok((baseline, report));
    }
    let mut best_exact: Option<Kp20TopologySplitTrial> = None;
    for strip_mode in Kp20TopologySplitStripMode::all() {
        for inner_mode in Kp20TopologySplitInnerMode::all() {
            for scale in [1.0, 0.75, 0.5] {
                report.evaluated += 1;
                let applied = apply_multipatch_candidate(
                    &baseline.macro_patches,
                    &baseline.summary,
                    &strips,
                    &interfaces,
                    &inner,
                    strip_mode.mode(),
                    inner_mode.mode(),
                    scale,
                );
                if applied.changed_piece_count == 0 {
                    report.rejected += 1;
                    continue;
                }
                let fast = fast_candidate_summary_from_changed_patches(
                    tensor_border,
                    applied.candidate.clone(),
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
                    applied.candidate.clone(),
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
                let trial = Kp20TopologySplitTrial {
                    macro_patches: applied.candidate,
                    exact,
                    changed_piece_count: applied.changed_piece_count,
                    strip_mode,
                    inner_mode,
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
    report.changed_piece_count = best.changed_piece_count;
    report.strip_mode = best.strip_mode.label();
    report.inner_mode = best.inner_mode.label();
    report.scale = best.scale;
    report.decline_reason = "none";
    let mut refined = baseline;
    refined.summary = best.exact.0;
    refined.validation = best.exact.1;
    refined.macro_patches = best.macro_patches;
    Ok((refined, report))
}

#[derive(Debug, Clone, Copy)]
struct Kp20DirectMultipatchMeasurement {
    patch_count: usize,
    fold_count: usize,
}

fn measure_direct_multipatch_surfaces(
    strips: &super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport,
    inner: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport,
) -> Kp20DirectMultipatchMeasurement {
    let mut measurement = Kp20DirectMultipatchMeasurement {
        patch_count: 0,
        fold_count: 0,
    };
    for strip in &strips.strip_surfaces {
        measure_direct_patch(strip.surface, &mut measurement);
    }
    for inner_patch in &inner.inner_patches {
        measure_direct_patch(inner_patch.surface, &mut measurement);
    }
    measurement
}

fn measure_direct_patch(patch: Kp20Bi6Patch, measurement: &mut Kp20DirectMultipatchMeasurement) {
    measurement.patch_count += 1;
    if gregory_from_bi6(patch)
        .self_intersection_report(SELF_INT_TOL)
        .intersects()
    {
        measurement.fold_count += 1;
    }
}

fn gregory_from_bi6(patch: Kp20Bi6Patch) -> GregoryPatch {
    GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
        grid: elevate_bi6_grid_to_bi9(patch.cps),
    }))
}

#[derive(Debug, Clone)]
struct Kp20TopologySplitTrial {
    macro_patches: Vec<Kp20MacroPatch>,
    exact: FoldSummary,
    changed_piece_count: usize,
    strip_mode: Kp20TopologySplitStripMode,
    inner_mode: Kp20TopologySplitInnerMode,
    scale: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kp20TopologySplitStripMode {
    FairBlend,
    WorstEdgeParameterExpand,
    WorstEdgeParameterContract,
}

impl Kp20TopologySplitStripMode {
    fn all() -> [Self; 3] {
        [
            Self::FairBlend,
            Self::WorstEdgeParameterExpand,
            Self::WorstEdgeParameterContract,
        ]
    }

    fn mode(self) -> Kp20StripInteriorOptimizerMode {
        match self {
            Self::FairBlend => Kp20StripInteriorOptimizerMode::FairBlend,
            Self::WorstEdgeParameterExpand => {
                Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
            }
            Self::WorstEdgeParameterContract => {
                Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract
            }
        }
    }

    fn label(self) -> &'static str {
        self.mode().label()
    }
}

#[derive(Debug, Clone, Copy)]
enum Kp20TopologySplitInnerMode {
    SupportPropagation,
    SideAlignedPropagatedInnerPatch,
    ExactFoldSetSharedSupportInnerExpand,
    ExactFoldSetSharedSupportInnerContract,
}

impl Kp20TopologySplitInnerMode {
    fn all() -> [Self; 4] {
        [
            Self::SupportPropagation,
            Self::SideAlignedPropagatedInnerPatch,
            Self::ExactFoldSetSharedSupportInnerExpand,
            Self::ExactFoldSetSharedSupportInnerContract,
        ]
    }

    fn mode(self) -> Kp20InnerConstructionMode {
        match self {
            Self::SupportPropagation => Kp20InnerConstructionMode::SupportPropagation,
            Self::SideAlignedPropagatedInnerPatch => {
                Kp20InnerConstructionMode::SideAlignedPropagatedInnerPatch
            }
            Self::ExactFoldSetSharedSupportInnerExpand => {
                Kp20InnerConstructionMode::ExactFoldSetSharedSupportInnerExpand
            }
            Self::ExactFoldSetSharedSupportInnerContract => {
                Kp20InnerConstructionMode::ExactFoldSetSharedSupportInnerContract
            }
        }
    }

    fn label(self) -> &'static str {
        self.mode().label()
    }
}

fn trial_is_better(candidate: &FoldSummary, best: &FoldSummary) -> bool {
    candidate.0.fold_count < best.0.fold_count
        || (candidate.0.fold_count == best.0.fold_count
            && candidate.0.min_jacobian_normal_cos > best.0.min_jacobian_normal_cos)
}
