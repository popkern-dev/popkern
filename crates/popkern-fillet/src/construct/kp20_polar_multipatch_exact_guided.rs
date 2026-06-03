//! Exact-guided KP20 multipatch retry selection.

use popkern_core::geometry::gregory::TensorBorderBi3;

use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_apply::{Kp20InnerConstructionMode, apply_multipatch_candidate};
use super::kp20_polar_multipatch_evaluate::{
    Kp20MultiPatchTrialCandidate, trial_is_better_by_exact,
};
use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport;
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;
use super::kp20_polar_multipatch_strip_optimize::Kp20StripInteriorOptimizerMode;
use popkern_core::geometry::gregory::{Kp20Bi6ValidationReport, Kp20MacroPatch};

const EXACT_GUIDED_DAMPING_SCALES: [f64; 3] = [2.0, 1.5, 1.0];

pub(super) struct Kp20ExactGuidedTrial {
    pub(super) trial: Kp20MultiPatchTrialCandidate,
    pub(super) exact: (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
}

pub(super) struct Kp20ExactGuidedSelection {
    pub(super) accepted: Option<Kp20ExactGuidedTrial>,
    pub(super) best_rejected: Option<Kp20ExactGuidedTrial>,
}

pub(super) fn select_exact_guided_trial(
    tensor_border: &TensorBorderBi3,
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    exact_guidance: &Kp20SamplerFoldSummary,
    strips: &Kp20MultiPatchStripBuildReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    trial_count: &mut usize,
) -> Result<Kp20ExactGuidedSelection, String> {
    let mut best_rejected = None;
    for strip_mode in Kp20StripInteriorOptimizerMode::exact_guided_all() {
        for inner_mode in Kp20InnerConstructionMode::exact_guided_all() {
            for scale in EXACT_GUIDED_DAMPING_SCALES {
                *trial_count += 1;
                let applied = apply_multipatch_candidate(
                    accepted_patches,
                    exact_guidance,
                    strips,
                    interfaces,
                    inner,
                    strip_mode,
                    inner_mode,
                    scale,
                );
                if applied.changed_piece_count == 0 {
                    continue;
                }
                let fast = fast_candidate_summary_from_changed_patches(
                    tensor_border,
                    applied.candidate.clone(),
                    accepted_summary,
                    &applied.changed_patches,
                    accepted_summary.fold_count,
                )?;
                if !(validation_is_clean(fast.1) && fast.0.fold_count < accepted_summary.fold_count)
                {
                    continue;
                }
                let exact = candidate_summary_from_changed_patches(
                    tensor_border,
                    applied.candidate.clone(),
                    accepted_summary,
                    &applied.changed_patches,
                    accepted_summary.fold_count,
                )?;
                let trial = Kp20ExactGuidedTrial {
                    trial: Kp20MultiPatchTrialCandidate {
                        applied,
                        strip_mode,
                        inner_mode,
                        scale,
                        fast,
                    },
                    exact,
                };
                if validation_is_clean(trial.exact.1)
                    && trial.exact.0.fold_count < accepted_summary.fold_count
                {
                    return Ok(Kp20ExactGuidedSelection {
                        accepted: Some(trial),
                        best_rejected,
                    });
                }
                if best_rejected
                    .as_ref()
                    .is_none_or(|best| trial_is_better_by_exact(&trial, best))
                {
                    best_rejected = Some(trial);
                }
            }
        }
    }
    Ok(Kp20ExactGuidedSelection {
        accepted: None,
        best_rejected,
    })
}
