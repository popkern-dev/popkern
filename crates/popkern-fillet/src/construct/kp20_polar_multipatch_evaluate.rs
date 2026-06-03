//! Candidate-only KP20 multipatch fold evaluation.

use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_apply::{
    Kp20AppliedMultiPatchCandidate, Kp20InnerConstructionMode, apply_multipatch_candidate,
};
use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport;
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;
use super::kp20_polar_multipatch_strip_optimize::Kp20StripInteriorOptimizerMode;
use super::kp20_polar_multipatch_validate::Kp20MultiPatchValidationReport;

const STRIP_INTERIOR_DAMPING_SCALES: [f64; 6] = [2.0, 1.5, 1.0, 0.75, 0.5, 0.25];
const EXACT_GUIDED_DAMPING_SCALES: [f64; 3] = [2.0, 1.5, 1.0];

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchEvaluationReport {
    pub(super) attempted_strip_patch_count: usize,
    pub(super) attempted_inner_patch_count: usize,
    pub(super) changed_strip_piece_count: usize,
    pub(super) changed_inner_piece_count: usize,
    pub(super) changed_piece_count: usize,
    pub(super) strip_optimizer_trial_count: usize,
    pub(super) strip_optimizer_mode: Option<&'static str>,
    pub(super) strip_optimizer_selected_scale: Option<f64>,
    pub(super) inner_construction_mode: Option<&'static str>,
    pub(super) fast_evaluated: bool,
    pub(super) exact_evaluated: bool,
    pub(super) fast_fold_count: Option<usize>,
    pub(super) exact_fold_count: Option<usize>,
    pub(super) fast_support_fold_count: Option<usize>,
    pub(super) fast_inner_fold_count: Option<usize>,
    pub(super) exact_support_fold_count: Option<usize>,
    pub(super) exact_inner_fold_count: Option<usize>,
    pub(super) fast_piece_fold_counts: Option<[usize; 4]>,
    pub(super) exact_piece_fold_counts: Option<[usize; 4]>,
    pub(super) fast_min_jacobian_normal_cos: Option<f64>,
    pub(super) exact_min_jacobian_normal_cos: Option<f64>,
    pub(super) fast_worst_fold_site: Option<Kp20FoldSite>,
    pub(super) exact_worst_fold_site: Option<Kp20FoldSite>,
    pub(super) fast_validation: Option<Kp20CandidateValidationResiduals>,
    pub(super) exact_validation: Option<Kp20CandidateValidationResiduals>,
    pub(super) validation_clean: bool,
    pub(super) accepted: bool,
    pub(super) decline_reason: Option<&'static str>,
}

pub(super) fn evaluate_multipatch_candidate(
    tensor_border: &popkern_core::geometry::gregory::TensorBorderBi3,
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    validation: &Kp20MultiPatchValidationReport,
    strips: &Kp20MultiPatchStripBuildReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
) -> Result<Kp20MultiPatchEvaluationReport, String> {
    if validation.decline_reason.is_some() {
        return Ok(declined_report(
            strips,
            inner,
            0,
            0,
            0,
            None,
            None,
            None,
            Some("multipatch_validation_declined"),
        ));
    }

    let mut best_trial: Option<Kp20MultiPatchTrialCandidate> = None;
    let mut trial_count = 0_usize;
    for strip_mode in Kp20StripInteriorOptimizerMode::all() {
        for inner_mode in Kp20InnerConstructionMode::all() {
            for scale in STRIP_INTERIOR_DAMPING_SCALES {
                trial_count += 1;
                let applied = apply_multipatch_candidate(
                    accepted_patches,
                    accepted_summary,
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
                let trial = Kp20MultiPatchTrialCandidate {
                    applied,
                    strip_mode,
                    inner_mode,
                    scale,
                    fast,
                };
                if best_trial
                    .as_ref()
                    .is_none_or(|best| trial_is_better(&trial, best))
                {
                    best_trial = Some(trial);
                }
            }
        }
    }
    let Some(best_trial) = best_trial else {
        return Ok(declined_report(
            strips,
            inner,
            0,
            0,
            trial_count,
            None,
            None,
            None,
            Some("multipatch_no_changed_pieces"),
        ));
    };

    let fast = best_trial.fast.clone();
    let fast_clean = validation_is_clean(fast.1);
    let fast_ownership = fold_ownership_counts(&fast.0);
    if !(fast_clean && fast.0.fold_count < accepted_summary.fold_count) {
        return Ok(Kp20MultiPatchEvaluationReport {
            attempted_strip_patch_count: strips.attempted_strip_patch_count,
            attempted_inner_patch_count: inner.attempted_inner_patch_count,
            changed_strip_piece_count: best_trial.applied.changed_strip_piece_count,
            changed_inner_piece_count: best_trial.applied.changed_inner_piece_count,
            changed_piece_count: best_trial.applied.changed_piece_count,
            strip_optimizer_trial_count: trial_count,
            strip_optimizer_mode: Some(best_trial.strip_mode.label()),
            strip_optimizer_selected_scale: Some(best_trial.scale),
            inner_construction_mode: Some(best_trial.inner_mode.label()),
            fast_evaluated: true,
            exact_evaluated: false,
            fast_fold_count: Some(fast.0.fold_count),
            exact_fold_count: None,
            fast_support_fold_count: Some(fast_ownership.support_fold_count),
            fast_inner_fold_count: Some(fast_ownership.inner_fold_count),
            exact_support_fold_count: None,
            exact_inner_fold_count: None,
            fast_piece_fold_counts: Some(fast.0.piece_fold_counts),
            exact_piece_fold_counts: None,
            fast_min_jacobian_normal_cos: Some(fast.0.min_jacobian_normal_cos),
            exact_min_jacobian_normal_cos: None,
            fast_worst_fold_site: fast.0.worst_fold_site(),
            exact_worst_fold_site: None,
            fast_validation: Some(Kp20CandidateValidationResiduals::from(fast.1)),
            exact_validation: None,
            validation_clean: fast_clean,
            accepted: false,
            decline_reason: Some(if fast_clean {
                "multipatch_no_fast_fold_reduction"
            } else {
                "multipatch_fast_validation_not_clean"
            }),
        });
    }

    let exact = candidate_summary_from_changed_patches(
        tensor_border,
        best_trial.applied.candidate,
        accepted_summary,
        &best_trial.applied.changed_patches,
        accepted_summary.fold_count,
    )?;
    let exact_clean = validation_is_clean(exact.1);
    let exact_ownership = fold_ownership_counts(&exact.0);
    let accepted = exact_clean && exact.0.fold_count < accepted_summary.fold_count;
    if !accepted
        && exact_clean
        && let Some(refined) = exact_guided_accepted_trial(
            tensor_border,
            accepted_patches,
            accepted_summary,
            &exact.0,
            strips,
            interfaces,
            inner,
            &mut trial_count,
        )?
    {
        let refined_fast_ownership = fold_ownership_counts(&refined.trial.fast.0);
        let refined_exact_ownership = fold_ownership_counts(&refined.exact.0);
        return Ok(Kp20MultiPatchEvaluationReport {
            attempted_strip_patch_count: strips.attempted_strip_patch_count,
            attempted_inner_patch_count: inner.attempted_inner_patch_count,
            changed_strip_piece_count: refined.trial.applied.changed_strip_piece_count,
            changed_inner_piece_count: refined.trial.applied.changed_inner_piece_count,
            changed_piece_count: refined.trial.applied.changed_piece_count,
            strip_optimizer_trial_count: trial_count,
            strip_optimizer_mode: Some(refined.trial.strip_mode.label()),
            strip_optimizer_selected_scale: Some(refined.trial.scale),
            inner_construction_mode: Some(refined.trial.inner_mode.label()),
            fast_evaluated: true,
            exact_evaluated: true,
            fast_fold_count: Some(refined.trial.fast.0.fold_count),
            exact_fold_count: Some(refined.exact.0.fold_count),
            fast_support_fold_count: Some(refined_fast_ownership.support_fold_count),
            fast_inner_fold_count: Some(refined_fast_ownership.inner_fold_count),
            exact_support_fold_count: Some(refined_exact_ownership.support_fold_count),
            exact_inner_fold_count: Some(refined_exact_ownership.inner_fold_count),
            fast_piece_fold_counts: Some(refined.trial.fast.0.piece_fold_counts),
            exact_piece_fold_counts: Some(refined.exact.0.piece_fold_counts),
            fast_min_jacobian_normal_cos: Some(refined.trial.fast.0.min_jacobian_normal_cos),
            exact_min_jacobian_normal_cos: Some(refined.exact.0.min_jacobian_normal_cos),
            fast_worst_fold_site: refined.trial.fast.0.worst_fold_site(),
            exact_worst_fold_site: refined.exact.0.worst_fold_site(),
            fast_validation: Some(Kp20CandidateValidationResiduals::from(refined.trial.fast.1)),
            exact_validation: Some(Kp20CandidateValidationResiduals::from(refined.exact.1)),
            validation_clean: true,
            accepted: true,
            decline_reason: None,
        });
    }
    Ok(Kp20MultiPatchEvaluationReport {
        attempted_strip_patch_count: strips.attempted_strip_patch_count,
        attempted_inner_patch_count: inner.attempted_inner_patch_count,
        changed_strip_piece_count: best_trial.applied.changed_strip_piece_count,
        changed_inner_piece_count: best_trial.applied.changed_inner_piece_count,
        changed_piece_count: best_trial.applied.changed_piece_count,
        strip_optimizer_trial_count: trial_count,
        strip_optimizer_mode: Some(best_trial.strip_mode.label()),
        strip_optimizer_selected_scale: Some(best_trial.scale),
        inner_construction_mode: Some(best_trial.inner_mode.label()),
        fast_evaluated: true,
        exact_evaluated: true,
        fast_fold_count: Some(fast.0.fold_count),
        exact_fold_count: Some(exact.0.fold_count),
        fast_support_fold_count: Some(fast_ownership.support_fold_count),
        fast_inner_fold_count: Some(fast_ownership.inner_fold_count),
        exact_support_fold_count: Some(exact_ownership.support_fold_count),
        exact_inner_fold_count: Some(exact_ownership.inner_fold_count),
        fast_piece_fold_counts: Some(fast.0.piece_fold_counts),
        exact_piece_fold_counts: Some(exact.0.piece_fold_counts),
        fast_min_jacobian_normal_cos: Some(fast.0.min_jacobian_normal_cos),
        exact_min_jacobian_normal_cos: Some(exact.0.min_jacobian_normal_cos),
        fast_worst_fold_site: fast.0.worst_fold_site(),
        exact_worst_fold_site: exact.0.worst_fold_site(),
        fast_validation: Some(Kp20CandidateValidationResiduals::from(fast.1)),
        exact_validation: Some(Kp20CandidateValidationResiduals::from(exact.1)),
        validation_clean: exact_clean,
        accepted,
        decline_reason: if accepted {
            None
        } else if exact_clean {
            Some("multipatch_no_exact_fold_reduction")
        } else {
            Some("multipatch_exact_validation_not_clean")
        },
    })
}

struct Kp20ExactGuidedAcceptedTrial {
    trial: Kp20MultiPatchTrialCandidate,
    exact: (
        Kp20SamplerFoldSummary,
        popkern_core::geometry::gregory::Kp20Bi6ValidationReport,
    ),
}

fn exact_guided_accepted_trial(
    tensor_border: &popkern_core::geometry::gregory::TensorBorderBi3,
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    exact_guidance: &Kp20SamplerFoldSummary,
    strips: &Kp20MultiPatchStripBuildReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    trial_count: &mut usize,
) -> Result<Option<Kp20ExactGuidedAcceptedTrial>, String> {
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
                if validation_is_clean(exact.1) && exact.0.fold_count < accepted_summary.fold_count
                {
                    return Ok(Some(Kp20ExactGuidedAcceptedTrial {
                        trial: Kp20MultiPatchTrialCandidate {
                            applied,
                            strip_mode,
                            inner_mode,
                            scale,
                            fast,
                        },
                        exact,
                    }));
                }
            }
        }
    }
    Ok(None)
}

fn declined_report(
    strips: &Kp20MultiPatchStripBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    changed_strip_piece_count: usize,
    changed_piece_count: usize,
    strip_optimizer_trial_count: usize,
    strip_optimizer_mode: Option<&'static str>,
    strip_optimizer_selected_scale: Option<f64>,
    inner_construction_mode: Option<&'static str>,
    decline_reason: Option<&'static str>,
) -> Kp20MultiPatchEvaluationReport {
    Kp20MultiPatchEvaluationReport {
        attempted_strip_patch_count: strips.attempted_strip_patch_count,
        attempted_inner_patch_count: inner.attempted_inner_patch_count,
        changed_strip_piece_count,
        changed_inner_piece_count: 0,
        changed_piece_count,
        strip_optimizer_trial_count,
        strip_optimizer_mode,
        strip_optimizer_selected_scale,
        inner_construction_mode,
        fast_evaluated: false,
        exact_evaluated: false,
        fast_fold_count: None,
        exact_fold_count: None,
        fast_support_fold_count: None,
        fast_inner_fold_count: None,
        exact_support_fold_count: None,
        exact_inner_fold_count: None,
        fast_piece_fold_counts: None,
        exact_piece_fold_counts: None,
        fast_min_jacobian_normal_cos: None,
        exact_min_jacobian_normal_cos: None,
        fast_worst_fold_site: None,
        exact_worst_fold_site: None,
        fast_validation: None,
        exact_validation: None,
        validation_clean: false,
        accepted: false,
        decline_reason,
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20CandidateValidationResiduals {
    pub(super) tensor_border: f64,
    pub(super) internal_c2_position: f64,
    pub(super) internal_c2_first: f64,
    pub(super) origin_g2_strip: f64,
}

impl From<popkern_core::geometry::gregory::Kp20Bi6ValidationReport>
    for Kp20CandidateValidationResiduals
{
    fn from(report: popkern_core::geometry::gregory::Kp20Bi6ValidationReport) -> Self {
        Self {
            tensor_border: report.max_tensor_border_residual,
            internal_c2_position: report.max_internal_c2_position_residual,
            internal_c2_first: report.max_internal_c2_first_residual,
            origin_g2_strip: report.max_origin_g2_strip_residual,
        }
    }
}

#[derive(Debug, Clone)]
struct Kp20MultiPatchTrialCandidate {
    applied: Kp20AppliedMultiPatchCandidate,
    strip_mode: Kp20StripInteriorOptimizerMode,
    inner_mode: Kp20InnerConstructionMode,
    scale: f64,
    fast: (
        Kp20SamplerFoldSummary,
        popkern_core::geometry::gregory::Kp20Bi6ValidationReport,
    ),
}

fn trial_is_better(
    candidate: &Kp20MultiPatchTrialCandidate,
    best: &Kp20MultiPatchTrialCandidate,
) -> bool {
    let candidate_clean = validation_is_clean(candidate.fast.1);
    let best_clean = validation_is_clean(best.fast.1);
    if candidate_clean != best_clean {
        return candidate_clean;
    }
    candidate.fast.0.fold_count < best.fast.0.fold_count
        || (candidate.fast.0.fold_count == best.fast.0.fold_count
            && candidate.fast.0.min_jacobian_normal_cos > best.fast.0.min_jacobian_normal_cos)
}

#[derive(Debug, Clone, Copy)]
struct Kp20FoldOwnershipCounts {
    support_fold_count: usize,
    inner_fold_count: usize,
}

fn fold_ownership_counts(summary: &Kp20SamplerFoldSummary) -> Kp20FoldOwnershipCounts {
    let support_fold_count = summary.piece_fold_counts[Kp20PieceId::NearMidpointLeft.index()]
        + summary.piece_fold_counts[Kp20PieceId::NearMidpointRight.index()];
    let inner_fold_count = summary.piece_fold_counts[Kp20PieceId::NearOriginLeft.index()]
        + summary.piece_fold_counts[Kp20PieceId::NearOriginRight.index()];
    Kp20FoldOwnershipCounts {
        support_fold_count,
        inner_fold_count,
    }
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::gregory::Kp20PieceId;

    use super::*;

    #[test]
    fn fold_ownership_counts_split_support_and_inner_pieces() {
        let mut piece_fold_counts = [0_usize; 4];
        piece_fold_counts[Kp20PieceId::NearMidpointLeft.index()] = 2;
        piece_fold_counts[Kp20PieceId::NearMidpointRight.index()] = 3;
        piece_fold_counts[Kp20PieceId::NearOriginLeft.index()] = 1;
        piece_fold_counts[Kp20PieceId::NearOriginRight.index()] = 1;
        let summary = Kp20SamplerFoldSummary {
            fold_count: 7,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts,
            fold_details: Vec::new(),
            fold_sites: Vec::new(),
        };

        let ownership = fold_ownership_counts(&summary);

        assert_eq!(ownership.support_fold_count, 5);
        assert_eq!(ownership.inner_fold_count, 2);
    }
}
