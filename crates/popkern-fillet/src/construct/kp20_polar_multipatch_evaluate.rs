//! Candidate-only KP20 multipatch fold evaluation.

use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;
use super::kp20_polar_multipatch_strip_optimize::{
    Kp20StripInteriorOptimizerMode, optimized_strip_surface,
};
use super::kp20_polar_multipatch_validate::Kp20MultiPatchValidationReport;
use super::kp20_polar_support_surface_split::{
    propagate_support_pole_interfaces, support_side_index,
};

const STRIP_INTERIOR_DAMPING_SCALES: [f64; 4] = [1.0, 0.75, 0.5, 0.25];

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
            Some("multipatch_validation_declined"),
        ));
    }

    let mut best_trial: Option<Kp20StripInteriorTrialCandidate> = None;
    let mut trial_count = 0_usize;
    for mode in Kp20StripInteriorOptimizerMode::all() {
        for scale in STRIP_INTERIOR_DAMPING_SCALES {
            trial_count += 1;
            let mut candidate = accepted_patches.to_vec();
            let mut changed_patches = vec![[false; 4]; candidate.len()];
            let strip_apply = apply_strip_patches_to_macro_candidate(
                accepted_patches,
                &mut candidate,
                &mut changed_patches,
                accepted_summary,
                strips,
                mode,
                scale,
            );
            let changed_inner_piece_count = propagate_strip_owned_inner_pieces(
                &mut candidate,
                &mut changed_patches,
                strip_apply.support_changed,
            );
            let changed_piece_count = changed_piece_count(&changed_patches);
            if changed_piece_count == 0 {
                continue;
            }
            let fast = fast_candidate_summary_from_changed_patches(
                tensor_border,
                candidate.clone(),
                accepted_summary,
                &changed_patches,
                accepted_summary.fold_count,
            )?;
            let trial = Kp20StripInteriorTrialCandidate {
                candidate,
                changed_patches,
                changed_strip_piece_count: strip_apply.changed_strip_piece_count,
                changed_inner_piece_count,
                changed_piece_count,
                mode,
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
    let Some(best_trial) = best_trial else {
        return Ok(declined_report(
            strips,
            inner,
            0,
            0,
            trial_count,
            None,
            None,
            Some("multipatch_no_changed_pieces"),
        ));
    };
    let changed_strip_piece_count = best_trial.changed_strip_piece_count;
    let changed_inner_piece_count = best_trial.changed_inner_piece_count;
    let changed_piece_count = best_trial.changed_piece_count;
    if changed_piece_count == 0 {
        return Ok(declined_report(
            strips,
            inner,
            changed_strip_piece_count,
            changed_piece_count,
            trial_count,
            None,
            None,
            Some("multipatch_no_changed_pieces"),
        ));
    }

    let fast = best_trial.fast.clone();
    let fast_clean = validation_is_clean(fast.1);
    let fast_better = fast_clean && fast.0.fold_count < accepted_summary.fold_count;
    if !fast_better {
        let fast_ownership = fold_ownership_counts(&fast.0);
        return Ok(Kp20MultiPatchEvaluationReport {
            attempted_strip_patch_count: strips.attempted_strip_patch_count,
            attempted_inner_patch_count: inner.attempted_inner_patch_count,
            changed_strip_piece_count,
            changed_inner_piece_count,
            changed_piece_count,
            strip_optimizer_trial_count: trial_count,
            strip_optimizer_mode: Some(best_trial.mode.label()),
            strip_optimizer_selected_scale: Some(best_trial.scale),
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
        best_trial.candidate,
        accepted_summary,
        &best_trial.changed_patches,
        accepted_summary.fold_count,
    )?;
    let exact_clean = validation_is_clean(exact.1);
    let accepted = exact_clean && exact.0.fold_count < accepted_summary.fold_count;
    let fast_ownership = fold_ownership_counts(&fast.0);
    let exact_ownership = fold_ownership_counts(&exact.0);
    Ok(Kp20MultiPatchEvaluationReport {
        attempted_strip_patch_count: strips.attempted_strip_patch_count,
        attempted_inner_patch_count: inner.attempted_inner_patch_count,
        changed_strip_piece_count,
        changed_inner_piece_count,
        changed_piece_count,
        strip_optimizer_trial_count: trial_count,
        strip_optimizer_mode: Some(best_trial.mode.label()),
        strip_optimizer_selected_scale: Some(best_trial.scale),
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

fn declined_report(
    strips: &Kp20MultiPatchStripBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
    changed_strip_piece_count: usize,
    changed_piece_count: usize,
    strip_optimizer_trial_count: usize,
    strip_optimizer_mode: Option<&'static str>,
    strip_optimizer_selected_scale: Option<f64>,
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
        validation_clean: false,
        accepted: false,
        decline_reason,
    }
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

#[derive(Debug, Clone)]
struct Kp20StripInteriorTrialCandidate {
    candidate: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
    changed_strip_piece_count: usize,
    changed_inner_piece_count: usize,
    changed_piece_count: usize,
    mode: Kp20StripInteriorOptimizerMode,
    scale: f64,
    fast: (
        Kp20SamplerFoldSummary,
        popkern_core::geometry::gregory::Kp20Bi6ValidationReport,
    ),
}

struct Kp20StripApplyResult {
    changed_strip_piece_count: usize,
    support_changed: Vec<[bool; 2]>,
}

fn trial_is_better(
    candidate: &Kp20StripInteriorTrialCandidate,
    best: &Kp20StripInteriorTrialCandidate,
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

fn apply_strip_patches_to_macro_candidate(
    accepted_patches: &[Kp20MacroPatch],
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    accepted_summary: &Kp20SamplerFoldSummary,
    strips: &Kp20MultiPatchStripBuildReport,
    mode: Kp20StripInteriorOptimizerMode,
    interior_scale: f64,
) -> Kp20StripApplyResult {
    let mut changed_count = 0_usize;
    let mut support_changed = vec![[false; 2]; macro_patches.len()];
    for strip in &strips.strip_surfaces {
        if strip.sector_index >= macro_patches.len() {
            continue;
        }
        *macro_patches[strip.sector_index].piece_mut(strip.strip_piece) = optimized_strip_surface(
            accepted_patches,
            accepted_summary,
            strip,
            mode,
            interior_scale,
        );
        if !changed_patches[strip.sector_index][strip.strip_piece.index()] {
            changed_count += 1;
        }
        changed_patches[strip.sector_index][strip.strip_piece.index()] = true;
        support_changed[strip.sector_index][support_side_index(strip.strip_piece)] = true;
    }
    Kp20StripApplyResult {
        changed_strip_piece_count: changed_count,
        support_changed,
    }
}

fn propagate_strip_owned_inner_pieces(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    support_changed: Vec<[bool; 2]>,
) -> usize {
    let before = changed_patches.to_vec();
    propagate_support_pole_interfaces(macro_patches, changed_patches, support_changed);
    changed_patches
        .iter()
        .zip(before)
        .map(|(after, before)| {
            [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight]
                .into_iter()
                .filter(|piece| after[piece.index()] && !before[piece.index()])
                .count()
        })
        .sum()
}

fn changed_piece_count(changed_patches: &[[bool; 4]]) -> usize {
    changed_patches
        .iter()
        .flat_map(|pieces| pieces.iter())
        .filter(|changed| **changed)
        .count()
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
    use crate::construct::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;

    #[test]
    fn adapter_marks_support_and_near_origin_pieces_for_candidate() {
        let inventory = multipatch_inventory_from_fold_summary(&Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.5,
                v: 0.5,
                min_jacobian_ratio: -0.1,
                min_jacobian_normal_cos: -0.1,
            }],
        });
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
        let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch]);
        let mut candidate = vec![Kp20MacroPatch::zero()];
        let accepted = candidate.clone();
        let mut changed = vec![[false; 4]; candidate.len()];

        let strip_apply = apply_strip_patches_to_macro_candidate(
            &accepted,
            &mut candidate,
            &mut changed,
            &empty_fold_summary(),
            &strips,
            Kp20StripInteriorOptimizerMode::FairBlend,
            1.0,
        );
        let changed_inner = propagate_strip_owned_inner_pieces(
            &mut candidate,
            &mut changed,
            strip_apply.support_changed,
        );

        assert_eq!(strip_apply.changed_strip_piece_count, 1);
        assert_eq!(changed_inner, 1);
        assert!(changed[0][Kp20PieceId::NearMidpointLeft.index()]);
        assert!(changed[0][Kp20PieceId::NearOriginLeft.index()]);
        assert!(!changed[0][Kp20PieceId::NearOriginRight.index()]);
        assert_eq!(changed_piece_count(&changed), 2);
        assert_eq!(
            candidate[0].piece(Kp20PieceId::NearMidpointLeft).cps[6][3],
            Point3::new(1.0, 2.0, 3.0)
        );
        assert_eq!(
            candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[0][3],
            Point3::new(1.0, 2.0, 3.0)
        );
    }

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

    fn empty_fold_summary() -> Kp20SamplerFoldSummary {
        Kp20SamplerFoldSummary {
            fold_count: 0,
            min_jacobian_normal_cos: 1.0,
            piece_fold_counts: [0; 4],
            fold_details: Vec::new(),
            fold_sites: Vec::new(),
        }
    }
}
