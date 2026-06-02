//! Candidate-only KP20 multipatch fold evaluation.

use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;
use super::kp20_polar_multipatch_validate::Kp20MultiPatchValidationReport;

const STRIP_INTERIOR_DAMPING_SCALES: [f64; 4] = [1.0, 0.75, 0.5, 0.25];

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchEvaluationReport {
    pub(super) attempted_strip_patch_count: usize,
    pub(super) attempted_inner_patch_count: usize,
    pub(super) changed_strip_piece_count: usize,
    pub(super) changed_inner_piece_count: usize,
    pub(super) changed_piece_count: usize,
    pub(super) strip_optimizer_trial_count: usize,
    pub(super) strip_optimizer_selected_scale: Option<f64>,
    pub(super) fast_evaluated: bool,
    pub(super) exact_evaluated: bool,
    pub(super) fast_fold_count: Option<usize>,
    pub(super) exact_fold_count: Option<usize>,
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
            Some("multipatch_validation_declined"),
        ));
    }

    let mut best_trial: Option<Kp20StripInteriorTrialCandidate> = None;
    let mut trial_count = 0_usize;
    for scale in STRIP_INTERIOR_DAMPING_SCALES {
        trial_count += 1;
        let mut candidate = accepted_patches.to_vec();
        let mut changed_patches = vec![[false; 4]; candidate.len()];
        let changed_strip_piece_count = apply_strip_patches_to_macro_candidate(
            accepted_patches,
            &mut candidate,
            &mut changed_patches,
            strips,
            scale,
        );
        let changed_inner_piece_count =
            apply_inner_patches_to_macro_candidate(&mut candidate, &mut changed_patches, inner);
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
            changed_strip_piece_count,
            changed_inner_piece_count,
            changed_piece_count,
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
    let Some(best_trial) = best_trial else {
        return Ok(declined_report(
            strips,
            inner,
            0,
            0,
            trial_count,
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
            Some("multipatch_no_changed_pieces"),
        ));
    }

    let fast = best_trial.fast.clone();
    let fast_clean = validation_is_clean(fast.1);
    let fast_better = fast_clean && fast.0.fold_count < accepted_summary.fold_count;
    if !fast_better {
        return Ok(Kp20MultiPatchEvaluationReport {
            attempted_strip_patch_count: strips.attempted_strip_patch_count,
            attempted_inner_patch_count: inner.attempted_inner_patch_count,
            changed_strip_piece_count,
            changed_inner_piece_count,
            changed_piece_count,
            strip_optimizer_trial_count: trial_count,
            strip_optimizer_selected_scale: Some(best_trial.scale),
            fast_evaluated: true,
            exact_evaluated: false,
            fast_fold_count: Some(fast.0.fold_count),
            exact_fold_count: None,
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
    Ok(Kp20MultiPatchEvaluationReport {
        attempted_strip_patch_count: strips.attempted_strip_patch_count,
        attempted_inner_patch_count: inner.attempted_inner_patch_count,
        changed_strip_piece_count,
        changed_inner_piece_count,
        changed_piece_count,
        strip_optimizer_trial_count: trial_count,
        strip_optimizer_selected_scale: Some(best_trial.scale),
        fast_evaluated: true,
        exact_evaluated: true,
        fast_fold_count: Some(fast.0.fold_count),
        exact_fold_count: Some(exact.0.fold_count),
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
        strip_optimizer_selected_scale,
        fast_evaluated: false,
        exact_evaluated: false,
        fast_fold_count: None,
        exact_fold_count: None,
        fast_min_jacobian_normal_cos: None,
        exact_min_jacobian_normal_cos: None,
        validation_clean: false,
        accepted: false,
        decline_reason,
    }
}

#[derive(Debug, Clone)]
struct Kp20StripInteriorTrialCandidate {
    candidate: Vec<Kp20MacroPatch>,
    changed_patches: Vec<[bool; 4]>,
    changed_strip_piece_count: usize,
    changed_inner_piece_count: usize,
    changed_piece_count: usize,
    scale: f64,
    fast: (
        Kp20SamplerFoldSummary,
        popkern_core::geometry::gregory::Kp20Bi6ValidationReport,
    ),
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
    strips: &Kp20MultiPatchStripBuildReport,
    interior_scale: f64,
) -> usize {
    let mut changed_count = 0_usize;
    for strip in &strips.strip_surfaces {
        if strip.sector_index >= macro_patches.len() {
            continue;
        }
        *macro_patches[strip.sector_index].piece_mut(strip.strip_piece) =
            scaled_strip_surface(accepted_patches, strip, interior_scale);
        if !changed_patches[strip.sector_index][strip.strip_piece.index()] {
            changed_count += 1;
        }
        changed_patches[strip.sector_index][strip.strip_piece.index()] = true;
    }
    changed_count
}

fn scaled_strip_surface(
    accepted_patches: &[Kp20MacroPatch],
    strip: &super::kp20_polar_multipatch_strip::Kp20MultiPatchStripSurface,
    interior_scale: f64,
) -> popkern_core::geometry::gregory::Kp20Bi6Patch {
    let mut surface = strip.surface;
    if let Some(accepted) = accepted_patches
        .get(strip.sector_index)
        .map(|patch| *patch.piece(strip.strip_piece))
    {
        for row in 2..=4 {
            for column in 0..7 {
                let base = accepted.cps[row][column];
                let target = strip.surface.cps[row][column];
                surface.cps[row][column] = base + (target - base) * interior_scale;
            }
        }
    }
    surface
}

fn apply_inner_patches_to_macro_candidate(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    inner: &Kp20MultiPatchInnerBuildReport,
) -> usize {
    let mut changed_count = 0_usize;
    for patch in &inner.inner_patches {
        if patch.sector_index >= macro_patches.len() {
            continue;
        }
        for piece in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            *macro_patches[patch.sector_index].piece_mut(piece) = patch.surface;
            if !changed_patches[patch.sector_index][piece.index()] {
                changed_count += 1;
            }
            changed_patches[patch.sector_index][piece.index()] = true;
        }
    }
    changed_count
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
    use crate::construct::kp20_polar_multipatch_inner::build_multipatch_inner_patches;
    use crate::construct::kp20_polar_multipatch_interface::build_multipatch_interfaces;
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
        let interfaces = build_multipatch_interfaces(&strips);
        let inner = build_multipatch_inner_patches(&inventory, &interfaces);
        let mut candidate = vec![Kp20MacroPatch::zero()];
        let accepted = candidate.clone();
        let mut changed = vec![[false; 4]; candidate.len()];

        let changed_strips = apply_strip_patches_to_macro_candidate(
            &accepted,
            &mut candidate,
            &mut changed,
            &strips,
            1.0,
        );
        let changed_inner =
            apply_inner_patches_to_macro_candidate(&mut candidate, &mut changed, &inner);

        assert_eq!(changed_strips, 1);
        assert_eq!(changed_inner, 2);
        assert!(changed[0][Kp20PieceId::NearMidpointLeft.index()]);
        assert!(changed[0][Kp20PieceId::NearOriginLeft.index()]);
        assert!(changed[0][Kp20PieceId::NearOriginRight.index()]);
        assert_eq!(changed_piece_count(&changed), 3);
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
    fn strip_damping_preserves_constrained_rows() {
        let mut accepted = Kp20MacroPatch::zero();
        let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
        piece.cps[1][2] = Point3::new(1.0, 0.0, 0.0);
        piece.cps[2][2] = Point3::new(2.0, 0.0, 0.0);
        piece.cps[5][2] = Point3::new(5.0, 0.0, 0.0);
        piece.cps[6][2] = Point3::new(6.0, 0.0, 0.0);
        let mut strip_surface = *accepted.piece(Kp20PieceId::NearMidpointLeft);
        strip_surface.cps[2][2] = Point3::new(10.0, 0.0, 0.0);
        let strip = super::super::kp20_polar_multipatch_strip::Kp20MultiPatchStripSurface {
            patch_name: "strip".to_string(),
            sector_index: 0,
            strip_piece: Kp20PieceId::NearMidpointLeft,
            surface: strip_surface,
        };

        let scaled = scaled_strip_surface(&[accepted], &strip, 0.25);

        assert_eq!(scaled.cps[1][2], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(scaled.cps[2][2], Point3::new(4.0, 0.0, 0.0));
        assert_eq!(scaled.cps[5][2], Point3::new(5.0, 0.0, 0.0));
        assert_eq!(scaled.cps[6][2], Point3::new(6.0, 0.0, 0.0));
    }
}
