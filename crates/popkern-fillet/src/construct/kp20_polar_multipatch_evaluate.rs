//! Candidate-only KP20 multipatch fold evaluation.

use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, validation_is_clean,
};
use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_validate::Kp20MultiPatchValidationReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchEvaluationReport {
    pub(super) attempted_inner_patch_count: usize,
    pub(super) changed_piece_count: usize,
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
    inner: &Kp20MultiPatchInnerBuildReport,
) -> Result<Kp20MultiPatchEvaluationReport, String> {
    if validation.decline_reason.is_some() {
        return Ok(declined_report(
            inner,
            0,
            Some("multipatch_validation_declined"),
        ));
    }

    let mut candidate = accepted_patches.to_vec();
    let changed_patches = apply_inner_patches_to_macro_candidate(&mut candidate, inner);
    let changed_piece_count = changed_piece_count(&changed_patches);
    if changed_piece_count == 0 {
        return Ok(declined_report(
            inner,
            changed_piece_count,
            Some("multipatch_no_changed_pieces"),
        ));
    }

    let fast = fast_candidate_summary_from_changed_patches(
        tensor_border,
        candidate.clone(),
        accepted_summary,
        &changed_patches,
        accepted_summary.fold_count,
    )?;
    let fast_clean = validation_is_clean(fast.1);
    let fast_better = fast_clean && fast.0.fold_count < accepted_summary.fold_count;
    if !fast_better {
        return Ok(Kp20MultiPatchEvaluationReport {
            attempted_inner_patch_count: inner.built_inner_patch_count,
            changed_piece_count,
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
        candidate,
        accepted_summary,
        &changed_patches,
        accepted_summary.fold_count,
    )?;
    let exact_clean = validation_is_clean(exact.1);
    let accepted = exact_clean && exact.0.fold_count < accepted_summary.fold_count;
    Ok(Kp20MultiPatchEvaluationReport {
        attempted_inner_patch_count: inner.built_inner_patch_count,
        changed_piece_count,
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
    inner: &Kp20MultiPatchInnerBuildReport,
    changed_piece_count: usize,
    decline_reason: Option<&'static str>,
) -> Kp20MultiPatchEvaluationReport {
    Kp20MultiPatchEvaluationReport {
        attempted_inner_patch_count: inner.built_inner_patch_count,
        changed_piece_count,
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

fn apply_inner_patches_to_macro_candidate(
    macro_patches: &mut [Kp20MacroPatch],
    inner: &Kp20MultiPatchInnerBuildReport,
) -> Vec<[bool; 4]> {
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    for patch in &inner.inner_patches {
        if patch.sector_index >= macro_patches.len() {
            continue;
        }
        for piece in [Kp20PieceId::NearOriginLeft, Kp20PieceId::NearOriginRight] {
            *macro_patches[patch.sector_index].piece_mut(piece) = patch.surface;
            changed_patches[patch.sector_index][piece.index()] = true;
        }
    }
    changed_patches
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
    fn adapter_marks_near_origin_pieces_for_inner_candidate() {
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

        let changed = apply_inner_patches_to_macro_candidate(&mut candidate, &inner);

        assert!(changed[0][Kp20PieceId::NearOriginLeft.index()]);
        assert!(changed[0][Kp20PieceId::NearOriginRight.index()]);
        assert_eq!(changed_piece_count(&changed), 2);
        assert_eq!(
            candidate[0].piece(Kp20PieceId::NearOriginLeft).cps[0][3],
            Point3::new(1.0, 2.0, 3.0)
        );
    }
}
