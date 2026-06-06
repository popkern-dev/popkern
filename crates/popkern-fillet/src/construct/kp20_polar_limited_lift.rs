//! Bounded KP20 lift candidate for high-valence dispatch.

use std::time::Instant;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, TensorBorderBi3, assemble_kp20_bi6_cap,
};

use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::super::kp20_polar_fold_summary::gregory_patches_from_macro_patches;
use super::super::kp20_polar_support_lift::{
    Kp20LiftSearchProfile, Kp20SupportNormalLiftReport,
    lift_candidate_from_macro_patches_with_profile,
};
use super::labels::{piece_counts_label, validation_label};

pub(super) struct Kp20LimitedLiftCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) report: Kp20SupportNormalLiftReport,
    elapsed_ms: u128,
}

impl Kp20LimitedLiftCandidate {
    pub(super) fn detail(
        &self,
        default: &Kp20SamplerFoldSummary,
        lifted: &Kp20SamplerFoldSummary,
    ) -> String {
        format!(
            "status=sampled;strategy=direct_support_normal_lift;default_folds={};folds={};piece_folds={};validation:{};optimizer_evaluated={};optimizer_accepted={};accepted_steps={};elapsed_ms={}",
            default.fold_count,
            lifted.fold_count,
            piece_counts_label(lifted.piece_fold_counts),
            validation_label(self.report.validation),
            self.report.optimizer_stats.evaluated_candidates,
            self.report.optimizer_stats.accepted_candidates,
            self.report.accepted_steps.len(),
            self.elapsed_ms,
        )
    }
}

pub(super) fn try_limited_nonoverlap_lift(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20LimitedLiftCandidate, String> {
    let cap = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let start = Instant::now();
    let report = lift_candidate_from_macro_patches_with_profile(
        star_point,
        tensor_border,
        cap.macro_patches,
        Kp20LiftSearchProfile::General,
        8,
    )?;
    let patches = gregory_patches_from_macro_patches(&report.macro_patches);
    Ok(Kp20LimitedLiftCandidate {
        patches,
        report,
        elapsed_ms: start.elapsed().as_millis(),
    })
}
