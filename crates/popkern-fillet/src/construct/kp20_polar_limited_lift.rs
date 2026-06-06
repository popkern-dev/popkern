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
use super::super::kp20_polar_support_reparameterization::Kp20SupportDomainReport;
use super::labels::{piece_counts_label, validation_label};
use super::support_domain::support_domain_basis_candidate;

pub(super) struct Kp20LimitedLiftCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) report: Kp20SupportNormalLiftReport,
    strategy: &'static str,
    support_domain: Option<Kp20SupportDomainReport>,
    elapsed_ms: u128,
}

impl Kp20LimitedLiftCandidate {
    pub(super) fn detail(
        &self,
        default: &Kp20SamplerFoldSummary,
        lifted: &Kp20SamplerFoldSummary,
    ) -> String {
        let support_domain = self
            .support_domain
            .as_ref()
            .map(|report| {
                format!(
                    ";support_domain_evaluated={};support_domain_accepted={};support_domain_rejected={};support_domain_start_folds={};support_domain_end_folds={};support_domain_decline_reason={}",
                    report.evaluated,
                    report.accepted,
                    report.rejected,
                    report.start_folds,
                    report.end_folds,
                    report.decline_reason,
                )
            })
            .unwrap_or_default();
        format!(
            "status=sampled;strategy={};default_folds={};folds={};piece_folds={};validation:{};optimizer_evaluated={};optimizer_accepted={};accepted_steps={};elapsed_ms={}{}",
            self.strategy,
            default.fold_count,
            lifted.fold_count,
            piece_counts_label(lifted.piece_fold_counts),
            validation_label(self.report.validation),
            self.report.optimizer_stats.evaluated_candidates,
            self.report.optimizer_stats.accepted_candidates,
            self.report.accepted_steps.len(),
            self.elapsed_ms,
            support_domain,
        )
    }
}

pub(super) fn try_limited_nonoverlap_lift(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20LimitedLiftCandidate, String> {
    match try_support_domain_limited_lift(star_point, tensor_border) {
        Ok(candidate) => return Ok(candidate),
        Err(support_domain_error) => default_limited_lift(star_point, tensor_border)
            .map_err(|direct_error| format!("{support_domain_error};{direct_error}")),
    }
}

fn try_support_domain_limited_lift(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20LimitedLiftCandidate, String> {
    let start = Instant::now();
    let (macro_patches, support_domain, _, _, _) =
        support_domain_basis_candidate(star_point, tensor_border)?;
    let report = lift_candidate_from_macro_patches_with_profile(
        star_point,
        tensor_border,
        macro_patches,
        Kp20LiftSearchProfile::General,
        8,
    )?;
    let patches = gregory_patches_from_macro_patches(&report.macro_patches);
    Ok(Kp20LimitedLiftCandidate {
        patches,
        report,
        strategy: "support_domain_basis_lift",
        support_domain: Some(support_domain),
        elapsed_ms: start.elapsed().as_millis(),
    })
}

fn default_limited_lift(
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
        strategy: "direct_support_normal_lift",
        support_domain: None,
        elapsed_ms: start.elapsed().as_millis(),
    })
}
