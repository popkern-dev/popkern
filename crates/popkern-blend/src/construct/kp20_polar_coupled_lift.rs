//! Full KP20 nonoverlap-coupled lift candidate assembly.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, TensorBorderBi3};

use super::super::kp20_polar_fold_summary::{
    gregory_patches_from_macro_patches, validation_is_clean,
};
use super::super::kp20_polar_support_image::refine_support_side_reparameterization;
use super::super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_macro_patches,
};
use super::super::kp20_polar_support_reparameterization::{
    Kp20CoupledDomainReport, Kp20SplitDomainReport, Kp20SupportDomainRefineReport,
    Kp20SupportDomainReport, Kp20SupportReparameterizationReport, Kp20SupportTransitionReport,
};
use super::labels::validation_label;
use super::origin_inner_row::{
    Kp20OriginInnerRowReport, Kp20PairedInnerRowReport, Kp20SupportInnerRowReport,
    refine_origin_inner_interior_row, refine_paired_inner_interior_row,
    refine_support_inner_interior_row,
};
use super::support_domain::support_domain_basis_candidate;
use super::support_transition::refine_support_transition_sheet;
use super::timing::{Kp20DispatchStageTiming, timed};

pub(super) struct Kp20LiftedDispatchCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) report: Kp20SupportNormalLiftReport,
    pub(super) origin_inner: Kp20OriginInnerRowReport,
    pub(super) paired_inner: Kp20PairedInnerRowReport,
    pub(super) support_inner: Kp20SupportInnerRowReport,
    pub(super) support_reparameterization: Kp20SupportReparameterizationReport,
    pub(super) support_transition: Kp20SupportTransitionReport,
    pub(super) support_domain: Kp20SupportDomainReport,
    pub(super) support_domain_refine: Kp20SupportDomainRefineReport,
    pub(super) coupled_domain: Kp20CoupledDomainReport,
    pub(super) split_domain: Kp20SplitDomainReport,
    pub(super) timing: Kp20DispatchStageTiming,
}

pub(super) fn try_nonoverlap_coupled_lift_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20LiftedDispatchCandidate, String> {
    let mut timing = Kp20DispatchStageTiming::default();
    let (support_domain_result, elapsed) =
        timed(|| support_domain_basis_candidate(star_point, tensor_border));
    timing.support_domain_ms = elapsed;
    let (macro_patches, support_domain, support_domain_refine, coupled_domain, split_domain) =
        support_domain_result.map_err(|e| timed_error(e, timing))?;
    let (transition_result, elapsed) =
        timed(|| refine_support_transition_sheet(tensor_border, macro_patches));
    timing.support_transition_ms = elapsed;
    let (macro_patches, support_transition) =
        transition_result.map_err(|e| timed_error(e, timing))?;
    let (lift_result, elapsed) =
        timed(|| lift_candidate_from_macro_patches(star_point, tensor_border, macro_patches));
    timing.support_normal_lift_ms = elapsed;
    let report = lift_result.map_err(|e| timed_error(e, timing))?;
    let (origin_inner_result, elapsed) =
        timed(|| refine_origin_inner_interior_row(star_point, tensor_border, report));
    timing.origin_inner_ms = elapsed;
    let (report, origin_inner) = origin_inner_result.map_err(|e| timed_error(e, timing))?;
    let (paired_inner_result, elapsed) =
        timed(|| refine_paired_inner_interior_row(star_point, tensor_border, report));
    timing.paired_inner_ms = elapsed;
    let (report, paired_inner) = paired_inner_result.map_err(|e| timed_error(e, timing))?;
    let (support_inner_result, elapsed) =
        timed(|| refine_support_inner_interior_row(star_point, tensor_border, report));
    timing.support_inner_ms = elapsed;
    let (report, support_inner) = support_inner_result.map_err(|e| timed_error(e, timing))?;
    let (support_reparameterization_result, elapsed) =
        timed(|| refine_support_side_reparameterization(star_point, tensor_border, report));
    timing.support_reparameterization_ms = elapsed;
    let (report, support_reparameterization) =
        support_reparameterization_result.map_err(|e| timed_error(e, timing))?;
    if !validation_is_clean(report.validation) {
        return Err(format!(
            "validation residuals exceed tolerance: {};{}",
            validation_label(report.validation),
            timing.label()
        ));
    }
    let patches = gregory_patches_from_macro_patches(&report.macro_patches);
    Ok(Kp20LiftedDispatchCandidate {
        patches,
        report,
        origin_inner,
        paired_inner,
        support_inner,
        support_reparameterization,
        support_transition,
        support_domain,
        support_domain_refine,
        coupled_domain,
        split_domain,
        timing,
    })
}

fn timed_error(error: String, timing: Kp20DispatchStageTiming) -> String {
    format!("{error};{}", timing.label())
}
