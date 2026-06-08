//! Bounded KP20 lift candidate for high-valence dispatch.

use std::time::Instant;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, TensorBorderBi3, assemble_kp20_bi6_cap,
};

use super::super::kp20_polar_expanded_support_diagnostics::expanded_support_topology_detail;
use super::super::kp20_polar_expanded_support_topology::{
    Kp20ExpandedSupportTopologyReport, measure_expanded_support_topology,
};
use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::super::kp20_polar_fold_summary::gregory_patches_from_macro_patches;
use super::super::kp20_polar_multipatch_topology_split::{
    Kp20TopologySplitReport, refine_multipatch_topology_split,
};
use super::super::kp20_polar_support_image::refine_support_side_reparameterization;
use super::super::kp20_polar_support_image_ownership::{
    refine_multi_support_ownership_reparameterization, refine_support_ownership_reparameterization,
};
use super::super::kp20_polar_support_image_separated::refine_fold_site_separated_reparameterization;
use super::super::kp20_polar_support_lift::{
    Kp20LiftSearchProfile, Kp20SupportNormalLiftReport,
    lift_candidate_from_macro_patches_with_profile,
};
use super::super::kp20_polar_support_macro_construction::{
    Kp20SupportMacroConstructionReport, refine_support_macro_construction,
};
use super::super::kp20_polar_support_reparameterization::{
    Kp20SupportDomainReport, Kp20SupportReparameterizationReport,
};
use super::super::kp20_polar_support_strip_reauthor::{
    Kp20SupportStripReauthorReport, refine_support_strip_reauthoring,
};
use super::super::kp20_polar_support_transition_strip::{
    Kp20SupportTransitionStripReport, refine_support_transition_strip,
};
use super::super::kp20_polar_support_two_lane::{
    Kp20SupportTwoLaneReport, refine_support_two_lane,
};
use super::labels::{fold_ownership_label, piece_counts_label, validation_label};
use super::origin_inner_row::{
    Kp20OriginInnerRowReport, Kp20PairedInnerRowReport, Kp20SupportInnerRowReport,
    refine_origin_inner_interior_row, refine_paired_inner_interior_row,
    refine_support_inner_interior_row,
};
use super::support_domain::support_domain_initial_basis_candidate;

pub(super) struct Kp20LimitedLiftCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) report: Kp20SupportNormalLiftReport,
    strategy: &'static str,
    support_domain: Option<Kp20SupportDomainReport>,
    origin_inner: Option<Kp20OriginInnerRowReport>,
    paired_inner: Option<Kp20PairedInnerRowReport>,
    support_inner: Option<Kp20SupportInnerRowReport>,
    support_reparameterization: Option<Kp20SupportReparameterizationReport>,
    support_ownership_reparameterization: Option<Kp20SupportReparameterizationReport>,
    multi_support_reparameterization: Option<Kp20SupportReparameterizationReport>,
    separated_support_reparameterization: Option<Kp20SupportReparameterizationReport>,
    support_transition_strip: Option<Kp20SupportTransitionStripReport>,
    support_macro_construction: Option<Kp20SupportMacroConstructionReport>,
    multipatch_topology_split: Option<Kp20TopologySplitReport>,
    support_two_lane: Option<Kp20SupportTwoLaneReport>,
    expanded_support_topology: Option<Kp20ExpandedSupportTopologyReport>,
    support_strip_reauthor: Option<Kp20SupportStripReauthorReport>,
    elapsed_ms: u128,
}

impl Kp20LimitedLiftCandidate {
    pub(super) fn accepted_expanded_support_patches(&self) -> Option<&[GregoryPatch]> {
        let report = self.expanded_support_topology.as_ref()?;
        (report.accepted > 0 && !report.emitted_patches.is_empty())
            .then_some(report.emitted_patches.as_slice())
    }

    pub(super) fn accepted_support_cell_patches(&self) -> Option<&[GregoryPatch]> {
        let report = self.expanded_support_topology.as_ref()?;
        let cells = &report.support_cell_primitive;
        (cells.folded_patch_count == 0 && !cells.emitted_patches.is_empty())
            .then_some(cells.emitted_patches.as_slice())
    }

    pub(super) fn support_cell_floor_patches(&self) -> Option<&[GregoryPatch]> {
        let report = self.expanded_support_topology.as_ref()?;
        let cells = &report.support_cell_primitive;
        (!cells.committed_floor_patches.is_empty())
            .then_some(cells.committed_floor_patches.as_slice())
    }

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
        let origin_inner = self
            .origin_inner
            .as_ref()
            .map(|report| {
                format!(
                    ";origin_inner_evaluated={};origin_inner_accepted={};origin_inner_rejected={};origin_inner_start_folds={};origin_inner_end_folds={};origin_inner_mode={};origin_inner_scale={:.3e}",
                    report.evaluated,
                    report.accepted,
                    report.rejected,
                    report.start_folds,
                    report.end_folds,
                    report.mode,
                    report.scale,
                )
            })
            .unwrap_or_default();
        let paired_inner = self
            .paired_inner
            .as_ref()
            .map(|report| inner_report_detail("paired_inner", report))
            .unwrap_or_default();
        let support_inner = self
            .support_inner
            .as_ref()
            .map(|report| inner_report_detail("support_inner", report))
            .unwrap_or_default();
        let support_reparameterization = self
            .support_reparameterization
            .as_ref()
            .map(|report| support_reparameterization_detail("support_reparam", report))
            .unwrap_or_default();
        let support_strip_reauthor = self
            .support_strip_reauthor
            .as_ref()
            .map(|report| support_strip_reauthor_detail("support_strip_reauthor", report))
            .unwrap_or_default();
        let support_ownership_reparameterization = self
            .support_ownership_reparameterization
            .as_ref()
            .map(|report| support_reparameterization_detail("support_owner_reparam", report))
            .unwrap_or_default();
        let multi_support_reparameterization = self
            .multi_support_reparameterization
            .as_ref()
            .map(|report| support_reparameterization_detail("support_multi_reparam", report))
            .unwrap_or_default();
        let separated_support_reparameterization = self
            .separated_support_reparameterization
            .as_ref()
            .map(|report| support_reparameterization_detail("support_separated_reparam", report))
            .unwrap_or_default();
        let support_transition_strip = self
            .support_transition_strip
            .as_ref()
            .map(|report| support_transition_strip_detail("support_transition_strip", report))
            .unwrap_or_default();
        let multipatch_topology_split = self
            .multipatch_topology_split
            .as_ref()
            .map(|report| multipatch_topology_split_detail("multipatch_topology_split", report))
            .unwrap_or_default();
        let support_macro_construction = self
            .support_macro_construction
            .as_ref()
            .map(|report| support_macro_construction_detail("support_macro_construction", report))
            .unwrap_or_default();
        let support_two_lane = self
            .support_two_lane
            .as_ref()
            .map(|report| support_two_lane_detail("support_two_lane", report))
            .unwrap_or_default();
        let expanded_support_topology = self
            .expanded_support_topology
            .as_ref()
            .map(|report| expanded_support_topology_detail("expanded_support_topology", report))
            .unwrap_or_default();
        format!(
            "status=sampled;strategy={};default_folds={};folds={};piece_folds={};{};validation:{};optimizer_evaluated={};optimizer_accepted={};accepted_steps={};elapsed_ms={}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
            self.strategy,
            default.fold_count,
            lifted.fold_count,
            piece_counts_label(lifted.piece_fold_counts),
            fold_ownership_label(lifted),
            validation_label(self.report.validation),
            self.report.optimizer_stats.evaluated_candidates,
            self.report.optimizer_stats.accepted_candidates,
            self.report.accepted_steps.len(),
            self.elapsed_ms,
            support_domain,
            origin_inner,
            paired_inner,
            support_inner,
            support_reparameterization,
            support_ownership_reparameterization,
            multi_support_reparameterization,
            separated_support_reparameterization,
            support_transition_strip,
            support_macro_construction,
            multipatch_topology_split,
            support_two_lane,
            expanded_support_topology,
            support_strip_reauthor,
        )
    }
}
fn inner_report_detail(label: &str, report: &Kp20OriginInnerRowReport) -> String {
    format!(
        ";{label}_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_mode={};{label}_scale={:.3e};{label}_pre_evaluated={};{label}_pre_accepted={};{label}_pre_rejected={};{label}_pre_mode={}",
        report.evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.mode,
        report.scale,
        report.pre_evaluated,
        report.pre_accepted,
        report.pre_rejected,
        report.pre_mode,
    )
}

fn support_reparameterization_detail(
    label: &str,
    report: &Kp20SupportReparameterizationReport,
) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_lift_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_target_piece={:?};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.lift_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.target_piece,
        report.decline_reason,
    )
}

fn support_strip_reauthor_detail(label: &str, report: &Kp20SupportStripReauthorReport) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_target_piece={:?};{label}_target_sector={};{label}_mode={};{label}_coupling={};{label}_coupling_family={};{label}_scale={:.3e};{label}_worst_before={};{label}_worst_after={};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.exact_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.target_piece,
        report
            .target_sector
            .map(|sector| sector.to_string())
            .unwrap_or_else(|| "none".to_string()),
        report.mode,
        report.coupling,
        report.coupling_family,
        report.scale,
        report.worst_before,
        report.worst_after,
        report.decline_reason,
    )
}

fn support_transition_strip_detail(
    label: &str,
    report: &Kp20SupportTransitionStripReport,
) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_target_piece={:?};{label}_target_sector={};{label}_band_width={};{label}_mode={};{label}_scale={:.3e};{label}_fold_u={:.3};{label}_fold_v={:.3};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.exact_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.target_piece,
        report
            .target_sector
            .map(|sector| sector.to_string())
            .unwrap_or_else(|| "none".to_string()),
        report.band_width,
        report.mode,
        report.scale,
        report.fold_u,
        report.fold_v,
        report.decline_reason,
    )
}

fn multipatch_topology_split_detail(label: &str, report: &Kp20TopologySplitReport) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_strip_count={};{label}_inner_count={};{label}_direct_patch_count={};{label}_direct_patch_folds={};{label}_changed_piece_count={};{label}_strip_mode={};{label}_inner_mode={};{label}_scale={:.3e};{label}_validation_decline_reason={};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.exact_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.strip_count,
        report.inner_count,
        report.direct_patch_count,
        report.direct_patch_folds,
        report.changed_piece_count,
        report.strip_mode,
        report.inner_mode,
        report.scale,
        report.validation_decline_reason,
        report.decline_reason,
    )
}

fn support_macro_construction_detail(
    label: &str,
    report: &Kp20SupportMacroConstructionReport,
) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_target_piece={:?};{label}_target_sector={};{label}_mode={};{label}_scale={:.3e};{label}_fold_u={:.3};{label}_fold_v={:.3};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.exact_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.target_piece,
        report
            .target_sector
            .map(|sector| sector.to_string())
            .unwrap_or_else(|| "none".to_string()),
        report.mode,
        report.scale,
        report.fold_u,
        report.fold_v,
        report.decline_reason,
    )
}

fn support_two_lane_detail(label: &str, report: &Kp20SupportTwoLaneReport) -> String {
    format!(
        ";{label}_evaluated={};{label}_fast_accepted={};{label}_fast_rejected={};{label}_exact_evaluated={};{label}_accepted={};{label}_rejected={};{label}_start_folds={};{label}_end_folds={};{label}_target_piece={:?};{label}_target_sector={};{label}_lane_count={};{label}_split_parameter={:.3};{label}_max_interface_residual={:.3e};{label}_mode={};{label}_scale={:.3e};{label}_decline_reason={}",
        report.evaluated,
        report.fast_accepted,
        report.fast_rejected,
        report.exact_evaluated,
        report.accepted,
        report.rejected,
        report.start_folds,
        report.end_folds,
        report.target_piece,
        report
            .target_sector
            .map(|sector| sector.to_string())
            .unwrap_or_else(|| "none".to_string()),
        report.lane_count,
        report.split_parameter,
        report.max_interface_residual,
        report.mode,
        report.scale,
        report.decline_reason,
    )
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
    let (macro_patches, support_domain) =
        support_domain_initial_basis_candidate(star_point, tensor_border)?;
    let report = lift_candidate_from_macro_patches_with_profile(
        star_point,
        tensor_border,
        macro_patches,
        Kp20LiftSearchProfile::General,
        8,
    )?;
    let (report, origin_inner) =
        refine_origin_inner_interior_row(star_point, tensor_border, report)?;
    let (report, paired_inner) =
        refine_paired_inner_interior_row(star_point, tensor_border, report)?;
    let (report, support_inner) =
        refine_support_inner_interior_row(star_point, tensor_border, report)?;
    let (report, support_reparameterization) =
        refine_support_side_reparameterization(star_point, tensor_border, report)?;
    let (report, support_ownership_reparameterization) =
        refine_support_ownership_reparameterization(star_point, tensor_border, report)?;
    let (report, multi_support_reparameterization) =
        refine_multi_support_ownership_reparameterization(star_point, tensor_border, report)?;
    let (report, separated_support_reparameterization) =
        refine_fold_site_separated_reparameterization(star_point, tensor_border, report)?;
    let (report, support_transition_strip) =
        refine_support_transition_strip(tensor_border, report)?;
    let (report, support_macro_construction) =
        refine_support_macro_construction(tensor_border, report)?;
    let (report, multipatch_topology_split) =
        refine_multipatch_topology_split(tensor_border, report)?;
    let (report, support_two_lane) = refine_support_two_lane(tensor_border, report)?;
    let (report, expanded_support_topology) =
        measure_expanded_support_topology(tensor_border, report)?;
    let (report, support_strip_reauthor) = refine_support_strip_reauthoring(tensor_border, report)?;
    let patches = gregory_patches_from_macro_patches(&report.macro_patches);
    Ok(Kp20LimitedLiftCandidate {
        patches,
        report,
        strategy: "support_domain_basis_lift",
        support_domain: Some(support_domain),
        origin_inner: Some(origin_inner),
        paired_inner: Some(paired_inner),
        support_inner: Some(support_inner),
        support_reparameterization: Some(support_reparameterization),
        support_ownership_reparameterization: Some(support_ownership_reparameterization),
        multi_support_reparameterization: Some(multi_support_reparameterization),
        separated_support_reparameterization: Some(separated_support_reparameterization),
        support_transition_strip: Some(support_transition_strip),
        support_macro_construction: Some(support_macro_construction),
        multipatch_topology_split: Some(multipatch_topology_split),
        support_two_lane: Some(support_two_lane),
        expanded_support_topology: Some(expanded_support_topology),
        support_strip_reauthor: Some(support_strip_reauthor),
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
        origin_inner: None,
        paired_inner: None,
        support_inner: None,
        support_reparameterization: None,
        support_ownership_reparameterization: None,
        multi_support_reparameterization: None,
        separated_support_reparameterization: None,
        support_transition_strip: None,
        support_macro_construction: None,
        multipatch_topology_split: None,
        support_two_lane: None,
        expanded_support_topology: None,
        support_strip_reauthor: None,
        elapsed_ms: start.elapsed().as_millis(),
    })
}
