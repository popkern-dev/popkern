//! Shared-trim quality facts for KP20 support-cell patchwork.

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_shared_trim_measure.rs"]
mod shared_trim_measure;

pub(super) use shared_trim_measure::measure_support_cell_shared_trim_quality;

const CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET: f64 = 5.0e-4;
const COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET: f64 = 1.0e-10;
const SUPPORT_CELL_SHARED_TRIM_REPAIR_MONOTONIC_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SupportCellSharedTrimQuality {
    pub(super) shared_trim_count: usize,
    pub(super) partial_shared_trim_count: usize,
    pub(super) pcurve_out_of_unit_count: usize,
    pub(super) degenerate_pcurve_count: usize,
    pub(super) max_boundary_row_shared_trim_count: usize,
    pub(super) worst_abs_uv: f64,
    pub(super) worst_g0_gap: f64,
    pub(super) worst_normal_disagreement: f64,
    pub(super) worst_trim: Option<Kp20SupportCellWorstTrim>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SupportCellWorstTrim {
    pub(super) face_a: FaceHandle,
    pub(super) row_a: Option<u8>,
    pub(super) face_b: FaceHandle,
    pub(super) row_b: Option<u8>,
    pub(super) sample_fraction: f64,
    pub(super) full_span: bool,
}

impl Kp20SupportCellSharedTrimQuality {
    pub(super) fn is_no_worse_than(self, before: Self) -> bool {
        self.pcurve_out_of_unit_count <= before.pcurve_out_of_unit_count
            && self.degenerate_pcurve_count <= before.degenerate_pcurve_count
            && is_monotonic_metric(self.worst_g0_gap, before.worst_g0_gap)
    }
}

fn is_monotonic_metric(after: f64, before: f64) -> bool {
    after <= before + SUPPORT_CELL_SHARED_TRIM_REPAIR_MONOTONIC_EPS
}

pub(super) fn support_cell_shared_trim_quality_fact(
    n_sides: usize,
    selected: &'static str,
    measurement_scope: &'static str,
    quality: Kp20SupportCellSharedTrimQuality,
) -> ApexConstructionDiagnostic {
    let classification = classify_support_cell_shared_trim_quality(quality);
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_shared_trim_quality",
        n_sides,
        Some(format!(
            "selected={selected};measurement_scope={measurement_scope};current_g0_budget={:.15e};commercial_g0_target={:.15e};quality_status={};shared_trim_count={};partial_shared_trim_count={};pcurve_out_of_unit_count={};degenerate_pcurve_count={};max_boundary_row_shared_trim_count={};worst_abs_uv={:.15e};worst_g0_gap={:.15e};worst_normal_disagreement={:.15e};{}",
            CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET,
            COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET,
            classification.status,
            quality.shared_trim_count,
            quality.partial_shared_trim_count,
            quality.pcurve_out_of_unit_count,
            quality.degenerate_pcurve_count,
            quality.max_boundary_row_shared_trim_count,
            quality.worst_abs_uv,
            quality.worst_g0_gap,
            quality.worst_normal_disagreement,
            worst_trim_detail(quality)
        )),
        None,
        None,
    )
}

fn worst_trim_detail(quality: Kp20SupportCellSharedTrimQuality) -> String {
    let Some(worst) = quality.worst_trim else {
        return "worst_trim=none".to_string();
    };
    format!(
        "worst_trim=present;worst_trim_span={};worst_face_a={:?};worst_row_a={};worst_face_b={:?};worst_row_b={};worst_sample_fraction={:.15e}",
        if worst.full_span { "full" } else { "partial" },
        worst.face_a,
        row_label(worst.row_a),
        worst.face_b,
        row_label(worst.row_b),
        worst.sample_fraction
    )
}

fn row_label(row: Option<u8>) -> String {
    row.map(|row| row.to_string())
        .unwrap_or_else(|| "none".to_string())
}

pub(super) fn support_cell_shared_trim_repair_fact(
    n_sides: usize,
    selected: &'static str,
    selected_candidate: &'static str,
    accepted: bool,
    before: Kp20SupportCellSharedTrimQuality,
    after: Kp20SupportCellSharedTrimQuality,
) -> ApexConstructionDiagnostic {
    let status = if accepted {
        "accepted_non_regressing_final_repair"
    } else {
        "declined_to_preserve_final_quality"
    };
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_shared_trim_final_repair",
        n_sides,
        Some(format!(
            "selected={selected};selected_candidate={selected_candidate};repair_status={status};before_shared_trim_count={};after_shared_trim_count={};before_partial_shared_trim_count={};after_partial_shared_trim_count={};before_pcurve_out_of_unit_count={};after_pcurve_out_of_unit_count={};before_degenerate_pcurve_count={};after_degenerate_pcurve_count={};before_max_boundary_row_shared_trim_count={};after_max_boundary_row_shared_trim_count={};before_worst_abs_uv={:.15e};after_worst_abs_uv={:.15e};before_worst_g0_gap={:.15e};after_worst_g0_gap={:.15e};before_worst_normal_disagreement={:.15e};after_worst_normal_disagreement={:.15e}",
            before.shared_trim_count,
            after.shared_trim_count,
            before.partial_shared_trim_count,
            after.partial_shared_trim_count,
            before.pcurve_out_of_unit_count,
            after.pcurve_out_of_unit_count,
            before.degenerate_pcurve_count,
            after.degenerate_pcurve_count,
            before.max_boundary_row_shared_trim_count,
            after.max_boundary_row_shared_trim_count,
            before.worst_abs_uv,
            after.worst_abs_uv,
            before.worst_g0_gap,
            after.worst_g0_gap,
            before.worst_normal_disagreement,
            after.worst_normal_disagreement
        )),
        None,
        None,
    )
}

pub(super) fn support_cell_shared_trim_candidate_fact(
    n_sides: usize,
    selected: &'static str,
    candidate: &'static str,
    candidate_status: &'static str,
    operation_status: &'static str,
    before: Kp20SupportCellSharedTrimQuality,
    after: Kp20SupportCellSharedTrimQuality,
) -> ApexConstructionDiagnostic {
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_shared_trim_candidate_score",
        n_sides,
        Some(format!(
            "selected={selected};candidate={candidate};candidate_status={candidate_status};operation_status={operation_status};before_pcurve_out_of_unit_count={};after_pcurve_out_of_unit_count={};before_degenerate_pcurve_count={};after_degenerate_pcurve_count={};before_worst_g0_gap={:.15e};after_worst_g0_gap={:.15e};before_worst_abs_uv={:.15e};after_worst_abs_uv={:.15e};before_worst_normal_disagreement={:.15e};after_worst_normal_disagreement={:.15e}",
            before.pcurve_out_of_unit_count,
            after.pcurve_out_of_unit_count,
            before.degenerate_pcurve_count,
            after.degenerate_pcurve_count,
            before.worst_g0_gap,
            after.worst_g0_gap,
            before.worst_abs_uv,
            after.worst_abs_uv,
            before.worst_normal_disagreement,
            after.worst_normal_disagreement
        )),
        None,
        None,
    )
}

struct Kp20SupportCellSharedTrimQualityClassification {
    status: &'static str,
}

fn classify_support_cell_shared_trim_quality(
    quality: Kp20SupportCellSharedTrimQuality,
) -> Kp20SupportCellSharedTrimQualityClassification {
    let domain_valid =
        quality.pcurve_out_of_unit_count == 0 && quality.degenerate_pcurve_count == 0;
    let current_g0_valid = quality.worst_g0_gap <= CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET;
    let commercial_g0_valid = quality.worst_g0_gap <= COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET;
    let status = if !domain_valid {
        "pcurve_domain_invalid"
    } else if !current_g0_valid {
        "current_g0_budget_exceeded"
    } else if !commercial_g0_valid {
        "current_budget_pass_commercial_g0_exceeded"
    } else {
        "commercial_g0_pass"
    };
    Kp20SupportCellSharedTrimQualityClassification { status }
}

pub(crate) fn append_final_support_cell_shared_trim_quality(
    brep: &BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let selected = support_cell_selected_label(diagnostics);
    for group in support_cell_face_groups(corner_patches) {
        let quality = measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        diagnostics.push(support_cell_shared_trim_quality_fact(
            group.n_sides,
            selected,
            "final_post_pipeline",
            quality,
        ));
    }
}

pub(super) struct SupportCellFaceGroup {
    pub(super) n_sides: usize,
    pub(super) faces: Vec<FaceHandle>,
}

pub(super) fn support_cell_face_groups(
    corner_patches: &[CornerPatchInfo],
) -> Vec<SupportCellFaceGroup> {
    let mut groups: Vec<(VertexHandle, SupportCellFaceGroup)> = Vec::new();
    for patch in corner_patches {
        if !patch
            .sub_patch_provenance
            .is_some_and(|provenance| matches!(provenance.role, CornerSubPatchRole::SupportCell))
        {
            continue;
        }
        if let Some((_, group)) = groups
            .iter_mut()
            .find(|(vertex, group)| *vertex == patch.vertex && group.n_sides == patch.n_sides)
        {
            group.faces.push(patch.face);
        } else {
            groups.push((
                patch.vertex,
                SupportCellFaceGroup {
                    n_sides: patch.n_sides,
                    faces: vec![patch.face],
                },
            ));
        }
    }
    groups
        .into_iter()
        .map(|(_, group)| group)
        .filter(|group| !group.faces.is_empty())
        .collect()
}

pub(super) fn support_cell_selected_label(
    diagnostics: &[ApexConstructionDiagnostic],
) -> &'static str {
    diagnostics
        .iter()
        .rev()
        .find_map(|diagnostic| {
            let ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
                return None;
            };
            (cap.stage == "kp20_polar_macro.dispatch.support_cell_assembled_quality")
                .then(|| cap.detail_code.as_deref())
                .flatten()
                .and_then(selected_label_from_detail)
        })
        .unwrap_or("unknown")
}

fn selected_label_from_detail(detail: &str) -> Option<&'static str> {
    if detail.contains("selected=promoted") {
        Some("promoted")
    } else if detail.contains("selected=fallback_four_cell_floor") {
        Some("fallback_four_cell_floor")
    } else {
        None
    }
}
