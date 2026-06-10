//! Shared-trim quality facts for KP20 support-cell patchwork.

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_shared_trim_measure.rs"]
mod shared_trim_measure;

pub(super) use shared_trim_measure::measure_support_cell_shared_trim_quality;

const CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET: f64 = 5.0e-4;
const COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET: f64 = 1.0e-10;
/// Commercial G1 target on `1 − |n_a · n_b|` across a shared
/// trim. Calibrated to the ring-gate angular threshold
/// `G1_GATE = 5e-2 (≈ 2.9°)` via `1 − cos θ ≈ θ²/2`, so both
/// gates accept the same worst-case normal mismatch.
const COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G1_TARGET: f64 = 1.25e-3;
/// Commercial G2 target on the scale-normalized curvature
/// residual `|κ_a − κ_b| / (1 + max(|κ_a|, |κ_b|))` — the same
/// form and value as the ring-gate `G2_GATE`.
const COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G2_TARGET: f64 = 5.0e-2;
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
    pub(super) worst_g2_residual: f64,
    pub(super) g2_unmeasured_sample_count: usize,
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

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SupportCellTrimOwnershipPlan {
    pub(super) shared_trim_count: usize,
    pub(super) full_row_owned_trim_count: usize,
    pub(super) segment_owned_trim_count: usize,
    pub(super) max_boundary_row_shared_trim_count: usize,
}

impl Kp20SupportCellSharedTrimQuality {
    /// Candidate-ladder non-regression: pcurve domain validity
    /// and G0 are monotonic. G1/G2 are intentionally *not*
    /// gated here — intermediate repairs may trade normal or
    /// curvature agreement on the way to a better G0 (the
    /// tangent-row pass recovers G1 afterwards, and the final
    /// quality classification reports all three continuity
    /// orders). Continuity-order trade-offs are rejected at the
    /// final verification gate, not per ladder step, because
    /// G1/G2 measured against not-yet-repaired pcurve domains
    /// are not comparable across a domain-changing repair.
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
    piecewise_boundary_committed: &'static str,
    quality: Kp20SupportCellSharedTrimQuality,
) -> ApexConstructionDiagnostic {
    let classification = classify_support_cell_shared_trim_quality(quality);
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_shared_trim_quality",
        n_sides,
        Some(format!(
            "selected={selected};measurement_scope={measurement_scope};piecewise_boundary_committed={piecewise_boundary_committed};current_g0_budget={:.15e};commercial_g0_target={:.15e};commercial_g1_target={:.15e};commercial_g2_target={:.15e};quality_status={};shared_trim_count={};partial_shared_trim_count={};pcurve_out_of_unit_count={};degenerate_pcurve_count={};max_boundary_row_shared_trim_count={};worst_abs_uv={:.15e};worst_g0_gap={:.15e};worst_normal_disagreement={:.15e};worst_g2_residual={:.15e};g2_unmeasured_sample_count={};{}",
            CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET,
            COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET,
            COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G1_TARGET,
            COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G2_TARGET,
            classification.status,
            quality.shared_trim_count,
            quality.partial_shared_trim_count,
            quality.pcurve_out_of_unit_count,
            quality.degenerate_pcurve_count,
            quality.max_boundary_row_shared_trim_count,
            quality.worst_abs_uv,
            quality.worst_g0_gap,
            quality.worst_normal_disagreement,
            quality.worst_g2_residual,
            quality.g2_unmeasured_sample_count,
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
            "selected={selected};selected_candidate={selected_candidate};repair_status={status};before_shared_trim_count={};after_shared_trim_count={};before_partial_shared_trim_count={};after_partial_shared_trim_count={};before_pcurve_out_of_unit_count={};after_pcurve_out_of_unit_count={};before_degenerate_pcurve_count={};after_degenerate_pcurve_count={};before_max_boundary_row_shared_trim_count={};after_max_boundary_row_shared_trim_count={};before_worst_abs_uv={:.15e};after_worst_abs_uv={:.15e};before_worst_g0_gap={:.15e};after_worst_g0_gap={:.15e};before_worst_normal_disagreement={:.15e};after_worst_normal_disagreement={:.15e};before_worst_g2_residual={:.15e};after_worst_g2_residual={:.15e};before_g2_unmeasured_sample_count={};after_g2_unmeasured_sample_count={}",
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
            after.worst_normal_disagreement,
            before.worst_g2_residual,
            after.worst_g2_residual,
            before.g2_unmeasured_sample_count,
            after.g2_unmeasured_sample_count
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
    operation_status: &str,
    before: Kp20SupportCellSharedTrimQuality,
    after: Kp20SupportCellSharedTrimQuality,
) -> ApexConstructionDiagnostic {
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_shared_trim_candidate_score",
        n_sides,
        Some(format!(
            "selected={selected};candidate={candidate};candidate_status={candidate_status};operation_status={operation_status};before_pcurve_out_of_unit_count={};after_pcurve_out_of_unit_count={};before_degenerate_pcurve_count={};after_degenerate_pcurve_count={};before_worst_g0_gap={:.15e};after_worst_g0_gap={:.15e};before_worst_abs_uv={:.15e};after_worst_abs_uv={:.15e};before_worst_normal_disagreement={:.15e};after_worst_normal_disagreement={:.15e};before_worst_g2_residual={:.15e};after_worst_g2_residual={:.15e};before_g2_unmeasured_sample_count={};after_g2_unmeasured_sample_count={}",
            before.pcurve_out_of_unit_count,
            after.pcurve_out_of_unit_count,
            before.degenerate_pcurve_count,
            after.degenerate_pcurve_count,
            before.worst_g0_gap,
            after.worst_g0_gap,
            before.worst_abs_uv,
            after.worst_abs_uv,
            before.worst_normal_disagreement,
            after.worst_normal_disagreement,
            before.worst_g2_residual,
            after.worst_g2_residual,
            before.g2_unmeasured_sample_count,
            after.g2_unmeasured_sample_count
        )),
        None,
        None,
    )
}

pub(super) fn support_cell_trim_ownership_plan_fact(
    n_sides: usize,
    selected: &'static str,
    measurement_scope: &'static str,
    piecewise_boundary_committed: &'static str,
    plan: Kp20SupportCellTrimOwnershipPlan,
) -> ApexConstructionDiagnostic {
    let ownership_status = if plan.segment_owned_trim_count > 0 {
        "piecewise_segment_ownership_required"
    } else {
        "full_row_ownership_sufficient"
    };
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_trim_ownership_plan",
        n_sides,
        Some(format!(
            "selected={selected};measurement_scope={measurement_scope};ownership_status={ownership_status};shared_trim_count={};full_row_owned_trim_count={};segment_owned_trim_count={};max_boundary_row_shared_trim_count={};piecewise_boundary_committed={piecewise_boundary_committed}",
            plan.shared_trim_count,
            plan.full_row_owned_trim_count,
            plan.segment_owned_trim_count,
            plan.max_boundary_row_shared_trim_count
        )),
        None,
        None,
    )
}

fn trim_ownership_plan(
    quality: Kp20SupportCellSharedTrimQuality,
) -> Kp20SupportCellTrimOwnershipPlan {
    Kp20SupportCellTrimOwnershipPlan {
        shared_trim_count: quality.shared_trim_count,
        full_row_owned_trim_count: quality
            .shared_trim_count
            .saturating_sub(quality.partial_shared_trim_count),
        segment_owned_trim_count: quality.partial_shared_trim_count,
        max_boundary_row_shared_trim_count: quality.max_boundary_row_shared_trim_count,
    }
}

struct Kp20SupportCellSharedTrimQualityClassification {
    status: &'static str,
}

/// Status lattice for shared-trim quality, ordered by the first
/// gate that fails: pcurve domain validity, then the interim G0
/// budget, then the commercial G0/G1/G2 targets in continuity
/// order. `commercial_g2_unmeasured` ranks below
/// `commercial_pass` — unmeasured curvature samples fail closed
/// rather than passing silently.
fn classify_support_cell_shared_trim_quality(
    quality: Kp20SupportCellSharedTrimQuality,
) -> Kp20SupportCellSharedTrimQualityClassification {
    let domain_valid =
        quality.pcurve_out_of_unit_count == 0 && quality.degenerate_pcurve_count == 0;
    let current_g0_valid = quality.worst_g0_gap <= CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET;
    let commercial_g0_valid = quality.worst_g0_gap <= COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET;
    let commercial_g1_valid =
        quality.worst_normal_disagreement <= COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G1_TARGET;
    let commercial_g2_valid =
        quality.worst_g2_residual <= COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G2_TARGET;
    let g2_measured = quality.g2_unmeasured_sample_count == 0;
    let status = if !domain_valid {
        "pcurve_domain_invalid"
    } else if !current_g0_valid {
        "current_g0_budget_exceeded"
    } else if !commercial_g0_valid {
        "current_budget_pass_commercial_g0_exceeded"
    } else if !commercial_g1_valid {
        "commercial_g1_exceeded"
    } else if !g2_measured {
        "commercial_g2_unmeasured"
    } else if !commercial_g2_valid {
        "commercial_g2_exceeded"
    } else {
        "commercial_pass"
    };
    Kp20SupportCellSharedTrimQualityClassification { status }
}

pub(crate) fn append_final_support_cell_shared_trim_quality(
    brep: &BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let selected = support_cell_selected_label(diagnostics);
    let piecewise_boundary_committed = support_cell_piecewise_boundary_committed_label(diagnostics);
    for group in support_cell_face_groups(corner_patches) {
        let quality = measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        diagnostics.push(support_cell_shared_trim_quality_fact(
            group.n_sides,
            selected,
            "final_post_pipeline",
            piecewise_boundary_committed,
            quality,
        ));
        diagnostics.push(support_cell_trim_ownership_plan_fact(
            group.n_sides,
            selected,
            "final_post_pipeline",
            piecewise_boundary_committed,
            trim_ownership_plan(quality),
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

fn support_cell_piecewise_boundary_committed_label(
    diagnostics: &[ApexConstructionDiagnostic],
) -> &'static str {
    let mut saw_split_decline = false;
    for diagnostic in diagnostics.iter().rev() {
        let ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
            continue;
        };
        if cap.stage != "kp20_polar_macro.dispatch.support_cell_shared_trim_candidate_score" {
            continue;
        }
        let Some(detail) = cap.detail_code.as_deref() else {
            continue;
        };
        if detail.contains("segment_topology_split_committed") {
            return "split_committed";
        }
        if detail.contains("segment_topology_already_segmented") {
            return "already_segmented";
        }
        saw_split_decline |= detail.contains("partial_boundary_segment_edge_split_rolled_back")
            || detail.contains("partial_boundary_segment_topology_address_invalid");
    }
    if saw_split_decline {
        "declined"
    } else {
        "false"
    }
}
