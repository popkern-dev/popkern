//! Shared-trim quality facts for KP20 support-cell patchwork.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

const CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET: f64 = 2.0e-3;
const COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET: f64 = 1.0e-10;
const SUPPORT_CELL_SHARED_TRIM_REPAIR_MONOTONIC_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SupportCellSharedTrimQuality {
    pub(super) shared_trim_count: usize,
    pub(super) pcurve_out_of_unit_count: usize,
    pub(super) degenerate_pcurve_count: usize,
    pub(super) max_boundary_row_shared_trim_count: usize,
    pub(super) worst_abs_uv: f64,
    pub(super) worst_g0_gap: f64,
    pub(super) worst_normal_disagreement: f64,
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

pub(super) fn measure_support_cell_shared_trim_quality(
    brep: &BRep,
    faces: &[FaceHandle],
    n_samples: usize,
) -> Kp20SupportCellSharedTrimQuality {
    let face_set = faces.iter().copied().collect::<HashSet<_>>();
    let mut visited_edges = HashSet::new();
    let mut boundary_row_counts = HashMap::new();
    let mut quality = Kp20SupportCellSharedTrimQuality {
        shared_trim_count: 0,
        pcurve_out_of_unit_count: 0,
        degenerate_pcurve_count: 0,
        max_boundary_row_shared_trim_count: 0,
        worst_abs_uv: 0.0,
        worst_g0_gap: 0.0,
        worst_normal_disagreement: 0.0,
    };
    for &face in faces {
        let Some(start) = brep
            .faces
            .get(face)
            .and_then(|face_data| brep.loops.get(face_data.outer_loop))
            .map(|loop_| loop_.half_edge)
        else {
            continue;
        };
        let mut cursor = start;
        for _ in 0..brep.half_edges.len().max(1) {
            let half_edge = &brep.half_edges[cursor];
            if visited_edges.insert(half_edge.edge) {
                let twin = &brep.half_edges[half_edge.twin];
                if let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_| loop_.face)
                    && face_set.contains(&twin_face)
                {
                    measure_shared_trim_pair(
                        brep,
                        face,
                        cursor,
                        twin_face,
                        half_edge.twin,
                        n_samples,
                        &mut quality,
                        &mut boundary_row_counts,
                    );
                }
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
    }
    quality
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
            "selected={selected};measurement_scope={measurement_scope};current_g0_budget={:.15e};commercial_g0_target={:.15e};quality_status={};shared_trim_count={};pcurve_out_of_unit_count={};degenerate_pcurve_count={};max_boundary_row_shared_trim_count={};worst_abs_uv={:.15e};worst_g0_gap={:.15e};worst_normal_disagreement={:.15e}",
            CURRENT_SUPPORT_CELL_SHARED_TRIM_G0_BUDGET,
            COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET,
            classification.status,
            quality.shared_trim_count,
            quality.pcurve_out_of_unit_count,
            quality.degenerate_pcurve_count,
            quality.max_boundary_row_shared_trim_count,
            quality.worst_abs_uv,
            quality.worst_g0_gap,
            quality.worst_normal_disagreement
        )),
        None,
        None,
    )
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
            "selected={selected};selected_candidate={selected_candidate};repair_status={status};before_shared_trim_count={};after_shared_trim_count={};before_pcurve_out_of_unit_count={};after_pcurve_out_of_unit_count={};before_degenerate_pcurve_count={};after_degenerate_pcurve_count={};before_max_boundary_row_shared_trim_count={};after_max_boundary_row_shared_trim_count={};before_worst_abs_uv={:.15e};after_worst_abs_uv={:.15e};before_worst_g0_gap={:.15e};after_worst_g0_gap={:.15e};before_worst_normal_disagreement={:.15e};after_worst_normal_disagreement={:.15e}",
            before.shared_trim_count,
            after.shared_trim_count,
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

fn measure_shared_trim_pair(
    brep: &BRep,
    face_a: FaceHandle,
    half_edge_a: popkern_core::topology::HalfEdgeHandle,
    face_b: FaceHandle,
    half_edge_b: popkern_core::topology::HalfEdgeHandle,
    n_samples: usize,
    quality: &mut Kp20SupportCellSharedTrimQuality,
    boundary_row_counts: &mut HashMap<(FaceHandle, u8), usize>,
) {
    let Some(surface_a) = brep.faces[face_a].surface.as_ref() else {
        return;
    };
    let Some(surface_b) = brep.faces[face_b].surface.as_ref() else {
        return;
    };
    if !matches!(surface_a.geometry, SurfaceGeometry::Gregory(_))
        || !matches!(surface_b.geometry, SurfaceGeometry::Gregory(_))
    {
        return;
    }
    let Some(pcurve_a) = brep.half_edges[half_edge_a].pcurve.as_ref() else {
        return;
    };
    let Some(pcurve_b) = brep.half_edges[half_edge_b].pcurve.as_ref() else {
        return;
    };
    quality.shared_trim_count += 1;
    record_boundary_row_count(quality, boundary_row_counts, face_a, pcurve_a);
    record_boundary_row_count(quality, boundary_row_counts, face_b, pcurve_b);
    record_pcurve_span(quality, pcurve_a, pcurve_b);
    for index in 0..n_samples {
        let fraction = if n_samples <= 1 {
            0.5
        } else {
            index as f64 / (n_samples - 1) as f64
        };
        let uv_a = pcurve_at_fraction(pcurve_a, fraction);
        let uv_b = pcurve_at_fraction(pcurve_b, 1.0 - fraction);
        record_uv(quality, uv_a);
        record_uv(quality, uv_b);
        let point_a = surface_a.geometry.point_at(uv_a.x, uv_a.y);
        let point_b = surface_b.geometry.point_at(uv_b.x, uv_b.y);
        quality.worst_g0_gap = quality.worst_g0_gap.max(point_a.distance_to(&point_b));
        let normal_a = oriented_normal(surface_a, uv_a.x, uv_a.y);
        let normal_b = oriented_normal(surface_b, uv_b.x, uv_b.y);
        if let (Some(a), Some(b)) = (normal_a.normalized(), normal_b.normalized()) {
            quality.worst_normal_disagreement =
                quality.worst_normal_disagreement.max(1.0 - a.dot(&b).abs());
        }
    }
}

fn record_boundary_row_count(
    quality: &mut Kp20SupportCellSharedTrimQuality,
    boundary_row_counts: &mut HashMap<(FaceHandle, u8), usize>,
    face: FaceHandle,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
) {
    let Some(row) = boundary_row(pcurve) else {
        return;
    };
    let count = boundary_row_counts
        .entry((face, row))
        .and_modify(|count| *count += 1)
        .or_insert(1);
    quality.max_boundary_row_shared_trim_count =
        quality.max_boundary_row_shared_trim_count.max(*count);
}

fn boundary_row(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> Option<u8> {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    if close(start.x, 0.0) && close(end.x, 0.0) {
        Some(0)
    } else if close(start.x, 1.0) && close(end.x, 1.0) {
        Some(1)
    } else if close(start.y, 0.0) && close(end.y, 0.0) {
        Some(2)
    } else if close(start.y, 1.0) && close(end.y, 1.0) {
        Some(3)
    } else {
        None
    }
}

fn record_pcurve_span(
    quality: &mut Kp20SupportCellSharedTrimQuality,
    pcurve_a: &popkern_core::geometry::pcurve::Pcurve,
    pcurve_b: &popkern_core::geometry::pcurve::Pcurve,
) {
    let span_a = pcurve_span(pcurve_a);
    let span_b = pcurve_span(pcurve_b);
    quality.degenerate_pcurve_count += usize::from(span_a <= 1.0e-12);
    quality.degenerate_pcurve_count += usize::from(span_b <= 1.0e-12);
}

fn pcurve_span(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> f64 {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    ((end.x - start.x).powi(2) + (end.y - start.y).powi(2)).sqrt()
}

fn record_uv(quality: &mut Kp20SupportCellSharedTrimQuality, uv: Vec3) {
    quality.worst_abs_uv = quality.worst_abs_uv.max(uv.x.abs()).max(uv.y.abs());
    if uv.x < -1.0e-9 || uv.x > 1.0 + 1.0e-9 || uv.y < -1.0e-9 || uv.y > 1.0 + 1.0e-9 {
        quality.pcurve_out_of_unit_count += 1;
    }
}

fn pcurve_at_fraction(pcurve: &popkern_core::geometry::pcurve::Pcurve, fraction: f64) -> Vec3 {
    let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
    let point = pcurve.point_at(t);
    Vec3::new(point.x, point.y, point.z)
}

fn oriented_normal(
    surface: &popkern_core::geometry::surface::SurfaceOnFace,
    u: f64,
    v: f64,
) -> Vec3 {
    let normal = surface.geometry.normal_at(u, v);
    if surface.orientation {
        normal
    } else {
        normal * -1.0
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
