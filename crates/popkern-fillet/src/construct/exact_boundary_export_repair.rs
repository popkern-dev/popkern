//! Exact-boundary export-readiness repair without broad support-face demotion.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::{exact_boundary_cap_trim_restore, pcurve_measured_repair, surface_driven_edge_repair};

pub(super) fn repair_exact_boundary_export_blockers(
    result: &mut BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    let exact_boundary_committed = repair_report.apex_routes.iter().any(|route| {
        route.exact_boundary_report.candidate_status
            == crate::types::ExactBoundaryCandidateOutcome::Accepted
    });
    if !exact_boundary_committed {
        return;
    }
    result.repair_line_edge_curve_endpoint_coherence(1.0e-6);
    if let Ok(edge_curves_completed) =
        popkern_core::topology::complete_missing_edge_curves_transactionally(
            result,
            exact_boundary_edge_curve_completion_options(),
        )
    {
        *result = edge_curves_completed.brep;
    }
    result.fix_stale_pcurves();
    pcurve_measured_repair::repair_measured_plane_cylinder_pcurves(result, 1.0e-6);
    surface_driven_edge_repair::repair_measured_edges_from_single_surface_pcurves(
        result, 1.0e-6, 0.06,
    );
    exact_boundary_cap_trim_restore::restore_exact_boundary_cap_pcurves(result, repair_report);
    record_exact_boundary_cap_outer_declines(result, repair_report, 5.0e-3);
    record_exact_boundary_support_scope_declines(result, repair_report, 1.0e-6);
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct CapOuterPairResidual {
    residual: f64,
    sample_parameter: f64,
}

fn cap_support_pair_residual_with_parameter(
    result: &BRep,
    cap_half_edge: HalfEdgeHandle,
    cap_pcurve: &Pcurve,
    support_pcurve: &Pcurve,
    cap_surface: &SurfaceGeometry,
    support_surface: &SurfaceGeometry,
) -> Option<CapOuterPairResidual> {
    let edge = result.half_edges.get(cap_half_edge)?.edge;
    let edge_curve = result.edges.get(edge)?.curve.as_ref()?;
    let (cap_start, cap_end) = cap_pcurve.parameter_range();
    let (support_start, support_end) = support_pcurve.parameter_range();
    let cap_forward = pcurve_follows_edge_curve(cap_surface, cap_pcurve, edge_curve);
    let support_forward = pcurve_follows_edge_curve(support_surface, support_pcurve, edge_curve);
    let samples = 100usize;
    let mut worst = 0.0_f64;
    let mut worst_parameter = 0.0_f64;
    for index in 0..samples {
        let fraction = 0.05 + 0.90 * index as f64 / (samples - 1) as f64;
        let cap_fraction = if cap_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let support_fraction = if support_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let cap_t = cap_start + (cap_end - cap_start) * cap_fraction;
        let support_t = support_start + (support_end - support_start) * support_fraction;
        let cap_uv = cap_pcurve.point_at(cap_t);
        let support_uv = support_pcurve.point_at(support_t);
        let cap_point = cap_surface.point_at(cap_uv.x, cap_uv.y);
        let support_point = support_surface.point_at(support_uv.x, support_uv.y);
        let gap = cap_point.distance_to(&support_point);
        if !gap.is_finite() {
            return None;
        }
        if gap > worst {
            worst = gap;
            worst_parameter = fraction;
        }
    }
    Some(CapOuterPairResidual {
        residual: worst,
        sample_parameter: worst_parameter,
    })
}

fn pcurve_follows_edge_curve(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    edge_curve: &popkern_core::geometry::curve::CurveOnEdge,
) -> bool {
    let (pc_start, pc_end) = pcurve.parameter_range();
    let pcurve_start = pcurve_surface_point(surface, pcurve, pc_start);
    let pcurve_end = pcurve_surface_point(surface, pcurve, pc_end);
    let curve_start = edge_curve.geometry.point_at(edge_curve.t_start);
    let curve_end = edge_curve.geometry.point_at(edge_curve.t_end);
    curve_start.distance_to(&pcurve_start) + curve_end.distance_to(&pcurve_end)
        <= curve_start.distance_to(&pcurve_end) + curve_end.distance_to(&pcurve_start)
}

fn pcurve_surface_point(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<SurfaceGeometry> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    brep.faces
        .get(face)?
        .surface
        .as_ref()
        .map(|surface| surface.geometry.clone())
}

/// Resolve the current (cap-side, support-side) half-edge pair for a
/// sector outer edge. The half-edge handles recorded at cap
/// construction are volatile — the V-free shell consolidation can
/// exchange a seam's halves between loops or re-author a rebuilt
/// support ring — while the EDGE and the support FACE keep their
/// identity. When the recorded support face owns one side of the
/// edge, that side is the support side and its partner the cap side;
/// otherwise the recorded handles stand and the existing structured
/// diagnostics (topology_missing / ownership_mismatch) report drift.
fn resolve_sector_outer_sides(
    result: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    recorded_cap: HalfEdgeHandle,
    recorded_support: HalfEdgeHandle,
    support_face: FaceHandle,
) -> (HalfEdgeHandle, HalfEdgeHandle) {
    let derived = result.edges.get(edge).and_then(|edge_data| {
        let a = edge_data.half_edge;
        let b = result.half_edges.get(a)?.twin;
        result.half_edges.get(b)?;
        let face_of = |handle: HalfEdgeHandle| -> Option<FaceHandle> {
            result
                .loops
                .get(result.half_edges.get(handle)?.loop_)
                .map(|loop_data| loop_data.face)
        };
        match (face_of(a)?, face_of(b)?) {
            (_, face_b) if face_b == support_face => Some((a, b)),
            (face_a, _) if face_a == support_face => Some((b, a)),
            _ => None,
        }
    });
    derived.unwrap_or((recorded_cap, recorded_support))
}

fn record_exact_boundary_cap_outer_declines(
    result: &BRep,
    repair_report: &mut crate::types::BlendRepairReport,
    tolerance: f64,
) {
    let mut apex_diagnostics = Vec::new();
    for route in &mut repair_report.apex_routes {
        if route.exact_boundary_report.candidate_status
            != crate::types::ExactBoundaryCandidateOutcome::Accepted
        {
            continue;
        }
        for sector in route.exact_boundary_report.sectors.clone() {
            for (side_index, edge, outer_side) in [
                (0usize, Some(sector.fillet_end_edge), "FilletEnd"),
                (1usize, sector.face_edge, "FaceEdge"),
            ] {
                // Welded shared corners have no bridge side.
                let Some(edge) = edge else {
                    continue;
                };
                let (Some(recorded_cap), Some(recorded_support), Some(side_support_face)) = (
                    sector.cap_half_edges[side_index],
                    sector.support_half_edges[side_index],
                    if side_index == 0 {
                        Some(sector.fillet_support_face)
                    } else {
                        sector.bridge_support_face
                    },
                ) else {
                    continue;
                };
                let (cap_half_edge, support_half_edge) = resolve_sector_outer_sides(
                    result,
                    edge,
                    recorded_cap,
                    recorded_support,
                    side_support_face,
                );
                let (code, residual, sample_parameter) = exact_boundary_cap_outer_decline(
                    result,
                    cap_half_edge,
                    support_half_edge,
                    tolerance,
                );
                let Some(code) = code else {
                    continue;
                };
                let diagnostic = crate::types::ExactBoundaryCapDiagnostic {
                    stage: "exact_boundary_cap.cap_outer".to_string(),
                    code,
                    detail_code: None,
                    sector_name: sector.sector_name.clone(),
                    edge: Some(edge),
                    support_half_edge: Some(support_half_edge),
                    support_face: Some(side_support_face),
                    residual,
                    outer_side: Some(outer_side.to_string()),
                    sample_parameter,
                    tolerance: Some(tolerance),
                    shared_fact_ids: Vec::new(),
                };
                route.exact_boundary_report.diagnostics.push(
                    crate::types::ExactBoundaryReportDiagnostic {
                        schema_version: 1,
                        stage: diagnostic.stage.clone(),
                        code: diagnostic.code.clone(),
                        detail_code: diagnostic.detail_code.clone(),
                        sector_name: diagnostic.sector_name.clone(),
                        edge: diagnostic.edge,
                        support_half_edge: diagnostic.support_half_edge,
                        support_face: diagnostic.support_face,
                        residual: diagnostic.residual,
                        outer_side: diagnostic.outer_side.clone(),
                        sample_parameter: diagnostic.sample_parameter,
                        tolerance: diagnostic.tolerance,
                        shared_fact_ids: Vec::new(),
                    },
                );
                apex_diagnostics.push(crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                    diagnostic,
                ));
            }
        }
    }
    repair_report.apex_diagnostics.extend(apex_diagnostics);
}

fn exact_boundary_cap_outer_decline(
    result: &BRep,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    tolerance: f64,
) -> (Option<String>, Option<f64>, Option<f64>) {
    let Some(cap_pcurve) = result
        .half_edges
        .get(cap_half_edge)
        .and_then(|half_edge| half_edge.pcurve.as_ref())
    else {
        return (
            Some("exact_boundary_cap.cap_outer.pcurve_missing".to_string()),
            None,
            None,
        );
    };
    let Some(support_pcurve) = result
        .half_edges
        .get(support_half_edge)
        .and_then(|half_edge| half_edge.pcurve.as_ref())
    else {
        return (
            Some("exact_boundary_cap.cap_outer.pcurve_missing".to_string()),
            None,
            None,
        );
    };
    let Some(cap_surface) = surface_for_half_edge(result, cap_half_edge) else {
        return (
            Some("exact_boundary_cap.cap_outer.surface_missing".to_string()),
            None,
            None,
        );
    };
    let Some(support_surface) = surface_for_half_edge(result, support_half_edge) else {
        return (
            Some("exact_boundary_cap.cap_outer.surface_missing".to_string()),
            None,
            None,
        );
    };
    let Some(report) = cap_support_pair_residual_with_parameter(
        result,
        cap_half_edge,
        cap_pcurve,
        support_pcurve,
        &cap_surface,
        &support_surface,
    ) else {
        return (
            Some("exact_boundary_cap.cap_outer.measurement_unavailable".to_string()),
            None,
            None,
        );
    };
    if report.residual > tolerance {
        return (
            Some("exact_boundary_cap.cap_outer.pair_residual_exceeded".to_string()),
            Some(report.residual),
            Some(report.sample_parameter),
        );
    }
    (None, None, None)
}

fn record_exact_boundary_support_scope_declines(
    result: &BRep,
    repair_report: &mut crate::types::BlendRepairReport,
    tolerance: f64,
) {
    let mut apex_diagnostics = Vec::new();
    for route in &mut repair_report.apex_routes {
        if route.exact_boundary_report.candidate_status
            != crate::types::ExactBoundaryCandidateOutcome::Accepted
        {
            continue;
        }
        for sector in route.exact_boundary_report.sectors.clone() {
            for (edge, support_half_edge, support_face) in [
                (
                    Some(sector.fillet_end_edge),
                    sector.support_half_edges[0],
                    Some(sector.fillet_support_face),
                ),
                (
                    sector.face_edge,
                    sector.support_half_edges[1],
                    sector.bridge_support_face,
                ),
            ] {
                // Welded shared corners have no bridge side.
                let (Some(edge), Some(support_half_edge), Some(support_face)) =
                    (edge, support_half_edge, support_face)
                else {
                    continue;
                };
                let (_, support_half_edge) = resolve_sector_outer_sides(
                    result,
                    edge,
                    support_half_edge,
                    support_half_edge,
                    support_face,
                );
                let (code, residual) = exact_boundary_support_scope_decline(
                    result,
                    support_half_edge,
                    support_face,
                    tolerance,
                );
                let Some(code) = code else {
                    continue;
                };
                let diagnostic = crate::types::ExactBoundaryCapDiagnostic {
                    stage: "exact_boundary_cap.support_scope".to_string(),
                    code,
                    detail_code: None,
                    sector_name: sector.sector_name.clone(),
                    edge: Some(edge),
                    support_half_edge: Some(support_half_edge),
                    support_face: Some(support_face),
                    residual,
                    outer_side: None,
                    sample_parameter: None,
                    tolerance: Some(tolerance),
                    shared_fact_ids: Vec::new(),
                };
                route.exact_boundary_report.diagnostics.push(
                    crate::types::ExactBoundaryReportDiagnostic {
                        schema_version: 1,
                        stage: diagnostic.stage.clone(),
                        code: diagnostic.code.clone(),
                        detail_code: diagnostic.detail_code.clone(),
                        sector_name: diagnostic.sector_name.clone(),
                        edge: diagnostic.edge,
                        support_half_edge: diagnostic.support_half_edge,
                        support_face: diagnostic.support_face,
                        residual: diagnostic.residual,
                        outer_side: diagnostic.outer_side.clone(),
                        sample_parameter: diagnostic.sample_parameter,
                        tolerance: Some(tolerance),
                        shared_fact_ids: Vec::new(),
                    },
                );
                apex_diagnostics.push(crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                    diagnostic,
                ));
            }
        }
    }
    repair_report.apex_diagnostics.extend(apex_diagnostics);
}

fn exact_boundary_support_scope_decline(
    result: &BRep,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    tolerance: f64,
) -> (Option<String>, Option<f64>) {
    let Some(half_edge) = result.half_edges.get(support_half_edge) else {
        return (
            Some("exact_boundary_cap.support_scope.topology_missing".to_string()),
            None,
        );
    };
    let Some(loop_) = result.loops.get(half_edge.loop_) else {
        return (
            Some("exact_boundary_cap.support_scope.topology_missing".to_string()),
            None,
        );
    };
    if loop_.face != support_face {
        return (
            Some("exact_boundary_cap.support_scope.ownership_mismatch".to_string()),
            None,
        );
    }
    if half_edge.pcurve.is_none() {
        return (
            Some("exact_boundary_cap.support_scope.pcurve_missing".to_string()),
            None,
        );
    }
    let Some(report) = result.measure_half_edge_pcurve_deviation(support_half_edge) else {
        return (
            Some("exact_boundary_cap.support_scope.measurement_unavailable".to_string()),
            None,
        );
    };
    if report.max_deviation > tolerance {
        return (
            Some("exact_boundary_cap.support_scope.residual_exceeded".to_string()),
            Some(report.max_deviation),
        );
    }
    (None, None)
}

fn exact_boundary_edge_curve_completion_options()
-> popkern_core::topology::BrepEdgeCurveCompletionOptions {
    popkern_core::topology::BrepEdgeCurveCompletionOptions {
        tolerance: 1e-6,
        samples: 9,
        allow_vertex_segments_without_surface: true,
        require_all_completed: false,
        validate_output: false,
    }
}
