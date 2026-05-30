//! Exact-boundary cap trim restoration after generic pcurve repair passes.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle};

pub(super) fn restore_exact_boundary_cap_pcurves(
    result: &mut BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    let sectors = repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.exact_boundary_report.candidate_status
                == crate::types::ExactBoundaryCandidateOutcome::Accepted
        })
        .flat_map(|route| route.exact_boundary_report.sectors.iter().cloned())
        .collect::<Vec<_>>();
    let mut resolved_cap_outer = Vec::new();
    for sector in sectors {
        for (side_index, outer_side) in [(0usize, "FilletEnd"), (1usize, "FaceEdge")] {
            let cap_half_edge = sector.cap_half_edges[side_index];
            let support_half_edge = sector.support_half_edges[side_index];
            let Some(surface) = surface_for_half_edge(result, cap_half_edge) else {
                continue;
            };
            let Some(candidate) =
                exact_boundary_cap_side_pcurve(result, cap_half_edge, support_half_edge, &surface)
            else {
                continue;
            };
            super::exact_boundary_cap_trim_surface::relock_exact_boundary_cap_side(
                result,
                cap_half_edge,
                support_half_edge,
                candidate.side,
                candidate.side_forward,
            );
            if let Some(half_edge) = result.half_edges.get_mut(cap_half_edge) {
                half_edge.pcurve = Some(candidate.pcurve);
            }
            if cap_support_pair_residual(result, cap_half_edge, support_half_edge)
                .is_some_and(|residual| residual <= 5.0e-3)
            {
                resolved_cap_outer.push((sector.sector_name.clone(), outer_side.to_string()));
            }
        }
    }
    prune_resolved_cap_outer_diagnostics(repair_report, &resolved_cap_outer);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExactBoundaryCapSurfaceSide {
    VMin,
    VMax,
    UMin,
    UMax,
}

struct ExactBoundaryCapPcurveCandidate {
    pcurve: Pcurve,
    side: ExactBoundaryCapSurfaceSide,
    side_forward: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct ExactBoundaryCapPcurveScore {
    g0_gap: f64,
    g1_deficit: f64,
}

fn exact_boundary_cap_side_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<ExactBoundaryCapPcurveCandidate> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
        return None;
    }
    let candidates = [
        (
            ExactBoundaryCapSurfaceSide::VMin,
            Pcurve::linear((u_min, v_min), (u_max, v_min), 0.0, 1.0),
        ),
        (
            ExactBoundaryCapSurfaceSide::VMax,
            Pcurve::linear((u_min, v_max), (u_max, v_max), 0.0, 1.0),
        ),
        (
            ExactBoundaryCapSurfaceSide::UMin,
            Pcurve::linear((u_min, v_min), (u_min, v_max), 0.0, 1.0),
        ),
        (
            ExactBoundaryCapSurfaceSide::UMax,
            Pcurve::linear((u_max, v_min), (u_max, v_max), 0.0, 1.0),
        ),
    ];
    candidates
        .into_iter()
        .flat_map(|(side, candidate)| {
            let reversed = reverse_linear_pcurve(&candidate);
            [(side, true, candidate), (side, false, reversed)]
        })
        .filter_map(|(side, side_forward, candidate)| {
            let score = exact_boundary_cap_pcurve_score(
                brep,
                half_edge,
                support_half_edge,
                surface,
                &candidate,
            )?;
            Some((
                score,
                ExactBoundaryCapPcurveCandidate {
                    pcurve: candidate,
                    side,
                    side_forward,
                },
            ))
        })
        .min_by(|(a, _), (b, _)| compare_cap_pcurve_scores(*a, *b))
        .map(|(_, candidate)| candidate)
}

fn exact_boundary_cap_pcurve_score(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> Option<ExactBoundaryCapPcurveScore> {
    let endpoint_gap = pcurve_endpoint_gap(brep, half_edge, surface, pcurve)?;
    let edge = brep.half_edges.get(half_edge)?.edge;
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let trace_gap = pcurve_curve_trace_gap(brep, half_edge, surface, pcurve, curve)?;
    Some(ExactBoundaryCapPcurveScore {
        g0_gap: endpoint_gap.max(trace_gap),
        g1_deficit: support_normal_deficit_for_candidate(brep, support_half_edge, surface, pcurve)
            .unwrap_or(f64::INFINITY),
    })
}

fn compare_cap_pcurve_scores(
    a: ExactBoundaryCapPcurveScore,
    b: ExactBoundaryCapPcurveScore,
) -> std::cmp::Ordering {
    if (a.g0_gap - b.g0_gap).abs() <= 1.0e-8 {
        a.g1_deficit
            .total_cmp(&b.g1_deficit)
            .then_with(|| a.g0_gap.total_cmp(&b.g0_gap))
    } else {
        a.g0_gap
            .total_cmp(&b.g0_gap)
            .then_with(|| a.g1_deficit.total_cmp(&b.g1_deficit))
    }
}

fn support_normal_deficit_for_candidate(
    brep: &BRep,
    support_half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    cap_pcurve: &Pcurve,
) -> Option<f64> {
    let support_pcurve = brep.half_edges.get(support_half_edge)?.pcurve.as_ref()?;
    let support_surface = surface_for_half_edge(brep, support_half_edge)?;
    let (cap_start_t, cap_end_t) = cap_pcurve.parameter_range();
    let (support_start_t, support_end_t) = support_pcurve.parameter_range();
    let cap_start = pcurve_surface_point(cap_surface, cap_pcurve, cap_start_t);
    let cap_end = pcurve_surface_point(cap_surface, cap_pcurve, cap_end_t);
    let support_start = pcurve_surface_point(&support_surface, support_pcurve, support_start_t);
    let support_end = pcurve_surface_point(&support_surface, support_pcurve, support_end_t);
    let support_forward = cap_start.distance_to(&support_start) + cap_end.distance_to(&support_end)
        <= cap_start.distance_to(&support_end) + cap_end.distance_to(&support_start);
    let mut worst = 0.0_f64;
    for index in 0..17 {
        let parameter = 0.05 + 0.90 * index as f64 / 16.0;
        let support_parameter = if support_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let cap_uv = cap_pcurve.point_at(cap_start_t + (cap_end_t - cap_start_t) * parameter);
        let support_uv = support_pcurve
            .point_at(support_start_t + (support_end_t - support_start_t) * support_parameter);
        let deficit = interior_normal_deficit_between_surfaces(
            cap_surface,
            cap_uv.x,
            cap_uv.y,
            &support_surface,
            support_uv.x,
            support_uv.y,
        )?;
        worst = worst.max(deficit);
    }
    Some(worst)
}

fn pcurve_surface_point(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

fn interior_normal_deficit_between_surfaces(
    cap_surface: &SurfaceGeometry,
    cap_u: f64,
    cap_v: f64,
    support_surface: &SurfaceGeometry,
    support_u: f64,
    support_v: f64,
) -> Option<f64> {
    let ((cap_u_min, cap_u_max), (cap_v_min, cap_v_max)) = cap_surface.parameter_range();
    let ((support_u_min, support_u_max), (support_v_min, support_v_max)) =
        support_surface.parameter_range();
    let cap_u = interior_parameter(cap_u, cap_u_min, cap_u_max);
    let cap_v = interior_parameter(cap_v, cap_v_min, cap_v_max);
    let support_u = interior_parameter(support_u, support_u_min, support_u_max);
    let support_v = interior_parameter(support_v, support_v_min, support_v_max);
    let cap_normal = cap_surface
        .partial_u(cap_u, cap_v)
        .cross(&cap_surface.partial_v(cap_u, cap_v));
    let support_normal = support_surface
        .partial_u(support_u, support_v)
        .cross(&support_surface.partial_v(support_u, support_v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length <= 1.0e-300 || support_length <= 1.0e-300 {
        return None;
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}

fn pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let destination = brep
        .vertices
        .get(half_edge_destination_vertex(brep, half_edge)?)?
        .position;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    Some(
        start
            .distance_to(&origin)
            .max(end.distance_to(&destination)),
    )
}

fn pcurve_curve_trace_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    curve: &CurveOnEdge,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    let forward = origin.distance_to(&curve_start) <= origin.distance_to(&curve_end);
    let (pc_start, pc_end) = pcurve.parameter_range();
    let mut worst = 0.0_f64;
    for index in 0..17 {
        let fraction = index as f64 / 16.0;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pc_t = pc_start + (pc_end - pc_start) * fraction;
        let uv = pcurve.point_at(pc_t);
        let surface_point = surface.point_at(uv.x, uv.y);
        let gap = surface_point.distance_to(&curve.geometry.point_at(curve_t));
        if !gap.is_finite() {
            return None;
        }
        worst = worst.max(gap);
    }
    Some(worst)
}

fn reverse_linear_pcurve(pcurve: &Pcurve) -> Pcurve {
    let (t_start, t_end) = pcurve.parameter_range();
    let start = pcurve.point_at(t_start);
    let end = pcurve.point_at(t_end);
    Pcurve::linear((end.x, end.y), (start.x, start.y), t_start, t_end)
}

fn cap_support_pair_residual(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
) -> Option<f64> {
    let cap_pcurve = brep.half_edges.get(cap_half_edge)?.pcurve.as_ref()?;
    let support_pcurve = brep.half_edges.get(support_half_edge)?.pcurve.as_ref()?;
    let cap_surface = surface_for_half_edge(brep, cap_half_edge)?;
    let support_surface = surface_for_half_edge(brep, support_half_edge)?;
    let (cap_start, cap_end) = cap_pcurve.parameter_range();
    let (support_start, support_end) = support_pcurve.parameter_range();
    let cap_start_point = pcurve_surface_point(&cap_surface, cap_pcurve, cap_start);
    let cap_end_point = pcurve_surface_point(&cap_surface, cap_pcurve, cap_end);
    let support_start_point = pcurve_surface_point(&support_surface, support_pcurve, support_start);
    let support_end_point = pcurve_surface_point(&support_surface, support_pcurve, support_end);
    let support_forward = cap_start_point.distance_to(&support_start_point)
        + cap_end_point.distance_to(&support_end_point)
        <= cap_start_point.distance_to(&support_end_point)
            + cap_end_point.distance_to(&support_start_point);
    let mut worst = 0.0_f64;
    for index in 0..17 {
        let fraction = index as f64 / 16.0;
        let support_fraction = if support_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let cap_t = cap_start + (cap_end - cap_start) * fraction;
        let support_t = support_start + (support_end - support_start) * support_fraction;
        let cap_uv = cap_pcurve.point_at(cap_t);
        let support_uv = support_pcurve.point_at(support_t);
        let cap_point = cap_surface.point_at(cap_uv.x, cap_uv.y);
        let support_point = support_surface.point_at(support_uv.x, support_uv.y);
        let gap = cap_point.distance_to(&support_point);
        if !gap.is_finite() {
            return None;
        }
        worst = worst.max(gap);
    }
    Some(worst)
}

fn prune_resolved_cap_outer_diagnostics(
    repair_report: &mut crate::types::BlendRepairReport,
    resolved: &[(String, String)],
) {
    if resolved.is_empty() {
        return;
    }
    for route in &mut repair_report.apex_routes {
        route
            .exact_boundary_report
            .diagnostics
            .retain(|diagnostic| !is_resolved_cap_outer_diagnostic(diagnostic, resolved));
    }
    repair_report.apex_diagnostics.retain(|diagnostic| {
        let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
            return true;
        };
        !resolved.iter().any(|(sector_name, outer_side)| {
            cap.stage == "exact_boundary_cap.cap_outer"
                && cap.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded"
                && &cap.sector_name == sector_name
                && cap.outer_side.as_ref() == Some(outer_side)
        })
    });
}

fn is_resolved_cap_outer_diagnostic(
    diagnostic: &crate::types::ExactBoundaryReportDiagnostic,
    resolved: &[(String, String)],
) -> bool {
    resolved.iter().any(|(sector_name, outer_side)| {
        diagnostic.stage == "exact_boundary_cap.cap_outer"
            && diagnostic.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded"
            && &diagnostic.sector_name == sector_name
            && diagnostic.outer_side.as_ref() == Some(outer_side)
    })
}

fn half_edge_destination_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<popkern_core::topology::VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin| twin.origin)
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<SurfaceGeometry> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    brep.faces
        .get(face)?
        .surface
        .as_ref()
        .map(|surface| surface.geometry.clone())
}
