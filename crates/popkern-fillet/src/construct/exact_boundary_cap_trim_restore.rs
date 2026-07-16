//! Exact-boundary cap trim restoration after generic pcurve repair passes.

#[path = "exact_boundary_cap_trim_restore_diagnostics.rs"]
mod diagnostics;
#[path = "exact_boundary_cap_trim_restore_measure.rs"]
mod measure;

use measure::{
    interior_normal_deficit_between_surfaces, pcurve_curve_trace_gap, pcurve_endpoint_gap,
    pcurve_surface_point, reverse_linear_pcurve, support_normal_deficit_for_candidate,
};

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle};

pub(crate) fn restore_exact_boundary_cap_pcurves(
    result: &mut BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    let routes = repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.exact_boundary_report.candidate_status
                == crate::types::ExactBoundaryCandidateOutcome::Accepted
        })
        .map(|route| route.exact_boundary_report.sectors.clone())
        .collect::<Vec<_>>();
    let mut resolved_cap_outer = Vec::new();
    for sectors in routes {
        restore_exact_boundary_route_sides(result, &sectors, &mut resolved_cap_outer);
    }
    diagnostics::prune_resolved_cap_outer_diagnostics(repair_report, &resolved_cap_outer);
}

/// A side already agreeing to these tolerances (the commercial cap-outer
/// gates) is left untouched: the committed exact-boundary rows are the
/// ground truth and a relock's fixed-knot trace fit can only approximate
/// them (measured: rewriting an 8e-6 FilletEnd row produced a 9.6e-3
/// boundary displacement — PK-124). Both position AND tangency gate the
/// skip: oblique-bridge routes (pyramid lateral faces) rely on the
/// adjacent-row relock for G1 even when G0 is already tight.
const RESTORE_G0_TOLERANCE: f64 = 5.0e-3;
const RESTORE_G1_TOLERANCE: f64 = 5.0e-2;

#[derive(Debug, Clone, Copy, PartialEq)]
struct CapOuterSideQuality {
    g0: f64,
    g1: f64,
}

impl CapOuterSideQuality {
    fn within_gates(self) -> bool {
        self.g0 <= RESTORE_G0_TOLERANCE && self.g1 <= RESTORE_G1_TOLERANCE
    }

    /// A side has not been materially damaged: each metric stays within
    /// its gate or within a small factor of its previous value.
    fn no_worse_than(self, other: CapOuterSideQuality) -> bool {
        self.g0 <= RESTORE_G0_TOLERANCE.max(other.g0 * 1.5)
            && self.g1 <= RESTORE_G1_TOLERANCE.max(other.g1 * 1.5)
    }

    fn better_than(self, other: CapOuterSideQuality) -> bool {
        // Position is never traded for tangency: a relock only counts as
        // an improvement while G0 stays within the gate (or at least does
        // not worsen). The relock's fixed-knot trace fit can displace a
        // committed boundary by ~1e-2 on partition-route caps, which would
        // otherwise be accepted for a G1 gain (measured: pyramid4_mixed
        // shell gap 2.26e-2 + 9.6e-3 seams — PK-124).
        if self.g0 > RESTORE_G0_TOLERANCE.max(other.g0) {
            return false;
        }
        let self_failing = !self.within_gates();
        let other_failing = !other.within_gates();
        match (self_failing, other_failing) {
            (false, true) => true,
            (true, false) => false,
            // Both pass or both fail: compare the dominant defect.
            _ => self.g0.max(self.g1) < other.g0.max(other.g1),
        }
    }
}

fn restore_exact_boundary_route_sides(
    result: &mut BRep,
    sectors: &[crate::types::ExactBoundaryReportSector],
    resolved_cap_outer: &mut Vec<(String, String)>,
) {
    // Every route side participates in the cross-side veto: a relock's
    // full adjacent-row write can displace transverse boundary rows of
    // the SAME surface, so acceptance re-measures all sides.
    let route_sides: Vec<(String, &'static str, HalfEdgeHandle, HalfEdgeHandle)> = sectors
        .iter()
        .flat_map(|sector| {
            [(0usize, "FilletEnd"), (1usize, "FaceEdge")].map(|(side_index, outer_side)| {
                (
                    sector.sector_name.clone(),
                    outer_side,
                    sector.cap_half_edges[side_index],
                    sector.support_half_edges[side_index],
                )
            })
        })
        .collect();
    let measure_all = |brep: &BRep| -> Vec<Option<CapOuterSideQuality>> {
        route_sides
            .iter()
            .map(|(_, _, cap, support)| cap_support_pair_quality(brep, *cap, *support))
            .collect()
    };
    for index in 0..route_sides.len() {
        let (sector_name, outer_side, cap_half_edge, support_half_edge) = (
            route_sides[index].0.clone(),
            route_sides[index].1,
            route_sides[index].2,
            route_sides[index].3,
        );
        let pre_all = measure_all(result);
        let pre_quality = pre_all[index];
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace exact-boundary/relock] {sector_name}/{outer_side}: pre {pre_quality:?}"
            );
        }
        if pre_quality.is_some_and(CapOuterSideQuality::within_gates) {
            resolved_cap_outer.push((sector_name, outer_side.to_string()));
            continue;
        }
        let Some(surface) = surface_for_half_edge(result, cap_half_edge) else {
            continue;
        };
        let Some(candidate) =
            exact_boundary_cap_side_pcurve(result, cap_half_edge, support_half_edge, &surface)
        else {
            continue;
        };
        // Transactional repair: keep the relocked side only when its own
        // measured pair quality improved AND no other route side
        // degraded past its gates (the full adjacent-row write shares
        // corner control points with the transverse boundary rows).
        let snapshot = result.clone();
        super::exact_boundary_cap_trim_surface::relock_exact_boundary_cap_side(
            result,
            cap_half_edge,
            support_half_edge,
            candidate.side,
            candidate.side_forward,
            sectors.len() >= 4,
        );
        if let Some(half_edge) = result.half_edges.get_mut(cap_half_edge) {
            half_edge.pcurve = Some(candidate.pcurve);
        }
        let post_all = measure_all(result);
        let post_quality = post_all[index];
        let own_improved = match (pre_quality, post_quality) {
            (_, None) => false,
            (None, Some(_)) => true,
            (Some(pre), Some(post)) => post.better_than(pre),
        };
        let cross_side_damage = own_improved
            && route_sides.iter().enumerate().any(|(other, _)| {
                if other == index {
                    return false;
                }
                match (pre_all[other], post_all[other]) {
                    (Some(pre), Some(post)) => !post.no_worse_than(pre),
                    (Some(_), None) => true,
                    _ => false,
                }
            });
        if !own_improved || cross_side_damage {
            *result = snapshot;
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace exact-boundary/relock] {sector_name}/{outer_side}: \
                     declined (own {pre_quality:?} -> {post_quality:?}, cross_damage={cross_side_damage})"
                );
            }
            continue;
        }
        if post_quality.is_some_and(|quality| quality.g0 <= RESTORE_G0_TOLERANCE) {
            resolved_cap_outer.push((sector_name, outer_side.to_string()));
        }
    }
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

fn cap_support_pair_quality(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
) -> Option<CapOuterSideQuality> {
    let cap_pcurve = brep.half_edges.get(cap_half_edge)?.pcurve.as_ref()?;
    let support_pcurve = brep.half_edges.get(support_half_edge)?.pcurve.as_ref()?;
    let cap_surface = surface_for_half_edge(brep, cap_half_edge)?;
    let support_surface = surface_for_half_edge(brep, support_half_edge)?;
    let edge = brep.half_edges.get(cap_half_edge)?.edge;
    let edge_curve = brep.edges.get(edge)?.curve.as_ref()?;
    let (cap_start, cap_end) = cap_pcurve.parameter_range();
    let (support_start, support_end) = support_pcurve.parameter_range();
    let cap_forward = pcurve_follows_edge_curve(&cap_surface, cap_pcurve, edge_curve);
    let support_forward = pcurve_follows_edge_curve(&support_surface, support_pcurve, edge_curve);
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    for index in 0..100 {
        let fraction = 0.05 + 0.90 * index as f64 / 99.0;
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
        let deficit = interior_normal_deficit_between_surfaces(
            &cap_surface,
            cap_uv.x,
            cap_uv.y,
            &support_surface,
            support_uv.x,
            support_uv.y,
        )?;
        if !gap.is_finite() || !deficit.is_finite() {
            return None;
        }
        worst_g0 = worst_g0.max(gap);
        worst_g1 = worst_g1.max(deficit);
    }
    Some(CapOuterSideQuality {
        g0: worst_g0,
        g1: worst_g1,
    })
}

fn pcurve_follows_edge_curve(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    edge_curve: &CurveOnEdge,
) -> bool {
    let (pc_start, pc_end) = pcurve.parameter_range();
    let pcurve_start = pcurve_surface_point(surface, pcurve, pc_start);
    let pcurve_end = pcurve_surface_point(surface, pcurve, pc_end);
    let curve_start = edge_curve.geometry.point_at(edge_curve.t_start);
    let curve_end = edge_curve.geometry.point_at(edge_curve.t_end);
    curve_start.distance_to(&pcurve_start) + curve_end.distance_to(&pcurve_end)
        <= curve_start.distance_to(&pcurve_end) + curve_end.distance_to(&pcurve_start)
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
