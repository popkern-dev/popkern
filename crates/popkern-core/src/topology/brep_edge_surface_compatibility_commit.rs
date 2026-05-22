//! Commit helpers for edge/surface compatibility repair.

use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::edge_triple::EdgeTriple;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::brep_curve_trim_reconstruction_helpers::project_edge_curve_to_surface_pcurve;
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::brep_edge_surface_compatibility_helpers::{
    build_reanchored_intersection_triple, edge_vertices, endpoint_gap, move_vertex_with_tolerance,
    proposed_vertex_incident_endpoint_gap, residuals_for_candidate,
};
use crate::topology::{
    BRep, EdgeHandle, HalfEdgeHandle, SharedTrimEdgeInput, commit_shared_trim_edge_transactionally,
    edge_segment,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn try_endpoint_reanchored_repair(
    brep: &mut BRep,
    edge: EdgeHandle,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    start: crate::geometry::Point3,
    end: crate::geometry::Point3,
    options: BrepEdgeSurfaceCompatibilityOptions,
    mut evidence: BrepEdgeSurfaceRepairEvidence,
) -> Option<BrepEdgeSurfaceCompatibilityRepairDecision> {
    evidence.attempted_anchor_repair = true;
    match build_reanchored_intersection_triple(surf_a, surf_b, start, end, brep, options.samples) {
        Ok((triple, motion)) if anchor_motion_allowed(motion, options) => {
            evidence.anchor_motion = motion;
            Some(commit_edge_triple(
                brep, edge, he_a, he_b, triple, surf_a, surf_b, options, evidence, false,
            ))
        }
        Ok((_, motion)) => {
            evidence.anchor_motion = motion;
            Some(declined(
                edge,
                BrepEdgeSurfaceCompatibilityDeclineReason::AnchorMotionExceeded,
                evidence,
            ))
        }
        Err(anchor_label) => {
            evidence.intersection_failure = Some(anchor_label);
            None
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn commit_edge_triple(
    brep: &mut BRep,
    edge: EdgeHandle,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    triple: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    options: BrepEdgeSurfaceCompatibilityOptions,
    mut evidence: BrepEdgeSurfaceRepairEvidence,
    repaired_from_boundary: bool,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    let Some((start, end)) = edge_segment(brep, edge) else {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
            evidence,
        );
    };
    let endpoint_gap = endpoint_gap(&triple.curve_3d, triple.t_start, triple.t_end, start, end);
    let (residual_deviation_a, residual_deviation_b) = residuals_for_candidate(
        &triple.curve_3d,
        (triple.t_start, triple.t_end),
        surf_a,
        surf_b,
        options.samples,
    );
    evidence.residual_deviation_a = residual_deviation_a;
    evidence.residual_deviation_b = residual_deviation_b;
    let residual_deviation = residual_deviation_a.max(residual_deviation_b);
    if residual_deviation > options.tolerance {
        return BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::ResidualDeviation,
            ),
            endpoint_gap,
            residual_deviation,
            evidence,
            shared_fact_ids: Vec::new(),
        };
    }
    let Some((start_vertex, end_vertex)) = edge_vertices(brep, edge) else {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
            evidence,
        );
    };
    let rollback_snapshot = brep.clone();
    let repaired_start = triple.curve_3d.point_at(triple.t_start);
    let repaired_end = triple.curve_3d.point_at(triple.t_end);
    if endpoint_gap > 0.0 {
        let incident_gap = proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            start_vertex,
            repaired_start,
            options.tolerance,
        )
        .max(proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            end_vertex,
            repaired_end,
            options.tolerance,
        ));
        if incident_gap > options.max_anchor_motion {
            let mut evidence = evidence;
            evidence.max_incident_endpoint_gap = incident_gap;
            return BrepEdgeSurfaceCompatibilityRepairDecision {
                edge,
                outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                    BrepEdgeSurfaceCompatibilityDeclineReason::IncidentEndpointConflict,
                ),
                endpoint_gap,
                residual_deviation,
                evidence,
                shared_fact_ids: Vec::new(),
            };
        }
        move_vertex_with_tolerance(brep, start_vertex, repaired_start, endpoint_gap);
        move_vertex_with_tolerance(brep, end_vertex, repaired_end, endpoint_gap);
    }
    let geometry = triple
        .curve_3d
        .canonical_line_in_range(
            (triple.t_start, triple.t_end),
            options.tolerance.max(1.0e-8),
        )
        .unwrap_or_else(|| triple.curve_3d.clone());
    let curve = CurveOnEdge {
        geometry,
        t_start: triple.t_start,
        t_end: triple.t_end,
    };
    let mut pcurve_source = brep.clone();
    pcurve_source.edges[edge].curve = Some(curve.clone());
    let pcurve_a = project_edge_curve_to_surface_pcurve(&pcurve_source, he_a, options.samples)
        .unwrap_or_else(|| orient_pcurve_to_half_edge(brep, he_a, surf_a, triple.pcurve_a));
    let pcurve_b = project_edge_curve_to_surface_pcurve(&pcurve_source, he_b, options.samples)
        .unwrap_or_else(|| orient_pcurve_to_half_edge(brep, he_b, surf_b, triple.pcurve_b));
    if let Err(failure) = commit_shared_trim_edge_transactionally(
        brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a: he_a,
            half_edge_b: he_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: options.tolerance,
            validate_brep: options.validate_output,
            provenance: "edge_surface_compatibility_repair",
        },
    ) {
        *brep = rollback_snapshot;
        evidence.shared_trim_failure = Some(failure.reason);
        evidence.shared_trim_curve_surface_deviation =
            failure.validation.max_curve_surface_deviation;
        evidence.shared_trim_side_agreement_deviation =
            failure.validation.max_side_agreement_deviation;
        evidence.shared_trim_endpoint_gap = failure.validation.max_endpoint_gap;
        return BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::SharedTrimValidationFailed,
            ),
            endpoint_gap,
            residual_deviation,
            evidence,
            shared_fact_ids: Vec::new(),
        };
    }
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: if repaired_from_boundary {
            BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
        } else {
            BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceIntersection
        },
        endpoint_gap,
        residual_deviation,
        evidence,
        shared_fact_ids: Vec::new(),
    }
}

fn orient_pcurve_to_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: Pcurve,
) -> Pcurve {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return pcurve;
    };
    let Some(next) = brep.half_edges.get(half_edge_data.next) else {
        return pcurve;
    };
    let Some(start) = brep
        .vertices
        .get(half_edge_data.origin)
        .map(|vertex| vertex.position)
    else {
        return pcurve;
    };
    let Some(end) = brep.vertices.get(next.origin).map(|vertex| vertex.position) else {
        return pcurve;
    };
    let p_start = surface_point_at_pcurve_endpoint(surface, &pcurve, pcurve.t_start);
    let p_end = surface_point_at_pcurve_endpoint(surface, &pcurve, pcurve.t_end);
    let forward_gap = p_start.distance_to(&start).max(p_end.distance_to(&end));
    let reversed_gap = p_start.distance_to(&end).max(p_end.distance_to(&start));
    if reversed_gap < forward_gap {
        pcurve.reversed()
    } else {
        pcurve
    }
}

fn surface_point_at_pcurve_endpoint(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
) -> crate::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

pub(super) fn anchor_motion_allowed(
    motion: f64,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> bool {
    motion <= options.max_anchor_motion
}

fn declined(
    edge: EdgeHandle,
    reason: BrepEdgeSurfaceCompatibilityDeclineReason,
    evidence: BrepEdgeSurfaceRepairEvidence,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(reason),
        endpoint_gap: 0.0,
        residual_deviation: 0.0,
        evidence,
        shared_fact_ids: Vec::new(),
    }
}
