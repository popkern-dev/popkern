//! Commit helpers for edge/surface compatibility repair.

use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::edge_triple::EdgeTriple;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::brep_edge_surface_compatibility_helpers::{
    build_reanchored_intersection_triple, edge_vertices, endpoint_gap, move_vertex_with_tolerance,
    proposed_vertex_incident_endpoint_gap, residual_for_candidate,
};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, edge_segment};

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
    evidence: BrepEdgeSurfaceRepairEvidence,
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
    let residual_deviation = residual_for_candidate(
        &triple.curve_3d,
        (triple.t_start, triple.t_end),
        surf_a,
        surf_b,
        options.samples,
    );
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
        if incident_gap > options.tolerance {
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
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry,
        t_start: triple.t_start,
        t_end: triple.t_end,
    });
    brep.half_edges[he_a].pcurve = Some(triple.pcurve_a);
    brep.half_edges[he_b].pcurve = Some(triple.pcurve_b);
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
