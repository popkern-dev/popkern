//! Per-edge orchestration for edge/surface compatibility repair.

use crate::geometry::edge_triple::SurfaceSide;
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_boundary::is_repairable_surface_pair;
use crate::topology::brep_edge_surface_compatibility_commit::{
    anchor_motion_allowed, commit_edge_triple, try_endpoint_reanchored_repair,
};
use crate::topology::brep_edge_surface_compatibility_decline::{
    declined, declined_without_evidence,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::{
    BrepEdgeSurfaceBoundaryDeclineEvidence, BrepEdgeSurfaceRepairEvidence,
};
use crate::topology::brep_edge_surface_compatibility_helpers::{
    adjacent_surfaced_pair, build_edge_triple, build_pcurve_scoped_reanchored_intersection_triple,
    build_reanchored_intersection_triple, endpoint_gap, infer_boundary_hint,
    single_surfaced_half_edge,
};
use crate::topology::brep_edge_surface_compatibility_pcurve_trace::try_commit_pcurve_trace_edge;
use crate::topology::brep_edge_surface_compatibility_shared_facts::{
    edge_repair_diagnostic_codes, edge_repair_fact_ids, with_shared_fact_ids,
};
use crate::topology::brep_edge_surface_compatibility_single_surface::commit_single_surface_edge;
use crate::topology::{BRep, BrepRepairKnowledgeGraph, EdgeHandle, HalfEdgeHandle, edge_segment};

pub(super) fn repair_edge(
    brep: &mut BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    edge: EdgeHandle,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    let shared_fact_ids = edge_repair_fact_ids(knowledge, brep, edge);
    let Some((he_a, he_b, face_a, face_b)) = adjacent_surfaced_pair(brep, edge) else {
        if let Some((half_edge, face)) = single_surfaced_half_edge(brep, edge)
            && let Some(surface) = brep
                .faces
                .get(face)
                .and_then(|face| face.surface.as_ref())
                .map(|surface| surface.geometry.clone())
        {
            let mut evidence = BrepEdgeSurfaceRepairEvidence::from_single_surface(&surface);
            evidence.repair_knowledge_diagnostics = edge_repair_diagnostic_codes(knowledge, edge);
            return with_shared_fact_ids(
                commit_single_surface_edge(brep, edge, half_edge, &surface, options, evidence),
                shared_fact_ids,
            );
        }
        return declined_without_evidence(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingAdjacentSurfaces,
            shared_fact_ids,
        );
    };
    let Some((start, end)) = edge_segment(brep, edge) else {
        return declined_without_evidence(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
            shared_fact_ids,
        );
    };
    if start.distance_to(&end) <= options.tolerance.max(1e-12) {
        return declined_without_evidence(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::DegenerateEdge,
            shared_fact_ids,
        );
    }
    let (Some(surface_a), Some(surface_b)) = (
        brep.faces
            .get(face_a)
            .and_then(|face| face.surface.as_ref()),
        brep.faces
            .get(face_b)
            .and_then(|face| face.surface.as_ref()),
    ) else {
        return declined_without_evidence(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingAdjacentSurfaces,
            shared_fact_ids,
        );
    };
    let surf_a = surface_a.geometry.clone();
    let surf_b = surface_b.geometry.clone();
    let boundary_hint = infer_boundary_hint(brep, he_a, he_b, &surf_a, &surf_b, options.tolerance);
    let mut evidence =
        BrepEdgeSurfaceRepairEvidence::from_surfaces(&surf_a, &surf_b, boundary_hint.as_ref());
    evidence.repair_knowledge_diagnostics = edge_repair_diagnostic_codes(knowledge, edge);
    evidence.max_allowed_anchor_motion = options.max_anchor_motion;
    add_boundary_decline_evidence(brep, he_a, SurfaceSide::A, &mut evidence);
    add_boundary_decline_evidence(brep, he_b, SurfaceSide::B, &mut evidence);
    if let Some(decision) = try_commit_pcurve_trace_edge(
        brep,
        edge,
        he_a,
        he_b,
        &surf_a,
        &surf_b,
        options,
        evidence.clone(),
    ) && !matches!(
        decision.outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(_)
    ) {
        return with_shared_fact_ids(decision, shared_fact_ids);
    }
    let boundary_attempt = if let Some(hint) = boundary_hint.clone() {
        evidence.attempted_boundary_reconstruction = true;
        match build_edge_triple(&surf_a, &surf_b, start, end, Some(hint), brep) {
            Ok(triple) => Some(triple),
            Err(label) => {
                evidence.boundary_failure = Some(label);
                None
            }
        }
    } else {
        None
    };
    let (triple, repaired_from_boundary) = if let Some(triple) = boundary_attempt {
        (triple, true)
    } else if is_repairable_surface_pair(&surf_a, &surf_b) {
        evidence.attempted_surface_intersection = true;
        match build_edge_triple(&surf_a, &surf_b, start, end, None, brep) {
            Ok(triple) => (triple, false),
            Err(label) => {
                evidence.intersection_failure = Some(label);
                evidence.attempted_anchor_repair = true;
                match build_pcurve_scoped_reanchored_intersection_triple(
                    &surf_a,
                    &surf_b,
                    start,
                    end,
                    brep,
                    he_a,
                    he_b,
                    options.samples,
                )
                .or_else(|_| {
                    build_reanchored_intersection_triple(
                        &surf_a,
                        &surf_b,
                        start,
                        end,
                        brep,
                        options.samples,
                    )
                }) {
                    Ok((triple, motion)) if anchor_motion_allowed(motion, options) => {
                        evidence.anchor_motion = motion;
                        (triple, false)
                    }
                    Ok((_, motion)) => {
                        evidence.anchor_motion = motion;
                        return declined(
                            edge,
                            BrepEdgeSurfaceCompatibilityDeclineReason::AnchorMotionExceeded,
                            evidence,
                            shared_fact_ids,
                        );
                    }
                    Err(anchor_label) => {
                        evidence.intersection_failure = Some(anchor_label);
                        return declined(
                            edge,
                            BrepEdgeSurfaceCompatibilityDeclineReason::IntersectionFailed,
                            evidence,
                            shared_fact_ids,
                        );
                    }
                }
            }
        }
    } else if boundary_hint.is_some() {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::SurfaceBoundaryReconstructionFailed,
            evidence,
            shared_fact_ids,
        );
    } else {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::UnsupportedSurfacePair,
            evidence,
            shared_fact_ids,
        );
    };
    let endpoint_gap = endpoint_gap(&triple.curve_3d, triple.t_start, triple.t_end, start, end);
    if endpoint_gap > options.tolerance && evidence.anchor_motion == 0.0 {
        if !repaired_from_boundary
            && is_repairable_surface_pair(&surf_a, &surf_b)
            && let Some(decision) = try_endpoint_reanchored_repair(
                brep,
                edge,
                he_a,
                he_b,
                &surf_a,
                &surf_b,
                start,
                end,
                options,
                evidence.clone(),
            )
        {
            return with_shared_fact_ids(decision, shared_fact_ids);
        }
        return BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::EndpointMismatch,
            ),
            endpoint_gap,
            residual_deviation: 0.0,
            evidence,
            shared_fact_ids,
        };
    }
    with_shared_fact_ids(
        commit_edge_triple(
            brep,
            edge,
            he_a,
            he_b,
            triple,
            &surf_a,
            &surf_b,
            options,
            evidence,
            repaired_from_boundary,
        ),
        shared_fact_ids,
    )
}

fn add_boundary_decline_evidence(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    side: SurfaceSide,
    evidence: &mut BrepEdgeSurfaceRepairEvidence,
) {
    let Some(decline) = brep.surface_boundary_decline(half_edge) else {
        return;
    };
    evidence
        .boundary_declines
        .push(BrepEdgeSurfaceBoundaryDeclineEvidence {
            side,
            reason: decline.reason.clone(),
            best_boundary: decline.best_boundary,
            best_boundary_deviation: decline.best_boundary_deviation,
            endpoint_gap: decline.endpoint_gap,
        });
}
