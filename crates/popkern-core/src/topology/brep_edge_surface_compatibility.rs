//! Edge-curve compatibility with adjacent face surfaces.

use crate::geometry::edge_triple::SurfaceSide;
use crate::healing::deviation::edge_surface_deviation;
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
    build_reanchored_intersection_triple, endpoint_gap, face_surface_for_half_edge,
    infer_boundary_hint, unique_edge_half_edges,
};
use crate::topology::brep_edge_surface_compatibility_shared_facts::{
    edge_repair_fact_ids, edge_surface_repair_knowledge, shared_fact_ids_from_report,
    with_shared_fact_ids,
};
use crate::topology::brep_edge_surface_compatibility_vertex_star::{
    BrepEdgeSurfaceVertexStarReport, classify_edge_surface_vertex_stars,
};
use crate::topology::brep_edge_surface_compatibility_vertex_star_repair::{
    edge_is_currently_incompatible, record_repair_decision, repair_coherent_vertex_stars,
};
use crate::topology::{
    BRep, BrepRepairFactId, BrepRepairKnowledgeGraph, EdgeHandle, FaceHandle, HalfEdgeHandle,
    edge_segment,
};

/// Options for direct edge-curve to face-surface compatibility checks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityOptions {
    /// Maximum accepted curve-to-surface distance.
    pub tolerance: f64,
    /// Maximum allowed movement when re-anchoring edge vertices to a verified intersection.
    pub max_anchor_motion: f64,
    /// Number of samples used for direct deviation checks.
    pub samples: usize,
    /// Whether every incompatible edge must be repaired for the transaction to commit.
    pub require_all_compatible: bool,
    /// Whether repaired output must pass full BRep validation.
    pub validate_output: bool,
}

/// Read-only compatibility report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityReport {
    /// Number of half-edges with both a 3D edge curve and adjacent face surface.
    pub inspected_half_edges: usize,
    /// Half-edge uses whose edge curve does not lie on the owning face surface.
    pub incompatible_half_edges: usize,
    /// Worst measured direct edge-curve to face-surface deviation.
    pub max_deviation: f64,
    /// Per-half-edge findings above tolerance.
    pub findings: Vec<BrepEdgeSurfaceCompatibilityFinding>,
}

/// Per-half-edge direct edge/surface compatibility finding.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityFinding {
    /// Shared topological edge.
    pub edge: EdgeHandle,
    /// Face use of the edge.
    pub half_edge: HalfEdgeHandle,
    /// Adjacent face owning the surface.
    pub face: FaceHandle,
    /// Maximum sampled edge-curve to surface distance.
    pub max_deviation: f64,
    /// Curve parameter where the maximum deviation was measured.
    pub max_curve_parameter: f64,
}

/// Successful transactional repair result.
#[derive(Debug, Clone)]
pub struct BrepEdgeSurfaceCompatibilityRepairResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Repair report.
    pub report: BrepEdgeSurfaceCompatibilityRepairReport,
}

/// Aggregate edge/surface compatibility repair report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityRepairReport {
    /// Direct compatibility report before repair.
    pub before: BrepEdgeSurfaceCompatibilityReport,
    /// Direct compatibility report after repair attempts.
    pub after: BrepEdgeSurfaceCompatibilityReport,
    /// Number of topological edges considered for repair.
    pub attempted_edges: usize,
    /// Number of edges whose 3D curve and pcurves were replaced.
    pub repaired_edges: usize,
    /// Number of edges declined.
    pub declined_edges: usize,
    /// Maximum endpoint gap among accepted candidates.
    pub max_endpoint_gap: f64,
    /// Read-only coordinated vertex-star repair classification.
    pub vertex_stars: BrepEdgeSurfaceVertexStarReport,
    /// Shared repair-knowledge facts consumed by this repair report.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
    /// Per-edge decisions.
    pub decisions: Vec<BrepEdgeSurfaceCompatibilityRepairDecision>,
}

/// Per-edge compatibility repair decision.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityRepairDecision {
    /// Edge considered for repair.
    pub edge: EdgeHandle,
    /// Decision outcome.
    pub outcome: BrepEdgeSurfaceCompatibilityRepairOutcome,
    /// Maximum endpoint gap for an accepted candidate.
    pub endpoint_gap: f64,
    /// Maximum residual direct curve-to-adjacent-surface deviation after repair.
    pub residual_deviation: f64,
    /// Structured evidence describing the attempted reconstruction path.
    pub evidence: BrepEdgeSurfaceRepairEvidence,
    /// Shared repair-knowledge facts consumed by this edge decision.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Compatibility repair outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceCompatibilityRepairOutcome {
    /// Edge curve and pcurves were rederived from the adjacent surfaces.
    RepairedFromSurfaceIntersection,
    /// Edge curve and pcurves were rederived from a validated owning surface boundary.
    RepairedFromSurfaceBoundary,
    /// Repair was declined.
    Declined(BrepEdgeSurfaceCompatibilityDeclineReason),
}

/// Machine-readable repair decline reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceCompatibilityDeclineReason {
    /// Required edge, half-edge, loop, face, or vertex topology was missing.
    MissingTopology,
    /// Edge was not bounded by exactly two distinct surfaced face uses.
    MissingAdjacentSurfaces,
    /// Edge endpoints are coincident at the active tolerance.
    DegenerateEdge,
    /// Current generic repair path only commits analytic surface-pair intersections.
    UnsupportedSurfacePair,
    /// A bounded-surface boundary hint was available but no boundary reconstruction applied.
    SurfaceBoundaryReconstructionFailed,
    /// Surface-intersection construction did not produce a candidate.
    IntersectionFailed,
    /// Candidate did not hit the existing topological endpoints within tolerance.
    EndpointMismatch,
    /// Candidate required vertex motion beyond the caller's local repair policy.
    AnchorMotionExceeded,
    /// Candidate vertex motion would break another incident edge endpoint.
    IncidentEndpointConflict,
    /// Candidate still deviated from an adjacent surface beyond tolerance.
    ResidualDeviation,
}

/// Transactional repair failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceCompatibilityRepairFailure {
    /// Machine-readable failure reason.
    pub reason: BrepEdgeSurfaceCompatibilityRepairFailureReason,
    /// Report accumulated before rollback.
    pub report: BrepEdgeSurfaceCompatibilityRepairReport,
}

/// Reason the transaction declined to commit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepEdgeSurfaceCompatibilityRepairFailureReason {
    /// One or more incompatible edge uses remain.
    UnresolvedIncompatibleEdges,
    /// Repaired clone failed full BRep validation.
    ValidationFailed(String),
}

/// Analyze whether existing 3D edge curves lie on all adjacent face surfaces.
pub fn analyze_brep_edge_surface_compatibility(
    brep: &BRep,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceCompatibilityReport {
    let mut report = BrepEdgeSurfaceCompatibilityReport {
        inspected_half_edges: 0,
        incompatible_half_edges: 0,
        max_deviation: 0.0,
        findings: Vec::new(),
    };
    for (edge, edge_data) in &brep.edges {
        let Some(curve) = &edge_data.curve else {
            continue;
        };
        for half_edge in unique_edge_half_edges(brep, edge_data.half_edge) {
            let Some((face, surface)) = face_surface_for_half_edge(brep, half_edge) else {
                continue;
            };
            let deviation = edge_surface_deviation(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                surface,
                options.samples,
            );
            report.inspected_half_edges += 1;
            report.max_deviation = report.max_deviation.max(deviation.max_deviation);
            if deviation.max_deviation > options.tolerance {
                report.incompatible_half_edges += 1;
                report.findings.push(BrepEdgeSurfaceCompatibilityFinding {
                    edge,
                    half_edge,
                    face,
                    max_deviation: deviation.max_deviation,
                    max_curve_parameter: deviation.max_location,
                });
            }
        }
    }
    report
}

/// Repair incompatible edge curves on a cloned BRep.
pub fn repair_brep_edge_surface_compatibility_transactionally(
    brep: &BRep,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> Result<BrepEdgeSurfaceCompatibilityRepairResult, Box<BrepEdgeSurfaceCompatibilityRepairFailure>>
{
    let before = analyze_brep_edge_surface_compatibility(brep, options);
    let knowledge = edge_surface_repair_knowledge(brep, &before);
    let vertex_stars = classify_edge_surface_vertex_stars(brep, &knowledge, &before, options);
    let mut repaired = brep.clone();
    let mut report = BrepEdgeSurfaceCompatibilityRepairReport {
        before,
        after: BrepEdgeSurfaceCompatibilityReport {
            inspected_half_edges: 0,
            incompatible_half_edges: 0,
            max_deviation: 0.0,
            findings: Vec::new(),
        },
        attempted_edges: 0,
        repaired_edges: 0,
        declined_edges: 0,
        max_endpoint_gap: 0.0,
        vertex_stars,
        shared_fact_ids: Vec::new(),
        decisions: Vec::new(),
    };
    let mut edges: Vec<_> = report
        .before
        .findings
        .iter()
        .map(|finding| finding.edge)
        .collect();
    edges.sort_by_key(|edge| format!("{edge:?}"));
    edges.dedup();

    repair_coherent_vertex_stars(&mut repaired, &knowledge, &mut report, options);

    for edge in edges {
        if !edge_is_currently_incompatible(&repaired, edge, options) {
            continue;
        }
        let decision = repair_edge(&mut repaired, &knowledge, edge, options);
        record_repair_decision(&mut report, decision);
    }
    report.after = analyze_brep_edge_surface_compatibility(&repaired, options);
    report.shared_fact_ids = shared_fact_ids_from_report(&report);

    if options.require_all_compatible && report.after.incompatible_half_edges > 0 {
        return Err(Box::new(BrepEdgeSurfaceCompatibilityRepairFailure {
            reason: BrepEdgeSurfaceCompatibilityRepairFailureReason::UnresolvedIncompatibleEdges,
            report,
        }));
    }
    if options.validate_output
        && let Err(err) = repaired.validate()
    {
        return Err(Box::new(BrepEdgeSurfaceCompatibilityRepairFailure {
            reason: BrepEdgeSurfaceCompatibilityRepairFailureReason::ValidationFailed(
                err.to_string(),
            ),
            report,
        }));
    }
    Ok(BrepEdgeSurfaceCompatibilityRepairResult {
        brep: repaired,
        report,
    })
}

pub(super) fn repair_edge(
    brep: &mut BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    edge: EdgeHandle,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    let shared_fact_ids = edge_repair_fact_ids(knowledge, brep, edge);
    let Some((he_a, he_b, face_a, face_b)) = adjacent_surfaced_pair(brep, edge) else {
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
    evidence.max_allowed_anchor_motion = options.max_anchor_motion;
    add_boundary_decline_evidence(brep, he_a, SurfaceSide::A, &mut evidence);
    add_boundary_decline_evidence(brep, he_b, SurfaceSide::B, &mut evidence);
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

#[cfg(test)]
#[path = "brep_edge_surface_compatibility_tests.rs"]
mod tests;
