//! Read-only vertex-star classification for coordinated edge/surface repair.

use std::collections::{BTreeSet, HashSet};

use crate::geometry::Point3;
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityOptions, BrepEdgeSurfaceCompatibilityReport,
};
use crate::topology::brep_edge_surface_compatibility_helpers::{
    edge_vertices, incident_edges, proposed_vertex_incident_endpoint_gap,
};
use crate::topology::brep_edge_surface_compatibility_vertex_candidate::{
    BrepEdgeSurfaceVertexCandidateDiagnostic, proposed_vertex_position_for_edge,
};
use crate::topology::{BRep, BrepRepairFactId, BrepRepairKnowledgeGraph, EdgeHandle, VertexHandle};

/// Read-only classification of coordinated vertex-star repair opportunities.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceVertexStarReport {
    /// Number of vertices incident to at least one incompatible edge.
    pub inspected_vertices: usize,
    /// Number of vertex stars with a coherent, policy-compliant candidate.
    pub coherent_vertices: usize,
    /// Number of vertex stars declined.
    pub declined_vertices: usize,
    /// Per-vertex decisions.
    pub decisions: Vec<BrepEdgeSurfaceVertexStarDecision>,
}

/// Per-vertex coordinated repair classification.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceVertexStarDecision {
    /// Vertex whose incident edge star was classified.
    pub vertex: VertexHandle,
    /// Incident incompatible edges considered for coordinated repair.
    pub incompatible_edges: Vec<EdgeHandle>,
    /// All live incident edges at the vertex.
    pub incident_edges: Vec<EdgeHandle>,
    /// Classification outcome.
    pub outcome: BrepEdgeSurfaceVertexStarOutcome,
    /// Largest proposed vertex motion among candidate incident edges.
    pub max_anchor_motion: f64,
    /// Largest disagreement between proposed repaired vertex positions.
    pub max_candidate_spread: f64,
    /// Largest gap a candidate move would introduce on preserved incident edges.
    pub max_preserved_incident_gap: f64,
    /// Candidate vertex position when the star has a coherent proposal.
    pub proposed_position: Option<Point3>,
    /// Shared repair-knowledge facts consumed while classifying this vertex star.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
    /// Per-edge, per-source candidate evidence collected without mutation.
    pub candidate_diagnostics: Vec<BrepEdgeSurfaceVertexCandidateDiagnostic>,
}

/// Vertex-star classification result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceVertexStarOutcome {
    /// The star has a coherent candidate, but this classifier does not mutate.
    CoherentCandidate,
    /// The star cannot be repaired by the current coordinated policy.
    Declined(BrepEdgeSurfaceVertexStarDeclineReason),
}

/// Machine-readable vertex-star decline reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceVertexStarDeclineReason {
    /// Required topology or geometry was missing.
    MissingTopology,
    /// No incident incompatible edge produced a reanchored candidate.
    NoRepairCandidate,
    /// Proposed repair would exceed the caller's anchor-motion policy.
    AnchorMotionExceeded,
    /// Incident incompatible edges proposed different repaired vertex positions.
    InconsistentCandidatePositions,
    /// Moving this vertex would break a preserved incident edge.
    PreservedIncidentConflict,
}

pub(super) fn classify_edge_surface_vertex_stars(
    brep: &BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    report: &BrepEdgeSurfaceCompatibilityReport,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceVertexStarReport {
    let incompatible_edges: HashSet<_> =
        report.findings.iter().map(|finding| finding.edge).collect();
    let mut vertices = Vec::new();
    for edge in &incompatible_edges {
        if let Some((start, end)) = edge_vertices(brep, *edge) {
            vertices.push(start);
            vertices.push(end);
        }
    }
    vertices.sort_by_key(|vertex| format!("{vertex:?}"));
    vertices.dedup();

    let mut decisions = Vec::with_capacity(vertices.len());
    let mut coherent_vertices = 0;
    let mut declined_vertices = 0;
    for vertex in vertices {
        let decision = classify_vertex_star(brep, knowledge, vertex, &incompatible_edges, options);
        match decision.outcome {
            BrepEdgeSurfaceVertexStarOutcome::CoherentCandidate => coherent_vertices += 1,
            BrepEdgeSurfaceVertexStarOutcome::Declined(_) => declined_vertices += 1,
        }
        decisions.push(decision);
    }
    BrepEdgeSurfaceVertexStarReport {
        inspected_vertices: decisions.len(),
        coherent_vertices,
        declined_vertices,
        decisions,
    }
}

fn classify_vertex_star(
    brep: &BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    vertex: VertexHandle,
    incompatible_edges: &HashSet<EdgeHandle>,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceVertexStarDecision {
    let incident_edges = incident_edges(brep, vertex);
    let incompatible_incident_edges: Vec<_> = incident_edges
        .iter()
        .copied()
        .filter(|edge| incompatible_edges.contains(edge))
        .collect();
    let Some(current) = brep.vertices.get(vertex).map(|vertex| vertex.position) else {
        return vertex_decision(
            vertex,
            incompatible_incident_edges,
            incident_edges,
            BrepEdgeSurfaceVertexStarOutcome::Declined(
                BrepEdgeSurfaceVertexStarDeclineReason::MissingTopology,
            ),
            0.0,
            0.0,
            0.0,
            None,
            Vec::new(),
            Vec::new(),
        );
    };
    let mut candidates = Vec::new();
    let mut candidate_diagnostics = Vec::new();
    for edge in &incompatible_incident_edges {
        let candidate = proposed_vertex_position_for_edge(brep, knowledge, *edge, vertex, options);
        candidate_diagnostics.extend(candidate.diagnostics);
        if let Some(point) = candidate.point {
            candidates.push((*edge, point));
        }
    }
    if candidates.len() != incompatible_incident_edges.len() {
        return vertex_decision(
            vertex,
            incompatible_incident_edges,
            incident_edges,
            BrepEdgeSurfaceVertexStarOutcome::Declined(
                BrepEdgeSurfaceVertexStarDeclineReason::NoRepairCandidate,
            ),
            0.0,
            0.0,
            0.0,
            None,
            shared_fact_ids(&candidate_diagnostics),
            candidate_diagnostics,
        );
    }
    let centroid = candidate_centroid(&candidates);
    let max_motion = candidates
        .iter()
        .map(|(_, point)| point.distance_to(&current))
        .fold(0.0, f64::max);
    let max_spread = candidates
        .iter()
        .map(|(_, point)| point.distance_to(&centroid))
        .fold(0.0, f64::max);
    let max_preserved_gap = incident_edges
        .iter()
        .filter(|edge| !incompatible_edges.contains(edge))
        .map(|edge| {
            proposed_vertex_incident_endpoint_gap(brep, *edge, vertex, centroid, options.tolerance)
        })
        .fold(0.0, f64::max);
    let outcome = if max_motion > options.max_anchor_motion {
        BrepEdgeSurfaceVertexStarOutcome::Declined(
            BrepEdgeSurfaceVertexStarDeclineReason::AnchorMotionExceeded,
        )
    } else if max_spread > options.tolerance {
        BrepEdgeSurfaceVertexStarOutcome::Declined(
            BrepEdgeSurfaceVertexStarDeclineReason::InconsistentCandidatePositions,
        )
    } else if max_preserved_gap > options.tolerance {
        BrepEdgeSurfaceVertexStarOutcome::Declined(
            BrepEdgeSurfaceVertexStarDeclineReason::PreservedIncidentConflict,
        )
    } else {
        BrepEdgeSurfaceVertexStarOutcome::CoherentCandidate
    };
    vertex_decision(
        vertex,
        incompatible_incident_edges,
        incident_edges,
        outcome,
        max_motion,
        max_spread,
        max_preserved_gap,
        matches!(outcome, BrepEdgeSurfaceVertexStarOutcome::CoherentCandidate).then_some(centroid),
        shared_fact_ids(&candidate_diagnostics),
        candidate_diagnostics,
    )
}

fn candidate_centroid(candidates: &[(EdgeHandle, Point3)]) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (_, point) in candidates {
        x += point.x;
        y += point.y;
        z += point.z;
    }
    let scale = 1.0 / candidates.len() as f64;
    Point3::new(x * scale, y * scale, z * scale)
}

#[allow(clippy::too_many_arguments)]
fn vertex_decision(
    vertex: VertexHandle,
    incompatible_edges: Vec<EdgeHandle>,
    incident_edges: Vec<EdgeHandle>,
    outcome: BrepEdgeSurfaceVertexStarOutcome,
    max_anchor_motion: f64,
    max_candidate_spread: f64,
    max_preserved_incident_gap: f64,
    proposed_position: Option<Point3>,
    shared_fact_ids: Vec<BrepRepairFactId>,
    candidate_diagnostics: Vec<BrepEdgeSurfaceVertexCandidateDiagnostic>,
) -> BrepEdgeSurfaceVertexStarDecision {
    BrepEdgeSurfaceVertexStarDecision {
        vertex,
        incompatible_edges,
        incident_edges,
        outcome,
        max_anchor_motion,
        max_candidate_spread,
        max_preserved_incident_gap,
        proposed_position,
        shared_fact_ids,
        candidate_diagnostics,
    }
}

fn shared_fact_ids(
    diagnostics: &[BrepEdgeSurfaceVertexCandidateDiagnostic],
) -> Vec<BrepRepairFactId> {
    let fact_ids = diagnostics
        .iter()
        .flat_map(|diagnostic| diagnostic.shared_fact_ids.iter().copied())
        .collect::<BTreeSet<_>>();
    fact_ids.into_iter().collect()
}
