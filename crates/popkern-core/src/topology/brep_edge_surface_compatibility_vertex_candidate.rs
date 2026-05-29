//! Candidate-source diagnostics for edge/surface vertex-star classification.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, SurfaceGeometry};
use crate::topology::brep_edge_surface_compatibility::BrepEdgeSurfaceCompatibilityOptions;
use crate::topology::brep_edge_surface_compatibility_boundary::{
    infer_surface_boundary_side_read_only, is_repairable_surface_pair, surface_boundary_curve,
};
use crate::topology::brep_edge_surface_compatibility_candidate_facts::{
    edge_surface_candidate_fact_ids, surface_boundary_fact_ids, trim_endpoint_fact_ids,
};
use crate::topology::brep_edge_surface_compatibility_helpers::{
    adjacent_surfaced_pair, build_pcurve_scoped_reanchored_intersection_triple,
    build_reanchored_intersection_triple, edge_vertices,
};
use crate::topology::brep_edge_surface_compatibility_trim_candidate::{
    TrimEndpointCandidate, trim_decline, trim_endpoint_candidates_for_side,
};
use crate::topology::brep_edge_surface_compatibility_vertex_candidate_geometry::{
    closest_boundary_parameter, point_centroid,
};
use crate::topology::{
    BRep, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, EdgeHandle,
    HalfEdgeHandle, VertexHandle,
};

/// Read-only candidate-source diagnostic for one edge incident to a vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeSurfaceVertexCandidateDiagnostic {
    /// Edge whose endpoint candidate was attempted.
    pub edge: EdgeHandle,
    /// Half-edge side for side-local sources.
    pub half_edge: Option<HalfEdgeHandle>,
    /// Candidate source that was attempted.
    pub source: BrepEdgeSurfaceVertexCandidateSource,
    /// Attempt outcome.
    pub outcome: BrepEdgeSurfaceVertexCandidateOutcome,
    /// Candidate point when the source accepted.
    pub point: Option<Point3>,
    /// Motion from the current topological vertex to the candidate point.
    pub anchor_motion: f64,
    /// Side-to-side disagreement for aggregate bounded-surface candidates.
    pub side_disagreement: f64,
    /// Shared repair-knowledge facts consumed by this diagnostic.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Source of vertex candidate evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceVertexCandidateSource {
    /// Reanchored adjacent-surface intersection.
    SurfaceIntersection,
    /// Canonical bounded-surface boundary curve.
    SurfaceBoundary,
    /// Existing half-edge trim endpoint evaluated on the owning surface.
    TrimEndpoint,
    /// Agreement check across side-local candidate sources.
    SideAgreement,
}

/// Machine-readable candidate-source outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeSurfaceVertexCandidateOutcome {
    /// Candidate source produced usable evidence.
    Accepted,
    /// Candidate was accepted by using the opposite trim endpoint because it matches this vertex.
    ReconciledEndpoint,
    /// Required topology or geometry was missing.
    MissingTopology,
    /// Source does not apply to this surface or half-edge.
    SourceUnavailable,
    /// Source evaluation failed numerically.
    EvaluationFailed,
    /// Candidate would exceed caller anchor-motion policy.
    AnchorMotionExceeded,
    /// Candidate sides disagreed beyond local tolerance.
    SideDisagreement,
}

pub(super) struct EdgeVertexCandidate {
    pub point: Option<Point3>,
    pub diagnostics: Vec<BrepEdgeSurfaceVertexCandidateDiagnostic>,
}

impl TrimEndpointCandidate {
    fn diagnostic(
        &self,
        edge: EdgeHandle,
        half_edge: HalfEdgeHandle,
        shared_fact_ids: Vec<BrepRepairFactId>,
    ) -> BrepEdgeSurfaceVertexCandidateDiagnostic {
        candidate_diagnostic(
            edge,
            Some(half_edge),
            CandidateDiagnosticInput {
                source: BrepEdgeSurfaceVertexCandidateSource::TrimEndpoint,
                outcome: self.outcome,
                point: self.point,
                anchor_motion: self.anchor_motion,
                side_disagreement: 0.0,
                shared_fact_ids,
            },
        )
    }
}

pub(super) fn proposed_vertex_position_for_edge(
    brep: &BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    edge: EdgeHandle,
    vertex: VertexHandle,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> EdgeVertexCandidate {
    let mut diagnostics = Vec::new();
    let Some((he_a, he_b, face_a, face_b)) = adjacent_surfaced_pair(brep, edge) else {
        diagnostics.push(missing_topology_diagnostic(edge));
        return no_candidate(diagnostics);
    };
    let Some(surf_a) = face_surface(brep, face_a) else {
        return missing_surface_candidate(edge, diagnostics);
    };
    let Some(surf_b) = face_surface(brep, face_b) else {
        return missing_surface_candidate(edge, diagnostics);
    };
    let Some((start, end)) = crate::topology::edge_segment(brep, edge) else {
        return missing_surface_candidate(edge, diagnostics);
    };
    let Some((start_vertex, end_vertex)) = edge_vertices(brep, edge) else {
        return missing_surface_candidate(edge, diagnostics);
    };
    if !has_boundary_owned_surface(surf_a, surf_b)
        && is_repairable_surface_pair(surf_a, surf_b)
        && let Ok((triple, _)) = build_pcurve_scoped_reanchored_intersection_triple(
            surf_a,
            surf_b,
            start,
            end,
            brep,
            he_a,
            he_b,
            options.samples,
        )
        .or_else(|_| {
            build_reanchored_intersection_triple(surf_a, surf_b, start, end, brep, options.samples)
        })
    {
        let point = if vertex == start_vertex {
            Some(triple.curve_3d.point_at(triple.t_start))
        } else if vertex == end_vertex {
            Some(triple.curve_3d.point_at(triple.t_end))
        } else {
            None
        };
        if let Some(point) = point {
            diagnostics.push(accepted_candidate_diagnostic(
                brep,
                edge,
                None,
                BrepEdgeSurfaceVertexCandidateSource::SurfaceIntersection,
                vertex,
                point,
            ));
            return EdgeVertexCandidate {
                point: Some(point),
                diagnostics,
            };
        }
    }
    diagnostics.push(candidate_diagnostic(
        edge,
        None,
        CandidateDiagnosticInput {
            source: BrepEdgeSurfaceVertexCandidateSource::SurfaceIntersection,
            outcome: BrepEdgeSurfaceVertexCandidateOutcome::SourceUnavailable,
            point: None,
            anchor_motion: 0.0,
            side_disagreement: 0.0,
            shared_fact_ids: Vec::new(),
        },
    ));
    let boundary = boundary_vertex_position_for_edge(
        brep,
        knowledge,
        edge,
        [(he_a, surf_a), (he_b, surf_b)],
        vertex,
        options,
    );
    diagnostics.extend(boundary.diagnostics);
    EdgeVertexCandidate {
        point: boundary.point,
        diagnostics,
    }
}

fn has_boundary_owned_surface(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    matches!(
        (a, b),
        (SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_), _)
            | (_, SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_))
    )
}

fn boundary_vertex_position_for_edge(
    brep: &BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    edge: EdgeHandle,
    half_edges: [(HalfEdgeHandle, &SurfaceGeometry); 2],
    vertex: VertexHandle,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> EdgeVertexCandidate {
    let Some(current) = brep.vertices.get(vertex).map(|vertex| vertex.position) else {
        return no_candidate(vec![missing_topology_diagnostic(edge)]);
    };
    let mut side_candidates = Vec::new();
    let mut diagnostics = Vec::new();
    let side_half_edges = [half_edges[0].0, half_edges[1].0];
    for (half_edge, surface) in half_edges {
        let surface_boundary_fact_ids =
            surface_boundary_fact_ids(knowledge, half_edge, BrepRepairFactKind::SurfaceBoundary);
        if let Some(point) =
            surface_boundary_endpoint_candidate(brep, half_edge, surface, current, options)
        {
            diagnostics.push(candidate_diagnostic(
                edge,
                Some(half_edge),
                CandidateDiagnosticInput {
                    source: BrepEdgeSurfaceVertexCandidateSource::SurfaceBoundary,
                    outcome: BrepEdgeSurfaceVertexCandidateOutcome::Accepted,
                    point: Some(point),
                    anchor_motion: point.distance_to(&current),
                    side_disagreement: 0.0,
                    shared_fact_ids: surface_boundary_fact_ids,
                },
            ));
            side_candidates.push(vec![TrimEndpointCandidate {
                point: Some(point),
                outcome: BrepEdgeSurfaceVertexCandidateOutcome::Accepted,
                anchor_motion: point.distance_to(&current),
            }]);
            continue;
        }
        diagnostics.push(candidate_diagnostic(
            edge,
            Some(half_edge),
            CandidateDiagnosticInput {
                source: BrepEdgeSurfaceVertexCandidateSource::SurfaceBoundary,
                outcome: BrepEdgeSurfaceVertexCandidateOutcome::SourceUnavailable,
                point: None,
                anchor_motion: 0.0,
                side_disagreement: 0.0,
                shared_fact_ids: surface_boundary_fact_ids,
            },
        ));
        let trim_candidates =
            trim_endpoint_candidates_for_side(brep, half_edge, surface, current, options);
        let shared_fact_ids = trim_endpoint_fact_ids(
            knowledge,
            half_edge,
            BrepRepairFactKind::TrimEndpointAgreement,
        );
        if trim_candidates.is_empty() {
            diagnostics.push(
                trim_decline(BrepEdgeSurfaceVertexCandidateOutcome::SourceUnavailable).diagnostic(
                    edge,
                    half_edge,
                    shared_fact_ids,
                ),
            );
        } else {
            diagnostics.extend(
                trim_candidates.iter().map(|candidate| {
                    candidate.diagnostic(edge, half_edge, shared_fact_ids.clone())
                }),
            );
            side_candidates.push(trim_candidates);
        }
    }
    if side_candidates.is_empty() {
        return no_candidate(diagnostics);
    }
    let Some((centroid, max_spread)) = best_side_candidate_agreement(&side_candidates) else {
        return no_candidate(diagnostics);
    };
    let accepted = max_spread <= options.tolerance;
    diagnostics.push(candidate_diagnostic(
        edge,
        None,
        CandidateDiagnosticInput {
            source: BrepEdgeSurfaceVertexCandidateSource::SideAgreement,
            outcome: if accepted {
                BrepEdgeSurfaceVertexCandidateOutcome::Accepted
            } else {
                BrepEdgeSurfaceVertexCandidateOutcome::SideDisagreement
            },
            point: accepted.then_some(centroid),
            anchor_motion: centroid.distance_to(&current),
            side_disagreement: max_spread,
            shared_fact_ids: edge_surface_candidate_fact_ids(knowledge, &side_half_edges),
        },
    ));
    EdgeVertexCandidate {
        point: accepted.then_some(centroid),
        diagnostics,
    }
}

fn face_surface(brep: &BRep, face: crate::topology::FaceHandle) -> Option<&SurfaceGeometry> {
    brep.faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
}

fn missing_surface_candidate(
    edge: EdgeHandle,
    mut diagnostics: Vec<BrepEdgeSurfaceVertexCandidateDiagnostic>,
) -> EdgeVertexCandidate {
    diagnostics.push(missing_topology_diagnostic(edge));
    no_candidate(diagnostics)
}

fn no_candidate(diagnostics: Vec<BrepEdgeSurfaceVertexCandidateDiagnostic>) -> EdgeVertexCandidate {
    EdgeVertexCandidate {
        point: None,
        diagnostics,
    }
}

fn surface_boundary_endpoint_candidate(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    current: Point3,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> Option<Point3> {
    let boundary =
        infer_surface_boundary_side_read_only(brep, half_edge, surface, options.tolerance)?;
    let curve = surface_boundary_curve(surface, boundary)?;
    boundary_endpoint_candidate(&curve, current, options.samples, options.tolerance)
}

fn boundary_endpoint_candidate(
    curve: &NurbsCurve,
    target: Point3,
    samples: usize,
    tolerance: f64,
) -> Option<Point3> {
    let (t0, t1) = curve.parameter_range();
    let t = closest_boundary_parameter(curve, target, t0, t1, samples.max(17))?;
    let point = curve.point_at(t);
    if point.distance_to(&target) <= tolerance {
        Some(target)
    } else {
        Some(point)
    }
}

fn best_side_candidate_agreement(
    candidates: &[Vec<TrimEndpointCandidate>],
) -> Option<(Point3, f64)> {
    if candidates.len() == 1 {
        let point = candidates[0].iter().find_map(|candidate| candidate.point)?;
        return Some((point, 0.0));
    }
    let mut best: Option<(Point3, f64)> = None;
    for first in &candidates[0] {
        let Some(first_point) = first.point else {
            continue;
        };
        for second in &candidates[1] {
            let Some(second_point) = second.point else {
                continue;
            };
            let centroid = point_centroid(&[first_point, second_point]);
            let spread = first_point
                .distance_to(&centroid)
                .max(second_point.distance_to(&centroid));
            if best
                .as_ref()
                .is_none_or(|(_, best_spread)| spread < *best_spread)
            {
                best = Some((centroid, spread));
            }
        }
    }
    best
}

fn accepted_candidate_diagnostic(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: Option<HalfEdgeHandle>,
    source: BrepEdgeSurfaceVertexCandidateSource,
    vertex: VertexHandle,
    point: Point3,
) -> BrepEdgeSurfaceVertexCandidateDiagnostic {
    let motion = brep
        .vertices
        .get(vertex)
        .map(|vertex| point.distance_to(&vertex.position))
        .unwrap_or(0.0);
    candidate_diagnostic(
        edge,
        half_edge,
        CandidateDiagnosticInput {
            source,
            outcome: BrepEdgeSurfaceVertexCandidateOutcome::Accepted,
            point: Some(point),
            anchor_motion: motion,
            side_disagreement: 0.0,
            shared_fact_ids: Vec::new(),
        },
    )
}

fn missing_topology_diagnostic(edge: EdgeHandle) -> BrepEdgeSurfaceVertexCandidateDiagnostic {
    candidate_diagnostic(
        edge,
        None,
        CandidateDiagnosticInput {
            source: BrepEdgeSurfaceVertexCandidateSource::SurfaceIntersection,
            outcome: BrepEdgeSurfaceVertexCandidateOutcome::MissingTopology,
            point: None,
            anchor_motion: 0.0,
            side_disagreement: 0.0,
            shared_fact_ids: Vec::new(),
        },
    )
}

struct CandidateDiagnosticInput {
    source: BrepEdgeSurfaceVertexCandidateSource,
    outcome: BrepEdgeSurfaceVertexCandidateOutcome,
    point: Option<Point3>,
    anchor_motion: f64,
    side_disagreement: f64,
    shared_fact_ids: Vec<BrepRepairFactId>,
}

fn candidate_diagnostic(
    edge: EdgeHandle,
    half_edge: Option<HalfEdgeHandle>,
    input: CandidateDiagnosticInput,
) -> BrepEdgeSurfaceVertexCandidateDiagnostic {
    BrepEdgeSurfaceVertexCandidateDiagnostic {
        edge,
        half_edge,
        source: input.source,
        outcome: input.outcome,
        point: input.point,
        anchor_motion: input.anchor_motion,
        side_disagreement: input.side_disagreement,
        shared_fact_ids: input.shared_fact_ids,
    }
}
