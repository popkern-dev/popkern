//! Transactional BRep sewing for caller-supplied clean edge pairs.

use std::collections::HashSet;

use crate::geometry::{Point3, SurfaceEval};
use crate::topology::brep_sewing_facts::{sewing_pair_fact_ids, sewing_pair_knowledge};
use crate::topology::{
    BRep, BrepEdgePairPolicyContext, BrepEdgePairPolicyDecision, BrepEdgePairPolicyDeclineReason,
    BrepEdgePairPolicyOutcome, BrepRepairFactId, BrepTopologyRepairPolicy, EdgeHandle,
    HalfEdgeHandle, VertexHandle, assess_edge_pair_topology_policy,
};

#[derive(Debug, Clone, PartialEq, Eq)]
/// Caller-provided labels for diagnosing where a sewing candidate came from.
pub struct BrepSewingProvenance {
    /// Label for the first edge's source.
    pub first: String,
    /// Label for the second edge's source.
    pub second: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// A caller-supplied pair of boundary edges to try to sew.
pub struct BrepSewingEdgePair {
    /// First boundary edge.
    pub first: EdgeHandle,
    /// Second boundary edge, expected to run opposite the first edge.
    pub second: EdgeHandle,
    /// Source labels retained in diagnostics.
    pub provenance: BrepSewingProvenance,
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Local policy for exact/clean boundary-edge sewing.
pub struct BrepSewingOptions {
    /// Maximum accepted endpoint and sampled surface disagreement.
    pub tolerance: f64,
    /// Number of samples used for pcurve/surface agreement checks.
    pub samples: usize,
    /// Whether both edge pcurves must sample to the same 3D curve.
    pub require_surface_agreement: bool,
    /// Whether to run `BRep::validate` before returning the sewn clone.
    pub validate_output: bool,
    /// Whether pairs with existing radial/twin topology may be accepted.
    pub allow_nonmanifold_edges: bool,
}

#[derive(Debug, Clone)]
/// Successful transactional sewing result.
pub struct BrepSewingResult {
    /// Sewn clone of the input BRep.
    pub brep: BRep,
    /// Structured acceptance and decline diagnostics.
    pub report: BrepSewingReport,
}

#[derive(Debug, Clone, PartialEq)]
/// Aggregate diagnostics from a sewing attempt.
pub struct BrepSewingReport {
    /// Number of caller-supplied candidate pairs.
    pub candidate_pairs: usize,
    /// Number of pairs accepted and welded.
    pub accepted_pairs: usize,
    /// Pair-level declines with explicit reasons.
    pub declined_pairs: Vec<BrepSewingDecline>,
    /// Accepted edge pairs that were sewn.
    pub accepted_pair_edges: Vec<(EdgeHandle, EdgeHandle)>,
    /// Boundary edge count after accepted pairs are sewn.
    pub remaining_boundary_edges: usize,
    /// Maximum reversed endpoint gap observed while assessing pairs.
    pub max_endpoint_gap: f64,
    /// Maximum sampled pcurve/surface disagreement observed.
    pub max_surface_gap: f64,
    /// Shared repair knowledge facts referenced by this sewing attempt.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

#[derive(Debug, Clone, PartialEq)]
/// Pair-level sewing decline.
pub struct BrepSewingDecline {
    /// Candidate pair that declined.
    pub pair: BrepSewingEdgePair,
    /// Reason the pair was not sewn.
    pub reason: BrepSewingDeclineReason,
    /// Detailed topology policy explanation when the decline came from policy.
    pub policy: Option<BrepEdgePairPolicyDecision>,
    /// Shared repair knowledge facts used to decline this pair.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

#[derive(Debug, Clone, PartialEq)]
/// Reason a candidate pair was not sewn.
pub enum BrepSewingDeclineReason {
    /// A pair referenced the same edge twice.
    SameEdge,
    /// At least one edge handle was missing or had incomplete topology.
    MissingEdge,
    /// One of the edges was already consumed by an earlier accepted pair.
    EdgeAlreadyUsed,
    /// At least one edge was not marked as a boundary edge.
    NotBoundaryEdge,
    /// Pair would touch existing twin/radial non-manifold topology.
    NonManifoldConflict,
    /// Pair touched an already joined edge while non-manifold repair was disabled.
    AlreadyJoined,
    /// Pair has internally inconsistent topology.
    AmbiguousTopology,
    /// Reversed endpoints were farther apart than the local tolerance.
    EndpointGapExceeded {
        /// Observed reversed endpoint gap.
        gap: f64,
        /// Local sewing tolerance.
        tolerance: f64,
    },
    /// Surface agreement was requested but pcurve/surface data was missing.
    MissingSurfaceSampleData,
    /// Sampled pcurve/surface points disagreed beyond tolerance.
    SurfaceGapExceeded {
        /// Observed sampled surface gap.
        gap: f64,
        /// Local sewing tolerance.
        tolerance: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
/// Transaction-level failure after pair assessment and mutation.
pub struct BrepSewingFailure {
    /// Failure reason.
    pub reason: BrepSewingFailureReason,
    /// Report accumulated before the transaction failed.
    pub report: BrepSewingReport,
}

#[derive(Debug, Clone, PartialEq)]
/// Reason a sewing transaction failed after mutation.
pub enum BrepSewingFailureReason {
    /// The cloned output failed `BRep::validate`.
    ValidationFailed(String),
}

/// Sew caller-supplied exact/clean boundary-edge pairs into a cloned BRep.
pub fn sew_clean_boundary_edge_pairs(
    brep: &BRep,
    pairs: &[BrepSewingEdgePair],
    options: BrepSewingOptions,
) -> Result<BrepSewingResult, BrepSewingFailure> {
    let knowledge = sewing_pair_knowledge(brep, pairs, options.tolerance, options.samples);
    if !knowledge.validate_fresh(brep) {
        return Err(BrepSewingFailure {
            reason: BrepSewingFailureReason::ValidationFailed(
                "stale sewing repair knowledge graph".to_string(),
            ),
            report: empty_report(brep, pairs),
        });
    }
    let mut candidate = brep.clone();
    let mut used = HashSet::new();
    let mut report = BrepSewingReport {
        candidate_pairs: pairs.len(),
        accepted_pairs: 0,
        declined_pairs: Vec::new(),
        accepted_pair_edges: Vec::new(),
        remaining_boundary_edges: candidate.boundary_edges.len(),
        max_endpoint_gap: 0.0,
        max_surface_gap: 0.0,
        shared_fact_ids: Vec::new(),
    };

    for pair in pairs {
        let shared_fact_ids = sewing_pair_fact_ids(&knowledge, brep, pair);
        match assess_pair(&candidate, pair, options, &used) {
            Ok(assessment) => {
                report.max_endpoint_gap = report.max_endpoint_gap.max(assessment.endpoint_gap);
                report.max_surface_gap = report.max_surface_gap.max(assessment.surface_gap);
                report.shared_fact_ids.extend(shared_fact_ids);
                weld_boundary_pair(&mut candidate, pair.first, pair.second);
                used.insert(pair.first);
                used.insert(pair.second);
                report.accepted_pairs += 1;
                report.accepted_pair_edges.push((pair.first, pair.second));
            }
            Err(assessment) => {
                let reason = assessment.reason;
                if let BrepSewingDeclineReason::EndpointGapExceeded { gap, .. } = reason {
                    report.max_endpoint_gap = report.max_endpoint_gap.max(gap);
                }
                if let BrepSewingDeclineReason::SurfaceGapExceeded { gap, .. } = reason {
                    report.max_surface_gap = report.max_surface_gap.max(gap);
                }
                report.declined_pairs.push(BrepSewingDecline {
                    pair: pair.clone(),
                    reason,
                    policy: assessment.policy,
                    shared_fact_ids: shared_fact_ids.clone(),
                });
                report.shared_fact_ids.extend(shared_fact_ids);
            }
        }
    }

    repair_vertex_outgoing(&mut candidate);
    report.remaining_boundary_edges = candidate.boundary_edges.len();
    report.shared_fact_ids.sort();
    report.shared_fact_ids.dedup();
    if options.validate_output
        && let Err(err) = candidate.validate()
    {
        return Err(BrepSewingFailure {
            reason: BrepSewingFailureReason::ValidationFailed(err.to_string()),
            report,
        });
    }
    Ok(BrepSewingResult {
        brep: candidate,
        report,
    })
}

fn empty_report(brep: &BRep, pairs: &[BrepSewingEdgePair]) -> BrepSewingReport {
    BrepSewingReport {
        candidate_pairs: pairs.len(),
        accepted_pairs: 0,
        declined_pairs: Vec::new(),
        accepted_pair_edges: Vec::new(),
        remaining_boundary_edges: brep.boundary_edges.len(),
        max_endpoint_gap: 0.0,
        max_surface_gap: 0.0,
        shared_fact_ids: Vec::new(),
    }
}

#[derive(Debug, Clone, Copy)]
struct PairAssessment {
    endpoint_gap: f64,
    surface_gap: f64,
}

#[derive(Debug, Clone)]
struct PairDecline {
    reason: BrepSewingDeclineReason,
    policy: Option<BrepEdgePairPolicyDecision>,
}

fn assess_pair(
    brep: &BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSewingOptions,
    used: &HashSet<EdgeHandle>,
) -> Result<PairAssessment, Box<PairDecline>> {
    let policy_decision = assess_edge_pair_topology_policy(
        brep,
        pair,
        BrepTopologyRepairPolicy {
            allow_nonmanifold_edges: options.allow_nonmanifold_edges,
        },
        BrepEdgePairPolicyContext { used_edges: used },
    );
    if !policy_decision.is_accepted() {
        return Err(Box::new(PairDecline {
            reason: map_policy_decline(&policy_decision),
            policy: Some(policy_decision),
        }));
    }
    let Some(first) = edge_segment(brep, pair.first) else {
        return Err(simple_decline(BrepSewingDeclineReason::MissingEdge));
    };
    let Some(second) = edge_segment(brep, pair.second) else {
        return Err(simple_decline(BrepSewingDeclineReason::MissingEdge));
    };

    let endpoint_gap = reversed_endpoint_gap(first, second);
    if endpoint_gap > options.tolerance {
        return Err(simple_decline(
            BrepSewingDeclineReason::EndpointGapExceeded {
                gap: endpoint_gap,
                tolerance: options.tolerance,
            },
        ));
    }
    let surface_gap = if options.require_surface_agreement {
        reversed_sample_gap(brep, pair.first, pair.second, options.samples)
            .ok_or_else(|| simple_decline(BrepSewingDeclineReason::MissingSurfaceSampleData))?
    } else {
        0.0
    };
    if surface_gap > options.tolerance {
        return Err(simple_decline(
            BrepSewingDeclineReason::SurfaceGapExceeded {
                gap: surface_gap,
                tolerance: options.tolerance,
            },
        ));
    }
    Ok(PairAssessment {
        endpoint_gap,
        surface_gap,
    })
}

fn simple_decline(reason: BrepSewingDeclineReason) -> Box<PairDecline> {
    Box::new(PairDecline {
        reason,
        policy: None,
    })
}

fn map_policy_decline(decision: &BrepEdgePairPolicyDecision) -> BrepSewingDeclineReason {
    let BrepEdgePairPolicyOutcome::Declined(reason) = &decision.outcome else {
        return BrepSewingDeclineReason::AmbiguousTopology;
    };
    match reason {
        BrepEdgePairPolicyDeclineReason::SameEdge => BrepSewingDeclineReason::SameEdge,
        BrepEdgePairPolicyDeclineReason::MissingEdge => BrepSewingDeclineReason::MissingEdge,
        BrepEdgePairPolicyDeclineReason::EdgeAlreadyUsed => {
            BrepSewingDeclineReason::EdgeAlreadyUsed
        }
        BrepEdgePairPolicyDeclineReason::NotBoundaryEdge => {
            BrepSewingDeclineReason::NotBoundaryEdge
        }
        BrepEdgePairPolicyDeclineReason::AlreadyJoined => BrepSewingDeclineReason::AlreadyJoined,
        BrepEdgePairPolicyDeclineReason::NonManifoldRadialConflict => {
            BrepSewingDeclineReason::NonManifoldConflict
        }
        BrepEdgePairPolicyDeclineReason::AmbiguousTopology => {
            BrepSewingDeclineReason::AmbiguousTopology
        }
    }
}

/// Return the 3D start and end positions of an edge's canonical half-edge.
pub fn edge_segment(brep: &BRep, edge: EdgeHandle) -> Option<(Point3, Point3)> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    let start = brep.vertices.get(half_edge.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    Some((start, end))
}

fn reversed_sample_gap(
    brep: &BRep,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
    samples: usize,
) -> Option<f64> {
    let a = sample_edge_on_face(brep, edge_a, samples)?;
    let b = sample_edge_on_face(brep, edge_b, samples)?;
    Some(
        a.iter()
            .zip(b.iter().rev())
            .map(|(pa, pb)| pa.distance_to(pb))
            .fold(0.0_f64, f64::max),
    )
}

fn sample_edge_on_face(brep: &BRep, edge: EdgeHandle, samples: usize) -> Option<Vec<Point3>> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let pcurve = half_edge.pcurve.as_ref()?;
    let loop_ = brep.loops.get(half_edge.loop_)?;
    let face = brep.faces.get(loop_.face)?;
    let surface = face.surface.as_ref()?;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for i in 0..samples {
        let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * i as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(t);
        points.push(surface.geometry.point_at(uv.x, uv.y));
    }
    Some(points)
}

fn reversed_endpoint_gap(a: (Point3, Point3), b: (Point3, Point3)) -> f64 {
    a.0.distance_to(&b.1).max(a.1.distance_to(&b.0))
}

fn weld_boundary_pair(brep: &mut BRep, edge_a: EdgeHandle, edge_b: EdgeHandle) {
    let he_a = brep.edges[edge_a].half_edge;
    let he_b = brep.edges[edge_b].half_edge;
    let a_start = brep.half_edges[he_a].origin;
    let a_end = brep.half_edges[brep.half_edges[he_a].next].origin;
    let b_start = brep.half_edges[he_b].origin;
    let b_end = brep.half_edges[brep.half_edges[he_b].next].origin;

    replace_vertex(brep, b_start, a_end);
    replace_vertex(brep, b_end, a_start);
    brep.half_edges[he_a].twin = he_b;
    brep.half_edges[he_b].twin = he_a;
    brep.half_edges[he_b].edge = edge_a;
    brep.boundary_edges.remove(&edge_a);
    brep.boundary_edges.remove(&edge_b);
    brep.edges.remove(edge_b);
}

fn replace_vertex(brep: &mut BRep, from: VertexHandle, to: VertexHandle) {
    if from == to || !brep.vertices.contains_key(from) || !brep.vertices.contains_key(to) {
        return;
    }
    for (_, half_edge) in &mut brep.half_edges {
        if half_edge.origin == from {
            half_edge.origin = to;
        }
    }
    if brep.vertices[to].outgoing == HalfEdgeHandle::default() {
        brep.vertices[to].outgoing = brep.vertices[from].outgoing;
    }
    brep.vertices.remove(from);
}

fn repair_vertex_outgoing(brep: &mut BRep) {
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for vertex in vertices {
        brep.vertices[vertex].outgoing = HalfEdgeHandle::default();
    }
    for (half_edge, data) in &brep.half_edges {
        if brep.vertices[data.origin].outgoing == HalfEdgeHandle::default() {
            brep.vertices[data.origin].outgoing = half_edge;
        }
    }
}

#[cfg(test)]
#[path = "brep_sewing_tests.rs"]
mod tests;
