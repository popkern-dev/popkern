//! Transactional completion of missing BRep edge curves.

use crate::topology::{BRep, EdgeHandle};

use super::brep_edge_curve_completion_helpers::{
    curve_from_surface_trim, curve_from_vertices, surface_trim_sources, validate_candidate,
    validated_anchored_surface_trim_curve, validated_vertex_curve_for_surfaced_edge,
};

/// Options controlling missing edge-curve completion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepEdgeCurveCompletionOptions {
    /// Maximum allowed endpoint and sampled geometry disagreement.
    pub tolerance: f64,
    /// Number of samples used when lifting trims and validating candidates.
    pub samples: usize,
    /// Whether non-surfaced edges may be completed as vertex-to-vertex segments.
    pub allow_vertex_segments_without_surface: bool,
    /// Whether any declined edge rolls back the whole completion transaction.
    pub require_all_completed: bool,
    /// Whether completed output must pass full BRep validation.
    pub validate_output: bool,
}

/// Successful transactional completion result.
#[derive(Debug, Clone)]
pub struct BrepEdgeCurveCompletionResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Completion report.
    pub report: BrepEdgeCurveCompletionReport,
}

/// Aggregate missing edge-curve completion report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeCurveCompletionReport {
    /// Number of edges inspected because they had no 3D curve.
    pub inspected_missing_edges: usize,
    /// Number of missing edge curves completed.
    pub completed_edges: usize,
    /// Number of missing edge curves declined.
    pub declined_edges: usize,
    /// Maximum endpoint gap observed while validating accepted candidates.
    pub max_endpoint_gap: f64,
    /// Maximum sampled pcurve/surface disagreement for accepted candidates.
    pub max_sampled_disagreement: f64,
    /// Per-edge completion decisions.
    pub decisions: Vec<BrepEdgeCurveCompletionDecision>,
}

impl BrepEdgeCurveCompletionReport {
    /// True when every inspected missing curve was completed.
    pub fn is_complete(&self) -> bool {
        self.declined_edges == 0
    }
}

/// Per-edge completion decision.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeCurveCompletionDecision {
    /// Edge being considered.
    pub edge: EdgeHandle,
    /// Decision outcome.
    pub outcome: BrepEdgeCurveCompletionOutcome,
    /// Number of trim/surface sources used to validate the candidate.
    pub source_trims: usize,
    /// Maximum endpoint gap observed for this edge.
    pub endpoint_gap: f64,
    /// Maximum sampled trim/surface disagreement for this edge.
    pub sampled_disagreement: f64,
}

/// Completion decision outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeCurveCompletionOutcome {
    /// A curve was completed from one or more pcurve/surface sources.
    CompletedFromSurfaceTrims,
    /// A non-surfaced edge was completed as a vertex-to-vertex segment.
    CompletedFromVertices,
    /// Completion was declined.
    Declined(BrepEdgeCurveCompletionDeclineReason),
}

/// Machine-readable decline reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepEdgeCurveCompletionDeclineReason {
    /// Edge or canonical half-edge topology was missing.
    MissingTopology,
    /// Edge endpoints were coincident or too short for a bounded curve.
    DegenerateEdge,
    /// A surfaced use exists but no pcurve/surface lift was available.
    MissingSurfaceTrimSource,
    /// Endpoint validation exceeded tolerance.
    EndpointMismatch,
    /// Candidate curve disagreed with at least one owning trim/surface.
    SampledDisagreement,
    /// Curve fitting failed for sampled trim/surface points.
    FittingFailed,
    /// Vertex-only fallback was disabled.
    VertexSegmentFallbackDisabled,
}

/// Transactional completion failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgeCurveCompletionFailure {
    /// Machine-readable failure reason.
    pub reason: BrepEdgeCurveCompletionFailureReason,
    /// Report accumulated before rollback.
    pub report: BrepEdgeCurveCompletionReport,
}

/// Reason the completion transaction failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepEdgeCurveCompletionFailureReason {
    /// At least one missing edge curve could not be completed.
    UnresolvedMissingEdgeCurves,
    /// Completed clone failed validation.
    ValidationFailed(String),
}

/// Analyze and complete missing 3D edge curves on a cloned BRep.
pub fn complete_missing_edge_curves_transactionally(
    brep: &BRep,
    options: BrepEdgeCurveCompletionOptions,
) -> Result<BrepEdgeCurveCompletionResult, Box<BrepEdgeCurveCompletionFailure>> {
    let mut repaired = brep.clone();
    let mut report = BrepEdgeCurveCompletionReport {
        inspected_missing_edges: 0,
        completed_edges: 0,
        declined_edges: 0,
        max_endpoint_gap: 0.0,
        max_sampled_disagreement: 0.0,
        decisions: Vec::new(),
    };

    let missing: Vec<_> = brep
        .edges
        .iter()
        .filter_map(|(edge, data)| data.curve.is_none().then_some(edge))
        .collect();

    for edge in missing {
        report.inspected_missing_edges += 1;
        let decision = complete_edge_curve(brep, &mut repaired, edge, options);
        if matches!(
            decision.outcome,
            BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims
                | BrepEdgeCurveCompletionOutcome::CompletedFromVertices
        ) {
            report.completed_edges += 1;
            report.max_endpoint_gap = report.max_endpoint_gap.max(decision.endpoint_gap);
            report.max_sampled_disagreement = report
                .max_sampled_disagreement
                .max(decision.sampled_disagreement);
        } else {
            report.declined_edges += 1;
        }
        report.decisions.push(decision);
    }

    if options.require_all_completed && !report.is_complete() {
        return Err(Box::new(BrepEdgeCurveCompletionFailure {
            reason: BrepEdgeCurveCompletionFailureReason::UnresolvedMissingEdgeCurves,
            report,
        }));
    }
    if options.validate_output
        && let Err(err) = repaired.validate()
    {
        return Err(Box::new(BrepEdgeCurveCompletionFailure {
            reason: BrepEdgeCurveCompletionFailureReason::ValidationFailed(err.to_string()),
            report,
        }));
    }
    Ok(BrepEdgeCurveCompletionResult {
        brep: repaired,
        report,
    })
}

fn complete_edge_curve(
    original: &BRep,
    repaired: &mut BRep,
    edge: EdgeHandle,
    options: BrepEdgeCurveCompletionOptions,
) -> BrepEdgeCurveCompletionDecision {
    let sources = surface_trim_sources(original, edge);
    let source_trims = sources.len();
    let candidate = if let Some(source) = sources.first() {
        curve_from_surface_trim(original, edge, source, options)
    } else {
        curve_from_vertices(original, edge, options)
    };
    let candidate = match candidate {
        Ok(curve) => curve,
        Err(reason) => return decision(edge, reason, source_trims, 0.0, 0.0),
    };
    let Some((endpoint_gap, sampled_disagreement)) =
        validate_candidate(original, edge, &candidate, &sources, options)
    else {
        return decision(
            edge,
            BrepEdgeCurveCompletionDeclineReason::MissingTopology,
            source_trims,
            0.0,
            0.0,
        );
    };
    if endpoint_gap > options.tolerance {
        if let Some((anchored_candidate, anchored_deviation)) =
            validated_anchored_surface_trim_curve(original, edge, &sources, options)
        {
            repaired.edges[edge].curve = Some(anchored_candidate);
            return BrepEdgeCurveCompletionDecision {
                edge,
                outcome: BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims,
                source_trims,
                endpoint_gap: 0.0,
                sampled_disagreement: anchored_deviation,
            };
        }
        if let Some(vertex_candidate) =
            validated_vertex_curve_for_surfaced_edge(original, edge, &sources, options)
        {
            repaired.edges[edge].curve = Some(vertex_candidate);
            return BrepEdgeCurveCompletionDecision {
                edge,
                outcome: BrepEdgeCurveCompletionOutcome::CompletedFromVertices,
                source_trims,
                endpoint_gap: 0.0,
                sampled_disagreement: 0.0,
            };
        }
        return decision(
            edge,
            BrepEdgeCurveCompletionDeclineReason::EndpointMismatch,
            source_trims,
            endpoint_gap,
            sampled_disagreement,
        );
    }
    if sampled_disagreement > options.tolerance {
        if let Some((anchored_candidate, anchored_deviation)) =
            validated_anchored_surface_trim_curve(original, edge, &sources, options)
        {
            repaired.edges[edge].curve = Some(anchored_candidate);
            return BrepEdgeCurveCompletionDecision {
                edge,
                outcome: BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims,
                source_trims,
                endpoint_gap: 0.0,
                sampled_disagreement: anchored_deviation,
            };
        }
        if let Some(vertex_candidate) =
            validated_vertex_curve_for_surfaced_edge(original, edge, &sources, options)
        {
            repaired.edges[edge].curve = Some(vertex_candidate);
            return BrepEdgeCurveCompletionDecision {
                edge,
                outcome: BrepEdgeCurveCompletionOutcome::CompletedFromVertices,
                source_trims,
                endpoint_gap: 0.0,
                sampled_disagreement: 0.0,
            };
        }
        return decision(
            edge,
            BrepEdgeCurveCompletionDeclineReason::SampledDisagreement,
            source_trims,
            endpoint_gap,
            sampled_disagreement,
        );
    }
    repaired.edges[edge].curve = Some(candidate);
    BrepEdgeCurveCompletionDecision {
        edge,
        outcome: if source_trims > 0 {
            BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims
        } else {
            BrepEdgeCurveCompletionOutcome::CompletedFromVertices
        },
        source_trims,
        endpoint_gap,
        sampled_disagreement,
    }
}

fn decision(
    edge: EdgeHandle,
    reason: BrepEdgeCurveCompletionDeclineReason,
    source_trims: usize,
    endpoint_gap: f64,
    sampled_disagreement: f64,
) -> BrepEdgeCurveCompletionDecision {
    BrepEdgeCurveCompletionDecision {
        edge,
        outcome: BrepEdgeCurveCompletionOutcome::Declined(reason),
        source_trims,
        endpoint_gap,
        sampled_disagreement,
    }
}

#[cfg(test)]
#[path = "brep_edge_curve_completion_tests.rs"]
mod tests;
