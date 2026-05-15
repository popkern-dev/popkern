//! Non-manifold policy checks and explanations for topology repair.

use std::collections::HashSet;

use crate::topology::{BRep, BrepSewingEdgePair, EdgeHandle, HalfEdgeHandle};

/// Policy flags for edge-pair topology repair decisions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepTopologyRepairPolicy {
    /// Whether a repair may touch existing non-manifold radial topology.
    pub allow_nonmanifold_edges: bool,
}

/// Read-only context for assessing one candidate edge pair.
#[derive(Debug, Clone, Copy)]
pub struct BrepEdgePairPolicyContext<'a> {
    /// Edges already consumed by earlier accepted pairs in the same transaction.
    pub used_edges: &'a HashSet<EdgeHandle>,
}

/// Result of a topology policy check for one edge pair.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdgePairPolicyDecision {
    /// Candidate pair that was assessed.
    pub pair: BrepSewingEdgePair,
    /// Machine-readable outcome.
    pub outcome: BrepEdgePairPolicyOutcome,
    /// Topology state for the first edge.
    pub first_edge: BrepEdgePolicyState,
    /// Topology state for the second edge.
    pub second_edge: BrepEdgePolicyState,
    /// User-facing explanation suitable for API or GUI diagnostics.
    pub explanation: String,
}

impl BrepEdgePairPolicyDecision {
    /// Returns true when the candidate may proceed to geometric checks.
    pub fn is_accepted(&self) -> bool {
        matches!(self.outcome, BrepEdgePairPolicyOutcome::Accepted)
    }
}

/// Machine-readable topology policy outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepEdgePairPolicyOutcome {
    /// Candidate satisfies topology policy and may proceed.
    Accepted,
    /// Candidate violates topology policy.
    Declined(BrepEdgePairPolicyDeclineReason),
}

/// Reason a topology policy declined a candidate pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepEdgePairPolicyDeclineReason {
    /// Both sides reference the same edge.
    SameEdge,
    /// At least one edge handle does not resolve to complete topology.
    MissingEdge,
    /// One of the edges was consumed earlier in the same transaction.
    EdgeAlreadyUsed,
    /// At least one edge is not a boundary edge.
    NotBoundaryEdge,
    /// At least one edge is already joined to a twin half-edge.
    AlreadyJoined,
    /// Existing non-manifold radial topology conflicts with the policy.
    NonManifoldRadialConflict,
    /// Edge state is internally inconsistent.
    AmbiguousTopology,
}

/// Read-only topology state for one candidate edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepEdgePolicyState {
    /// Edge handle being inspected.
    pub edge: EdgeHandle,
    /// Canonical half-edge when the edge exists.
    pub half_edge: Option<HalfEdgeHandle>,
    /// Whether the edge is currently marked as a boundary edge.
    pub is_boundary: bool,
    /// Whether the canonical half-edge has an existing valid twin.
    pub has_existing_twin: bool,
    /// Whether the canonical half-edge participates in radial topology.
    pub has_radial_cycle: bool,
    /// Whether the edge points at incomplete or contradictory topology.
    pub ambiguous: bool,
}

/// Apply generic topology policy to one edge-pair repair candidate.
pub fn assess_edge_pair_topology_policy(
    brep: &BRep,
    pair: &BrepSewingEdgePair,
    policy: BrepTopologyRepairPolicy,
    context: BrepEdgePairPolicyContext<'_>,
) -> BrepEdgePairPolicyDecision {
    let first_edge = edge_policy_state(brep, pair.first);
    let second_edge = edge_policy_state(brep, pair.second);
    let outcome = policy_outcome(pair, policy, context, first_edge, second_edge);
    let explanation = explain_policy_outcome(&outcome, pair, first_edge, second_edge);
    BrepEdgePairPolicyDecision {
        pair: pair.clone(),
        outcome,
        first_edge,
        second_edge,
        explanation,
    }
}

fn policy_outcome(
    pair: &BrepSewingEdgePair,
    policy: BrepTopologyRepairPolicy,
    context: BrepEdgePairPolicyContext<'_>,
    first: BrepEdgePolicyState,
    second: BrepEdgePolicyState,
) -> BrepEdgePairPolicyOutcome {
    use BrepEdgePairPolicyDeclineReason as Reason;
    if pair.first == pair.second {
        return BrepEdgePairPolicyOutcome::Declined(Reason::SameEdge);
    }
    if context.used_edges.contains(&pair.first) || context.used_edges.contains(&pair.second) {
        return BrepEdgePairPolicyOutcome::Declined(Reason::EdgeAlreadyUsed);
    }
    if first.half_edge.is_none() || second.half_edge.is_none() {
        return BrepEdgePairPolicyOutcome::Declined(Reason::MissingEdge);
    }
    if first.ambiguous || second.ambiguous {
        return BrepEdgePairPolicyOutcome::Declined(Reason::AmbiguousTopology);
    }
    if !first.is_boundary || !second.is_boundary {
        return BrepEdgePairPolicyOutcome::Declined(Reason::NotBoundaryEdge);
    }
    if !policy.allow_nonmanifold_edges && (first.has_existing_twin || second.has_existing_twin) {
        return BrepEdgePairPolicyOutcome::Declined(Reason::AlreadyJoined);
    }
    if !policy.allow_nonmanifold_edges && (first.has_radial_cycle || second.has_radial_cycle) {
        return BrepEdgePairPolicyOutcome::Declined(Reason::NonManifoldRadialConflict);
    }
    BrepEdgePairPolicyOutcome::Accepted
}

fn edge_policy_state(brep: &BRep, edge: EdgeHandle) -> BrepEdgePolicyState {
    let edge_data = brep.edges.get(edge);
    let half_edge = edge_data.map(|data| data.half_edge);
    let half_edge_data = half_edge.and_then(|handle| brep.half_edges.get(handle));
    let has_existing_twin =
        half_edge_data.is_some_and(|data| brep.half_edges.contains_key(data.twin));
    let has_radial_cycle =
        half_edge.is_some_and(|handle| brep.nonmanifold_radial.contains_key(&handle));
    let ambiguous = edge_data.is_some()
        && (half_edge_data.is_none() || half_edge_data.is_some_and(|data| data.edge != edge));
    BrepEdgePolicyState {
        edge,
        half_edge,
        is_boundary: brep.boundary_edges.contains(&edge),
        has_existing_twin,
        has_radial_cycle,
        ambiguous,
    }
}

fn explain_policy_outcome(
    outcome: &BrepEdgePairPolicyOutcome,
    pair: &BrepSewingEdgePair,
    first: BrepEdgePolicyState,
    second: BrepEdgePolicyState,
) -> String {
    use BrepEdgePairPolicyDeclineReason as Reason;
    match outcome {
        BrepEdgePairPolicyOutcome::Accepted => {
            "Candidate edge pair satisfies topology policy and may proceed to geometric checks."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::SameEdge) => {
            "Candidate references the same edge on both sides; sewing an edge to itself would corrupt topology."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::MissingEdge) => {
            "Candidate references an edge whose canonical half-edge is missing, so repair cannot preserve topology ownership."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::EdgeAlreadyUsed) => {
            "Candidate reuses an edge already accepted earlier in this transaction."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::NotBoundaryEdge) => {
            format!(
                "Candidate requires boundary edges, but boundary state is first={} second={}.",
                first.is_boundary, second.is_boundary
            )
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::AlreadyJoined) => {
            "Candidate touches an edge that is already joined to another face-use; policy forbids creating or extending non-manifold topology."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::NonManifoldRadialConflict) => {
            "Candidate touches existing radial non-manifold topology; policy requires an explicit non-manifold repair mode before mutation."
                .to_string()
        }
        BrepEdgePairPolicyOutcome::Declined(Reason::AmbiguousTopology) => {
            format!(
                "Candidate has inconsistent edge ownership for pair {:?}/{:?}; repair declined before mutation.",
                pair.first, pair.second
            )
        }
    }
}
