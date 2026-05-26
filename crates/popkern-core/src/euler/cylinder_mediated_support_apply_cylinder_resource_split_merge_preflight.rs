//! Preflight gates for split cylinder-resource merge candidates.

use std::collections::BTreeSet;

use super::{SplitResourceEndpointShape, SplitResourceMergeOperand, split_resource_endpoint_shape};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SplitResourceMergePreflight {
    status: SplitResourceMergePreflightStatus,
    operands: usize,
    cylinder_runs: usize,
    boundary_ready_runs: usize,
    live_twin_blocked_runs: usize,
    not_boundary_blocked_runs: usize,
    mixed_endpoint_operands: usize,
    disconnected_endpoint_operands: usize,
    outside_endpoint_matches: usize,
    connector_authoring_links: usize,
    connector_authoring_operands: usize,
}

impl SplitResourceMergePreflight {
    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "split_resource_merge_preflight=status={} operands={} cylinder_runs={} boundary_ready_runs={} live_twin_blocked_runs={} not_boundary_blocked_runs={} mixed_endpoint_operands={} disconnected_endpoint_operands={} outside_endpoint_matches={} connector_authoring_links={} connector_authoring_operands={}",
            self.status.diagnostic_label(),
            self.operands,
            self.cylinder_runs,
            self.boundary_ready_runs,
            self.live_twin_blocked_runs,
            self.not_boundary_blocked_runs,
            self.mixed_endpoint_operands,
            self.disconnected_endpoint_operands,
            self.outside_endpoint_matches,
            self.connector_authoring_links,
            self.connector_authoring_operands
        )
    }

    pub(super) fn requires_connector_authoring(self) -> bool {
        self.disconnected_endpoint_operands > 0
            && self.live_twin_blocked_runs == 0
            && self.not_boundary_blocked_runs == 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SplitResourceMergePreflightStatus {
    NotApplicable,
    Ready,
    BlockedByDisconnectedEndpointTopology,
    BlockedByMixedEndpointTopology,
    BlockedByLiveTwinOwner,
    BlockedByNonBoundarySide,
}

impl SplitResourceMergePreflightStatus {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::NotApplicable => "not-applicable",
            Self::Ready => "ready",
            Self::BlockedByDisconnectedEndpointTopology => {
                "blocked-by-disconnected-endpoint-topology"
            }
            Self::BlockedByMixedEndpointTopology => "blocked-by-mixed-endpoint-topology",
            Self::BlockedByLiveTwinOwner => "blocked-by-live-twin-owner",
            Self::BlockedByNonBoundarySide => "blocked-by-non-boundary-side",
        }
    }
}

pub(super) fn split_resource_merge_preflight(
    operands: &[SplitResourceMergeOperand],
) -> SplitResourceMergePreflight {
    let mut preflight = SplitResourceMergePreflight {
        status: SplitResourceMergePreflightStatus::NotApplicable,
        operands: operands.len(),
        cylinder_runs: 0,
        boundary_ready_runs: 0,
        live_twin_blocked_runs: 0,
        not_boundary_blocked_runs: 0,
        mixed_endpoint_operands: 0,
        disconnected_endpoint_operands: 0,
        outside_endpoint_matches: 0,
        connector_authoring_links: 0,
        connector_authoring_operands: 0,
    };
    if operands.is_empty() {
        return preflight;
    }
    for operand in operands {
        if split_resource_endpoint_shape(operand) == SplitResourceEndpointShape::Mixed {
            preflight.mixed_endpoint_operands += 1;
        }
        let endpoint_topology = split_resource_endpoint_topology(operand);
        if endpoint_topology.disconnected {
            preflight.disconnected_endpoint_operands += 1;
        }
        preflight.outside_endpoint_matches += endpoint_topology.outside_endpoint_matches;
        let connector_plan = split_resource_connector_authoring_plan(operand, endpoint_topology);
        preflight.connector_authoring_links += connector_plan.links_required;
        if connector_plan.links_required > 0 {
            preflight.connector_authoring_operands += 1;
        }
        for run in &operand.cylinder_runs {
            preflight.cylinder_runs += 1;
            if run.live_twin_half_edges > 0 {
                preflight.live_twin_blocked_runs += 1;
            }
            if run.not_boundary_half_edges > 0 {
                preflight.not_boundary_blocked_runs += 1;
            }
            if run.boundary_ready_half_edges == run.half_edges.len() {
                preflight.boundary_ready_runs += 1;
            }
        }
    }
    preflight.status = if preflight.not_boundary_blocked_runs > 0 {
        SplitResourceMergePreflightStatus::BlockedByNonBoundarySide
    } else if preflight.live_twin_blocked_runs > 0 {
        SplitResourceMergePreflightStatus::BlockedByLiveTwinOwner
    } else if preflight.disconnected_endpoint_operands > 0 {
        SplitResourceMergePreflightStatus::BlockedByDisconnectedEndpointTopology
    } else if preflight.mixed_endpoint_operands > 0 {
        SplitResourceMergePreflightStatus::BlockedByMixedEndpointTopology
    } else {
        SplitResourceMergePreflightStatus::Ready
    };
    preflight
}

pub(super) fn split_resource_connector_authoring_summary(
    operands: &[SplitResourceMergeOperand],
) -> String {
    let plans = operands
        .iter()
        .filter_map(split_resource_connector_authoring_plan_label)
        .collect::<Vec<_>>();
    let edges = operands
        .iter()
        .map(|operand| {
            split_resource_connector_authoring_plan(
                operand,
                split_resource_endpoint_topology(operand),
            )
            .links_required
        })
        .sum::<usize>();
    format!(
        "split_resource_connector_authoring_plan=operands={} links={} ordering={} total_connector_length={:.6e} max_connector_length={:.6e} connectors=[{}]",
        plans.len(),
        edges,
        connector_ordering_validation(operands).diagnostic_label(),
        connector_total_length(operands),
        connector_max_length(operands),
        plans.join(";")
    )
}

pub(super) fn split_resource_connector_transaction_summary(
    operands: &[SplitResourceMergeOperand],
) -> String {
    let links = connector_authoring_links(operands);
    let multi_run_faces = operands
        .iter()
        .filter(|operand| operand.cylinder_runs.len() > 1)
        .count();
    let status = if links > 0 {
        "ready-for-topology-mutation"
    } else {
        "not-applicable"
    };
    format!(
        "split_resource_connector_transaction=status={} primitive=multi-run-bridge-face-authoring multi_run_faces={} resource_links={} existing_boundary_spans={} expected_delta=computed-at-mutation-boundary rollback_boundary=candidate-brep-before-commit postconditions=half-edge-pairing+cycles+manifold+delta",
        status,
        multi_run_faces,
        links,
        existing_boundary_connector_spans(operands)
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConnectorOrderingValidation {
    NotRequired,
    PartitionOrderValidated,
    GeometricOrderCandidate,
}

impl ConnectorOrderingValidation {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::NotRequired => "not-required",
            Self::PartitionOrderValidated => "partition-order-validated",
            Self::GeometricOrderCandidate => "geometric-order-candidate",
        }
    }
}

fn connector_ordering_validation(
    operands: &[SplitResourceMergeOperand],
) -> ConnectorOrderingValidation {
    let mut needs_connector = false;
    for operand in operands {
        let endpoint_topology = split_resource_endpoint_topology(operand);
        let plan = split_resource_connector_authoring_plan(operand, endpoint_topology);
        if plan.links_required == 0 {
            continue;
        }
        needs_connector = true;
        if endpoint_topology.outside_endpoint_matches == 0
            && geometric_order_candidate(operand).is_some()
        {
            return ConnectorOrderingValidation::GeometricOrderCandidate;
        }
    }
    if needs_connector {
        ConnectorOrderingValidation::PartitionOrderValidated
    } else {
        ConnectorOrderingValidation::NotRequired
    }
}

fn connector_total_length(operands: &[SplitResourceMergeOperand]) -> f64 {
    operands
        .iter()
        .flat_map(connector_lengths)
        .fold(0.0, |total, length| total + length)
}

fn connector_authoring_links(operands: &[SplitResourceMergeOperand]) -> usize {
    operands
        .iter()
        .map(|operand| {
            split_resource_connector_authoring_plan(
                operand,
                split_resource_endpoint_topology(operand),
            )
            .links_required
        })
        .sum()
}

fn existing_boundary_connector_spans(operands: &[SplitResourceMergeOperand]) -> usize {
    operands
        .iter()
        .flat_map(existing_boundary_connector_paths)
        .count()
}

fn connector_max_length(operands: &[SplitResourceMergeOperand]) -> f64 {
    operands
        .iter()
        .flat_map(connector_lengths)
        .fold(0.0, f64::max)
}

fn connector_lengths(operand: &SplitResourceMergeOperand) -> Vec<f64> {
    operand
        .cylinder_runs
        .windows(2)
        .map(|pair| pair[0].end_position.distance_to(&pair[1].start_position))
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SplitResourceConnectorAuthoringPlan {
    links_required: usize,
}

fn split_resource_connector_authoring_plan(
    operand: &SplitResourceMergeOperand,
    endpoint_topology: SplitResourceEndpointTopology,
) -> SplitResourceConnectorAuthoringPlan {
    let links_required = if endpoint_topology.disconnected && !operand.cylinder_runs.is_empty() {
        operand.cylinder_runs.len() - 1
    } else {
        0
    };
    SplitResourceConnectorAuthoringPlan { links_required }
}

fn split_resource_connector_authoring_plan_label(
    operand: &SplitResourceMergeOperand,
) -> Option<String> {
    let endpoint_topology = split_resource_endpoint_topology(operand);
    let plan = split_resource_connector_authoring_plan(operand, endpoint_topology);
    if plan.links_required == 0 {
        return None;
    }
    let connectors = operand
        .cylinder_runs
        .windows(2)
        .enumerate()
        .map(|(index, pair)| {
            let path = directly_adjacent_connector_path(&pair[0], &pair[1]);
            format!(
                "{}:{}:{:?}->{:?}:direct_span={}",
                pair[0].cylinder_segment_index,
                pair[1].cylinder_segment_index,
                pair[0].end,
                pair[1].start,
                path.as_ref()
                    .map(|_| format!("{index}:direct"))
                    .unwrap_or_else(|| "none".to_string())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    Some(format!(
        "outside={} requests={} links={} endpoint_matches={} order={} connectors={}",
        operand.outside_segment_index,
        request_set_label(&operand.request_indices),
        plan.links_required,
        endpoint_topology.outside_endpoint_matches,
        connector_order_label(operand),
        connectors
    ))
}

fn existing_boundary_connector_paths(operand: &SplitResourceMergeOperand) -> Vec<()> {
    operand
        .cylinder_runs
        .windows(2)
        .filter_map(|pair| directly_adjacent_connector_path(&pair[0], &pair[1]))
        .collect()
}

fn directly_adjacent_connector_path(
    first: &super::SplitResourceCylinderRun,
    second: &super::SplitResourceCylinderRun,
) -> Option<()> {
    (first.loop_ == second.loop_ && Some(first.after_last) == second.half_edges.first().copied())
        .then_some(())
}

fn connector_order_label(operand: &SplitResourceMergeOperand) -> String {
    geometric_order_candidate(operand)
        .map(|candidate| {
            candidate
                .order
                .iter()
                .map(|index| {
                    operand.cylinder_runs[*index]
                        .cylinder_segment_index
                        .to_string()
                })
                .collect::<Vec<_>>()
                .join(">")
        })
        .unwrap_or_else(|| "unresolved".to_string())
}

#[derive(Debug, Clone, PartialEq)]
struct GeometricOrderCandidate {
    order: Vec<usize>,
    length: f64,
}

fn geometric_order_candidate(
    operand: &SplitResourceMergeOperand,
) -> Option<GeometricOrderCandidate> {
    if operand.cylinder_runs.len() < 2 {
        return None;
    }
    let mut best = None;
    let mut order = Vec::with_capacity(operand.cylinder_runs.len());
    let mut remaining = (0..operand.cylinder_runs.len()).collect::<Vec<_>>();
    search_geometric_order(operand, &mut order, &mut remaining, &mut best);
    best
}

fn search_geometric_order(
    operand: &SplitResourceMergeOperand,
    order: &mut Vec<usize>,
    remaining: &mut Vec<usize>,
    best: &mut Option<GeometricOrderCandidate>,
) {
    if remaining.is_empty() {
        let length = geometric_order_length(operand, order);
        if best
            .as_ref()
            .is_none_or(|candidate| length < candidate.length)
        {
            *best = Some(GeometricOrderCandidate {
                order: order.clone(),
                length,
            });
        }
        return;
    }
    for index in 0..remaining.len() {
        let run_index = remaining.remove(index);
        order.push(run_index);
        search_geometric_order(operand, order, remaining, best);
        order.pop();
        remaining.insert(index, run_index);
    }
}

fn geometric_order_length(operand: &SplitResourceMergeOperand, order: &[usize]) -> f64 {
    order
        .windows(2)
        .map(|pair| {
            operand.cylinder_runs[pair[0]]
                .end_position
                .distance_to(&operand.cylinder_runs[pair[1]].start_position)
        })
        .sum()
}

fn request_set_label(set: &[usize]) -> String {
    set.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("+")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SplitResourceEndpointTopology {
    disconnected: bool,
    outside_endpoint_matches: usize,
}

fn split_resource_endpoint_topology(
    operand: &SplitResourceMergeOperand,
) -> SplitResourceEndpointTopology {
    let starts = operand
        .cylinder_runs
        .iter()
        .map(|run| run.start)
        .collect::<BTreeSet<_>>();
    let ends = operand
        .cylinder_runs
        .iter()
        .map(|run| run.end)
        .collect::<BTreeSet<_>>();
    let endpoints = starts.union(&ends).copied().collect::<BTreeSet<_>>();
    let outside_endpoint_matches = [operand.outside_start, operand.outside_end]
        .into_iter()
        .filter(|endpoint| endpoints.contains(endpoint))
        .count();
    let disconnected = operand.cylinder_runs.len() > 1
        && starts.iter().filter(|start| ends.contains(start)).count()
            + ends.iter().filter(|end| starts.contains(end)).count()
            == 0;
    SplitResourceEndpointTopology {
        disconnected,
        outside_endpoint_matches,
    }
}
