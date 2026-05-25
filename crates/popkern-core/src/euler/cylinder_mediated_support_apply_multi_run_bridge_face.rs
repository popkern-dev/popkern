//! Typed plan for bridge faces backed by multiple disconnected cylinder runs.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::BRep;

use super::SplitResourceMergeOperand;

#[path = "cylinder_mediated_support_apply_multi_run_bridge_face_mutation.rs"]
mod mutation;

pub(in crate::euler::cylinder_mediated_support_apply) use mutation::MultiRunBridgeFaceMutationResult;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MultiRunBridgeFacePlan {
    outside_segment_index: usize,
    request_indices: Vec<usize>,
    outside_half_edges: usize,
    cylinder_runs: usize,
    resource_links: usize,
    existing_boundary_spans: usize,
    ordered_run_indices: Vec<usize>,
    total_connector_length: f64,
    max_connector_length: f64,
}

impl MultiRunBridgeFacePlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "outside={} requests={} outside_half_edges={} cylinder_runs={} resource_links={} existing_boundary_spans={} ordering=geometric-min-length order={} total_connector_length={:.6e} max_connector_length={:.6e}",
            self.outside_segment_index,
            request_set_label(&self.request_indices),
            self.outside_half_edges,
            self.cylinder_runs,
            self.resource_links,
            self.existing_boundary_spans,
            run_order_label(self),
            self.total_connector_length,
            self.max_connector_length
        )
    }
}

pub(super) fn multi_run_bridge_face_plan_summary(operands: &[SplitResourceMergeOperand]) -> String {
    let plans = operands
        .iter()
        .filter_map(multi_run_bridge_face_plan)
        .collect::<Vec<_>>();
    let resource_links = plans.iter().map(|plan| plan.resource_links).sum::<usize>();
    let existing_boundary_spans = plans
        .iter()
        .map(|plan| plan.existing_boundary_spans)
        .sum::<usize>();
    let status = if plans.is_empty() {
        "not-applicable"
    } else if existing_boundary_spans == resource_links {
        "ready-existing-boundary-spans"
    } else {
        "requires-custom-topology-authoring"
    };
    format!(
        "multi_run_bridge_face_plan=plans={} status={} resource_links={} existing_boundary_spans={} topology_authoring={} postconditions=half-edge-pairing+cycles+manifold+delta plans=[{}]",
        plans.len(),
        status,
        resource_links,
        existing_boundary_spans,
        if status == "requires-custom-topology-authoring" {
            "required"
        } else {
            "not-required"
        },
        plans
            .iter()
            .map(MultiRunBridgeFacePlan::diagnostic_label)
            .collect::<Vec<_>>()
            .join(";")
    )
}

pub(super) fn apply_multi_run_bridge_face_plans(
    brep: &mut BRep,
    operands: &[SplitResourceMergeOperand],
) -> Result<Option<MultiRunBridgeFaceMutationResult>, CoreError> {
    mutation::apply_multi_run_bridge_face_plans(brep, operands)
}

fn multi_run_bridge_face_plan(
    operand: &SplitResourceMergeOperand,
) -> Option<MultiRunBridgeFacePlan> {
    if operand.cylinder_runs.len() < 2 {
        return None;
    }
    let order = geometric_order_candidate(operand)?;
    let resource_links = order.len().saturating_sub(1);
    Some(MultiRunBridgeFacePlan {
        outside_segment_index: operand.outside_segment_index,
        request_indices: operand.request_indices.clone(),
        outside_half_edges: operand.outside_half_edges,
        cylinder_runs: operand.cylinder_runs.len(),
        resource_links,
        existing_boundary_spans: direct_boundary_span_count(operand, &order),
        total_connector_length: geometric_order_length(operand, &order),
        max_connector_length: geometric_order_max_link(operand, &order),
        ordered_run_indices: order,
    })
}

fn direct_boundary_span_count(operand: &SplitResourceMergeOperand, order: &[usize]) -> usize {
    order
        .windows(2)
        .filter(|pair| {
            let first = &operand.cylinder_runs[pair[0]];
            let second = &operand.cylinder_runs[pair[1]];
            first.loop_ == second.loop_
                && Some(first.after_last) == second.half_edges.first().copied()
        })
        .count()
}

fn geometric_order_candidate(operand: &SplitResourceMergeOperand) -> Option<Vec<usize>> {
    let run_count = operand.cylinder_runs.len();
    if run_count < 2 {
        return None;
    }
    let mut best_order = None;
    let mut order = Vec::with_capacity(run_count);
    let mut remaining = (0..run_count).collect::<Vec<_>>();
    search_geometric_order(operand, &mut order, &mut remaining, &mut best_order);
    best_order.map(|candidate: GeometricOrderCandidate| candidate.order)
}

#[derive(Debug, Clone, PartialEq)]
struct GeometricOrderCandidate {
    order: Vec<usize>,
    length: f64,
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
    connector_lengths(operand, order).sum()
}

fn geometric_order_max_link(operand: &SplitResourceMergeOperand, order: &[usize]) -> f64 {
    connector_lengths(operand, order).fold(0.0, f64::max)
}

fn connector_lengths<'a>(
    operand: &'a SplitResourceMergeOperand,
    order: &'a [usize],
) -> impl Iterator<Item = f64> + 'a {
    order.windows(2).map(|pair| {
        operand.cylinder_runs[pair[0]]
            .end_position
            .distance_to(&operand.cylinder_runs[pair[1]].start_position)
    })
}

fn run_order_label(plan: &MultiRunBridgeFacePlan) -> String {
    plan.ordered_run_indices
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(">")
}

fn request_set_label(set: &[usize]) -> String {
    set.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("+")
}

#[allow(dead_code)]
fn _endpoint_set_size(operand: &SplitResourceMergeOperand) -> usize {
    operand
        .cylinder_runs
        .iter()
        .flat_map(|run| [run.start, run.end])
        .collect::<BTreeSet<_>>()
        .len()
}
