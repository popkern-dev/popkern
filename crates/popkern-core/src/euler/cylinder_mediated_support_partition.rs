//! Canonical outside/cap interval partitioning for cylinder-mediated support.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::CoreError;
use crate::topology::{BRep, HalfEdgeHandle, LoopHandle};

use super::cylinder_mediated_support::{
    CylinderMediatedOutsideIntervalInteraction, CylinderMediatedSupportContinuationAggregateDelta,
    CylinderMediatedSupportContinuationPlan,
};

/// Canonical partition of selected outside/cap intervals by loop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedOutsideIntervalPartition {
    /// Per-loop partitions in deterministic loop-handle order.
    pub loops: Vec<CylinderMediatedOutsideLoopPartition>,
}

/// Canonical partition of selected outside/cap intervals on one loop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedOutsideLoopPartition {
    /// Outside loop being partitioned.
    pub outside_loop: LoopHandle,
    /// Maximal directed selected spans with identical request ownership.
    pub segments: Vec<CylinderMediatedOutsideIntervalPartitionSegment>,
}

/// A maximal directed outside/cap span selected by the same request set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedOutsideIntervalPartitionSegment {
    /// Half-open selected outside/cap half-edges in loop direction.
    pub half_edges: Vec<HalfEdgeHandle>,
    /// Batch request indices that select every half-edge in this segment.
    pub request_indices: Vec<usize>,
}

pub(super) fn partition_outside_intervals(
    brep: &BRep,
    plans: &[CylinderMediatedSupportContinuationPlan],
) -> Result<CylinderMediatedOutsideIntervalPartition, CoreError> {
    let mut selected_by_loop: BTreeMap<LoopHandle, BTreeMap<HalfEdgeHandle, BTreeSet<usize>>> =
        BTreeMap::new();
    for (request_index, plan) in plans.iter().enumerate() {
        let loop_membership = selected_by_loop
            .entry(plan.request.outside_loop)
            .or_default();
        for &half_edge in half_open_interval(&plan.outside_interval)? {
            let data = brep.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "outside interval partition references stale half-edge {half_edge:?}"
                ))
            })?;
            if data.loop_ != plan.request.outside_loop {
                return Err(CoreError::PreconditionFailed(format!(
                    "outside interval partition half-edge {half_edge:?} belongs to {:?}, expected {:?}",
                    data.loop_, plan.request.outside_loop
                )));
            }
            loop_membership
                .entry(half_edge)
                .or_default()
                .insert(request_index);
        }
    }

    let loops = selected_by_loop
        .into_iter()
        .map(|(outside_loop, membership)| {
            partition_loop(brep, outside_loop, &membership).map(|segments| {
                CylinderMediatedOutsideLoopPartition {
                    outside_loop,
                    segments,
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(CylinderMediatedOutsideIntervalPartition { loops })
}

pub(super) fn outside_interval_interaction(
    plans: &[CylinderMediatedSupportContinuationPlan],
) -> Result<CylinderMediatedOutsideIntervalInteraction, CoreError> {
    let mut seen = BTreeSet::new();
    let mut shared = BTreeSet::new();
    for plan in plans {
        for half_edge in half_open_interval(&plan.outside_interval)? {
            if !seen.insert(*half_edge) {
                shared.insert(*half_edge);
            }
        }
    }
    if shared.is_empty() {
        Ok(CylinderMediatedOutsideIntervalInteraction::Disjoint)
    } else {
        Ok(CylinderMediatedOutsideIntervalInteraction::Overlapping {
            shared_half_edges: shared.into_iter().collect(),
        })
    }
}

pub(super) fn aggregate_delta_from_partition(
    partition: &CylinderMediatedOutsideIntervalPartition,
) -> CylinderMediatedSupportContinuationAggregateDelta {
    let segments = partition
        .loops
        .iter()
        .map(|loop_partition| loop_partition.segments.len() as i64)
        .sum::<i64>();
    CylinderMediatedSupportContinuationAggregateDelta::Known {
        vertices: 0,
        edges: segments * 2,
        faces: segments,
        loops: segments,
    }
}

fn partition_loop(
    brep: &BRep,
    outside_loop: LoopHandle,
    membership: &BTreeMap<HalfEdgeHandle, BTreeSet<usize>>,
) -> Result<Vec<CylinderMediatedOutsideIntervalPartitionSegment>, CoreError> {
    let mut ordered = Vec::new();
    let start = brep.loops.get(outside_loop).ok_or_else(|| {
        CoreError::InvalidHandle(format!("outside loop {outside_loop:?} is not live"))
    })?;
    let mut cursor = start.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let data = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "outside interval partition references stale loop half-edge {cursor:?}"
            ))
        })?;
        if data.loop_ != outside_loop {
            return Err(CoreError::PreconditionFailed(format!(
                "outside interval partition traversal left loop {outside_loop:?} at {cursor:?}"
            )));
        }
        ordered.push(cursor);
        cursor = data.next;
        if cursor == start.half_edge {
            break;
        }
    }
    if cursor != start.half_edge {
        return Err(CoreError::PreconditionFailed(format!(
            "outside interval partition could not close loop {outside_loop:?}"
        )));
    }

    let mut segments: Vec<CylinderMediatedOutsideIntervalPartitionSegment> = Vec::new();
    let mut current_half_edges = Vec::new();
    let mut current_request_indices = Vec::new();
    for half_edge in ordered {
        let Some(request_indices) = membership.get(&half_edge) else {
            push_segment(
                &mut segments,
                &mut current_half_edges,
                &mut current_request_indices,
            );
            continue;
        };
        let request_indices = request_indices.iter().copied().collect::<Vec<_>>();
        if current_half_edges.is_empty() || current_request_indices == request_indices {
            current_half_edges.push(half_edge);
            current_request_indices = request_indices;
        } else {
            push_segment(
                &mut segments,
                &mut current_half_edges,
                &mut current_request_indices,
            );
            current_half_edges.push(half_edge);
            current_request_indices = request_indices;
        }
    }
    push_segment(
        &mut segments,
        &mut current_half_edges,
        &mut current_request_indices,
    );
    merge_wrapped_segment(&mut segments);
    Ok(segments)
}

fn push_segment(
    segments: &mut Vec<CylinderMediatedOutsideIntervalPartitionSegment>,
    current_half_edges: &mut Vec<HalfEdgeHandle>,
    current_request_indices: &mut Vec<usize>,
) {
    if current_half_edges.is_empty() {
        return;
    }
    segments.push(CylinderMediatedOutsideIntervalPartitionSegment {
        half_edges: std::mem::take(current_half_edges),
        request_indices: std::mem::take(current_request_indices),
    });
}

fn merge_wrapped_segment(segments: &mut Vec<CylinderMediatedOutsideIntervalPartitionSegment>) {
    if segments.len() < 2 {
        return;
    }
    let first_request_indices = segments[0].request_indices.clone();
    let last_index = segments.len() - 1;
    if segments[last_index].request_indices != first_request_indices {
        return;
    }
    let mut first = segments.remove(0);
    segments[last_index - 1]
        .half_edges
        .append(&mut first.half_edges);
}

fn half_open_interval(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.is_empty() {
        return Err(CoreError::PreconditionFailed(
            "outside endpoint interval is empty".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}
