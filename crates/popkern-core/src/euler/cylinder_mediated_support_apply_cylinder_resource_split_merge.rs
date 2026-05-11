//! Split cylinder-resource merge operands for mediated support continuation.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, VertexHandle};

use super::super::super::apply_bridge_cycles::SegmentBridgeCycleOperandSet;
use super::super::CylinderBridgeResourceSegment;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SplitResourceMergeOperand {
    outside_segment_index: usize,
    request_indices: Vec<usize>,
    outside_start: VertexHandle,
    outside_end: VertexHandle,
    outside_half_edges: usize,
    cylinder_runs: Vec<SplitResourceCylinderRun>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SplitResourceCylinderRun {
    cylinder_segment_index: usize,
    request_indices: Vec<usize>,
    start: VertexHandle,
    end: VertexHandle,
    half_edges: usize,
}

pub(super) fn split_resource_merge_operands(
    brep: &BRep,
    bridge_cycles: &SegmentBridgeCycleOperandSet,
    missing_outside_segments: &[(usize, Vec<usize>)],
    cylinder_segments: &[CylinderBridgeResourceSegment],
) -> Result<Vec<SplitResourceMergeOperand>, CoreError> {
    let mut operands = Vec::new();
    let cylinder_sets = cylinder_segments
        .iter()
        .map(|segment| segment.request_indices.clone())
        .collect::<Vec<_>>();
    for (outside_segment_index, missing) in missing_outside_segments {
        if !cylinder_request_sets_cover_missing_set(missing, &cylinder_sets) {
            continue;
        }
        let cylinder_runs = cylinder_segments
            .iter()
            .enumerate()
            .filter(|(_, segment)| {
                cylinder_request_set_participates_in_missing_set(missing, &segment.request_indices)
            })
            .map(|(cylinder_segment_index, segment)| {
                split_resource_cylinder_run(brep, cylinder_segment_index, segment)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if cylinder_runs.len() < 2 {
            continue;
        }
        let Some((outside_start, outside_end)) =
            bridge_cycles.boundary_endpoint_pair(*outside_segment_index)
        else {
            continue;
        };
        operands.push(SplitResourceMergeOperand {
            outside_segment_index: *outside_segment_index,
            request_indices: missing.clone(),
            outside_start,
            outside_end,
            outside_half_edges: bridge_cycles
                .boundary_half_edge_count(*outside_segment_index)
                .unwrap_or_default(),
            cylinder_runs,
        });
    }
    Ok(operands)
}

fn split_resource_cylinder_run(
    brep: &BRep,
    cylinder_segment_index: usize,
    segment: &CylinderBridgeResourceSegment,
) -> Result<SplitResourceCylinderRun, CoreError> {
    let first = segment.half_edges.first().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "split cylinder resource run {cylinder_segment_index} is empty"
        ))
    })?;
    let last = segment.half_edges.last().copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "split cylinder resource run {cylinder_segment_index} is empty"
        ))
    })?;
    let start = brep
        .half_edges
        .get(first)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "split cylinder resource run {cylinder_segment_index} has stale first half-edge {first:?}"
            ))
        })?;
    let after_last = brep
        .half_edges
        .get(last)
        .map(|half_edge| half_edge.next)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "split cylinder resource run {cylinder_segment_index} has stale last half-edge {last:?}"
            ))
        })?;
    let end = brep
        .half_edges
        .get(after_last)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "split cylinder resource run {cylinder_segment_index} has stale after-last half-edge {after_last:?}"
            ))
        })?;
    Ok(SplitResourceCylinderRun {
        cylinder_segment_index,
        request_indices: segment.request_indices.clone(),
        start,
        end,
        half_edges: segment.half_edges.len(),
    })
}

pub(super) fn cylinder_request_sets_cover_missing_set(
    missing: &[usize],
    cylinder_sets: &[Vec<usize>],
) -> bool {
    let covered = cylinder_sets
        .iter()
        .filter(|set| cylinder_request_set_participates_in_missing_set(missing, set))
        .flat_map(|set| set.iter().copied())
        .collect::<BTreeSet<_>>();
    missing.iter().all(|index| covered.contains(index))
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct SplitResourceMergeSummary {
    missing_sets: usize,
    coverable_sets: usize,
    candidate_cylinder_runs: usize,
}

impl SplitResourceMergeSummary {
    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "split_resource_merge_candidates=missing_sets={} coverable_sets={} cylinder_runs={}",
            self.missing_sets, self.coverable_sets, self.candidate_cylinder_runs
        )
    }
}

pub(super) fn split_resource_merge_summary(
    operands: &[SplitResourceMergeOperand],
) -> SplitResourceMergeSummary {
    let mut summary = SplitResourceMergeSummary {
        missing_sets: operands.len(),
        ..SplitResourceMergeSummary::default()
    };
    for operand in operands {
        summary.coverable_sets += 1;
        summary.candidate_cylinder_runs += operand.cylinder_runs.len();
    }
    summary
}

pub(super) fn split_resource_merge_operand_summary(
    operands: &[SplitResourceMergeOperand],
) -> String {
    format!(
        "{} split_resource_merge_operands=[{}]",
        split_resource_endpoint_shape_summary(operands).diagnostic_label(),
        operands
            .iter()
            .map(split_resource_merge_operand_label)
            .collect::<Vec<_>>()
            .join(";")
    )
}

pub(super) fn split_resource_operands_are_shared_start(
    operands: &[SplitResourceMergeOperand],
) -> bool {
    !operands.is_empty()
        && operands
            .iter()
            .all(|operand| split_resource_endpoint_shape(operand).is_shared_start())
}

fn split_resource_merge_operand_label(operand: &SplitResourceMergeOperand) -> String {
    format!(
        "outside={} requests={} endpoint_shape={} outside_endpoints={:?}->{:?} outside_half_edges={} cylinder_runs=[{}]",
        operand.outside_segment_index,
        request_set_label(&operand.request_indices),
        split_resource_endpoint_shape(operand).diagnostic_label(),
        operand.outside_start,
        operand.outside_end,
        operand.outside_half_edges,
        operand
            .cylinder_runs
            .iter()
            .map(split_resource_cylinder_run_label)
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SplitResourceEndpointShape {
    Aligned,
    SharedStart,
    SharedEnd,
    Mixed,
}

impl SplitResourceEndpointShape {
    fn is_shared_start(self) -> bool {
        self == Self::SharedStart
    }

    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::Aligned => "aligned",
            Self::SharedStart => "shared-start",
            Self::SharedEnd => "shared-end",
            Self::Mixed => "mixed",
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct SplitResourceEndpointShapeSummary {
    shared_start: usize,
    shared_end: usize,
    aligned: usize,
    mixed: usize,
}

impl SplitResourceEndpointShapeSummary {
    fn diagnostic_label(self) -> String {
        format!(
            "split_resource_endpoint_shapes=shared_start:{} shared_end:{} aligned:{} mixed:{}",
            self.shared_start, self.shared_end, self.aligned, self.mixed
        )
    }
}

fn split_resource_endpoint_shape_summary(
    operands: &[SplitResourceMergeOperand],
) -> SplitResourceEndpointShapeSummary {
    let mut summary = SplitResourceEndpointShapeSummary::default();
    for operand in operands {
        match split_resource_endpoint_shape(operand) {
            SplitResourceEndpointShape::Aligned => summary.aligned += 1,
            SplitResourceEndpointShape::SharedStart => summary.shared_start += 1,
            SplitResourceEndpointShape::SharedEnd => summary.shared_end += 1,
            SplitResourceEndpointShape::Mixed => summary.mixed += 1,
        }
    }
    summary
}

fn split_resource_endpoint_shape(
    operand: &SplitResourceMergeOperand,
) -> SplitResourceEndpointShape {
    let unique_starts = operand
        .cylinder_runs
        .iter()
        .map(|run| run.start)
        .collect::<BTreeSet<_>>()
        .len();
    let unique_ends = operand
        .cylinder_runs
        .iter()
        .map(|run| run.end)
        .collect::<BTreeSet<_>>()
        .len();
    match (unique_starts, unique_ends) {
        (1, 1) => SplitResourceEndpointShape::Aligned,
        (1, _) => SplitResourceEndpointShape::SharedStart,
        (_, 1) => SplitResourceEndpointShape::SharedEnd,
        _ => SplitResourceEndpointShape::Mixed,
    }
}

fn split_resource_cylinder_run_label(run: &SplitResourceCylinderRun) -> String {
    format!(
        "{}:{}:{:?}->{:?}:{}",
        run.cylinder_segment_index,
        request_set_label(&run.request_indices),
        run.start,
        run.end,
        run.half_edges
    )
}

fn cylinder_request_set_participates_in_missing_set(missing: &[usize], set: &[usize]) -> bool {
    set.iter().all(|index| missing.contains(index))
        || missing.iter().all(|index| set.contains(index))
}

fn request_set_label(set: &[usize]) -> String {
    set.iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("+")
}
