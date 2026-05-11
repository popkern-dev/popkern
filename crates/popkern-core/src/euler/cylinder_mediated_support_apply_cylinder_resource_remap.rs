//! Remap pre-exposure cylinder resources to post-exposure boundary sides.

use std::collections::BTreeMap;

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::CylinderMediatedSupportContinuationBatchPlan;
use super::apply_cylinder_resources::CylinderBridgeResourcePartition;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourceRemap {
    segments: Vec<CylinderBridgeResourceRemappedSegment>,
    unavailable_edges: Vec<CylinderBridgeUnavailableEdge>,
    owner_roles: CylinderBridgeUnavailableOwnerRoles,
    blocking_owner_roles: CylinderBridgeUnavailableOwnerRoles,
    blocking_spans: CylinderBridgeBlockingSpanSummary,
    non_contiguous_segments: usize,
}

impl CylinderBridgeResourceRemap {
    pub(super) fn diagnostic_label(&self) -> String {
        let half_edges = self
            .segments
            .iter()
            .map(|segment| segment.half_edges.len())
            .sum::<usize>();
        let reversed_segments = self
            .segments
            .iter()
            .filter(|segment| segment.orientation == CylinderBridgeResourceOrientation::Reversed)
            .count();
        format!(
            "cylinder_resource_remap=segments={} half_edges={} reversed_segments={} unavailable_edges={} unavailable_not_boundary={} unavailable_live_twin={} unavailable_desired_live={} unavailable_desired_retired={} unavailable_owner_roles={} unavailable_blocking_roles={} {} unavailable_owner_faces=[{}] non_contiguous_segments={}",
            self.segments.len(),
            half_edges,
            reversed_segments,
            self.unavailable_edges.len(),
            self.unavailable_not_boundary_count(),
            self.unavailable_live_twin_count(),
            self.unavailable_desired_live_count(),
            self.unavailable_desired_retired_count(),
            self.owner_roles.diagnostic_label(),
            self.blocking_owner_roles.diagnostic_label(),
            self.blocking_spans.diagnostic_label(),
            unavailable_owner_face_summary(&self.unavailable_edges),
            self.non_contiguous_segments
        )
    }

    fn unavailable_not_boundary_count(&self) -> usize {
        self.unavailable_edges
            .iter()
            .filter(|edge| edge.not_boundary)
            .count()
    }

    fn unavailable_live_twin_count(&self) -> usize {
        self.unavailable_edges
            .iter()
            .filter(|edge| edge.live_twin)
            .count()
    }

    fn unavailable_desired_live_count(&self) -> usize {
        self.unavailable_edges
            .iter()
            .filter(|edge| edge.desired_side_live)
            .count()
    }

    fn unavailable_desired_retired_count(&self) -> usize {
        self.unavailable_edges
            .iter()
            .filter(|edge| !edge.desired_side_live)
            .count()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CylinderBridgeResourceRemappedSegment {
    half_edges: Vec<HalfEdgeHandle>,
    orientation: CylinderBridgeResourceOrientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CylinderBridgeResourceOrientation {
    Forward,
    Reversed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CylinderBridgeUnavailableEdge {
    not_boundary: bool,
    live_twin: bool,
    desired_side_live: bool,
    owner_faces: Vec<FaceHandle>,
    blocking_faces: Vec<FaceHandle>,
    blocking_half_edges: Vec<HalfEdgeHandle>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct CylinderBridgeUnavailableOwnerRoles {
    cylinder: usize,
    support: usize,
    context: usize,
}

impl CylinderBridgeUnavailableOwnerRoles {
    fn diagnostic_label(self) -> String {
        format!(
            "cylinder:{} support:{} context:{}",
            self.cylinder, self.support, self.context
        )
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct CylinderBridgeBlockingSpanSummary {
    spans: usize,
    half_edges: usize,
    non_contiguous_spans: usize,
}

impl CylinderBridgeBlockingSpanSummary {
    fn diagnostic_label(self) -> String {
        format!(
            "blocking_spans={} blocking_half_edges={} non_contiguous_blocking_spans={}",
            self.spans, self.half_edges, self.non_contiguous_spans
        )
    }
}

pub(super) fn remap_cylinder_bridge_resources(
    brep: &BRep,
    partition: &CylinderBridgeResourcePartition,
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> Result<CylinderBridgeResourceRemap, CoreError> {
    let mut segments = Vec::new();
    let mut unavailable_edges = Vec::new();
    let mut non_contiguous_segments = 0usize;
    for segment in partition.segments() {
        let mut half_edges = Vec::with_capacity(segment.edges.len());
        for (&edge, &desired_half_edge) in segment.edges.iter().zip(&segment.half_edges) {
            let remap = remapped_boundary_half_edge(brep, edge, desired_half_edge)?;
            let Some(half_edge) = remap.half_edge else {
                unavailable_edges.push(remap.unavailable);
                continue;
            };
            half_edges.push(half_edge);
        }
        if !half_edges.is_empty() {
            match orient_remapped_segment(brep, half_edges)? {
                Some(segment) => segments.push(segment),
                None => non_contiguous_segments += 1,
            }
        }
    }
    let blocking_spans = summarize_blocking_spans(brep, &unavailable_edges)?;
    Ok(CylinderBridgeResourceRemap {
        segments,
        owner_roles: classify_unavailable_owner_roles(&unavailable_edges, plan),
        blocking_owner_roles: classify_unavailable_blocking_roles(&unavailable_edges, plan),
        blocking_spans,
        unavailable_edges,
        non_contiguous_segments,
    })
}

struct BoundaryHalfEdgeRemap {
    half_edge: Option<HalfEdgeHandle>,
    unavailable: CylinderBridgeUnavailableEdge,
}

fn remapped_boundary_half_edge(
    brep: &BRep,
    edge: EdgeHandle,
    desired_half_edge: HalfEdgeHandle,
) -> Result<BoundaryHalfEdgeRemap, CoreError> {
    let half_edge = brep
        .edges
        .get(edge)
        .map(|edge_data| edge_data.half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "cylinder resource remap references stale edge {edge:?}"
            ))
        })?;
    let half_edge_data = brep.half_edges.get(half_edge).ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "cylinder resource remap edge {edge:?} has stale half-edge"
        ))
    })?;
    let not_boundary = !brep.boundary_edges.contains(&edge);
    let live_twin = half_edge_data.twin != HalfEdgeHandle::default();
    let unavailable = CylinderBridgeUnavailableEdge {
        not_boundary,
        live_twin,
        desired_side_live: brep.half_edges.contains_key(desired_half_edge),
        owner_faces: edge_owner_faces(brep, half_edge),
        blocking_faces: desired_side_blocking_faces(brep, desired_half_edge),
        blocking_half_edges: desired_side_blocking_half_edges(brep, desired_half_edge),
    };
    if not_boundary || live_twin {
        return Ok(BoundaryHalfEdgeRemap {
            half_edge: None,
            unavailable,
        });
    }
    Ok(BoundaryHalfEdgeRemap {
        half_edge: Some(half_edge),
        unavailable,
    })
}

fn desired_side_blocking_faces(brep: &BRep, desired_half_edge: HalfEdgeHandle) -> Vec<FaceHandle> {
    let Some(desired) = brep.half_edges.get(desired_half_edge) else {
        return Vec::new();
    };
    let twin = desired.twin;
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return Vec::new();
    };
    brep.loops
        .get(twin_data.loop_)
        .map(|loop_data| vec![loop_data.face])
        .unwrap_or_default()
}

fn desired_side_blocking_half_edges(
    brep: &BRep,
    desired_half_edge: HalfEdgeHandle,
) -> Vec<HalfEdgeHandle> {
    let Some(desired) = brep.half_edges.get(desired_half_edge) else {
        return Vec::new();
    };
    if brep.half_edges.contains_key(desired.twin) {
        vec![desired.twin]
    } else {
        Vec::new()
    }
}

fn summarize_blocking_spans(
    brep: &BRep,
    unavailable_edges: &[CylinderBridgeUnavailableEdge],
) -> Result<CylinderBridgeBlockingSpanSummary, CoreError> {
    let mut spans = Vec::<Vec<HalfEdgeHandle>>::new();
    let mut current_loop: Option<LoopHandle> = None;
    let mut current_span = Vec::new();
    for blocking_half_edge in unavailable_edges
        .iter()
        .flat_map(|edge| edge.blocking_half_edges.iter().copied())
    {
        let loop_ = brep
            .half_edges
            .get(blocking_half_edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "cylinder resource blocker {blocking_half_edge:?} is stale"
                ))
            })?
            .loop_;
        if current_loop == Some(loop_) {
            current_span.push(blocking_half_edge);
            continue;
        }
        if !current_span.is_empty() {
            spans.push(std::mem::take(&mut current_span));
        }
        current_loop = Some(loop_);
        current_span.push(blocking_half_edge);
    }
    if !current_span.is_empty() {
        spans.push(current_span);
    }
    let half_edges = spans.iter().map(Vec::len).sum();
    let mut non_contiguous_spans = 0usize;
    for span in &mut spans {
        span.reverse();
        if !is_directed_path(brep, span)? {
            non_contiguous_spans += 1;
        }
    }
    Ok(CylinderBridgeBlockingSpanSummary {
        spans: spans.len(),
        half_edges,
        non_contiguous_spans,
    })
}

fn edge_owner_faces(brep: &BRep, half_edge: HalfEdgeHandle) -> Vec<FaceHandle> {
    [
        Some(half_edge),
        brep.half_edges.get(half_edge).map(|data| data.twin),
    ]
    .into_iter()
    .flatten()
    .filter_map(|candidate| {
        let loop_ = brep.half_edges.get(candidate)?.loop_;
        brep.loops.get(loop_).map(|loop_data| loop_data.face)
    })
    .fold(Vec::new(), |mut faces, face| {
        if !faces.contains(&face) {
            faces.push(face);
        }
        faces
    })
}

fn classify_unavailable_owner_roles(
    unavailable_edges: &[CylinderBridgeUnavailableEdge],
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> CylinderBridgeUnavailableOwnerRoles {
    let mut roles = CylinderBridgeUnavailableOwnerRoles::default();
    for unavailable in unavailable_edges {
        for &face in &unavailable.owner_faces {
            if plan
                .plans
                .iter()
                .any(|request| request.request.cylinder_face == face)
            {
                roles.cylinder += 1;
            } else if plan
                .plans
                .iter()
                .any(|request| request.request.support_face == face)
            {
                roles.support += 1;
            } else {
                roles.context += 1;
            }
        }
    }
    roles
}

fn classify_unavailable_blocking_roles(
    unavailable_edges: &[CylinderBridgeUnavailableEdge],
    plan: &CylinderMediatedSupportContinuationBatchPlan,
) -> CylinderBridgeUnavailableOwnerRoles {
    let mut roles = CylinderBridgeUnavailableOwnerRoles::default();
    for unavailable in unavailable_edges {
        for &face in &unavailable.blocking_faces {
            if plan
                .plans
                .iter()
                .any(|request| request.request.cylinder_face == face)
            {
                roles.cylinder += 1;
            } else if plan
                .plans
                .iter()
                .any(|request| request.request.support_face == face)
            {
                roles.support += 1;
            } else {
                roles.context += 1;
            }
        }
    }
    roles
}

fn orient_remapped_segment(
    brep: &BRep,
    mut half_edges: Vec<HalfEdgeHandle>,
) -> Result<Option<CylinderBridgeResourceRemappedSegment>, CoreError> {
    if is_directed_path(brep, &half_edges)? {
        return Ok(Some(CylinderBridgeResourceRemappedSegment {
            half_edges,
            orientation: CylinderBridgeResourceOrientation::Forward,
        }));
    }
    half_edges.reverse();
    if is_directed_path(brep, &half_edges)? {
        return Ok(Some(CylinderBridgeResourceRemappedSegment {
            half_edges,
            orientation: CylinderBridgeResourceOrientation::Reversed,
        }));
    }
    Ok(None)
}

fn is_directed_path(brep: &BRep, half_edges: &[HalfEdgeHandle]) -> Result<bool, CoreError> {
    for pair in half_edges.windows(2) {
        let first = brep.half_edges.get(pair[0]).ok_or_else(|| {
            CoreError::InvalidHandle("cylinder resource remap path half-edge is stale".into())
        })?;
        if first.next != pair[1] {
            return Ok(false);
        }
    }
    Ok(true)
}

fn unavailable_owner_face_summary(unavailable_edges: &[CylinderBridgeUnavailableEdge]) -> String {
    let mut counts = BTreeMap::<FaceHandle, usize>::new();
    for unavailable in unavailable_edges {
        for &face in &unavailable.owner_faces {
            *counts.entry(face).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .map(|(face, count)| format!("{face:?}:{count}"))
        .collect::<Vec<_>>()
        .join(",")
}
