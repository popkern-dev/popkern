//! Remap pre-exposure cylinder resources to post-exposure boundary sides.

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

use super::apply_cylinder_resources::CylinderBridgeResourcePartition;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CylinderBridgeResourceRemap {
    segments: Vec<CylinderBridgeResourceRemappedSegment>,
    unavailable_edges: Vec<EdgeHandle>,
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
            "cylinder_resource_remap=segments={} half_edges={} reversed_segments={} unavailable_edges={} non_contiguous_segments={}",
            self.segments.len(),
            half_edges,
            reversed_segments,
            self.unavailable_edges.len(),
            self.non_contiguous_segments
        )
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

pub(super) fn remap_cylinder_bridge_resources(
    brep: &BRep,
    partition: &CylinderBridgeResourcePartition,
) -> Result<CylinderBridgeResourceRemap, CoreError> {
    let mut segments = Vec::new();
    let mut unavailable_edges = Vec::new();
    let mut non_contiguous_segments = 0usize;
    for segment in partition.segments() {
        let mut half_edges = Vec::with_capacity(segment.edges.len());
        for &edge in &segment.edges {
            let Some(half_edge) = remapped_boundary_half_edge(brep, edge)? else {
                unavailable_edges.push(edge);
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
    Ok(CylinderBridgeResourceRemap {
        segments,
        unavailable_edges,
        non_contiguous_segments,
    })
}

fn remapped_boundary_half_edge(
    brep: &BRep,
    edge: EdgeHandle,
) -> Result<Option<HalfEdgeHandle>, CoreError> {
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
    if !brep.boundary_edges.contains(&edge) || half_edge_data.twin != HalfEdgeHandle::default() {
        return Ok(None);
    }
    Ok(Some(half_edge))
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
