//! Topology preflight for KP20 support-cell segment-owned boundaries.

use std::collections::HashSet;

use popkern_core::topology::HalfEdgeHandle;

use super::trim_ownership::SegmentOwnedBoundaryTrim;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SegmentBoundaryTopologyStatus {
    AlreadySegmentTopology,
    EdgeSplitRequired,
    InvalidTopologyAddress,
    DuplicateSegmentEdge,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SegmentBoundaryTopologyPlan {
    pub(super) segment_trim_count: usize,
    pub(super) already_segment_topology_count: usize,
    pub(super) edge_split_required_count: usize,
    pub(super) edge_split_parameter_count: usize,
    pub(super) invalid_topology_address_count: usize,
    pub(super) duplicate_segment_edge_count: usize,
    pub(super) status: SegmentBoundaryTopologyStatus,
}

impl SegmentBoundaryTopologyPlan {
    pub(super) fn is_already_segmented(self) -> bool {
        self.status == SegmentBoundaryTopologyStatus::AlreadySegmentTopology
            && self.segment_trim_count == self.already_segment_topology_count
            && self.edge_split_required_count == 0
            && self.edge_split_parameter_count == 0
            && self.invalid_topology_address_count == 0
            && self.duplicate_segment_edge_count == 0
    }
}

pub(super) fn segment_boundary_topology_plan(
    records: &[SegmentOwnedBoundaryTrim],
) -> SegmentBoundaryTopologyPlan {
    let mut seen_edges = HashSet::new();
    let mut seen_half_edge_pairs: HashSet<(String, String)> = HashSet::new();
    let mut plan = SegmentBoundaryTopologyPlan {
        segment_trim_count: records.len(),
        already_segment_topology_count: 0,
        edge_split_required_count: 0,
        edge_split_parameter_count: 0,
        invalid_topology_address_count: 0,
        duplicate_segment_edge_count: 0,
        status: SegmentBoundaryTopologyStatus::AlreadySegmentTopology,
    };
    for &record in records {
        if !record.has_topology_address() || !record.has_paired_half_edges() {
            plan.invalid_topology_address_count += 1;
            continue;
        }
        if !seen_edges.insert(record.edge)
            || !seen_half_edge_pairs.insert(normalized_half_edge_pair(
                record.half_edge_a,
                record.half_edge_b,
            ))
        {
            plan.duplicate_segment_edge_count += 1;
            continue;
        }
        if record.requires_edge_split() {
            plan.edge_split_required_count += 1;
            plan.edge_split_parameter_count += record.split_parameters().len();
        } else {
            plan.already_segment_topology_count += 1;
        }
    }
    plan.status = segment_boundary_topology_status(plan);
    plan
}

fn segment_boundary_topology_status(
    plan: SegmentBoundaryTopologyPlan,
) -> SegmentBoundaryTopologyStatus {
    if plan.invalid_topology_address_count > 0 {
        SegmentBoundaryTopologyStatus::InvalidTopologyAddress
    } else if plan.duplicate_segment_edge_count > 0 {
        SegmentBoundaryTopologyStatus::DuplicateSegmentEdge
    } else if plan.edge_split_required_count > 0 {
        SegmentBoundaryTopologyStatus::EdgeSplitRequired
    } else {
        SegmentBoundaryTopologyStatus::AlreadySegmentTopology
    }
}

fn normalized_half_edge_pair(a: HalfEdgeHandle, b: HalfEdgeHandle) -> (String, String) {
    if handle_key(a) <= handle_key(b) {
        (handle_key(a), handle_key(b))
    } else {
        (handle_key(b), handle_key(a))
    }
}

fn handle_key(handle: HalfEdgeHandle) -> String {
    format!("{handle:?}")
}
