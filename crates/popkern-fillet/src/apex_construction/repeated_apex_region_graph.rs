//! Loop-level region graph views for repeated-apex side-support planning.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_facts::{EndpointRegionMembership, SideSupportRegionFacts};
use super::repeated_apex_topology_knowledge::LoopKind;

#[derive(Clone, Debug)]
pub(super) struct RegionGraph {
    pub(super) nodes: Vec<RegionGraphNode>,
    pub(super) trim_endpoint_pairs: Vec<RegionTrimEndpointPair>,
}

#[derive(Clone, Debug)]
pub(super) struct RegionGraphNode {
    pub(super) face: FaceHandle,
    pub(super) loop_handle: LoopHandle,
    pub(super) kind: LoopKind,
    pub(super) boundary_half_edges: Vec<HalfEdgeHandle>,
    pub(super) boundary_vertices: Vec<VertexHandle>,
}

#[derive(Clone, Debug)]
pub(super) struct RegionTrimEndpointPair {
    pub(super) fillet_index: usize,
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) start_membership: EndpointRegionMembership,
    pub(super) end_membership: EndpointRegionMembership,
}

impl RegionGraph {
    pub(super) fn build(facts: &SideSupportRegionFacts) -> Self {
        let nodes = facts
            .support_families
            .iter()
            .flat_map(|family| {
                family.loops.iter().map(|loop_facts| RegionGraphNode {
                    face: loop_facts.face,
                    loop_handle: loop_facts.loop_handle,
                    kind: loop_facts.kind,
                    boundary_half_edges: loop_facts.boundary_half_edges.clone(),
                    boundary_vertices: loop_facts.boundary_vertices.clone(),
                })
            })
            .collect::<Vec<_>>();

        let trim_endpoint_pairs = facts
            .trim_endpoint_pairs
            .iter()
            .map(|trim| RegionTrimEndpointPair {
                fillet_index: trim.fillet_index,
                cylinder_face: trim.cylinder_face,
                support_face: trim.support_face,
                start: trim.start,
                end: trim.end,
                start_membership: trim.start_membership.clone(),
                end_membership: trim.end_membership.clone(),
            })
            .collect();

        Self {
            nodes,
            trim_endpoint_pairs,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "region_graph nodes=[{}] generator_pairs=[{}]",
            self.nodes
                .iter()
                .map(RegionGraphNode::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.trim_endpoint_pairs
                .iter()
                .map(RegionTrimEndpointPair::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl RegionGraphNode {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?} half_edges={:?} vertices={:?}",
            self.face,
            self.loop_handle,
            self.kind,
            self.boundary_half_edges,
            self.boundary_vertices
        )
    }
}

impl RegionTrimEndpointPair {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} support={:?} {:?}->{:?} start={} end={}",
            self.fillet_index,
            self.cylinder_face,
            self.support_face,
            self.start,
            self.end,
            self.start_membership.diagnostic_label(),
            self.end_membership.diagnostic_label()
        )
    }
}
