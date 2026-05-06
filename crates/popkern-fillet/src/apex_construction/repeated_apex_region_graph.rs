//! Loop-level region graph views for repeated-apex side-support planning.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_facts::{EndpointRegionMembership, SideSupportRegionFacts};
use super::repeated_apex_topology_knowledge::LoopKind;

#[derive(Clone, Debug)]
pub(super) struct RegionGraph {
    nodes: Vec<RegionGraphNode>,
    trim_endpoint_pairs: Vec<RegionTrimEndpointPair>,
}

#[derive(Clone, Debug)]
pub(super) struct RegionGraphNode {
    face: FaceHandle,
    loop_handle: LoopHandle,
    kind: LoopKind,
    boundary_half_edges: Vec<HalfEdgeHandle>,
    boundary_vertices: Vec<VertexHandle>,
}

#[derive(Clone, Debug)]
pub(super) struct RegionTrimEndpointPair {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    start_membership: EndpointRegionMembership,
    end_membership: EndpointRegionMembership,
}

impl RegionGraph {
    pub(super) fn build(facts: &SideSupportRegionFacts) -> Self {
        let nodes = facts
            .support_families()
            .flat_map(|family| {
                family.loops().iter().map(|loop_facts| RegionGraphNode {
                    face: loop_facts.face(),
                    loop_handle: loop_facts.loop_handle(),
                    kind: loop_facts.kind(),
                    boundary_half_edges: loop_facts.boundary_half_edges().to_vec(),
                    boundary_vertices: loop_facts.boundary_vertices().to_vec(),
                })
            })
            .collect::<Vec<_>>();

        let trim_endpoint_pairs = facts
            .trim_endpoint_pairs()
            .map(|trim| RegionTrimEndpointPair {
                fillet_index: trim.fillet_index(),
                cylinder_face: trim.cylinder_face(),
                support_face: trim.support_face(),
                start: trim.start(),
                end: trim.end(),
                start_membership: trim.start_membership().clone(),
                end_membership: trim.end_membership().clone(),
            })
            .collect();

        Self {
            nodes,
            trim_endpoint_pairs,
        }
    }

    pub(super) fn trim_endpoint_pairs(&self) -> impl Iterator<Item = &RegionTrimEndpointPair> {
        self.trim_endpoint_pairs.iter()
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
    pub(super) fn fillet_index(&self) -> usize {
        self.fillet_index
    }

    pub(super) fn cylinder_face(&self) -> FaceHandle {
        self.cylinder_face
    }

    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn start(&self) -> VertexHandle {
        self.start
    }

    pub(super) fn end(&self) -> VertexHandle {
        self.end
    }

    pub(super) fn start_membership(&self) -> &EndpointRegionMembership {
        &self.start_membership
    }

    pub(super) fn end_membership(&self) -> &EndpointRegionMembership {
        &self.end_membership
    }

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
