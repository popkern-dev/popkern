//! Classified side-support facts derived from the local knowledge graph.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_topology_knowledge::{LoopKind, TopologyKnowledgeGraph};
use super::repeated_apex_trim_subcomplex::TrimEdgeConstruction;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionFacts {
    pub(super) support_families: Vec<SideSupportFamilyFacts>,
    pub(super) trim_endpoint_pairs: Vec<SideSupportTrimEndpointFacts>,
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportFamilyFacts {
    pub(super) support_face: FaceHandle,
    pub(super) loops: Vec<SideSupportLoopFacts>,
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportLoopFacts {
    pub(super) face: FaceHandle,
    pub(super) loop_handle: LoopHandle,
    pub(super) kind: LoopKind,
    pub(super) boundary_half_edges: Vec<HalfEdgeHandle>,
    pub(super) boundary_vertices: Vec<VertexHandle>,
}

#[derive(Clone, Debug)]
pub(super) struct SideSupportTrimEndpointFacts {
    pub(super) fillet_index: usize,
    pub(super) support_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) start_membership: EndpointRegionMembership,
    pub(super) end_membership: EndpointRegionMembership,
}

#[derive(Clone, Debug)]
pub(super) enum EndpointRegionMembership {
    SameSupportLoop {
        face: FaceHandle,
        loop_handle: LoopHandle,
        half_edge: HalfEdgeHandle,
    },
    OutsideSupportFamily {
        face: FaceHandle,
        loop_handle: LoopHandle,
        half_edge: HalfEdgeHandle,
    },
    Ambiguous {
        containing_loop_count: usize,
    },
    Unknown,
}

impl SideSupportRegionFacts {
    pub(super) fn classify(knowledge: &TopologyKnowledgeGraph) -> Self {
        let support_families = knowledge
            .support_families
            .iter()
            .map(|family| SideSupportFamilyFacts {
                support_face: family.support_face,
                loops: family
                    .loops
                    .iter()
                    .map(|loop_facts| SideSupportLoopFacts {
                        face: loop_facts.face,
                        loop_handle: loop_facts.loop_handle,
                        kind: loop_facts.kind,
                        boundary_half_edges: loop_facts.boundary_half_edges.clone(),
                        boundary_vertices: loop_facts.boundary_vertices.clone(),
                    })
                    .collect(),
            })
            .collect();

        let trim_endpoint_pairs = knowledge
            .trim_intents
            .iter()
            .filter(|trim| trim.construction == TrimEdgeConstruction::NewExactGeneratorTrim)
            .filter_map(|trim| {
                let support_face = trim.owner_pair.support_face()?;
                Some(SideSupportTrimEndpointFacts {
                    fillet_index: trim.fillet_index,
                    support_face,
                    start: trim.start,
                    end: trim.end,
                    start_membership: endpoint_region_membership(
                        knowledge,
                        support_face,
                        trim.start,
                    ),
                    end_membership: endpoint_region_membership(knowledge, support_face, trim.end),
                })
            })
            .collect();

        Self {
            support_families,
            trim_endpoint_pairs,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "region_facts support_families=[{}] endpoint_pairs=[{}]",
            self.support_families
                .iter()
                .map(SideSupportFamilyFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.trim_endpoint_pairs
                .iter()
                .map(SideSupportTrimEndpointFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl SideSupportFamilyFacts {
    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} loops=[{}]",
            self.support_face,
            self.loops
                .iter()
                .map(SideSupportLoopFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

impl SideSupportLoopFacts {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?} vertices={:?}",
            self.face, self.loop_handle, self.kind, self.boundary_vertices
        )
    }
}

impl SideSupportTrimEndpointFacts {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] support={:?} {:?}->{:?} start={} end={}",
            self.fillet_index,
            self.support_face,
            self.start,
            self.end,
            self.start_membership.diagnostic_label(),
            self.end_membership.diagnostic_label()
        )
    }
}

impl EndpointRegionMembership {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            EndpointRegionMembership::SameSupportLoop {
                face,
                loop_handle,
                half_edge,
            } => {
                format!("support-loop({face:?}/{loop_handle:?}@{half_edge:?})")
            }
            EndpointRegionMembership::OutsideSupportFamily {
                face,
                loop_handle,
                half_edge,
            } => {
                format!("outside-support-family({face:?}/{loop_handle:?}@{half_edge:?})")
            }
            EndpointRegionMembership::Ambiguous {
                containing_loop_count,
            } => {
                format!("ambiguous({containing_loop_count} loops)")
            }
            EndpointRegionMembership::Unknown => "unknown".to_string(),
        }
    }
}

fn endpoint_region_membership(
    knowledge: &TopologyKnowledgeGraph,
    support_face: FaceHandle,
    vertex: VertexHandle,
) -> EndpointRegionMembership {
    let Some(membership) = knowledge.endpoint_membership(vertex) else {
        return EndpointRegionMembership::Unknown;
    };
    if let Some(loop_membership) = membership
        .containing_loops
        .iter()
        .find(|loop_membership| loop_membership.face == support_face)
    {
        return EndpointRegionMembership::SameSupportLoop {
            face: loop_membership.face,
            loop_handle: loop_membership.loop_handle,
            half_edge: loop_membership.half_edge,
        };
    }
    if membership.containing_loops.is_empty() {
        EndpointRegionMembership::Unknown
    } else if membership.containing_loops.len() == 1 {
        let loop_membership = &membership.containing_loops[0];
        EndpointRegionMembership::OutsideSupportFamily {
            face: loop_membership.face,
            loop_handle: loop_membership.loop_handle,
            half_edge: loop_membership.half_edge,
        }
    } else {
        EndpointRegionMembership::Ambiguous {
            containing_loop_count: membership.containing_loops.len(),
        }
    }
}
