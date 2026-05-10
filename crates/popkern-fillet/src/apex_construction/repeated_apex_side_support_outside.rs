//! Outside/outside support-continuation primitive derivation.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportContinuationPlan, PlannedTopologyOperation, SelectedSideSupportRegionPlan,
    SideSupportRegionPlanSelection,
};
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportPrimitiveDerivation {
    regions: usize,
    ready_regions: usize,
    details: Vec<OutsideOutsideSupportRegionDerivation>,
}

#[derive(Clone, Debug)]
struct OutsideOutsideSupportRegionDerivation {
    region: OutsideOutsideSupportContinuationPlan,
    topology_class: Option<CylinderMediatedOutsideOutsideClass>,
    reason: OutsideOutsideSupportRefusal,
}

#[derive(Clone, Debug)]
enum OutsideOutsideSupportRefusal {
    NoOutsideCylinderAdjacency,
    NoSupportCylinderAdjacency,
    TopologyOperatorNotImplemented,
}

#[derive(Clone, Debug)]
struct CylinderMediatedOutsideOutsideClass {
    cylinder_face: FaceHandle,
    outside_adjacencies: Vec<LoopAdjacencyOperand>,
    support_adjacencies: Vec<LoopAdjacencyOperand>,
}

#[derive(Clone, Debug)]
struct LoopAdjacencyOperand {
    source_loop: LoopHandle,
    source_half_edge: HalfEdgeHandle,
    cylinder_loop: LoopHandle,
    cylinder_half_edge: HalfEdgeHandle,
}

impl OutsideOutsideSupportPrimitiveDerivation {
    pub(super) fn build(
        selection: &SideSupportRegionPlanSelection,
        knowledge: &TopologyKnowledgeGraph,
    ) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => {
                Self::from_selected_plan(plan, knowledge)
            }
            SideSupportRegionPlanSelection::Refused(_) => Self {
                regions: 0,
                ready_regions: 0,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "outside-outside-support primitive_ready={} regions={} ready_regions={} details=[{}]",
            self.is_ready(),
            self.regions,
            self.ready_regions,
            self.details
                .iter()
                .map(OutsideOutsideSupportRegionDerivation::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }

    fn is_ready(&self) -> bool {
        self.regions > 0 && self.ready_regions == self.regions
    }

    fn from_selected_plan(
        plan: &SelectedSideSupportRegionPlan,
        knowledge: &TopologyKnowledgeGraph,
    ) -> Self {
        let details = plan
            .operations()
            .filter_map(|operation| match operation {
                PlannedTopologyOperation::OutsideOutsideSupportContinuation(region) => Some(
                    OutsideOutsideSupportRegionDerivation::derive(region, knowledge),
                ),
                PlannedTopologyOperation::SplitFaceWithEdge { .. }
                | PlannedTopologyOperation::BridgeSupportRegion(_) => None,
            })
            .collect::<Vec<_>>();
        Self {
            regions: details.len(),
            ready_regions: 0,
            details,
        }
    }
}

impl OutsideOutsideSupportRegionDerivation {
    fn derive(
        region: &OutsideOutsideSupportContinuationPlan,
        knowledge: &TopologyKnowledgeGraph,
    ) -> Self {
        let mut outside_adjacencies = Vec::new();
        let mut support_adjacencies = Vec::new();
        let mut cylinder_face = None;
        for trim in region.trims() {
            cylinder_face = Some(trim.cylinder_face());
            outside_adjacencies.extend(outside_cylinder_adjacencies(
                knowledge,
                trim.start_face(),
                trim.start_loop(),
                trim.cylinder_face(),
            ));
            outside_adjacencies.extend(outside_cylinder_adjacencies(
                knowledge,
                trim.end_face(),
                trim.end_loop(),
                trim.cylinder_face(),
            ));
            support_adjacencies.extend(support_cylinder_adjacencies(
                knowledge,
                trim.support_face(),
                trim.cylinder_face(),
            ));
        }
        outside_adjacencies.sort_by_key(LoopAdjacencyOperand::sort_key);
        outside_adjacencies.dedup_by_key(|operand| operand.sort_key());
        support_adjacencies.sort_by_key(LoopAdjacencyOperand::sort_key);
        support_adjacencies.dedup_by_key(|operand| operand.sort_key());

        let topology_class =
            cylinder_face.map(|cylinder_face| CylinderMediatedOutsideOutsideClass {
                cylinder_face,
                outside_adjacencies: outside_adjacencies.clone(),
                support_adjacencies: support_adjacencies.clone(),
            });
        let reason = if outside_adjacencies.is_empty() {
            OutsideOutsideSupportRefusal::NoOutsideCylinderAdjacency
        } else if support_adjacencies.is_empty() {
            OutsideOutsideSupportRefusal::NoSupportCylinderAdjacency
        } else {
            OutsideOutsideSupportRefusal::TopologyOperatorNotImplemented
        };
        Self {
            region: region.clone(),
            topology_class,
            reason,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} trims={} {} {}",
            self.region.support_face(),
            self.region.trims().len(),
            self.topology_class
                .as_ref()
                .map(CylinderMediatedOutsideOutsideClass::diagnostic_label)
                .unwrap_or_else(|| "class=unknown".to_string()),
            self.reason.diagnostic_label()
        )
    }
}

impl OutsideOutsideSupportRefusal {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            OutsideOutsideSupportRefusal::NoOutsideCylinderAdjacency => {
                "missing-outside-cylinder-adjacency"
            }
            OutsideOutsideSupportRefusal::NoSupportCylinderAdjacency => {
                "missing-support-cylinder-adjacency"
            }
            OutsideOutsideSupportRefusal::TopologyOperatorNotImplemented => {
                "missing-topology-operator"
            }
        }
    }
}

impl CylinderMediatedOutsideOutsideClass {
    fn diagnostic_label(&self) -> String {
        format!(
            "class=cylinder-mediated cylinder={:?} outside_adjacencies=[{}] support_adjacencies=[{}]",
            self.cylinder_face,
            self.outside_adjacencies
                .iter()
                .map(LoopAdjacencyOperand::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", "),
            self.support_adjacencies
                .iter()
                .map(LoopAdjacencyOperand::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

impl LoopAdjacencyOperand {
    fn sort_key(&self) -> (LoopHandle, HalfEdgeHandle, LoopHandle, HalfEdgeHandle) {
        (
            self.source_loop,
            self.source_half_edge,
            self.cylinder_loop,
            self.cylinder_half_edge,
        )
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}@{:?}<->{:?}@{:?}",
            self.source_loop, self.source_half_edge, self.cylinder_loop, self.cylinder_half_edge
        )
    }
}

fn outside_cylinder_adjacencies(
    knowledge: &TopologyKnowledgeGraph,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    cylinder_face: FaceHandle,
) -> Vec<LoopAdjacencyOperand> {
    knowledge
        .loop_adjacencies()
        .filter_map(|adjacency| {
            if adjacency.face() == outside_face
                && adjacency.loop_handle() == outside_loop
                && adjacency.adjacent_face() == cylinder_face
            {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.loop_handle(),
                    source_half_edge: adjacency.half_edge(),
                    cylinder_loop: adjacency.adjacent_loop(),
                    cylinder_half_edge: adjacency.adjacent_half_edge(),
                });
            }
            if adjacency.adjacent_face() == outside_face
                && adjacency.adjacent_loop() == outside_loop
                && adjacency.face() == cylinder_face
            {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.adjacent_loop(),
                    source_half_edge: adjacency.adjacent_half_edge(),
                    cylinder_loop: adjacency.loop_handle(),
                    cylinder_half_edge: adjacency.half_edge(),
                });
            }
            None
        })
        .collect()
}

fn support_cylinder_adjacencies(
    knowledge: &TopologyKnowledgeGraph,
    support_face: FaceHandle,
    cylinder_face: FaceHandle,
) -> Vec<LoopAdjacencyOperand> {
    knowledge
        .loop_adjacencies()
        .filter_map(|adjacency| {
            if adjacency.face() == support_face && adjacency.adjacent_face() == cylinder_face {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.loop_handle(),
                    source_half_edge: adjacency.half_edge(),
                    cylinder_loop: adjacency.adjacent_loop(),
                    cylinder_half_edge: adjacency.adjacent_half_edge(),
                });
            }
            if adjacency.adjacent_face() == support_face && adjacency.face() == cylinder_face {
                return Some(LoopAdjacencyOperand {
                    source_loop: adjacency.adjacent_loop(),
                    source_half_edge: adjacency.adjacent_half_edge(),
                    cylinder_loop: adjacency.loop_handle(),
                    cylinder_half_edge: adjacency.half_edge(),
                });
            }
            None
        })
        .collect()
}
