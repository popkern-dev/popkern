//! Accessors and diagnostics for side-support selected plans.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPreconditions, BridgeSupportTrimPlan, ExpectedTopologyDelta,
    OutsideOutsideSupportPreconditions, OutsideOutsideSupportTrimPlan, PlannerStage,
    SideSupportRegionRefusal, SideSupportRegionRefusalReason,
};

impl OutsideOutsideSupportTrimPlan {
    pub(super) fn fillet_index(&self) -> usize {
        self.fillet_index
    }

    pub(super) fn cylinder_face(&self) -> FaceHandle {
        self.cylinder_face
    }

    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn start_face(&self) -> FaceHandle {
        self.start_face
    }

    pub(super) fn start_loop(&self) -> LoopHandle {
        self.start_loop
    }

    pub(super) fn start_half_edge(&self) -> HalfEdgeHandle {
        self.start_half_edge
    }

    pub(super) fn end_face(&self) -> FaceHandle {
        self.end_face
    }

    pub(super) fn end_loop(&self) -> LoopHandle {
        self.end_loop
    }

    pub(super) fn end_half_edge(&self) -> HalfEdgeHandle {
        self.end_half_edge
    }

    pub(super) fn start(&self) -> VertexHandle {
        self.start
    }

    pub(super) fn end(&self) -> VertexHandle {
        self.end
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} support={:?} start={:?}/{:?}@{:?} end={:?}/{:?}@{:?} {:?}->{:?}",
            self.fillet_index,
            self.cylinder_face,
            self.support_face,
            self.start_face,
            self.start_loop,
            self.start_half_edge,
            self.end_face,
            self.end_loop,
            self.end_half_edge,
            self.start,
            self.end
        )
    }
}

impl BridgeSupportTrimPlan {
    pub(super) fn fillet_index(&self) -> usize {
        self.fillet_index
    }

    pub(super) fn cylinder_face(&self) -> FaceHandle {
        self.cylinder_face
    }

    pub(super) fn support_loop(&self) -> LoopHandle {
        self.support_loop
    }

    pub(super) fn support_half_edge(&self) -> HalfEdgeHandle {
        self.support_half_edge
    }

    pub(super) fn outside_face(&self) -> FaceHandle {
        self.outside_face
    }

    pub(super) fn outside_loop(&self) -> LoopHandle {
        self.outside_loop
    }

    pub(super) fn outside_half_edge(&self) -> HalfEdgeHandle {
        self.outside_half_edge
    }

    pub(super) fn start(&self) -> VertexHandle {
        self.start
    }

    pub(super) fn end(&self) -> VertexHandle {
        self.end
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder={:?} support_loop={:?}@{:?} outside={:?}/{:?}@{:?} {:?}->{:?}",
            self.fillet_index,
            self.cylinder_face,
            self.support_loop,
            self.support_half_edge,
            self.outside_face,
            self.outside_loop,
            self.outside_half_edge,
            self.start,
            self.end
        )
    }
}

impl ExpectedTopologyDelta {
    pub(super) fn zero() -> Self {
        Self::Known {
            vertices: 0,
            edges: 0,
            faces: 0,
            loops: 0,
        }
    }

    pub(super) fn plus(self, other: Self) -> Self {
        match (self, other) {
            (
                Self::Known {
                    vertices,
                    edges,
                    faces,
                    loops,
                },
                Self::Known {
                    vertices: other_vertices,
                    edges: other_edges,
                    faces: other_faces,
                    loops: other_loops,
                },
            ) => Self::Known {
                vertices: vertices + other_vertices,
                edges: edges + other_edges,
                faces: faces + other_faces,
                loops: loops + other_loops,
            },
            (Self::Unknown { reason }, _) | (_, Self::Unknown { reason }) => {
                Self::Unknown { reason }
            }
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            Self::Known {
                vertices,
                edges,
                faces,
                loops,
            } => format!("dV={vertices} dE={edges} dF={faces} dL={loops}"),
            Self::Unknown { reason } => format!("unknown reason={reason}"),
        }
    }
}

impl BridgeSupportRegionPreconditions {
    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "one_endpoint_on_support_loop={} other_endpoint_known={} unique_support_owner={} face_family_bridge_operands_known={}",
            self.exactly_one_endpoint_on_support_loop,
            self.other_endpoint_has_known_membership,
            self.support_owner_is_unique,
            self.face_family_bridge_operands_known
        )
    }
}

impl OutsideOutsideSupportPreconditions {
    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "both_endpoints_outside_support_family={} endpoint_membership_known={} unique_support_owner={} topology_operator_required={}",
            self.both_endpoints_outside_support_family,
            self.endpoint_membership_known,
            self.support_owner_is_unique,
            self.topology_operator_required
        )
    }
}

impl SideSupportRegionRefusal {
    pub(super) fn stage(&self) -> PlannerStage {
        self.stage
    }

    pub(super) fn reason(&self) -> SideSupportRegionRefusalReason {
        self.reason
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!("refusal stage={:?} reason={:?}", self.stage, self.reason)
    }
}
