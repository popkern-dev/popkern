//! Live BRep operand preflight for selected side-support plans.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPlan, BridgeSupportTrimPlan, PlannedTopologyOperation,
    SelectedSideSupportRegionPlan, SideSupportRegionPlanSelection,
};

#[derive(Clone, Debug)]
pub(super) enum SideSupportRegionLivePreconditionFailure {
    EndpointHalfEdgeMismatch {
        support_half_edge: HalfEdgeHandle,
        outside_half_edge: HalfEdgeHandle,
        start: VertexHandle,
        end: VertexHandle,
    },
    MissingFace(FaceHandle),
    MissingHalfEdge(HalfEdgeHandle),
    MissingLoop(LoopHandle),
    MissingVertex(VertexHandle),
    HalfEdgeLoopMismatch {
        half_edge: HalfEdgeHandle,
        expected_loop: LoopHandle,
        actual_loop: LoopHandle,
    },
}

impl SideSupportRegionLivePreconditionFailure {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionLivePreconditionFailure::EndpointHalfEdgeMismatch {
                support_half_edge,
                outside_half_edge,
                start,
                end,
            } => {
                format!(
                    "endpoint-half-edge-mismatch support={support_half_edge:?} outside={outside_half_edge:?} endpoints={start:?}->{end:?}"
                )
            }
            SideSupportRegionLivePreconditionFailure::MissingFace(face) => {
                format!("missing-face {face:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingHalfEdge(half_edge) => {
                format!("missing-half-edge {half_edge:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingLoop(loop_handle) => {
                format!("missing-loop {loop_handle:?}")
            }
            SideSupportRegionLivePreconditionFailure::MissingVertex(vertex) => {
                format!("missing-vertex {vertex:?}")
            }
            SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                half_edge,
                expected_loop,
                actual_loop,
            } => {
                format!(
                    "half-edge-loop-mismatch {half_edge:?} expected={expected_loop:?} actual={actual_loop:?}"
                )
            }
        }
    }
}

impl SideSupportRegionPlanSelection {
    pub(super) fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        match self {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => plan.validate_live_operands(brep),
            SideSupportRegionPlanSelection::Refused(_) => Ok(()),
        }
    }
}

impl SelectedSideSupportRegionPlan {
    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        for operation in self.operations() {
            operation.validate_live_operands(brep)?;
        }
        Ok(())
    }
}

impl PlannedTopologyOperation {
    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        match self {
            PlannedTopologyOperation::SplitFaceWithEdge { .. } => Ok(()),
            PlannedTopologyOperation::BridgeSupportRegion(plan) => {
                plan.validate_live_operands(brep)
            }
        }
    }
}

impl BridgeSupportRegionPlan {
    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        if !brep.faces.contains_key(self.support_face()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.support_face(),
            ));
        }
        for trim in self.trims() {
            trim.validate_live_operands(brep)?;
        }
        Ok(())
    }
}

impl BridgeSupportTrimPlan {
    fn validate_live_operands(
        &self,
        brep: &BRep,
    ) -> Result<(), SideSupportRegionLivePreconditionFailure> {
        if !brep.loops.contains_key(self.support_loop()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.support_loop(),
            ));
        }
        let Some(support_half_edge) = brep.half_edges.get(self.support_half_edge()) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.support_half_edge(),
            ));
        };
        if support_half_edge.loop_ != self.support_loop() {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.support_half_edge(),
                    expected_loop: self.support_loop(),
                    actual_loop: support_half_edge.loop_,
                },
            );
        }
        if !brep.faces.contains_key(self.outside_face()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.outside_face(),
            ));
        }
        if !brep.loops.contains_key(self.outside_loop()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.outside_loop(),
            ));
        }
        let Some(outside_half_edge) = brep.half_edges.get(self.outside_half_edge()) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.outside_half_edge(),
            ));
        };
        if outside_half_edge.loop_ != self.outside_loop() {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.outside_half_edge(),
                    expected_loop: self.outside_loop(),
                    actual_loop: outside_half_edge.loop_,
                },
            );
        }
        if !brep.vertices.contains_key(self.start()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.start(),
            ));
        }
        if !brep.vertices.contains_key(self.end()) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.end(),
            ));
        }
        let support_origin = support_half_edge.origin;
        let outside_origin = outside_half_edge.origin;
        if !((support_origin == self.start() && outside_origin == self.end())
            || (support_origin == self.end() && outside_origin == self.start()))
        {
            return Err(
                SideSupportRegionLivePreconditionFailure::EndpointHalfEdgeMismatch {
                    support_half_edge: self.support_half_edge(),
                    outside_half_edge: self.outside_half_edge(),
                    start: self.start(),
                    end: self.end(),
                },
            );
        }
        Ok(())
    }
}
