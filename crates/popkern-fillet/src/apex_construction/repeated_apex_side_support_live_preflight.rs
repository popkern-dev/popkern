//! Live BRep operand preflight for selected side-support plans.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPlan, BridgeSupportTrimPlan, PlannedTopologyOperation,
    SelectedSideSupportRegionPlan, SideSupportRegionPlanSelection,
};

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeOperandSummary {
    bridge_regions: usize,
    ready_regions: usize,
    details: Vec<FaceFamilyBridgeRegionSummary>,
}

#[derive(Clone, Debug)]
struct FaceFamilyBridgeRegionSummary {
    support_face: FaceHandle,
    support_path_len: Option<usize>,
    outside_path_len: Option<usize>,
    reason: Option<FaceFamilyBridgeOperandGap>,
}

#[derive(Clone, Debug)]
enum FaceFamilyBridgeOperandGap {
    UnsupportedTrimCount,
    MixedSupportLoops,
    MixedOutsideLoops,
    MissingSupportPath,
    MissingOutsidePath,
}

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

impl FaceFamilyBridgeOperandSummary {
    pub(super) fn build(selection: &SideSupportRegionPlanSelection, brep: &BRep) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => {
                Self::from_selected_plan(plan, brep)
            }
            SideSupportRegionPlanSelection::Refused(_) => Self {
                bridge_regions: 0,
                ready_regions: 0,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "face-family-bridge operands bridge_regions={} ready_regions={} details=[{}]",
            self.bridge_regions,
            self.ready_regions,
            self.details
                .iter()
                .map(FaceFamilyBridgeRegionSummary::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }

    fn from_selected_plan(plan: &SelectedSideSupportRegionPlan, brep: &BRep) -> Self {
        let details = plan
            .operations
            .iter()
            .filter_map(|operation| match operation {
                PlannedTopologyOperation::BridgeSupportRegion(region) => {
                    Some(FaceFamilyBridgeRegionSummary::build(region, brep))
                }
                PlannedTopologyOperation::SplitFaceWithEdge { .. } => None,
            })
            .collect::<Vec<_>>();
        let bridge_regions = details.len();
        let ready_regions = details
            .iter()
            .filter(|detail| detail.reason.is_none())
            .count();
        Self {
            bridge_regions,
            ready_regions,
            details,
        }
    }
}

impl FaceFamilyBridgeRegionSummary {
    fn build(region: &BridgeSupportRegionPlan, brep: &BRep) -> Self {
        let Some([first, second]) = two_bridge_trims(&region.trims) else {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::UnsupportedTrimCount,
            );
        };
        if first.support_loop != second.support_loop {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::MixedSupportLoops,
            );
        }
        if first.outside_face != second.outside_face || first.outside_loop != second.outside_loop {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::MixedOutsideLoops,
            );
        }
        let support_path = loop_path_between_origins(
            brep,
            first.support_loop,
            first.support_half_edge,
            second.support_half_edge,
        );
        let outside_path = loop_path_between_origins(
            brep,
            first.outside_loop,
            first.outside_half_edge,
            second.outside_half_edge,
        );
        let reason = if support_path.is_none() {
            Some(FaceFamilyBridgeOperandGap::MissingSupportPath)
        } else if outside_path.is_none() {
            Some(FaceFamilyBridgeOperandGap::MissingOutsidePath)
        } else {
            None
        };
        Self {
            support_face: region.support_face,
            support_path_len: support_path.as_ref().map(Vec::len),
            outside_path_len: outside_path.as_ref().map(Vec::len),
            reason,
        }
    }

    fn refused(support_face: FaceHandle, reason: FaceFamilyBridgeOperandGap) -> Self {
        Self {
            support_face,
            support_path_len: None,
            outside_path_len: None,
            reason: Some(reason),
        }
    }

    fn diagnostic_label(&self) -> String {
        let reason = self
            .reason
            .as_ref()
            .map(FaceFamilyBridgeOperandGap::diagnostic_label)
            .unwrap_or("ready");
        format!(
            "support={:?} support_path_len={:?} outside_path_len={:?} {reason}",
            self.support_face, self.support_path_len, self.outside_path_len
        )
    }
}

impl FaceFamilyBridgeOperandGap {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            FaceFamilyBridgeOperandGap::UnsupportedTrimCount => "unsupported-trim-count",
            FaceFamilyBridgeOperandGap::MixedSupportLoops => "mixed-support-loops",
            FaceFamilyBridgeOperandGap::MixedOutsideLoops => "mixed-outside-loops",
            FaceFamilyBridgeOperandGap::MissingSupportPath => "missing-support-path",
            FaceFamilyBridgeOperandGap::MissingOutsidePath => "missing-outside-path",
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
        for operation in &self.operations {
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
        if !brep.faces.contains_key(self.support_face) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.support_face,
            ));
        }
        for trim in &self.trims {
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
        if !brep.loops.contains_key(self.support_loop) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.support_loop,
            ));
        }
        let Some(support_half_edge) = brep.half_edges.get(self.support_half_edge) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.support_half_edge,
            ));
        };
        if support_half_edge.loop_ != self.support_loop {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.support_half_edge,
                    expected_loop: self.support_loop,
                    actual_loop: support_half_edge.loop_,
                },
            );
        }
        if !brep.faces.contains_key(self.outside_face) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingFace(
                self.outside_face,
            ));
        }
        if !brep.loops.contains_key(self.outside_loop) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingLoop(
                self.outside_loop,
            ));
        }
        let Some(outside_half_edge) = brep.half_edges.get(self.outside_half_edge) else {
            return Err(SideSupportRegionLivePreconditionFailure::MissingHalfEdge(
                self.outside_half_edge,
            ));
        };
        if outside_half_edge.loop_ != self.outside_loop {
            return Err(
                SideSupportRegionLivePreconditionFailure::HalfEdgeLoopMismatch {
                    half_edge: self.outside_half_edge,
                    expected_loop: self.outside_loop,
                    actual_loop: outside_half_edge.loop_,
                },
            );
        }
        if !brep.vertices.contains_key(self.start) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.start,
            ));
        }
        if !brep.vertices.contains_key(self.end) {
            return Err(SideSupportRegionLivePreconditionFailure::MissingVertex(
                self.end,
            ));
        }
        let support_origin = support_half_edge.origin;
        let outside_origin = outside_half_edge.origin;
        if !((support_origin == self.start && outside_origin == self.end)
            || (support_origin == self.end && outside_origin == self.start))
        {
            return Err(
                SideSupportRegionLivePreconditionFailure::EndpointHalfEdgeMismatch {
                    support_half_edge: self.support_half_edge,
                    outside_half_edge: self.outside_half_edge,
                    start: self.start,
                    end: self.end,
                },
            );
        }
        Ok(())
    }
}

fn two_bridge_trims(trims: &[BridgeSupportTrimPlan]) -> Option<[&BridgeSupportTrimPlan; 2]> {
    match trims {
        [first, second] => Some([first, second]),
        _ => None,
    }
}

fn loop_path_between_origins(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
    end: HalfEdgeHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let loop_data = brep.loops.get(loop_handle)?;
    let mut cursor = start;
    let mut path = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        path.push(cursor);
        if cursor == end {
            return Some(path);
        }
        cursor = half_edge.next;
        if cursor == start || cursor == loop_data.half_edge && path.len() > brep.half_edges.len() {
            break;
        }
    }
    None
}
