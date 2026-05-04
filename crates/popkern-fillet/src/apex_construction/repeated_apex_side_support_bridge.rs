//! Face-family bridge primitive operands for repeated-apex side-support plans.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPlan, BridgeSupportTrimPlan, PlannedTopologyOperation,
    SelectedSideSupportRegionPlan, SideSupportRegionPlanSelection,
};

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgePrimitiveDerivation {
    bridge_regions: usize,
    ready_regions: usize,
    primitive: Option<FaceFamilyBridgePrimitivePlan>,
    details: Vec<FaceFamilyBridgeRegionDerivation>,
}

#[derive(Clone, Debug)]
struct FaceFamilyBridgePrimitivePlan {
    regions: Vec<FaceFamilyBridgeRegionPrimitivePlan>,
}

#[derive(Clone, Debug)]
struct FaceFamilyBridgeRegionPrimitivePlan {
    support_face: FaceHandle,
    support_loop: LoopHandle,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    trims: [BridgeSupportTrimPlan; 2],
    support_path: Vec<HalfEdgeHandle>,
    outside_path: Vec<HalfEdgeHandle>,
    expected_edges: usize,
    expected_faces: usize,
    expected_loops: usize,
}

#[derive(Clone, Debug)]
struct FaceFamilyBridgeRegionDerivation {
    support_face: FaceHandle,
    primitive: Option<FaceFamilyBridgeRegionPrimitivePlan>,
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

impl FaceFamilyBridgePrimitiveDerivation {
    pub(super) fn build(selection: &SideSupportRegionPlanSelection, brep: &BRep) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => {
                Self::from_selected_plan(plan, brep)
            }
            SideSupportRegionPlanSelection::Refused(_) => Self {
                bridge_regions: 0,
                ready_regions: 0,
                primitive: None,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn is_complete(&self) -> bool {
        self.bridge_regions > 0 && self.ready_regions == self.bridge_regions
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "face-family-bridge primitive_ready={} bridge_regions={} ready_regions={} totals={} details=[{}]",
            self.is_complete(),
            self.bridge_regions,
            self.ready_regions,
            self.primitive
                .as_ref()
                .map(FaceFamilyBridgePrimitivePlan::diagnostic_label)
                .unwrap_or_else(|| "none".to_string()),
            self.details
                .iter()
                .map(FaceFamilyBridgeRegionDerivation::diagnostic_label)
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
                    Some(FaceFamilyBridgeRegionDerivation::build(region, brep))
                }
                PlannedTopologyOperation::SplitFaceWithEdge { .. } => None,
            })
            .collect::<Vec<_>>();
        let bridge_regions = details.len();
        let ready_regions = details
            .iter()
            .filter(|detail| detail.reason.is_none())
            .count();
        let primitive = if bridge_regions > 0 && ready_regions == bridge_regions {
            Some(FaceFamilyBridgePrimitivePlan {
                regions: details
                    .iter()
                    .filter_map(|detail| detail.primitive.clone())
                    .collect(),
            })
        } else {
            None
        };
        Self {
            bridge_regions,
            ready_regions,
            primitive,
            details,
        }
    }
}

impl FaceFamilyBridgePrimitivePlan {
    fn diagnostic_label(&self) -> String {
        let edges = self
            .regions
            .iter()
            .map(|region| region.expected_edges)
            .sum::<usize>();
        let faces = self
            .regions
            .iter()
            .map(|region| region.expected_faces)
            .sum::<usize>();
        let loops = self
            .regions
            .iter()
            .map(|region| region.expected_loops)
            .sum::<usize>();
        format!(
            "regions={} expected_delta=dV=0 dE={edges} dF={faces} dL={loops}",
            self.regions.len()
        )
    }
}

impl FaceFamilyBridgeRegionDerivation {
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
        let Some(support_path) = loop_path_between_origins(
            brep,
            first.support_loop,
            first.support_half_edge,
            second.support_half_edge,
        ) else {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::MissingSupportPath,
            );
        };
        let Some(outside_path) = loop_path_between_origins(
            brep,
            first.outside_loop,
            first.outside_half_edge,
            second.outside_half_edge,
        ) else {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::MissingOutsidePath,
            );
        };
        let primitive = FaceFamilyBridgeRegionPrimitivePlan {
            support_face: region.support_face,
            support_loop: first.support_loop,
            outside_face: first.outside_face,
            outside_loop: first.outside_loop,
            trims: [first.clone(), second.clone()],
            support_path,
            outside_path,
            expected_edges: 2,
            expected_faces: 2,
            expected_loops: 2,
        };
        Self {
            support_face: region.support_face,
            primitive: Some(primitive),
            reason: None,
        }
    }

    fn refused(support_face: FaceHandle, reason: FaceFamilyBridgeOperandGap) -> Self {
        Self {
            support_face,
            primitive: None,
            reason: Some(reason),
        }
    }

    fn diagnostic_label(&self) -> String {
        match (&self.primitive, &self.reason) {
            (Some(primitive), None) => primitive.diagnostic_label(),
            (None, Some(reason)) => format!(
                "support={:?} {}",
                self.support_face,
                reason.diagnostic_label()
            ),
            _ => format!(
                "support={:?} inconsistent-bridge-derivation",
                self.support_face
            ),
        }
    }
}

impl FaceFamilyBridgeRegionPrimitivePlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} support_loop={:?} outside={:?}/{:?} trims={} support_path_len={} outside_path_len={} delta=dV=0 dE={} dF={} dL={}",
            self.support_face,
            self.support_loop,
            self.outside_face,
            self.outside_loop,
            self.trims.len(),
            self.support_path.len(),
            self.outside_path.len(),
            self.expected_edges,
            self.expected_faces,
            self.expected_loops
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
        if cursor == start {
            break;
        }
    }
    None
}
