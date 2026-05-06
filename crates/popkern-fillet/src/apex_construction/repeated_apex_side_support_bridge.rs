//! Face-family bridge primitive operands for repeated-apex side-support plans.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_boundary::{
    FaceFamilyBridgeBoundaryExposurePlan, FaceFamilyBridgeBoundaryReadiness,
};
use super::repeated_apex_side_support_bridge_operator::{
    FaceFamilyBridgeTopologyOperatorReport, apply_topology_operator, topology_operator_report,
};
use super::repeated_apex_side_support_plan::{
    BridgeSupportRegionPlan, BridgeSupportTrimPlan, PlannedTopologyOperation,
    SelectedSideSupportRegionPlan, SideSupportRegionPlanSelection,
};
#[path = "repeated_apex_side_support_bridge_aggregate.rs"]
mod bridge_aggregate;

use bridge_aggregate::select_aggregate_region_cycles;

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgePrimitiveDerivation {
    bridge_regions: usize,
    ready_regions: usize,
    primitive: Option<FaceFamilyBridgePrimitivePlan>,
    details: Vec<FaceFamilyBridgeRegionDerivation>,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgePrimitivePlan {
    pub(super) regions: Vec<FaceFamilyBridgeRegionPrimitivePlan>,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeRegionPrimitivePlan {
    pub(super) support_face: FaceHandle,
    pub(super) support_loop: LoopHandle,
    pub(super) outside_face: FaceHandle,
    pub(super) outside_loop: LoopHandle,
    pub(super) trims: [BridgeSupportTrimPlan; 2],
    support_paths: FaceFamilyBridgeLoopIntervals,
    outside_paths: FaceFamilyBridgeLoopIntervals,
    pub(super) boundary_cycles: Vec<FaceFamilyBridgeBoundaryCycle>,
    pub(super) boundary_exposure: FaceFamilyBridgeBoundaryExposurePlan,
    pub(super) boundary_readiness: FaceFamilyBridgeBoundaryReadiness,
    pub(super) expected_edges: usize,
    pub(super) expected_faces: usize,
    pub(super) expected_loops: usize,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeBoundaryCycle {
    pub(super) support_interval: Vec<HalfEdgeHandle>,
    pub(super) outside_interval: Vec<HalfEdgeHandle>,
}

#[derive(Clone, Debug)]
struct FaceFamilyBridgeLoopIntervals {
    first_to_second: Vec<HalfEdgeHandle>,
    second_to_first: Vec<HalfEdgeHandle>,
    loop_len: usize,
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
    NonComplementarySupportIntervals,
    NonComplementaryOutsideIntervals,
    SharedBridgeOperandEdges,
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

    pub(super) fn topology_operator_report(
        &self,
        brep: &BRep,
    ) -> FaceFamilyBridgeTopologyOperatorReport {
        topology_operator_report(self.primitive.as_ref(), brep)
    }

    pub(super) fn apply_topology_operator(&self, brep: &mut BRep) -> Result<(), String> {
        let primitive = self
            .primitive
            .as_ref()
            .ok_or_else(|| "face-family bridge primitive is incomplete".to_string())?;
        apply_topology_operator(primitive, brep)
    }

    fn from_selected_plan(plan: &SelectedSideSupportRegionPlan, brep: &BRep) -> Self {
        let details = plan
            .operations()
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
            let regions = details
                .iter()
                .filter_map(|detail| detail.primitive.clone())
                .collect::<Vec<_>>();
            Some(FaceFamilyBridgePrimitivePlan {
                regions: select_aggregate_region_cycles(brep, regions),
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
        let Some(support_paths) = loop_intervals_between_origins(
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
        if !support_paths.are_complementary() {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::NonComplementarySupportIntervals,
            );
        }
        let Some(outside_paths) = loop_intervals_between_origins(
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
        if !outside_paths.are_complementary() {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::NonComplementaryOutsideIntervals,
            );
        }
        if loops_share_edges(brep, first.support_loop, first.outside_loop) {
            return Self::refused(
                region.support_face,
                FaceFamilyBridgeOperandGap::SharedBridgeOperandEdges,
            );
        }
        let boundary_cycles = vec![
            FaceFamilyBridgeBoundaryCycle {
                support_interval: support_paths.first_to_second.clone(),
                outside_interval: outside_paths.second_to_first.clone(),
            },
            FaceFamilyBridgeBoundaryCycle {
                support_interval: support_paths.second_to_first.clone(),
                outside_interval: outside_paths.first_to_second.clone(),
            },
        ];
        let boundary_exposure = FaceFamilyBridgeBoundaryExposurePlan::build(brep, &boundary_cycles);
        let boundary_readiness = boundary_exposure.readiness.clone();
        let expected_edges = boundary_cycles.len() * 2
            + boundary_exposure.fragment_retire_count()
            + boundary_exposure.self_owned_exposure_count();
        let expected_faces = boundary_cycles.len();
        let primitive = FaceFamilyBridgeRegionPrimitivePlan {
            support_face: region.support_face,
            support_loop: first.support_loop,
            outside_face: first.outside_face,
            outside_loop: first.outside_loop,
            trims: [first.clone(), second.clone()],
            support_paths,
            outside_paths,
            boundary_cycles,
            boundary_exposure,
            boundary_readiness,
            expected_edges,
            expected_faces,
            expected_loops: expected_faces,
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
    fn with_boundary_cycles(
        mut self,
        brep: &BRep,
        boundary_cycles: Vec<FaceFamilyBridgeBoundaryCycle>,
    ) -> Self {
        let boundary_exposure = FaceFamilyBridgeBoundaryExposurePlan::build(brep, &boundary_cycles);
        self.boundary_readiness = boundary_exposure.readiness.clone();
        self.expected_edges = boundary_cycles.len() * 2
            + boundary_exposure.fragment_retire_count()
            + boundary_exposure.self_owned_exposure_count();
        self.expected_faces = boundary_cycles.len();
        self.expected_loops = boundary_cycles.len();
        self.boundary_cycles = boundary_cycles;
        self.boundary_exposure = boundary_exposure;
        self
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} support_loop={:?} outside={:?}/{:?} trims={} support_paths={} outside_paths={} boundary_readiness={} boundary_exposure={} cycles=[{}] delta=dV=0 dE={} dF={} dL={}",
            self.support_face,
            self.support_loop,
            self.outside_face,
            self.outside_loop,
            self.trims.len(),
            self.support_paths.diagnostic_label(),
            self.outside_paths.diagnostic_label(),
            self.boundary_readiness.diagnostic_label(),
            self.boundary_exposure.diagnostic_label(),
            self.boundary_cycles
                .iter()
                .map(FaceFamilyBridgeBoundaryCycle::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.expected_edges,
            self.expected_faces,
            self.expected_loops
        )
    }
}

impl FaceFamilyBridgeBoundaryCycle {
    fn diagnostic_label(&self) -> String {
        format!(
            "support_interval={} outside_interval={}",
            self.support_interval.len(),
            self.outside_interval.len()
        )
    }
}

impl FaceFamilyBridgeLoopIntervals {
    fn diagnostic_label(&self) -> String {
        format!(
            "loop_len={} first_to_second={} second_to_first={}",
            self.loop_len,
            self.first_to_second.len(),
            self.second_to_first.len()
        )
    }

    fn are_complementary(&self) -> bool {
        self.first_to_second.len() + self.second_to_first.len() == self.loop_len + 2
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
            FaceFamilyBridgeOperandGap::NonComplementarySupportIntervals => {
                "non-complementary-support-intervals"
            }
            FaceFamilyBridgeOperandGap::NonComplementaryOutsideIntervals => {
                "non-complementary-outside-intervals"
            }
            FaceFamilyBridgeOperandGap::SharedBridgeOperandEdges => "shared-bridge-operand-edges",
        }
    }
}

fn two_bridge_trims(trims: &[BridgeSupportTrimPlan]) -> Option<[&BridgeSupportTrimPlan; 2]> {
    match trims {
        [first, second] => Some([first, second]),
        _ => None,
    }
}

fn loop_intervals_between_origins(
    brep: &BRep,
    loop_handle: LoopHandle,
    first: HalfEdgeHandle,
    second: HalfEdgeHandle,
) -> Option<FaceFamilyBridgeLoopIntervals> {
    let loop_len = loop_half_edge_count(brep, loop_handle, first)?;
    Some(FaceFamilyBridgeLoopIntervals {
        first_to_second: loop_path_between_origins(brep, loop_handle, first, second)?,
        second_to_first: loop_path_between_origins(brep, loop_handle, second, first)?,
        loop_len,
    })
}

fn loop_half_edge_count(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
) -> Option<usize> {
    let mut cursor = start;
    let mut count = 0;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        count += 1;
        cursor = half_edge.next;
        if cursor == start {
            return Some(count);
        }
    }
    None
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

fn loops_share_edges(brep: &BRep, first: LoopHandle, second: LoopHandle) -> bool {
    let mut first_edges = std::collections::BTreeSet::new();
    let mut cursor = brep.loops[first].half_edge;
    loop {
        first_edges.insert(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[first].half_edge {
            break;
        }
    }
    let mut cursor = brep.loops[second].half_edge;
    loop {
        if first_edges.contains(&brep.half_edges[cursor].edge) {
            return true;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[second].half_edge {
            break;
        }
    }
    false
}
