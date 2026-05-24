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
use super::repeated_apex_side_support_result::SideSupportRegionContinuationResult;
#[path = "repeated_apex_side_support_bridge_aggregate.rs"]
mod bridge_aggregate;
#[path = "repeated_apex_side_support_bridge_loop_topology.rs"]
mod bridge_loop_topology;
#[path = "repeated_apex_side_support_bridge_plan_access.rs"]
mod bridge_plan_access;

use bridge_aggregate::select_aggregate_region_cycles;
use bridge_loop_topology::{loop_intervals_between_origins, loops_share_edges};
use bridge_plan_access::{full_loop_interval_from, shorter_interval, two_bridge_trims};

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgePrimitiveDerivation {
    bridge_regions: usize,
    ready_regions: usize,
    primitive: Option<FaceFamilyBridgePrimitivePlan>,
    details: Vec<FaceFamilyBridgeRegionDerivation>,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgePrimitivePlan {
    regions: Vec<FaceFamilyBridgeRegionPrimitivePlan>,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeRegionPrimitivePlan {
    support_face: FaceHandle,
    support_loop: LoopHandle,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    trims: [BridgeSupportTrimPlan; 2],
    support_paths: FaceFamilyBridgeLoopIntervals,
    outside_paths: FaceFamilyBridgeLoopIntervals,
    boundary_cycles: Vec<FaceFamilyBridgeBoundaryCycle>,
    boundary_exposure: FaceFamilyBridgeBoundaryExposurePlan,
    boundary_readiness: FaceFamilyBridgeBoundaryReadiness,
    expected_edges: usize,
    expected_faces: usize,
    expected_loops: usize,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeBoundaryCycle {
    support_interval: Vec<HalfEdgeHandle>,
    outside_interval: Vec<HalfEdgeHandle>,
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeLoopIntervals {
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

    pub(super) fn apply_topology_operator(
        &self,
        brep: &mut BRep,
    ) -> Result<SideSupportRegionContinuationResult, String> {
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
                PlannedTopologyOperation::SplitFaceWithEdge { .. }
                | PlannedTopologyOperation::OutsideOutsideSupportContinuation(_) => None,
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
    pub(super) fn regions(&self) -> &[FaceFamilyBridgeRegionPrimitivePlan] {
        &self.regions
    }

    fn diagnostic_label(&self) -> String {
        let edges = self
            .regions()
            .iter()
            .map(|region| region.expected_edges)
            .sum::<usize>();
        let faces = self
            .regions()
            .iter()
            .map(|region| region.expected_faces)
            .sum::<usize>();
        let loops = self
            .regions()
            .iter()
            .map(|region| region.expected_loops)
            .sum::<usize>();
        format!(
            "regions={} expected_delta=dV=0 dE={edges} dF={faces} dL={loops}",
            self.regions().len()
        )
    }
}

impl FaceFamilyBridgeRegionDerivation {
    fn build(region: &BridgeSupportRegionPlan, brep: &BRep) -> Self {
        let Some([first, second]) = two_bridge_trims(region.trims()) else {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::UnsupportedTrimCount,
            );
        };
        if first.support_loop() != second.support_loop() {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::MixedSupportLoops,
            );
        }
        if first.outside_face() != second.outside_face()
            || first.outside_loop() != second.outside_loop()
        {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::MixedOutsideLoops,
            );
        }
        let Some(support_paths) = loop_intervals_between_origins(
            brep,
            first.support_loop(),
            first.support_half_edge(),
            second.support_half_edge(),
        ) else {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::MissingSupportPath,
            );
        };
        let Some(outside_paths) = loop_intervals_between_origins(
            brep,
            first.outside_loop(),
            first.outside_half_edge(),
            second.outside_half_edge(),
        ) else {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::MissingOutsidePath,
            );
        };
        if first.support_half_edge() == second.support_half_edge() {
            if !outside_paths.are_complementary() {
                return Self::refused(
                    region.support_face(),
                    FaceFamilyBridgeOperandGap::NonComplementaryOutsideIntervals,
                );
            }
            let Some(support_interval) =
                full_loop_interval_from(brep, first.support_loop(), first.support_half_edge())
            else {
                return Self::refused(
                    region.support_face(),
                    FaceFamilyBridgeOperandGap::MissingSupportPath,
                );
            };
            let outside_interval = shorter_interval(
                outside_paths.first_to_second(),
                outside_paths.second_to_first(),
            );
            if outside_interval.len() < 2 {
                return Self::refused(
                    region.support_face(),
                    FaceFamilyBridgeOperandGap::MissingOutsidePath,
                );
            }
            let boundary_cycles = vec![FaceFamilyBridgeBoundaryCycle {
                support_interval,
                outside_interval,
            }];
            let boundary_exposure =
                FaceFamilyBridgeBoundaryExposurePlan::build(brep, &boundary_cycles);
            let boundary_readiness = boundary_exposure.readiness().clone();
            let expected_edges = boundary_cycles.len() * 2
                + boundary_exposure.fragment_retire_count()
                + boundary_exposure.self_owned_exposure_count();
            let expected_faces = boundary_cycles.len();
            let primitive = FaceFamilyBridgeRegionPrimitivePlan {
                support_face: region.support_face(),
                support_loop: first.support_loop(),
                outside_face: first.outside_face(),
                outside_loop: first.outside_loop(),
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
            return Self {
                support_face: region.support_face(),
                primitive: Some(primitive),
                reason: None,
            };
        }
        if !support_paths.are_complementary() {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::NonComplementarySupportIntervals,
            );
        }
        if !outside_paths.are_complementary() {
            return Self::refused(
                region.support_face(),
                FaceFamilyBridgeOperandGap::NonComplementaryOutsideIntervals,
            );
        }
        if loops_share_edges(brep, first.support_loop(), first.outside_loop()) {
            return Self::refused(
                region.support_face(),
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
        let boundary_readiness = boundary_exposure.readiness().clone();
        let expected_edges = boundary_cycles.len() * 2
            + boundary_exposure.fragment_retire_count()
            + boundary_exposure.self_owned_exposure_count();
        let expected_faces = boundary_cycles.len();
        let primitive = FaceFamilyBridgeRegionPrimitivePlan {
            support_face: region.support_face(),
            support_loop: first.support_loop(),
            outside_face: first.outside_face(),
            outside_loop: first.outside_loop(),
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
            support_face: region.support_face(),
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
