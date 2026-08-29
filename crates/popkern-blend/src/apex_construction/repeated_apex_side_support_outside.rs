//! Outside/outside support-continuation primitive derivation.

use popkern_core::error::CoreError;
use popkern_core::euler::{
    CylinderMediatedSupportContinuationPostconditions, CylinderMediatedSupportContinuationRequest,
};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_outside_core::{
    cylinder_interval_candidates, outside_cylinder_adjacencies, support_cylinder_adjacencies,
    validate_core_preflights,
};
use super::repeated_apex_side_support_outside_operator::{
    CylinderMediatedSupportContinuationTopologyOperatorReport, topology_operator_report,
};
use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportContinuationPlan, PlannedTopologyOperation, SelectedSideSupportRegionPlan,
    SideSupportRegionPlanSelection,
};
use super::repeated_apex_side_support_result::SideSupportRegionContinuationResult;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportPrimitiveDerivation {
    regions: usize,
    ready_regions: usize,
    primitive: Option<OutsideOutsideSupportPrimitivePlan>,
    details: Vec<OutsideOutsideSupportRegionDerivation>,
}

#[derive(Clone, Debug)]
pub(super) struct OutsideOutsideSupportPrimitivePlan {
    requests: Vec<CylinderMediatedSupportContinuationRequest>,
}

#[derive(Clone, Debug)]
struct OutsideOutsideSupportRegionDerivation {
    region: OutsideOutsideSupportContinuationPlan,
    topology_class: Option<CylinderMediatedOutsideOutsideClass>,
    core_preflights: Vec<CylinderMediatedCorePreflight>,
    reason: OutsideOutsideSupportRefusal,
}

#[derive(Clone, Debug)]
enum OutsideOutsideSupportRefusal {
    NoOutsideCylinderAdjacency,
    NoSupportCylinderAdjacency,
    NoDirectedCylinderInterval,
    CrossOutsideOwnersRequirePrimitive,
    CorePreflightRejected,
    TopologyOperatorNotImplemented,
}

#[derive(Clone, Debug)]
pub(super) struct CylinderMediatedOutsideOutsideClass {
    pub(super) cylinder_face: FaceHandle,
    pub(super) outside_adjacencies: Vec<LoopAdjacencyOperand>,
    pub(super) support_adjacencies: Vec<LoopAdjacencyOperand>,
    pub(super) cylinder_intervals: Vec<CylinderIntervalOperand>,
}

#[derive(Clone, Debug)]
pub(super) struct LoopAdjacencyOperand {
    pub(super) source_loop: LoopHandle,
    pub(super) source_half_edge: HalfEdgeHandle,
    pub(super) cylinder_loop: LoopHandle,
    pub(super) cylinder_half_edge: HalfEdgeHandle,
}

#[derive(Clone, Debug)]
pub(super) struct CylinderIntervalOperand {
    pub(super) outside_loop: LoopHandle,
    pub(super) outside_cylinder_source_half_edge: HalfEdgeHandle,
    pub(super) outside_cylinder_half_edge: HalfEdgeHandle,
    pub(super) support_loop: LoopHandle,
    pub(super) support_cylinder_source_half_edge: HalfEdgeHandle,
    pub(super) support_cylinder_half_edge: HalfEdgeHandle,
    pub(super) cylinder_loop: LoopHandle,
    pub(super) interval: Vec<HalfEdgeHandle>,
}

#[derive(Clone, Debug)]
pub(super) enum CylinderMediatedCorePreflight {
    Accepted {
        fillet_index: usize,
        outside_loop: LoopHandle,
        support_loop: LoopHandle,
        outside_interval_len: usize,
        outside_complement_interval_len: usize,
        request: CylinderMediatedSupportContinuationRequest,
        reason: &'static str,
    },
    Rejected {
        fillet_index: usize,
        reason: String,
    },
    CrossOutsideOwners {
        fillet_index: usize,
        start_face: FaceHandle,
        start_loop: LoopHandle,
        end_face: FaceHandle,
        end_loop: LoopHandle,
        start_half_edge: HalfEdgeHandle,
        end_half_edge: HalfEdgeHandle,
    },
}

impl OutsideOutsideSupportPrimitiveDerivation {
    pub(super) fn build(
        selection: &SideSupportRegionPlanSelection,
        knowledge: &TopologyKnowledgeGraph,
        brep: &BRep,
    ) -> Self {
        match selection {
            SideSupportRegionPlanSelection::Selected(plan)
            | SideSupportRegionPlanSelection::Blocked(plan, _) => {
                Self::from_selected_plan(plan, knowledge, brep)
            }
            SideSupportRegionPlanSelection::Refused(_) => Self {
                regions: 0,
                ready_regions: 0,
                primitive: None,
                details: Vec::new(),
            },
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "outside-outside-support primitive_ready={} regions={} ready_regions={} totals={} details=[{}]",
            self.is_ready(),
            self.regions,
            self.ready_regions,
            self.primitive
                .as_ref()
                .map(OutsideOutsideSupportPrimitivePlan::diagnostic_label)
                .unwrap_or_else(|| "none".to_string()),
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

    pub(super) fn is_complete(&self) -> bool {
        self.is_ready() && self.primitive.is_some()
    }

    pub(super) fn topology_operator_report(
        &self,
        brep: &BRep,
    ) -> CylinderMediatedSupportContinuationTopologyOperatorReport {
        topology_operator_report(self.primitive.as_ref(), brep)
    }

    pub(super) fn apply_topology_operator(
        &self,
        brep: &mut BRep,
    ) -> Result<SideSupportRegionContinuationResult, CoreError> {
        let primitive = self.primitive.as_ref().ok_or_else(|| {
            CoreError::PreconditionFailed(
                "outside/outside support continuation has no primitive plan".into(),
            )
        })?;
        let batch = brep
            .preflight_cylinder_mediated_support_continuation_batch(primitive.requests.clone())?;
        let result = brep.cylinder_mediated_support_continuation_batch(&batch)?;
        match result.postconditions {
            CylinderMediatedSupportContinuationPostconditions::Satisfied => Ok(
                SideSupportRegionContinuationResult::from_outside_outside_cylinder_payload(
                    result.bridge_faces,
                    result.bridge_loops,
                ),
            ),
            CylinderMediatedSupportContinuationPostconditions::Refused { reason } => {
                Err(CoreError::PreconditionFailed(format!(
                    "outside/outside support continuation refused after mutation boundary: {reason}"
                )))
            }
        }
    }

    fn from_selected_plan(
        plan: &SelectedSideSupportRegionPlan,
        knowledge: &TopologyKnowledgeGraph,
        brep: &BRep,
    ) -> Self {
        let details = plan
            .operations()
            .filter_map(|operation| match operation {
                PlannedTopologyOperation::OutsideOutsideSupportContinuation(region) => Some(
                    OutsideOutsideSupportRegionDerivation::derive(region, knowledge, brep),
                ),
                PlannedTopologyOperation::SplitFaceWithEdge { .. }
                | PlannedTopologyOperation::BridgeSupportRegion(_) => None,
            })
            .collect::<Vec<_>>();
        let regions = details.len();
        let ready_regions = details
            .iter()
            .filter(|detail| detail.is_core_preflight_ready())
            .count();
        let primitive = if regions > 0 && ready_regions == regions {
            let requests = details
                .iter()
                .flat_map(OutsideOutsideSupportRegionDerivation::accepted_core_requests)
                .cloned()
                .collect::<Vec<_>>();
            Some(OutsideOutsideSupportPrimitivePlan { requests })
        } else {
            None
        };
        Self {
            regions,
            ready_regions,
            primitive,
            details,
        }
    }
}

impl OutsideOutsideSupportPrimitivePlan {
    pub(super) fn requests(&self) -> &[CylinderMediatedSupportContinuationRequest] {
        &self.requests
    }

    fn diagnostic_label(&self) -> String {
        format!("requests={}", self.requests.len())
    }
}

impl OutsideOutsideSupportRegionDerivation {
    fn derive(
        region: &OutsideOutsideSupportContinuationPlan,
        knowledge: &TopologyKnowledgeGraph,
        brep: &BRep,
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
                cylinder_intervals: cylinder_interval_candidates(
                    brep,
                    &outside_adjacencies,
                    &support_adjacencies,
                ),
            });
        let core_preflights = topology_class
            .as_ref()
            .map(|class| validate_core_preflights(region, class, knowledge, brep))
            .unwrap_or_default();
        let all_core_preflights_accepted = !region.trims().is_empty()
            && core_preflights.len() == region.trims().len()
            && core_preflights
                .iter()
                .all(CylinderMediatedCorePreflight::is_accepted);
        let reason = if outside_adjacencies.is_empty() {
            OutsideOutsideSupportRefusal::NoOutsideCylinderAdjacency
        } else if support_adjacencies.is_empty() {
            OutsideOutsideSupportRefusal::NoSupportCylinderAdjacency
        } else if topology_class
            .as_ref()
            .is_none_or(|class| class.cylinder_intervals.is_empty())
        {
            OutsideOutsideSupportRefusal::NoDirectedCylinderInterval
        } else if core_preflights
            .iter()
            .any(CylinderMediatedCorePreflight::is_cross_outside_owner)
        {
            OutsideOutsideSupportRefusal::CrossOutsideOwnersRequirePrimitive
        } else if !all_core_preflights_accepted {
            OutsideOutsideSupportRefusal::CorePreflightRejected
        } else {
            OutsideOutsideSupportRefusal::TopologyOperatorNotImplemented
        };
        Self {
            region: region.clone(),
            topology_class,
            core_preflights,
            reason,
        }
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} trims={} {} core_preflight=[{}] {}",
            self.region.support_face(),
            self.region.trims().len(),
            self.topology_class
                .as_ref()
                .map(CylinderMediatedOutsideOutsideClass::diagnostic_label)
                .unwrap_or_else(|| "class=unknown".to_string()),
            self.core_preflights
                .iter()
                .map(CylinderMediatedCorePreflight::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", "),
            self.reason.diagnostic_label()
        )
    }

    fn is_core_preflight_ready(&self) -> bool {
        !self.region.trims().is_empty()
            && self.core_preflights.len() == self.region.trims().len()
            && self
                .core_preflights
                .iter()
                .all(CylinderMediatedCorePreflight::is_accepted)
    }

    fn accepted_core_requests(
        &self,
    ) -> impl Iterator<Item = &CylinderMediatedSupportContinuationRequest> {
        self.core_preflights.iter().filter_map(|preflight| {
            if let CylinderMediatedCorePreflight::Accepted { request, .. } = preflight {
                Some(request)
            } else {
                None
            }
        })
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
            OutsideOutsideSupportRefusal::NoDirectedCylinderInterval => {
                "missing-directed-cylinder-interval"
            }
            OutsideOutsideSupportRefusal::CrossOutsideOwnersRequirePrimitive => {
                "cross-outside-owner-primitive-required"
            }
            OutsideOutsideSupportRefusal::CorePreflightRejected => "core-preflight-rejected",
            OutsideOutsideSupportRefusal::TopologyOperatorNotImplemented => {
                "missing-topology-operator"
            }
        }
    }
}

impl CylinderMediatedOutsideOutsideClass {
    fn diagnostic_label(&self) -> String {
        format!(
            "class=cylinder-mediated cylinder={:?} outside_adjacencies=[{}] support_adjacencies=[{}] cylinder_intervals=[{}]",
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
                .join(", "),
            self.cylinder_intervals
                .iter()
                .map(CylinderIntervalOperand::diagnostic_label)
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

impl CylinderIntervalOperand {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}->{:?} cylinder={:?}->{:?} loop={:?} len={} path={:?}",
            self.outside_cylinder_source_half_edge,
            self.support_cylinder_source_half_edge,
            self.outside_cylinder_half_edge,
            self.support_cylinder_half_edge,
            self.cylinder_loop,
            self.interval.len(),
            self.interval
        )
    }
}

impl CylinderMediatedCorePreflight {
    fn is_accepted(&self) -> bool {
        matches!(self, CylinderMediatedCorePreflight::Accepted { .. })
    }

    fn is_cross_outside_owner(&self) -> bool {
        matches!(
            self,
            CylinderMediatedCorePreflight::CrossOutsideOwners { .. }
        )
    }

    fn diagnostic_label(&self) -> String {
        match self {
            CylinderMediatedCorePreflight::Accepted {
                fillet_index,
                outside_loop,
                support_loop,
                outside_interval_len,
                outside_complement_interval_len,
                request: _,
                reason,
            } => {
                format!(
                    "FilletEnd[{fillet_index}] accepted outside_loop={outside_loop:?} support_loop={support_loop:?} outside_interval={outside_interval_len} outside_complement={outside_complement_interval_len} delta=unknown reason={reason}"
                )
            }
            CylinderMediatedCorePreflight::Rejected {
                fillet_index,
                reason,
            } => {
                format!("FilletEnd[{fillet_index}] rejected reason={reason}")
            }
            CylinderMediatedCorePreflight::CrossOutsideOwners {
                fillet_index,
                start_face,
                start_loop,
                end_face,
                end_loop,
                start_half_edge,
                end_half_edge,
            } => {
                format!(
                    "FilletEnd[{fillet_index}] cross-outside-owner start={start_face:?}/{start_loop:?}@{start_half_edge:?} end={end_face:?}/{end_loop:?}@{end_half_edge:?} primitive=multi-outside-owner-support-continuation-required"
                )
            }
        }
    }
}
