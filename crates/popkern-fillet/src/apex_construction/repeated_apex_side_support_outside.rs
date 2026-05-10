//! Outside/outside support-continuation primitive derivation.

use popkern_core::euler::{
    CylinderMediatedSupportContinuationPlan, CylinderMediatedSupportContinuationRequest,
};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::repeated_apex_side_support_outside_operator::{
    CylinderMediatedSupportContinuationTopologyOperatorReport, topology_operator_report,
};
use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportContinuationPlan, PlannedTopologyOperation, SelectedSideSupportRegionPlan,
    SideSupportRegionPlanSelection,
};
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
    CorePreflightRejected,
    TopologyOperatorNotImplemented,
}

#[derive(Clone, Debug)]
struct CylinderMediatedOutsideOutsideClass {
    cylinder_face: FaceHandle,
    outside_adjacencies: Vec<LoopAdjacencyOperand>,
    support_adjacencies: Vec<LoopAdjacencyOperand>,
    cylinder_intervals: Vec<CylinderIntervalOperand>,
}

#[derive(Clone, Debug)]
struct LoopAdjacencyOperand {
    source_loop: LoopHandle,
    source_half_edge: HalfEdgeHandle,
    cylinder_loop: LoopHandle,
    cylinder_half_edge: HalfEdgeHandle,
}

#[derive(Clone, Debug)]
struct CylinderIntervalOperand {
    outside_loop: LoopHandle,
    outside_cylinder_source_half_edge: HalfEdgeHandle,
    outside_cylinder_half_edge: HalfEdgeHandle,
    support_loop: LoopHandle,
    support_cylinder_source_half_edge: HalfEdgeHandle,
    support_cylinder_half_edge: HalfEdgeHandle,
    cylinder_loop: LoopHandle,
    interval: Vec<HalfEdgeHandle>,
}

#[derive(Clone, Debug)]
enum CylinderMediatedCorePreflight {
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
    ) -> CylinderMediatedSupportContinuationTopologyOperatorReport {
        topology_operator_report(self.primitive.as_ref())
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
            .map(|class| validate_core_preflights(region, class, brep))
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
    fn sort_key(&self) -> (usize, HalfEdgeHandle, HalfEdgeHandle) {
        (
            self.interval.len(),
            self.outside_cylinder_source_half_edge,
            self.support_cylinder_source_half_edge,
        )
    }

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
        }
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

fn cylinder_interval_candidates(
    brep: &BRep,
    outside_adjacencies: &[LoopAdjacencyOperand],
    support_adjacencies: &[LoopAdjacencyOperand],
) -> Vec<CylinderIntervalOperand> {
    let mut intervals = outside_adjacencies
        .iter()
        .flat_map(|outside| {
            support_adjacencies
                .iter()
                .filter(move |support| outside.cylinder_loop == support.cylinder_loop)
                .filter_map(move |support| {
                    directed_loop_path(
                        brep,
                        outside.cylinder_loop,
                        outside.cylinder_half_edge,
                        support.cylinder_half_edge,
                    )
                    .map(|interval| CylinderIntervalOperand {
                        outside_loop: outside.source_loop,
                        outside_cylinder_source_half_edge: outside.source_half_edge,
                        outside_cylinder_half_edge: outside.cylinder_half_edge,
                        support_loop: support.source_loop,
                        support_cylinder_source_half_edge: support.source_half_edge,
                        support_cylinder_half_edge: support.cylinder_half_edge,
                        cylinder_loop: outside.cylinder_loop,
                        interval,
                    })
                })
        })
        .collect::<Vec<_>>();
    intervals.sort_by_key(CylinderIntervalOperand::sort_key);
    intervals
}

fn directed_loop_path(
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

fn validate_core_preflights(
    region: &OutsideOutsideSupportContinuationPlan,
    topology_class: &CylinderMediatedOutsideOutsideClass,
    brep: &BRep,
) -> Vec<CylinderMediatedCorePreflight> {
    region
        .trims()
        .iter()
        .map(|trim| validate_trim_core_preflight(trim, topology_class, brep))
        .collect()
}

fn validate_trim_core_preflight(
    trim: &super::repeated_apex_side_support_plan::OutsideOutsideSupportTrimPlan,
    topology_class: &CylinderMediatedOutsideOutsideClass,
    brep: &BRep,
) -> CylinderMediatedCorePreflight {
    if trim.start_face() != trim.end_face() || trim.start_loop() != trim.end_loop() {
        return CylinderMediatedCorePreflight::Rejected {
            fillet_index: trim.fillet_index(),
            reason: format!(
                "outside endpoints span multiple outside owners: start={:?}/{:?} end={:?}/{:?}",
                trim.start_face(),
                trim.start_loop(),
                trim.end_face(),
                trim.end_loop()
            ),
        };
    }
    let Some(interval) = topology_class
        .cylinder_intervals
        .iter()
        .find(|interval| interval.outside_loop == trim.start_loop())
    else {
        return CylinderMediatedCorePreflight::Rejected {
            fillet_index: trim.fillet_index(),
            reason: format!(
                "missing cylinder interval for outside loop {:?}",
                trim.start_loop()
            ),
        };
    };
    let request = CylinderMediatedSupportContinuationRequest {
        support_face: trim.support_face(),
        outside_face: trim.start_face(),
        cylinder_face: trim.cylinder_face(),
        start: trim.start(),
        end: trim.end(),
        outside_loop: trim.start_loop(),
        outside_start_half_edge: trim.start_half_edge(),
        outside_end_half_edge: trim.end_half_edge(),
        support_loop: interval.support_loop,
        support_cylinder_half_edge: interval.support_cylinder_source_half_edge,
        outside_cylinder_half_edge: interval.outside_cylinder_source_half_edge,
        cylinder_interval: interval.interval.clone(),
    };
    match brep.preflight_cylinder_mediated_support_continuation(request) {
        Ok(plan) => CylinderMediatedCorePreflight::Accepted {
            fillet_index: trim.fillet_index(),
            outside_loop: trim.start_loop(),
            support_loop: interval.support_loop,
            outside_interval_len: plan.outside_interval.len(),
            outside_complement_interval_len: plan.outside_complement_interval.len(),
            reason: preflight_delta_reason(&plan),
            request: plan.request,
        },
        Err(err) => CylinderMediatedCorePreflight::Rejected {
            fillet_index: trim.fillet_index(),
            reason: err.to_string(),
        },
    }
}

fn preflight_delta_reason(plan: &CylinderMediatedSupportContinuationPlan) -> &'static str {
    match plan.expected_delta {
        popkern_core::euler::CylinderMediatedSupportContinuationDelta::Unknown { reason } => reason,
    }
}
