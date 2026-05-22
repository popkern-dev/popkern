//! Execution preflight for face-family bridge topology operators.

use popkern_core::euler::{
    FaceFamilyBridgeCycleRequest, FaceFamilyBridgeResult, TwoLoopBridgeRequest,
};
use popkern_core::topology::BRep;

use super::repeated_apex_side_support_boundary::BridgeBoundaryExposureStatus;
use super::repeated_apex_side_support_bridge::{
    FaceFamilyBridgePrimitivePlan, FaceFamilyBridgeRegionPrimitivePlan,
};
use super::repeated_apex_side_support_result::SideSupportRegionContinuationResult;
#[path = "repeated_apex_side_support_bridge_operator_exposure.rs"]
mod bridge_operator_exposure;
#[path = "repeated_apex_side_support_bridge_operator_exposure_walk.rs"]
mod bridge_operator_exposure_walk;

use bridge_operator_exposure::{expose_aggregate_bridge_boundaries, expose_bridge_boundaries};

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeTopologyOperatorReport {
    regions: usize,
    boundary_cycles: usize,
    expected_edges: usize,
    expected_faces: usize,
    expected_loops: usize,
    status: FaceFamilyBridgeTopologyOperatorStatus,
}

#[derive(Clone, Debug)]
enum FaceFamilyBridgeTopologyOperatorStatus {
    Ready,
    Refused(FaceFamilyBridgeTopologyOperatorGap),
    RefusedWithDetail(FaceFamilyBridgeTopologyOperatorGap, String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceFamilyBridgeTopologyOperatorGap {
    MissingPrimitivePlan,
    EmptyPrimitivePlan,
    MissingBoundaryCycle,
    MissingSelfOwnedBridgeOperandExposureOperator,
    MissingPartialFaceBoundaryExposureOperator,
    BridgeBoundaryExposureRejected,
    CorePreflightRejected,
    CoreBridgeOperatorRejected,
    CoreAggregateBridgeOperatorRejected,
}

pub(super) fn topology_operator_report(
    primitive: Option<&FaceFamilyBridgePrimitivePlan>,
    brep: &BRep,
) -> FaceFamilyBridgeTopologyOperatorReport {
    let Some(primitive) = primitive else {
        return FaceFamilyBridgeTopologyOperatorReport::refused(
            0,
            0,
            0,
            0,
            0,
            FaceFamilyBridgeTopologyOperatorGap::MissingPrimitivePlan,
        );
    };
    let regions = primitive.regions().len();
    let boundary_cycles = primitive
        .regions()
        .iter()
        .map(|region| region.boundary_cycles().len())
        .sum::<usize>();
    let expected_edges = primitive
        .regions()
        .iter()
        .map(|region| region.expected_edges())
        .sum::<usize>();
    let expected_faces = primitive
        .regions()
        .iter()
        .map(|region| region.expected_faces())
        .sum::<usize>();
    let expected_loops = primitive
        .regions()
        .iter()
        .map(|region| region.expected_loops())
        .sum::<usize>();
    if regions == 0 {
        return FaceFamilyBridgeTopologyOperatorReport::refused(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            FaceFamilyBridgeTopologyOperatorGap::EmptyPrimitivePlan,
        );
    }
    if primitive.regions().iter().any(|region| {
        region
            .boundary_cycles()
            .iter()
            .any(|cycle| cycle.is_empty())
    }) {
        return FaceFamilyBridgeTopologyOperatorReport::refused(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            FaceFamilyBridgeTopologyOperatorGap::MissingBoundaryCycle,
        );
    }
    if let Some(not_ready) = primitive
        .regions()
        .iter()
        .find_map(bridge_boundary_readiness_failure)
    {
        return FaceFamilyBridgeTopologyOperatorReport::refused_with_detail(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            not_ready.0,
            not_ready.1,
        );
    }
    if let [region] = primitive.regions()
        && let Some(core_operator_failure) = core_bridge_operator_failure(region, brep)
    {
        return FaceFamilyBridgeTopologyOperatorReport::refused_with_detail(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            core_operator_failure.0,
            core_operator_failure.1,
        );
    }
    if let Err(detail) = apply_topology_operator(primitive, &mut brep.clone()) {
        return FaceFamilyBridgeTopologyOperatorReport::refused_with_detail(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            FaceFamilyBridgeTopologyOperatorGap::CoreAggregateBridgeOperatorRejected,
            detail,
        );
    }
    FaceFamilyBridgeTopologyOperatorReport::ready(
        regions,
        boundary_cycles,
        expected_edges,
        expected_faces,
        expected_loops,
    )
}

pub(super) fn apply_topology_operator(
    primitive: &FaceFamilyBridgePrimitivePlan,
    brep: &mut BRep,
) -> Result<SideSupportRegionContinuationResult, String> {
    let mut candidate = brep.clone();
    expose_aggregate_bridge_boundaries(primitive, &mut candidate)?;
    if let [region] = primitive.regions() {
        let bridge_result =
            apply_region_bridge(region, &mut candidate).map_err(|(_, detail)| detail)?;
        let result = SideSupportRegionContinuationResult::from_bridge_payload(
            vec![
                bridge_result.first_bridge_face,
                bridge_result.second_bridge_face,
            ],
            vec![
                bridge_result.first_bridge_loop,
                bridge_result.second_bridge_loop,
            ],
        );
        *brep = candidate;
        return Ok(result);
    }
    let aggregate_cycles = aggregate_bridge_cycles(primitive);
    let aggregate_plan = candidate
        .preflight_face_family_bridge(aggregate_cycles)
        .map_err(|err| format!("core aggregate face-family bridge preflight rejected: {err}"))?;
    let bridge_result = candidate
        .face_family_bridge(&aggregate_plan)
        .map_err(|err| format!("core aggregate face-family bridge rejected: {err}"))?;
    assign_aggregate_bridge_support_surfaces(
        &aggregate_plan.cycles,
        &bridge_result,
        &mut candidate,
    )?;
    let result = SideSupportRegionContinuationResult::from_bridge_payload(
        bridge_result
            .cycles
            .iter()
            .map(|cycle| cycle.bridge_face)
            .collect(),
        bridge_result
            .cycles
            .iter()
            .map(|cycle| cycle.bridge_loop)
            .collect(),
    );
    *brep = candidate;
    Ok(result)
}

impl FaceFamilyBridgeTopologyOperatorReport {
    fn ready(
        regions: usize,
        boundary_cycles: usize,
        expected_edges: usize,
        expected_faces: usize,
        expected_loops: usize,
    ) -> Self {
        Self {
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            status: FaceFamilyBridgeTopologyOperatorStatus::Ready,
        }
    }

    fn refused(
        regions: usize,
        boundary_cycles: usize,
        expected_edges: usize,
        expected_faces: usize,
        expected_loops: usize,
        reason: FaceFamilyBridgeTopologyOperatorGap,
    ) -> Self {
        Self {
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            status: FaceFamilyBridgeTopologyOperatorStatus::Refused(reason),
        }
    }

    fn refused_with_detail(
        regions: usize,
        boundary_cycles: usize,
        expected_edges: usize,
        expected_faces: usize,
        expected_loops: usize,
        reason: FaceFamilyBridgeTopologyOperatorGap,
        detail: String,
    ) -> Self {
        Self {
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            status: FaceFamilyBridgeTopologyOperatorStatus::RefusedWithDetail(reason, detail),
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "topology_operator_ready={} regions={} boundary_cycles={} expected_delta=dV=0 dE={} dF={} dL={} reason={}",
            self.status.is_ready(),
            self.regions,
            self.boundary_cycles,
            self.expected_edges,
            self.expected_faces,
            self.expected_loops,
            self.status.diagnostic_label()
        )
    }

    pub(super) fn is_ready(&self) -> bool {
        self.status.is_ready()
    }
}

impl FaceFamilyBridgeTopologyOperatorStatus {
    fn is_ready(&self) -> bool {
        match self {
            FaceFamilyBridgeTopologyOperatorStatus::Ready => true,
            FaceFamilyBridgeTopologyOperatorStatus::Refused(_) => false,
            FaceFamilyBridgeTopologyOperatorStatus::RefusedWithDetail(_, _) => false,
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
            FaceFamilyBridgeTopologyOperatorStatus::Ready => "ready".to_string(),
            FaceFamilyBridgeTopologyOperatorStatus::Refused(reason) => {
                reason.diagnostic_label().to_string()
            }
            FaceFamilyBridgeTopologyOperatorStatus::RefusedWithDetail(reason, detail) => {
                format!("{} detail={detail}", reason.diagnostic_label())
            }
        }
    }
}

impl FaceFamilyBridgeTopologyOperatorGap {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            FaceFamilyBridgeTopologyOperatorGap::MissingPrimitivePlan => "missing-primitive-plan",
            FaceFamilyBridgeTopologyOperatorGap::EmptyPrimitivePlan => "empty-primitive-plan",
            FaceFamilyBridgeTopologyOperatorGap::MissingBoundaryCycle => "missing-boundary-cycle",
            FaceFamilyBridgeTopologyOperatorGap::MissingSelfOwnedBridgeOperandExposureOperator => {
                "missing-self-owned-bridge-operand-exposure-operator"
            }
            FaceFamilyBridgeTopologyOperatorGap::MissingPartialFaceBoundaryExposureOperator => {
                "missing-partial-face-boundary-exposure-operator"
            }
            FaceFamilyBridgeTopologyOperatorGap::BridgeBoundaryExposureRejected => {
                "bridge-boundary-exposure-rejected"
            }
            FaceFamilyBridgeTopologyOperatorGap::CorePreflightRejected => "core-preflight-rejected",
            FaceFamilyBridgeTopologyOperatorGap::CoreBridgeOperatorRejected => {
                "core-bridge-operator-rejected"
            }
            FaceFamilyBridgeTopologyOperatorGap::CoreAggregateBridgeOperatorRejected => {
                "core-aggregate-face-family-bridge-operator-rejected"
            }
        }
    }
}

fn aggregate_bridge_cycles(
    primitive: &FaceFamilyBridgePrimitivePlan,
) -> Vec<FaceFamilyBridgeCycleRequest> {
    primitive
        .regions()
        .iter()
        .flat_map(|region| {
            region
                .boundary_cycles()
                .iter()
                .map(|cycle| FaceFamilyBridgeCycleRequest {
                    first_loop: region.support_loop(),
                    first_interval: cycle.support_interval().to_vec(),
                    second_loop: region.outside_loop(),
                    second_interval: cycle.outside_interval().to_vec(),
                })
        })
        .collect()
}

fn assign_aggregate_bridge_support_surfaces(
    cycles: &[FaceFamilyBridgeCycleRequest],
    result: &FaceFamilyBridgeResult,
    brep: &mut BRep,
) -> Result<(), String> {
    if cycles.len() != result.cycles.len() {
        return Err(format!(
            "aggregate bridge result count mismatch: {} request cycle(s), {} committed cycle(s)",
            cycles.len(),
            result.cycles.len()
        ));
    }
    for (request, committed) in cycles.iter().zip(&result.cycles) {
        let source_face = brep
            .loops
            .get(request.first_loop)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| {
                format!(
                    "aggregate bridge support surface assignment references stale source loop {:?}",
                    request.first_loop
                )
            })?;
        let source_surface = brep
            .faces
            .get(source_face)
            .and_then(|face| face.surface.clone())
            .ok_or_else(|| {
                format!(
                    "aggregate bridge support surface assignment requires surfaced support face {:?}",
                    source_face
                )
            })?;
        let bridge_face = brep.faces.get_mut(committed.bridge_face).ok_or_else(|| {
            format!(
                "aggregate bridge support surface assignment references stale bridge face {:?}",
                committed.bridge_face
            )
        })?;
        bridge_face.surface = Some(source_surface);
    }
    Ok(())
}

fn bridge_boundary_readiness_failure(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
) -> Option<(FaceFamilyBridgeTopologyOperatorGap, String)> {
    if region.boundary_readiness().is_boundary_ready() {
        return None;
    }
    let gap = match region.boundary_exposure().status() {
        BridgeBoundaryExposureStatus::AlreadyBoundaryReady
        | BridgeBoundaryExposureStatus::RetireWholeFaces => return None,
        BridgeBoundaryExposureStatus::PartialFaceOwnership => {
            if region
                .boundary_exposure()
                .is_partial_exposure_operator_ready()
            {
                return None;
            }
            if region
                .boundary_exposure()
                .has_self_owned_protected_operand_pair()
            {
                FaceFamilyBridgeTopologyOperatorGap::MissingSelfOwnedBridgeOperandExposureOperator
            } else {
                FaceFamilyBridgeTopologyOperatorGap::MissingPartialFaceBoundaryExposureOperator
            }
        }
    };
    Some((gap, region.boundary_exposure().diagnostic_label()))
}

fn core_bridge_operator_failure(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
    brep: &BRep,
) -> Option<(FaceFamilyBridgeTopologyOperatorGap, String)> {
    let mut candidate = brep.clone();
    apply_region_topology_operator(region, &mut candidate).err()
}

fn apply_region_topology_operator(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
    candidate: &mut BRep,
) -> Result<(), (FaceFamilyBridgeTopologyOperatorGap, String)> {
    if let Some(err) = expose_bridge_boundaries(region, candidate) {
        return Err((
            FaceFamilyBridgeTopologyOperatorGap::BridgeBoundaryExposureRejected,
            err,
        ));
    }
    apply_region_bridge(region, candidate).map(|_| ())
}

fn apply_region_bridge(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
    candidate: &mut BRep,
) -> Result<popkern_core::euler::TwoLoopBridgeResult, (FaceFamilyBridgeTopologyOperatorGap, String)>
{
    let plan = candidate
        .preflight_two_loop_bridge(TwoLoopBridgeRequest {
            first_loop: region.support_loop(),
            first_start: region.support_trims()[0].support_half_edge(),
            first_end: region.support_trims()[1].support_half_edge(),
            second_loop: region.outside_loop(),
            second_start: region.support_trims()[0].outside_half_edge(),
            second_end: region.support_trims()[1].outside_half_edge(),
        })
        .map_err(|err| {
            (
                FaceFamilyBridgeTopologyOperatorGap::CorePreflightRejected,
                format!("{err}"),
            )
        })?;
    candidate.two_loop_bridge(&plan).map_err(|err| {
        (
            FaceFamilyBridgeTopologyOperatorGap::CoreBridgeOperatorRejected,
            format!("{err}"),
        )
    })
}
