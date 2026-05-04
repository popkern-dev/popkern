//! Execution preflight for face-family bridge topology operators.

use popkern_core::euler::TwoLoopBridgeRequest;
use popkern_core::topology::BRep;

use super::repeated_apex_side_support_bridge::{
    FaceFamilyBridgePrimitivePlan, FaceFamilyBridgeRegionPrimitivePlan,
};

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
    Refused(FaceFamilyBridgeTopologyOperatorGap),
    RefusedWithDetail(FaceFamilyBridgeTopologyOperatorGap, String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceFamilyBridgeTopologyOperatorGap {
    MissingPrimitivePlan,
    EmptyPrimitivePlan,
    MissingBoundaryCycle,
    CorePreflightRejected,
    MissingCrossFaceBridgeEulerOperator,
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
    let regions = primitive.regions.len();
    let boundary_cycles = primitive
        .regions
        .iter()
        .map(|region| region.boundary_cycles.len())
        .sum::<usize>();
    let expected_edges = primitive
        .regions
        .iter()
        .map(|region| region.expected_edges)
        .sum::<usize>();
    let expected_faces = primitive
        .regions
        .iter()
        .map(|region| region.expected_faces)
        .sum::<usize>();
    let expected_loops = primitive
        .regions
        .iter()
        .map(|region| region.expected_loops)
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
    if primitive.regions.iter().any(|region| {
        region
            .boundary_cycles
            .iter()
            .any(|cycle| cycle.support_interval.is_empty() || cycle.outside_interval.is_empty())
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
    if let Some(core_preflight_failure) = primitive
        .regions
        .iter()
        .find_map(|region| core_bridge_preflight_failure(region, brep))
    {
        return FaceFamilyBridgeTopologyOperatorReport::refused_with_detail(
            regions,
            boundary_cycles,
            expected_edges,
            expected_faces,
            expected_loops,
            FaceFamilyBridgeTopologyOperatorGap::CorePreflightRejected,
            core_preflight_failure,
        );
    }
    FaceFamilyBridgeTopologyOperatorReport::refused(
        regions,
        boundary_cycles,
        expected_edges,
        expected_faces,
        expected_loops,
        FaceFamilyBridgeTopologyOperatorGap::MissingCrossFaceBridgeEulerOperator,
    )
}

impl FaceFamilyBridgeTopologyOperatorReport {
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
}

impl FaceFamilyBridgeTopologyOperatorStatus {
    fn is_ready(&self) -> bool {
        match self {
            FaceFamilyBridgeTopologyOperatorStatus::Refused(_) => false,
            FaceFamilyBridgeTopologyOperatorStatus::RefusedWithDetail(_, _) => false,
        }
    }

    fn diagnostic_label(&self) -> String {
        match self {
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
            FaceFamilyBridgeTopologyOperatorGap::CorePreflightRejected => "core-preflight-rejected",
            FaceFamilyBridgeTopologyOperatorGap::MissingCrossFaceBridgeEulerOperator => {
                "missing-cross-face-bridge-euler-operator"
            }
        }
    }
}

fn core_bridge_preflight_failure(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
    brep: &BRep,
) -> Option<String> {
    brep.preflight_two_loop_bridge(TwoLoopBridgeRequest {
        first_loop: region.support_loop,
        first_start: region.trims[0].support_half_edge,
        first_end: region.trims[1].support_half_edge,
        second_loop: region.outside_loop,
        second_start: region.trims[0].outside_half_edge,
        second_end: region.trims[1].outside_half_edge,
    })
    .err()
    .map(|err| format!("{err}"))
}
