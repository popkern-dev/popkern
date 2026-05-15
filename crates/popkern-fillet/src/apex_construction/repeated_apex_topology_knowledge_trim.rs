//! Trim endpoint agreement facts for repeated-apex topology knowledge.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::repeated_apex_topology_knowledge::LoopAdjacencyFacts;

#[derive(Clone, Debug)]
pub(super) struct TrimEndpointAgreementFacts {
    half_edge: HalfEdgeHandle,
    adjacent_half_edge: HalfEdgeHandle,
    status: TrimEndpointAgreementStatus,
    direct_disagreement: f64,
    reconciled_disagreement: f64,
    selected_disagreement: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TrimEndpointAgreementStatus {
    AcceptedDirect,
    AcceptedReconciled,
    SideDisagreement,
    MissingSurface,
    MissingPcurve,
    MissingTopology,
    EvaluationFailed,
}

pub(super) fn collect_trim_endpoint_agreements(
    brep: &BRep,
    adjacencies: &[LoopAdjacencyFacts],
) -> Vec<TrimEndpointAgreementFacts> {
    adjacencies
        .iter()
        .map(|adjacency| trim_endpoint_agreement(brep, adjacency))
        .collect()
}

impl TrimEndpointAgreementFacts {
    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{:?}<->{:?} status={:?} direct={:.3e} reconciled={:.3e} selected={:.3e}",
            self.half_edge,
            self.adjacent_half_edge,
            self.status,
            self.direct_disagreement,
            self.reconciled_disagreement,
            self.selected_disagreement
        )
    }
}

fn trim_endpoint_agreement(
    brep: &BRep,
    adjacency: &LoopAdjacencyFacts,
) -> TrimEndpointAgreementFacts {
    let result = measure_trim_endpoint_agreement(brep, adjacency);
    TrimEndpointAgreementFacts {
        half_edge: adjacency.half_edge(),
        adjacent_half_edge: adjacency.adjacent_half_edge(),
        status: result.status,
        direct_disagreement: result.direct_disagreement,
        reconciled_disagreement: result.reconciled_disagreement,
        selected_disagreement: result
            .direct_disagreement
            .min(result.reconciled_disagreement),
    }
}

struct TrimEndpointAgreementMeasurement {
    status: TrimEndpointAgreementStatus,
    direct_disagreement: f64,
    reconciled_disagreement: f64,
}

fn measure_trim_endpoint_agreement(
    brep: &BRep,
    adjacency: &LoopAdjacencyFacts,
) -> TrimEndpointAgreementMeasurement {
    let points_a = match trim_endpoint_points(brep, adjacency.half_edge()) {
        Ok(points) => points,
        Err(status) => return unavailable(status),
    };
    let points_b = match trim_endpoint_points(brep, adjacency.adjacent_half_edge()) {
        Ok(points) => points,
        Err(status) => return unavailable(status),
    };
    let direct = max_pair_gap(points_a, points_b);
    let reconciled = max_pair_gap(points_a, [points_b[1], points_b[0]]);
    let tolerance = brep.tolerance.linear();
    let status = if !direct.is_finite() || !reconciled.is_finite() {
        TrimEndpointAgreementStatus::EvaluationFailed
    } else if direct <= tolerance {
        TrimEndpointAgreementStatus::AcceptedDirect
    } else if reconciled <= tolerance {
        TrimEndpointAgreementStatus::AcceptedReconciled
    } else {
        TrimEndpointAgreementStatus::SideDisagreement
    };
    TrimEndpointAgreementMeasurement {
        status,
        direct_disagreement: direct,
        reconciled_disagreement: reconciled,
    }
}

fn trim_endpoint_points(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Result<[Point3; 2], TrimEndpointAgreementStatus> {
    let half_edge_data = brep
        .half_edges
        .get(half_edge)
        .ok_or(TrimEndpointAgreementStatus::MissingTopology)?;
    let loop_data = brep
        .loops
        .get(half_edge_data.loop_)
        .ok_or(TrimEndpointAgreementStatus::MissingTopology)?;
    let surface = &brep
        .faces
        .get(loop_data.face)
        .ok_or(TrimEndpointAgreementStatus::MissingTopology)?
        .surface
        .as_ref()
        .ok_or(TrimEndpointAgreementStatus::MissingSurface)?
        .geometry;
    let pcurve = half_edge_data
        .pcurve
        .as_ref()
        .ok_or(TrimEndpointAgreementStatus::MissingPcurve)?;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    (point_is_finite(start) && point_is_finite(end))
        .then_some([start, end])
        .ok_or(TrimEndpointAgreementStatus::EvaluationFailed)
}

fn max_pair_gap(a: [Point3; 2], b: [Point3; 2]) -> f64 {
    a[0].distance_to(&b[0]).max(a[1].distance_to(&b[1]))
}

fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn unavailable(status: TrimEndpointAgreementStatus) -> TrimEndpointAgreementMeasurement {
    TrimEndpointAgreementMeasurement {
        status,
        direct_disagreement: f64::INFINITY,
        reconciled_disagreement: f64::INFINITY,
    }
}
