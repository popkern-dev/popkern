//! Geometry feasibility facts for bridge diagnostic candidates.

use crate::euler::{
    BridgeDiagnosticOperationRequest, ConstructionFamily, DiagnosticTolerance,
    DiagnosticToleranceSource, GeometryFeasibility, GeometryFeasibilityFact,
};
use crate::topology::{BRep, HalfEdgeHandle};

pub(super) fn geometry_facts_for_candidate(
    brep: &BRep,
    operation: &BridgeDiagnosticOperationRequest,
    family: ConstructionFamily,
) -> Vec<GeometryFeasibilityFact> {
    if !candidate_needs_bridge_span_geometry(family) {
        return vec![not_applicable_fact()];
    }
    let tolerance = DiagnosticTolerance {
        source: DiagnosticToleranceSource::LocalTopology,
        nanometers: Some(nanometers(brep.tolerance.linear())),
    };
    match shortest_interval_edge_length(
        brep,
        [&operation.first_interval, &operation.second_interval],
    ) {
        Some(length) if length <= brep.tolerance.linear() => vec![GeometryFeasibilityFact {
            feasibility: GeometryFeasibility::RepairRequiredBeforeConstruction,
            reason: "bridge-interval-has-short-or-degenerate-edge",
            tolerance,
        }],
        Some(_) => vec![GeometryFeasibilityFact {
            feasibility: GeometryFeasibility::GeometryReady,
            reason: "bridge-interval-edge-lengths-clear-local-tolerance",
            tolerance,
        }],
        None => vec![GeometryFeasibilityFact {
            feasibility: GeometryFeasibility::Unknown,
            reason: "bridge-interval-geometry-could-not-be-measured",
            tolerance,
        }],
    }
}

fn candidate_needs_bridge_span_geometry(family: ConstructionFamily) -> bool {
    matches!(
        family,
        ConstructionFamily::DisjointCrossLoopBridge
            | ConstructionFamily::DisjointSameLoopBridge
            | ConstructionFamily::BoundarySideSeparatorBridge
    )
}

fn not_applicable_fact() -> GeometryFeasibilityFact {
    GeometryFeasibilityFact {
        feasibility: GeometryFeasibility::NotApplicable,
        reason: "topology-only-candidate-no-geometry-fact-required",
        tolerance: DiagnosticTolerance::NOT_APPLICABLE,
    }
}

fn shortest_interval_edge_length<'a>(
    brep: &BRep,
    intervals: impl IntoIterator<Item = &'a Vec<HalfEdgeHandle>>,
) -> Option<f64> {
    intervals
        .into_iter()
        .filter(|interval| interval.len() >= 2)
        .flat_map(|interval| &interval[..interval.len() - 1])
        .filter_map(|half_edge| edge_length(brep, *half_edge))
        .min_by(|lhs, rhs| lhs.total_cmp(rhs))
}

fn edge_length(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<f64> {
    let data = brep.half_edges.get(half_edge)?;
    let origin = brep.vertices.get(data.origin)?.position;
    let next = brep.half_edges.get(data.next)?;
    let destination = brep.vertices.get(next.origin)?.position;
    Some(origin.distance_to(&destination))
}

fn nanometers(value: f64) -> i64 {
    (value * 1.0e9).round() as i64
}
