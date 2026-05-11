//! Bridge diagnostic construction-candidate classification.

use std::collections::BTreeSet;

use super::{
    ConstructionCandidateDraft, ConstructionFamily, ConstructionFeasibility,
    ConstructionRequirement, DiagnosticHandle, DiagnosticOperationId, DiagnosticOperatorFamily,
    DiagnosticResourceKind, DiagnosticTopologyDelta, OperatorReadinessStatus, ResourceConflictId,
    ResourceConflictKind, TopologyDiagnosticGraph, geometry,
};
use crate::euler::BridgeDiagnosticRequest;
use crate::topology::BRep;

pub(super) fn classify_construction_candidates(
    brep: &BRep,
    request: &BridgeDiagnosticRequest,
    mut graph: TopologyDiagnosticGraph,
) -> TopologyDiagnosticGraph {
    let readiness_entries = graph.readiness.clone();
    for readiness in readiness_entries {
        if let Some((family, feasibility, reason)) = candidate_for_readiness(
            &graph,
            readiness.operation,
            readiness.operator,
            readiness.status,
        ) {
            let postcondition_preview =
                postcondition_preview_for_candidate(&mut graph, readiness.operation, family);
            graph.push_candidate(ConstructionCandidateDraft {
                family,
                feasibility,
                feasibility_reason: reason,
                operations: vec![readiness.operation],
                conflicts: conflicts_for_operation(&graph, readiness.operation),
                requirements: vec![
                    ConstructionRequirement::LiveTopologyRevalidation,
                    ConstructionRequirement::PostconditionVerification,
                ],
                postcondition_preview,
                geometry_facts: request
                    .operations
                    .get(readiness.operation.0)
                    .map(|operation| {
                        geometry::geometry_facts_for_candidate(brep, operation, family)
                    })
                    .unwrap_or_default(),
            });
        }
    }
    graph
}

fn candidate_for_readiness(
    graph: &TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
    operator: DiagnosticOperatorFamily,
    status: OperatorReadinessStatus,
) -> Option<(ConstructionFamily, ConstructionFeasibility, &'static str)> {
    match (operator, status) {
        (DiagnosticOperatorFamily::FaceFamilyBridge, OperatorReadinessStatus::Ready) => Some((
            ConstructionFamily::DisjointCrossLoopBridge,
            ConstructionFeasibility::ExecutableByExistingOperator,
            "face-family-bridge-existing-operator-ready",
        )),
        (DiagnosticOperatorFamily::SameLoopFaceBridge, OperatorReadinessStatus::Ready) => Some((
            ConstructionFamily::DisjointSameLoopBridge,
            ConstructionFeasibility::ExecutableByExistingOperator,
            "same-loop-face-bridge-existing-operator-ready",
        )),
        (
            DiagnosticOperatorFamily::SameLoopBoundaryExposure,
            OperatorReadinessStatus::ReadyAfterNormalization,
        ) => Some((
            ConstructionFamily::PairedSameLoopTwinExposure,
            ConstructionFeasibility::ExecutableAfterKnownNormalization,
            "same-loop-boundary-exposure-normalization-ready",
        )),
        (
            DiagnosticOperatorFamily::BoundarySideAllocationOrSeparator,
            OperatorReadinessStatus::MissingOperator,
        ) => Some(boundary_side_separator_candidate(graph, operation)),
        _ => None,
    }
}

fn boundary_side_separator_candidate(
    graph: &TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
) -> (ConstructionFamily, ConstructionFeasibility, &'static str) {
    let already_boundary = conflict_count_for_operation(
        graph,
        operation,
        ResourceConflictKind::AlreadyBoundaryOverlap,
    );
    let duplicate_boundary = conflict_count_for_operation(
        graph,
        operation,
        ResourceConflictKind::DuplicateBoundaryEdgeConsumption,
    );
    let duplicate_half_edge = conflict_count_for_operation(
        graph,
        operation,
        ResourceConflictKind::DuplicateHalfEdgeConsumption,
    );
    let already_boundary_edges = unique_conflict_edge_count_for_operation(
        graph,
        operation,
        ResourceConflictKind::AlreadyBoundaryOverlap,
    );
    let (feasibility, reason) = if already_boundary == 1
        && already_boundary_edges == 1
        && duplicate_boundary == 1
        && duplicate_half_edge == 1
    {
        (
            ConstructionFeasibility::RequiresSeparatorSubcomplex,
            "missing-boundary-side-separator-subcomplex",
        )
    } else if already_boundary > 1 || already_boundary_edges > 1 {
        (
            ConstructionFeasibility::RequiresCycleRegrouping,
            "multiple-boundary-side-overlaps-require-cycle-regrouping",
        )
    } else if duplicate_half_edge > 1 || duplicate_boundary > 1 {
        (
            ConstructionFeasibility::RequiresCustomAggregateAuthoring,
            "boundary-side-overlap-requires-aggregate-authoring",
        )
    } else {
        (
            ConstructionFeasibility::InvalidNonManifold,
            "boundary-side-separator-missing-required-duplicate-boundary-side",
        )
    };
    (
        ConstructionFamily::BoundarySideSeparatorBridge,
        feasibility,
        reason,
    )
}

fn unique_conflict_edge_count_for_operation(
    graph: &TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
    kind: ResourceConflictKind,
) -> usize {
    graph
        .conflicts
        .iter()
        .filter(|conflict| conflict.kind == kind && conflict.operations.contains(&operation))
        .flat_map(|conflict| conflict.resources.iter())
        .filter_map(|resource_id| graph.resources.get(resource_id.0))
        .filter(|resource| resource.kind == DiagnosticResourceKind::EdgeIdentity)
        .filter_map(|resource| match resource.handle {
            DiagnosticHandle::Edge(edge) => Some(edge),
            _ => None,
        })
        .collect::<BTreeSet<_>>()
        .len()
}

fn conflict_count_for_operation(
    graph: &TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
    kind: ResourceConflictKind,
) -> usize {
    graph
        .conflicts
        .iter()
        .filter(|conflict| conflict.kind == kind && conflict.operations.contains(&operation))
        .count()
}

fn postcondition_preview_for_candidate(
    graph: &mut TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
    family: ConstructionFamily,
) -> Option<super::PostconditionPredictionId> {
    let (expected_delta, label) = match family {
        ConstructionFamily::DisjointCrossLoopBridge => (
            DiagnosticTopologyDelta {
                vertices: 0,
                edges: 2,
                faces: 1,
                loops: 1,
            },
            "face-family-bridge-postconditions",
        ),
        ConstructionFamily::DisjointSameLoopBridge => (
            DiagnosticTopologyDelta {
                vertices: 0,
                edges: 2,
                faces: 1,
                loops: 1,
            },
            "same-loop-face-bridge-postconditions",
        ),
        ConstructionFamily::PairedSameLoopTwinExposure => (
            DiagnosticTopologyDelta {
                vertices: 0,
                edges: 1,
                faces: 0,
                loops: 0,
            },
            "same-loop-boundary-exposure-postconditions",
        ),
        ConstructionFamily::BoundarySideSeparatorBridge => (
            DiagnosticTopologyDelta {
                vertices: 0,
                edges: 3,
                faces: 1,
                loops: 1,
            },
            "boundary-side-separator-bridge-postconditions",
        ),
        ConstructionFamily::Unknown => return None,
    };
    Some(graph.push_postcondition(operation, expected_delta, label))
}

fn conflicts_for_operation(
    graph: &TopologyDiagnosticGraph,
    operation: DiagnosticOperationId,
) -> Vec<ResourceConflictId> {
    graph
        .conflicts
        .iter()
        .filter(|conflict| conflict.operations.contains(&operation))
        .map(|conflict| conflict.id)
        .collect()
}
#[cfg(test)]
#[path = "bridge_diagnostics_candidates_tests.rs"]
mod tests;
