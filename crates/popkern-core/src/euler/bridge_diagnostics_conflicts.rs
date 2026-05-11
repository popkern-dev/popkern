//! Bridge diagnostic resource conflict classification.

use std::collections::{BTreeMap, BTreeSet};

use crate::euler::{
    DiagnosticHandle, DiagnosticOperationId, DiagnosticResourceId, DiagnosticResourceKind,
    ResourceConflictKind, ResourceConflictStatus, TopologyDiagnosticGraph,
};

pub(super) fn classify_duplicate_consumption(
    mut graph: TopologyDiagnosticGraph,
) -> TopologyDiagnosticGraph {
    emit_duplicate_conflicts(
        &mut graph,
        DiagnosticResourceKind::HalfEdgeSide,
        ResourceConflictKind::DuplicateHalfEdgeConsumption,
    );
    emit_duplicate_conflicts(
        &mut graph,
        DiagnosticResourceKind::EdgeIdentity,
        ResourceConflictKind::DuplicateBoundaryEdgeConsumption,
    );
    graph
}

fn emit_duplicate_conflicts(
    graph: &mut TopologyDiagnosticGraph,
    resource_kind: DiagnosticResourceKind,
    conflict_kind: ResourceConflictKind,
) {
    let mut resources_by_handle = BTreeMap::<DiagnosticHandle, Vec<DiagnosticResourceId>>::new();
    for resource in &graph.resources {
        if resource.kind == resource_kind {
            resources_by_handle
                .entry(resource.handle)
                .or_default()
                .push(resource.id);
        }
    }
    for resource_ids in resources_by_handle.into_values() {
        if resource_ids.len() < 2 {
            continue;
        }
        let operation_ids = operations_for_resources(graph, &resource_ids);
        graph.push_conflict(
            conflict_kind,
            resource_ids,
            operation_ids,
            ResourceConflictStatus::Unknown,
        );
    }
}

fn operations_for_resources(
    graph: &TopologyDiagnosticGraph,
    resource_ids: &[DiagnosticResourceId],
) -> Vec<DiagnosticOperationId> {
    let resource_set = resource_ids.iter().copied().collect::<BTreeSet<_>>();
    let mut operation_ids = BTreeSet::new();
    for operation in &graph.operations {
        if operation
            .consumes
            .iter()
            .any(|resource_id| resource_set.contains(resource_id))
        {
            operation_ids.insert(operation.id);
        }
    }
    operation_ids.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use crate::euler::{
        BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationPurpose,
        ResourceConflictKind,
    };
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn duplicate_consumption_detects_reused_half_edge_side() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let loop_ = brep.loops.keys().next().expect("cube has a loop");
        let first = brep.loops[loop_].half_edge;
        let second = brep.half_edges[first].next;

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![first, second],
                    second_loop: loop_,
                    second_interval: vec![first, second],
                }],
            })
            .expect("report builds");

        assert!(report.graph.conflicts.iter().any(|conflict| {
            conflict.kind == ResourceConflictKind::DuplicateHalfEdgeConsumption
        }));
    }

    #[test]
    fn duplicate_consumption_detects_reused_edge_identity() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let edge = brep.edges.keys().next().expect("cube has edges");
        let first = brep.edges[edge].half_edge;
        let second = brep.half_edges[first].twin;
        let first_loop = brep.half_edges[first].loop_;
        let second_loop = brep.half_edges[second].loop_;
        let first_after = brep.half_edges[first].next;
        let second_after = brep.half_edges[second].next;

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop,
                    first_interval: vec![first, first_after],
                    second_loop,
                    second_interval: vec![second, second_after],
                }],
            })
            .expect("report builds");

        assert!(report.graph.conflicts.iter().any(|conflict| {
            conflict.kind == ResourceConflictKind::DuplicateBoundaryEdgeConsumption
        }));
        assert!(!report.graph.conflicts.iter().any(|conflict| {
            conflict.kind == ResourceConflictKind::DuplicateHalfEdgeConsumption
        }));
    }
}
