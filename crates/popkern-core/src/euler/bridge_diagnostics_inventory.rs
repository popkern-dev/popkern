//! Bridge diagnostic resource inventory.

use crate::euler::{
    BridgeDiagnosticRequest, DiagnosticHandle, DiagnosticOwner, DiagnosticResourceId,
    DiagnosticResourceKind, DiagnosticResourceStatus, TopologyDiagnosticGraph,
};
use crate::topology::{BRep, HalfEdgeHandle, LoopHandle};

pub(super) fn populate_bridge_inventory(
    brep: &BRep,
    request: &BridgeDiagnosticRequest,
    mut graph: TopologyDiagnosticGraph,
) -> TopologyDiagnosticGraph {
    for operation in &request.operations {
        let mut consumes = Vec::new();
        consumes.extend(inventory_interval(
            brep,
            operation.first_loop,
            &operation.first_interval,
            &mut graph,
        ));
        consumes.extend(inventory_interval(
            brep,
            operation.second_loop,
            &operation.second_interval,
            &mut graph,
        ));
        graph.push_operation(operation.purpose, consumes);
    }
    graph
}

fn inventory_interval(
    brep: &BRep,
    loop_: LoopHandle,
    interval: &[HalfEdgeHandle],
    graph: &mut TopologyDiagnosticGraph,
) -> Vec<DiagnosticResourceId> {
    let mut consumes = Vec::new();
    let interval_status = if directed_interval_path(brep, loop_, interval).is_some() {
        DiagnosticResourceStatus::Live
    } else {
        DiagnosticResourceStatus::Invalid
    };
    consumes.push(graph.push_resource(
        DiagnosticResourceKind::LoopInterval,
        DiagnosticHandle::Loop(loop_),
        DiagnosticOwner::Loop(loop_),
        interval_status,
    ));
    let Some(path) = directed_interval_path(brep, loop_, interval) else {
        return consumes;
    };
    for &half_edge in path {
        let Some(data) = brep.half_edges.get(half_edge) else {
            consumes.push(graph.push_resource(
                DiagnosticResourceKind::HalfEdgeSide,
                DiagnosticHandle::HalfEdge(half_edge),
                DiagnosticOwner::Unknown,
                DiagnosticResourceStatus::Invalid,
            ));
            continue;
        };
        consumes.push(graph.push_resource(
            DiagnosticResourceKind::HalfEdgeSide,
            DiagnosticHandle::HalfEdge(half_edge),
            DiagnosticOwner::Loop(data.loop_),
            DiagnosticResourceStatus::Live,
        ));
        consumes.push(graph.push_resource(
            DiagnosticResourceKind::EdgeIdentity,
            DiagnosticHandle::Edge(data.edge),
            DiagnosticOwner::Unknown,
            DiagnosticResourceStatus::Live,
        ));
        consumes.push(graph.push_resource(
            DiagnosticResourceKind::Vertex,
            DiagnosticHandle::Vertex(data.origin),
            DiagnosticOwner::Unknown,
            DiagnosticResourceStatus::Live,
        ));
    }
    consumes
}

fn directed_interval_path<'a>(
    brep: &BRep,
    loop_: LoopHandle,
    interval: &'a [HalfEdgeHandle],
) -> Option<&'a [HalfEdgeHandle]> {
    if interval.len() < 2 {
        return None;
    }
    for half_edge in interval {
        let data = brep.half_edges.get(*half_edge)?;
        if data.loop_ != loop_ {
            return None;
        }
    }
    for pair in interval.windows(2) {
        if brep.half_edges.get(pair[0])?.next != pair[1] {
            return None;
        }
    }
    Some(&interval[..interval.len() - 1])
}

#[cfg(test)]
mod tests {
    use crate::euler::{
        BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationPurpose,
        DiagnosticResourceKind, DiagnosticResourceStatus,
    };
    use crate::tolerance::Tolerance;
    use crate::topology::{BRep, HalfEdgeHandle};

    #[test]
    fn bridge_inventory_records_directed_interval_resources() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let loop_ = brep.loops.keys().next().expect("cube has a loop");
        let first = brep.loops[loop_].half_edge;
        let second = brep.half_edges[first].next;
        let third = brep.half_edges[second].next;

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![first, second],
                    second_loop: loop_,
                    second_interval: vec![second, third],
                }],
            })
            .expect("report builds");

        assert_eq!(report.summary.operations, 1);
        assert_eq!(report.summary.conflicts, 0);
        assert_eq!(
            report
                .graph
                .resources
                .iter()
                .filter(|resource| resource.kind == DiagnosticResourceKind::LoopInterval)
                .count(),
            2
        );
        assert_eq!(
            report
                .graph
                .resources
                .iter()
                .filter(|resource| resource.kind == DiagnosticResourceKind::HalfEdgeSide)
                .count(),
            2
        );
        assert!(
            report
                .graph
                .resources
                .iter()
                .all(|resource| resource.status == DiagnosticResourceStatus::Live)
        );
    }

    #[test]
    fn bridge_inventory_marks_invalid_interval_without_failing_report() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let loop_ = brep.loops.keys().next().expect("cube has a loop");
        let first = brep.loops[loop_].half_edge;
        let stale = HalfEdgeHandle::default();

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![first, stale],
                    second_loop: loop_,
                    second_interval: vec![first, stale],
                }],
            })
            .expect("invalid inventory still reports");

        assert_eq!(report.summary.operations, 1);
        assert_eq!(
            report
                .graph
                .resources
                .iter()
                .filter(|resource| resource.status == DiagnosticResourceStatus::Invalid)
                .count(),
            2
        );
    }
}
