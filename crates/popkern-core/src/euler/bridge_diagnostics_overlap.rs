//! Bridge diagnostic overlap-shape classification.

use std::collections::{BTreeMap, BTreeSet};

use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticHandle,
    DiagnosticOperationId, DiagnosticResourceId, DiagnosticResourceKind, ResourceConflictKind,
    ResourceConflictStatus, TopologyDiagnosticGraph,
};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

pub(super) fn classify_overlap_shapes(
    brep: &BRep,
    request: &BridgeDiagnosticRequest,
    mut graph: TopologyDiagnosticGraph,
) -> TopologyDiagnosticGraph {
    for (index, operation) in request.operations.iter().enumerate() {
        let operation_id = DiagnosticOperationId(index);
        classify_operation_overlap(brep, operation_id, operation, &mut graph);
    }
    graph
}

fn classify_operation_overlap(
    brep: &BRep,
    operation_id: DiagnosticOperationId,
    operation: &BridgeDiagnosticOperationRequest,
    graph: &mut TopologyDiagnosticGraph,
) {
    if operation.first_loop != operation.second_loop {
        return;
    }
    let Some(first) = bridge_path_edges(brep, &operation.first_interval) else {
        return;
    };
    let Some(second) = bridge_path_edges(brep, &operation.second_interval) else {
        return;
    };
    let overlap_edges = first
        .set
        .intersection(&second.set)
        .copied()
        .collect::<Vec<_>>();
    if overlap_edges.is_empty() {
        return;
    }
    let kind = overlap_kind(brep, &first, &second, &overlap_edges);
    let resources = edge_resource_ids(graph, &overlap_edges);
    graph.push_conflict(
        kind,
        resources,
        vec![operation_id],
        ResourceConflictStatus::Unknown,
    );
}

fn overlap_kind(
    brep: &BRep,
    first: &BridgePathEdges,
    second: &BridgePathEdges,
    overlap_edges: &[EdgeHandle],
) -> ResourceConflictKind {
    if overlap_edges
        .iter()
        .any(|edge| brep.boundary_edges.contains(edge))
    {
        return ResourceConflictKind::AlreadyBoundaryOverlap;
    }
    if overlap_edges
        .iter()
        .any(|edge| edge_has_same_loop_twin_pair(brep, *edge))
    {
        return ResourceConflictKind::PairedSideNeedsExposure;
    }
    if first.repeated_edges || second.repeated_edges {
        return ResourceConflictKind::RepeatedEdgeOccurrence;
    }
    if overlap_edges.len() > 1 && (first.first == second.first || first.last == second.last) {
        ResourceConflictKind::AlignedOverlap
    } else if first.first == second.first {
        ResourceConflictKind::SharedStart
    } else if first.last == second.last {
        ResourceConflictKind::SharedEnd
    } else if first.first == second.last || first.last == second.first {
        ResourceConflictKind::CrossedEndpoint
    } else if one_path_contains_the_other(first, second) {
        ResourceConflictKind::NestedOverlap
    } else {
        ResourceConflictKind::InteriorOverlap
    }
}

fn edge_has_same_loop_twin_pair(brep: &BRep, edge: EdgeHandle) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(half_edge.twin) else {
        return false;
    };
    half_edge.edge == edge && twin.edge == edge && half_edge.loop_ == twin.loop_
}

fn one_path_contains_the_other(first: &BridgePathEdges, second: &BridgePathEdges) -> bool {
    (first.set.len() < second.set.len() && first.set.is_subset(&second.set))
        || (second.set.len() < first.set.len() && second.set.is_subset(&first.set))
}

fn edge_resource_ids(
    graph: &TopologyDiagnosticGraph,
    overlap_edges: &[EdgeHandle],
) -> Vec<DiagnosticResourceId> {
    let overlap_edges = overlap_edges.iter().copied().collect::<BTreeSet<_>>();
    graph
        .resources
        .iter()
        .filter_map(|resource| {
            if resource.kind != DiagnosticResourceKind::EdgeIdentity {
                return None;
            }
            let DiagnosticHandle::Edge(edge) = resource.handle else {
                return None;
            };
            overlap_edges.contains(&edge).then_some(resource.id)
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BridgePathEdges {
    set: BTreeSet<EdgeHandle>,
    first: Option<EdgeHandle>,
    last: Option<EdgeHandle>,
    repeated_edges: bool,
}

fn bridge_path_edges(brep: &BRep, interval: &[HalfEdgeHandle]) -> Option<BridgePathEdges> {
    if interval.len() < 2 {
        return None;
    }
    let mut set = BTreeSet::new();
    let mut counts = BTreeMap::<EdgeHandle, usize>::new();
    let mut first = None;
    let mut last = None;
    for &half_edge in &interval[..interval.len() - 1] {
        let edge = brep.half_edges.get(half_edge)?.edge;
        if first.is_none() {
            first = Some(edge);
        }
        last = Some(edge);
        set.insert(edge);
        *counts.entry(edge).or_default() += 1;
    }
    Some(BridgePathEdges {
        set,
        first,
        last,
        repeated_edges: counts.values().any(|count| *count > 1),
    })
}

#[cfg(test)]
mod tests {
    use crate::euler::{
        BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationPurpose,
        ResourceConflictKind,
    };
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn overlap_classification_detects_shared_start() {
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
                    second_interval: vec![first, second, third],
                }],
            })
            .expect("report builds");

        assert!(
            report
                .graph
                .conflicts
                .iter()
                .any(|conflict| conflict.kind == ResourceConflictKind::SharedStart)
        );
    }

    #[test]
    fn overlap_classification_detects_already_boundary_overlap() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let loop_ = brep.loops.keys().next().expect("cube has a loop");
        let first = brep.loops[loop_].half_edge;
        let second = brep.half_edges[first].next;
        let edge = brep.half_edges[first].edge;
        brep.boundary_edges.insert(edge);

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

        assert!(
            report
                .graph
                .conflicts
                .iter()
                .any(|conflict| conflict.kind == ResourceConflictKind::AlreadyBoundaryOverlap)
        );
    }

    #[test]
    fn overlap_classification_distinguishes_paired_side_exposure() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let outgoing = brep.vertices[vertex].outgoing;
        let spur = brep
            .mev(outgoing, Point3::new(0.5, 0.5, 0.5))
            .expect("same-loop spur is created");
        let loop_ = brep.half_edges[spur.half_edge].loop_;
        let first_after = brep.half_edges[spur.half_edge].next;
        let second_after = brep.half_edges[spur.half_edge_twin].next;

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![spur.half_edge, first_after],
                    second_loop: loop_,
                    second_interval: vec![spur.half_edge_twin, second_after],
                }],
            })
            .expect("report builds");

        assert!(
            report
                .graph
                .conflicts
                .iter()
                .any(|conflict| conflict.kind == ResourceConflictKind::PairedSideNeedsExposure)
        );
        assert!(
            !report
                .graph
                .conflicts
                .iter()
                .any(|conflict| conflict.kind == ResourceConflictKind::AlreadyBoundaryOverlap)
        );
    }
}
