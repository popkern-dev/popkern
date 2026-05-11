//! Bridge diagnostic operator readiness.

use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationId,
    DiagnosticOperatorFamily, FaceFamilyBridgeCycleRequest, OperatorReadinessStatus,
    ResourceConflictKind, SameLoopFaceBridgeRequest, TopologyDiagnosticGraph,
};
use crate::topology::BRep;

pub(super) fn classify_operator_readiness(
    brep: &BRep,
    request: &BridgeDiagnosticRequest,
    mut graph: TopologyDiagnosticGraph,
) -> TopologyDiagnosticGraph {
    for (index, operation) in request.operations.iter().enumerate() {
        let operation_id = DiagnosticOperationId(index);
        classify_operation_readiness(brep, operation_id, operation, &mut graph);
    }
    graph
}

fn classify_operation_readiness(
    brep: &BRep,
    operation_id: DiagnosticOperationId,
    operation: &BridgeDiagnosticOperationRequest,
    graph: &mut TopologyDiagnosticGraph,
) {
    if operation_has_conflict(
        graph,
        operation_id,
        ResourceConflictKind::AlreadyBoundaryOverlap,
    ) {
        graph.push_readiness(
            operation_id,
            DiagnosticOperatorFamily::BoundarySideAllocationOrSeparator,
            OperatorReadinessStatus::MissingOperator,
            "already-boundary-side-selected-by-multiple-bridge-cycles",
        );
        return;
    }
    if operation_has_conflict(
        graph,
        operation_id,
        ResourceConflictKind::PairedSideNeedsExposure,
    ) {
        graph.push_readiness(
            operation_id,
            DiagnosticOperatorFamily::SameLoopBoundaryExposure,
            OperatorReadinessStatus::ReadyAfterNormalization,
            "paired-same-loop-side-needs-boundary-exposure",
        );
        return;
    }
    if operation.first_loop == operation.second_loop {
        let status = match brep.preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_: operation.first_loop,
            first_interval: operation.first_interval.clone(),
            second_interval: operation.second_interval.clone(),
        }) {
            Ok(_) => OperatorReadinessStatus::Ready,
            Err(_) => OperatorReadinessStatus::InvalidInput,
        };
        graph.push_readiness(
            operation_id,
            DiagnosticOperatorFamily::SameLoopFaceBridge,
            status,
            same_loop_reason(status),
        );
        return;
    }
    let status = match brep.preflight_face_family_bridge(vec![FaceFamilyBridgeCycleRequest {
        first_loop: operation.first_loop,
        first_interval: operation.first_interval.clone(),
        second_loop: operation.second_loop,
        second_interval: operation.second_interval.clone(),
    }]) {
        Ok(_) => OperatorReadinessStatus::Ready,
        Err(_) => OperatorReadinessStatus::InvalidInput,
    };
    graph.push_readiness(
        operation_id,
        DiagnosticOperatorFamily::FaceFamilyBridge,
        status,
        face_family_reason(status),
    );
}

fn operation_has_conflict(
    graph: &TopologyDiagnosticGraph,
    operation_id: DiagnosticOperationId,
    kind: ResourceConflictKind,
) -> bool {
    graph
        .conflicts
        .iter()
        .any(|conflict| conflict.kind == kind && conflict.operations.contains(&operation_id))
}

fn same_loop_reason(status: OperatorReadinessStatus) -> &'static str {
    match status {
        OperatorReadinessStatus::Ready => "same-loop-face-bridge-ready",
        OperatorReadinessStatus::InvalidInput => "same-loop-face-bridge-invalid-input",
        _ => "same-loop-face-bridge-not-evaluated",
    }
}

fn face_family_reason(status: OperatorReadinessStatus) -> &'static str {
    match status {
        OperatorReadinessStatus::Ready => "face-family-bridge-ready",
        OperatorReadinessStatus::InvalidInput => "face-family-bridge-invalid-input",
        _ => "face-family-bridge-not-evaluated",
    }
}

#[cfg(test)]
mod tests {
    use crate::euler::{
        BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, DiagnosticOperationPurpose,
        DiagnosticOperatorFamily, OperatorReadinessStatus,
    };
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

    #[test]
    fn readiness_reports_same_loop_bridge_ready() {
        let brep = same_loop_boundary_fixture();
        let face = brep.faces.keys().next().expect("fixture has a face");
        let loop_ = brep.faces[face].outer_loop;
        let first = brep.loops[loop_].half_edge;
        let first_after = brep.half_edges[first].next;
        let second = nth_half_edge(&brep, first, 3);
        let second_after = brep.half_edges[second].next;

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![first, first_after],
                    second_loop: loop_,
                    second_interval: vec![second, second_after],
                }],
            })
            .expect("report builds");

        assert!(report.graph.readiness.iter().any(|readiness| {
            readiness.operator == DiagnosticOperatorFamily::SameLoopFaceBridge
                && readiness.status == OperatorReadinessStatus::Ready
        }));
    }

    #[test]
    fn readiness_reports_already_boundary_overlap_missing_operator() {
        let mut brep = same_loop_boundary_fixture();
        let face = brep.faces.keys().next().expect("fixture has a face");
        let loop_ = brep.faces[face].outer_loop;
        let first = brep.loops[loop_].half_edge;
        let first_after = brep.half_edges[first].next;
        let edge = brep.half_edges[first].edge;
        brep.boundary_edges.insert(edge);

        let report = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![first, first_after],
                    second_loop: loop_,
                    second_interval: vec![first, first_after],
                }],
            })
            .expect("report builds");

        assert!(report.graph.readiness.iter().any(|readiness| {
            readiness.operator == DiagnosticOperatorFamily::BoundarySideAllocationOrSeparator
                && readiness.status == OperatorReadinessStatus::MissingOperator
        }));
    }

    #[test]
    fn readiness_reports_cross_loop_bridge_ready() {
        let brep = two_boundary_loop_fixture();
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let first_loop = brep.faces[faces[0]].outer_loop;
        let second_loop = brep.faces[faces[1]].outer_loop;
        let first = brep.loops[first_loop].half_edge;
        let second = brep.loops[second_loop].half_edge;
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

        assert!(report.graph.readiness.iter().any(|readiness| {
            readiness.operator == DiagnosticOperatorFamily::FaceFamilyBridge
                && readiness.status == OperatorReadinessStatus::Ready
        }));
    }

    fn two_boundary_loop_fixture() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let solid = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell = brep.shells.insert(Shell {
            face: Default::default(),
            solid,
        });
        brep.solids[solid].shell = shell;
        let first = add_boundary_polygon(
            &mut brep,
            shell,
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
        );
        add_boundary_polygon(
            &mut brep,
            shell,
            &[
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
                Point3::new(3.0, 1.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
            ],
        );
        brep.shells[shell].face = first;
        brep
    }

    fn same_loop_boundary_fixture() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let solid = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell = brep.shells.insert(Shell {
            face: Default::default(),
            solid,
        });
        brep.solids[solid].shell = shell;
        let face = add_boundary_polygon(
            &mut brep,
            shell,
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.0),
                Point3::new(5.0, 0.0, 0.0),
            ],
        );
        brep.shells[shell].face = face;
        brep
    }

    fn add_boundary_polygon(
        brep: &mut BRep,
        shell: crate::topology::ShellHandle,
        positions: &[Point3],
    ) -> crate::topology::FaceHandle {
        let face = brep.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let loop_ = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });
        brep.faces[face].outer_loop = loop_;
        let mut half_edges = Vec::new();
        for position in positions {
            let vertex = brep.vertices.insert(Vertex {
                position: *position,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            });
            let edge = brep.edges.insert(Edge {
                half_edge: Default::default(),
                curve: None,
                name: None,
                tolerance: None,
                canonical: Default::default(),
                role: Default::default(),
            });
            let half_edge = brep.half_edges.insert(HalfEdge {
                origin: vertex,
                twin: HalfEdgeHandle::default(),
                next: Default::default(),
                prev: Default::default(),
                edge,
                loop_,
                pcurve: None,
            });
            brep.edges[edge].half_edge = half_edge;
            brep.vertices[vertex].outgoing = half_edge;
            brep.boundary_edges.insert(edge);
            half_edges.push(half_edge);
        }
        for index in 0..half_edges.len() {
            brep.half_edges[half_edges[index]].next = half_edges[(index + 1) % half_edges.len()];
            brep.half_edges[half_edges[index]].prev =
                half_edges[(index + half_edges.len() - 1) % half_edges.len()];
        }
        brep.loops[loop_].half_edge = half_edges[0];
        face
    }

    fn nth_half_edge(brep: &BRep, start: HalfEdgeHandle, offset: usize) -> HalfEdgeHandle {
        let mut cursor = start;
        for _ in 0..offset {
            cursor = brep.half_edges[cursor].next;
        }
        cursor
    }
}
