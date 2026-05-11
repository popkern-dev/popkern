use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, ConstructionFamily,
    ConstructionFeasibility, ConstructionRequirement, DiagnosticOperationPurpose,
    GeometryFeasibility, OperatorReadinessStatus,
};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

#[test]
fn candidates_report_disjoint_cross_loop_bridge() {
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

    assert!(report.graph.candidates.iter().any(|candidate| {
        candidate.family == ConstructionFamily::DisjointCrossLoopBridge
            && candidate.feasibility == ConstructionFeasibility::ExecutableByExistingOperator
            && candidate.feasibility_reason == "face-family-bridge-existing-operator-ready"
            && candidate.postcondition_preview.is_some()
            && candidate
                .geometry_facts
                .iter()
                .any(|fact| fact.feasibility == GeometryFeasibility::NotApplicable)
            && candidate
                .requirements
                .contains(&ConstructionRequirement::LiveTopologyRevalidation)
            && candidate
                .requirements
                .contains(&ConstructionRequirement::PostconditionVerification)
    }));
}

#[test]
fn candidates_report_disjoint_same_loop_bridge() {
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

    assert!(report.graph.candidates.iter().any(|candidate| {
        candidate.family == ConstructionFamily::DisjointSameLoopBridge
            && candidate.feasibility == ConstructionFeasibility::ExecutableByExistingOperator
            && candidate.feasibility_reason == "same-loop-face-bridge-existing-operator-ready"
            && candidate.postcondition_preview.is_some()
    }));
    assert_eq!(report.summary.disjoint_same_loop_bridge_candidates, 1);
    assert_eq!(report.summary.executable_existing_operator_candidates, 1);
    assert_eq!(report.summary.candidates_with_postcondition_preview, 1);
    assert_eq!(report.summary.predicted_delta_edges, 2);
    assert_eq!(report.summary.predicted_delta_faces, 1);
    assert_eq!(report.summary.predicted_delta_loops, 1);
}

#[test]
fn candidates_report_boundary_side_separator_family() {
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

    assert!(report.graph.candidates.iter().any(|candidate| {
        candidate.family == ConstructionFamily::BoundarySideSeparatorBridge
            && candidate.feasibility == ConstructionFeasibility::RequiresSeparatorSubcomplex
            && candidate.feasibility_reason == "missing-boundary-side-separator-subcomplex"
            && candidate.postcondition_preview.is_some()
            && !candidate.conflicts.is_empty()
    }));
    assert_eq!(report.summary.boundary_side_separator_bridge_candidates, 1);
    assert_eq!(report.summary.requires_separator_subcomplex_candidates, 1);
    assert_eq!(report.summary.predicted_delta_edges, 3);
    assert_eq!(report.summary.predicted_delta_faces, 1);
    assert_eq!(report.summary.predicted_delta_loops, 1);
}

#[test]
fn candidates_report_multiple_boundary_side_overlaps_as_regrouping() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let first = brep.loops[loop_].half_edge;
    let second = brep.half_edges[first].next;
    let third = brep.half_edges[second].next;

    let report = brep
        .bridge_diagnostic_report(BridgeDiagnosticRequest {
            operations: vec![BridgeDiagnosticOperationRequest {
                purpose: DiagnosticOperationPurpose::BridgeCycle,
                first_loop: loop_,
                first_interval: vec![first, second, third],
                second_loop: loop_,
                second_interval: vec![first, second, third],
            }],
        })
        .expect("report builds");

    assert!(report.graph.candidates.iter().any(|candidate| {
        candidate.family == ConstructionFamily::BoundarySideSeparatorBridge
            && candidate.feasibility == ConstructionFeasibility::RequiresCycleRegrouping
            && candidate.feasibility_reason
                == "multiple-boundary-side-overlaps-require-cycle-regrouping"
    }));
    assert_eq!(report.summary.requires_cycle_regrouping_candidates, 1);
}

#[test]
fn matrix_rejects_supported_readiness_without_candidate() {
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

    for readiness in &report.graph.readiness {
        if matches!(
            readiness.status,
            OperatorReadinessStatus::Ready
                | OperatorReadinessStatus::ReadyAfterNormalization
                | OperatorReadinessStatus::MissingOperator
        ) {
            assert!(
                report.graph.candidates.iter().any(|candidate| {
                    candidate.operations.contains(&readiness.operation)
                        && candidate.postcondition_preview.is_some()
                }),
                "readiness entry {readiness:?} has no construction candidate"
            );
        }
    }
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
            Point3::new(4.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    brep.shells[shell].face = face;
    brep
}

fn add_boundary_polygon(
    brep: &mut BRep,
    shell: crate::topology::ShellHandle,
    points: &[Point3],
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
    for point in points {
        let vertex = brep.vertices.insert(Vertex {
            position: *point,
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
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });
        brep.vertices[vertex].outgoing = half_edge;
        brep.edges[edge].half_edge = half_edge;
        brep.boundary_edges.insert(edge);
        half_edges.push(half_edge);
    }
    for index in 0..half_edges.len() {
        let current = half_edges[index];
        let next = half_edges[(index + 1) % half_edges.len()];
        let prev = half_edges[(index + half_edges.len() - 1) % half_edges.len()];
        brep.half_edges[current].next = next;
        brep.half_edges[current].prev = prev;
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}

fn nth_half_edge(brep: &BRep, start: HalfEdgeHandle, n: usize) -> HalfEdgeHandle {
    let mut current = start;
    for _ in 0..n {
        current = brep.half_edges[current].next;
    }
    current
}
