use crate::euler::{
    BoundarySideSeparatorBridgeRequest, BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest,
    ConstructionFamily, FaceFamilyBridgeCycleRequest,
};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

#[test]
fn boundary_side_separator_preflight_accepts_one_already_boundary_overlap() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let first = brep.loops[loop_].half_edge;
    let after_first = brep.half_edges[first].next;

    let plan = brep
        .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: vec![FaceFamilyBridgeCycleRequest {
                first_loop: loop_,
                first_interval: vec![first, after_first],
                second_loop: loop_,
                second_interval: vec![first, after_first],
            }],
        })
        .expect("one already-boundary overlap preflights");

    assert_eq!(plan.separator_source_edge, brep.half_edges[first].edge);
    assert_eq!(plan.separator_source_half_edge, first);
    assert_eq!(
        (
            plan.expected_delta.vertices,
            plan.expected_delta.edges,
            plan.expected_delta.faces,
            plan.expected_delta.loops,
        ),
        (0, 3, 1, 1)
    );
}

#[test]
fn boundary_side_separator_diagnostic_candidate_is_gate_complete() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let first = brep.loops[loop_].half_edge;
    let after_first = brep.half_edges[first].next;

    let report = brep
        .bridge_diagnostic_report(BridgeDiagnosticRequest {
            operations: vec![BridgeDiagnosticOperationRequest {
                purpose: crate::euler::DiagnosticOperationPurpose::BridgeCycle,
                first_loop: loop_,
                first_interval: vec![first, after_first],
                second_loop: loop_,
                second_interval: vec![first, after_first],
            }],
        })
        .expect("diagnostic report builds");
    let candidate = report
        .graph
        .candidates
        .iter()
        .find(|candidate| candidate.family == ConstructionFamily::BoundarySideSeparatorBridge)
        .expect("diagnostics include separator candidate");

    candidate
        .validate_mutation_gate()
        .expect("separator candidate has complete gate facts");
    assert_eq!(report.summary.missing_operators, 1);
    assert_eq!(
        report.summary.missing_primitive_operation_dependency_facts,
        1
    );
    assert_eq!(report.summary.downstream_blocked_missing_primitive_facts, 1);
}

#[test]
fn boundary_side_separator_preflight_accepts_live_shared_start_shape() {
    let brep = separator_bridge_fixture();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let outside_loop = brep.faces[faces[2]].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let second_start = brep.loops[second_loop].half_edge;
    let outside_start = brep.loops[outside_loop].half_edge;
    let outside_second = brep.half_edges[outside_start].next;
    let outside_third = brep.half_edges[outside_second].next;

    let plan = brep
        .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: vec![
                FaceFamilyBridgeCycleRequest {
                    first_loop,
                    first_interval: vec![first_start, brep.half_edges[first_start].next],
                    second_loop: outside_loop,
                    second_interval: vec![outside_start, outside_second],
                },
                FaceFamilyBridgeCycleRequest {
                    first_loop: second_loop,
                    first_interval: vec![second_start, brep.half_edges[second_start].next],
                    second_loop: outside_loop,
                    second_interval: vec![outside_second, outside_third],
                },
                FaceFamilyBridgeCycleRequest {
                    first_loop: outside_loop,
                    first_interval: vec![outside_start, outside_second],
                    second_loop: outside_loop,
                    second_interval: vec![outside_start, outside_second],
                },
            ],
        })
        .expect("live shared-start separator shape preflights");

    assert_eq!(plan.cycles.len(), 3);
    assert_eq!(
        plan.separator_source_edge,
        brep.half_edges[outside_start].edge
    );
    assert_eq!(
        (
            plan.expected_delta.vertices,
            plan.expected_delta.edges,
            plan.expected_delta.faces,
            plan.expected_delta.loops,
        ),
        (0, 7, 3, 3)
    );
}

#[test]
fn boundary_side_separator_bridge_commits_live_shared_start_shape() {
    let mut brep = separator_bridge_fixture();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let outside_loop = brep.faces[faces[2]].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let second_start = brep.loops[second_loop].half_edge;
    let outside_start = brep.loops[outside_loop].half_edge;
    let outside_second = brep.half_edges[outside_start].next;
    let outside_third = brep.half_edges[outside_second].next;
    let before = topology_counts(&brep);
    let plan = brep
        .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: vec![
                FaceFamilyBridgeCycleRequest {
                    first_loop,
                    first_interval: vec![first_start, brep.half_edges[first_start].next],
                    second_loop: outside_loop,
                    second_interval: vec![outside_start, outside_second],
                },
                FaceFamilyBridgeCycleRequest {
                    first_loop: second_loop,
                    first_interval: vec![second_start, brep.half_edges[second_start].next],
                    second_loop: outside_loop,
                    second_interval: vec![outside_second, outside_third],
                },
                FaceFamilyBridgeCycleRequest {
                    first_loop: outside_loop,
                    first_interval: vec![outside_start, outside_second],
                    second_loop: outside_loop,
                    second_interval: vec![outside_start, outside_second],
                },
            ],
        })
        .expect("live shared-start separator shape preflights");

    let result = brep
        .boundary_side_separator_bridge(&plan)
        .expect("separator bridge commits transactionally");

    assert_eq!(result.bridge_faces.len(), 3);
    assert_eq!(result.bridge_loops.len(), 3);
    assert!(!brep.boundary_edges.contains(&result.separator_edge));
    assert_eq!(
        topology_delta(&brep, before),
        (
            plan.expected_delta.vertices,
            plan.expected_delta.edges,
            plan.expected_delta.faces,
            plan.expected_delta.loops,
        )
    );
    brep.validate_half_edge_pairing().unwrap();
    brep.validate_half_edge_cycles().unwrap();
    brep.validate_manifold_edges().unwrap();
}

#[test]
fn boundary_side_separator_rejects_multiple_already_boundary_overlaps() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let first = brep.loops[loop_].half_edge;
    let second = brep.half_edges[first].next;
    let third = brep.half_edges[second].next;

    let err = brep
        .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: vec![
                FaceFamilyBridgeCycleRequest {
                    first_loop: loop_,
                    first_interval: vec![first, second],
                    second_loop: loop_,
                    second_interval: vec![first, second],
                },
                FaceFamilyBridgeCycleRequest {
                    first_loop: loop_,
                    first_interval: vec![second, third],
                    second_loop: loop_,
                    second_interval: vec![second, third],
                },
            ],
        })
        .expect_err("multiple already-boundary overlaps are out of scope");

    assert!(
        err.to_string()
            .contains("requires exactly one already-boundary duplicate side")
    );
}

#[test]
fn boundary_side_separator_rejects_non_boundary_overlap() {
    let mut brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let first = brep.loops[loop_].half_edge;
    let after_first = brep.half_edges[first].next;
    let edge = brep.half_edges[first].edge;
    brep.boundary_edges.remove(&edge);

    let err = brep
        .preflight_boundary_side_separator_bridge(BoundarySideSeparatorBridgeRequest {
            cycles: vec![FaceFamilyBridgeCycleRequest {
                first_loop: loop_,
                first_interval: vec![first, after_first],
                second_loop: loop_,
                second_interval: vec![first, after_first],
            }],
        })
        .expect_err("non-boundary overlap is not a separator case");

    assert!(
        err.to_string()
            .contains("requires a typed construction candidate")
    );
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
        ],
    );
    brep.shells[shell].face = face;
    brep
}

fn separator_bridge_fixture() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    add_boundary_polygon(
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
    let outside = add_boundary_polygon(
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
    brep.shells[shell].face = outside;
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

fn topology_counts(brep: &BRep) -> (usize, usize, usize, usize) {
    (
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
        brep.loops.len(),
    )
}

fn topology_delta(brep: &BRep, before: (usize, usize, usize, usize)) -> (i64, i64, i64, i64) {
    (
        brep.vertices.len() as i64 - before.0 as i64,
        brep.edges.len() as i64 - before.1 as i64,
        brep.faces.len() as i64 - before.2 as i64,
        brep.loops.len() as i64 - before.3 as i64,
    )
}
