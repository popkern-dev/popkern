use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticRequest, CandidateMissingFacet,
    ConstructionFamily, DiagnosticOperationPurpose,
};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

#[test]
fn complete_candidate_passes_mutation_gate_contract() {
    let report = cross_loop_report();
    let candidate = report
        .graph
        .candidates
        .iter()
        .find(|candidate| candidate.family == ConstructionFamily::DisjointCrossLoopBridge)
        .expect("fixture emits cross-loop candidate");

    candidate
        .validate_mutation_gate()
        .expect("complete candidate satisfies mutation gate");
}

#[test]
fn mutation_gate_reports_missing_geometry_facet() {
    let mut candidate = cross_loop_candidate();
    candidate.geometry_facts.clear();

    let error = candidate
        .validate_mutation_gate()
        .expect_err("missing geometry facts reject mutation");

    assert_eq!(error.candidate, candidate.id);
    assert_eq!(
        error.missing_facet,
        CandidateMissingFacet::GeometryFeasibility
    );
    assert_eq!(
        error.missing_facet.diagnostic_label(),
        "geometry-feasibility"
    );
}

#[test]
fn mutation_gate_reports_missing_semantic_history_facet() {
    let mut candidate = cross_loop_candidate();
    candidate
        .operation_dependencies
        .retain(|dependency| !dependency.primitive.contains("history"));

    let error = candidate
        .validate_mutation_gate()
        .expect_err("missing history dependencies reject mutation");

    assert_eq!(
        error.missing_facet,
        CandidateMissingFacet::SemanticHistoryDependency
    );
}

#[test]
fn mutation_gate_reports_missing_downstream_readiness_facet() {
    let mut candidate = cross_loop_candidate();
    candidate.downstream_facts.clear();

    let error = candidate
        .validate_mutation_gate()
        .expect_err("missing downstream facts reject mutation");

    assert_eq!(
        error.missing_facet,
        CandidateMissingFacet::DownstreamReadiness
    );
}

fn cross_loop_candidate() -> crate::euler::ConstructionCandidate {
    cross_loop_report()
        .graph
        .candidates
        .into_iter()
        .find(|candidate| candidate.family == ConstructionFamily::DisjointCrossLoopBridge)
        .expect("fixture emits cross-loop candidate")
}

fn cross_loop_report() -> crate::euler::BridgeDiagnosticReport {
    let brep = two_boundary_loop_fixture();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let first = brep.loops[first_loop].half_edge;
    let second = brep.loops[second_loop].half_edge;
    let first_after = brep.half_edges[first].next;
    let second_after = brep.half_edges[second].next;

    brep.bridge_diagnostic_report(BridgeDiagnosticRequest {
        operations: vec![BridgeDiagnosticOperationRequest {
            purpose: DiagnosticOperationPurpose::BridgeCycle,
            first_loop,
            first_interval: vec![first, first_after],
            second_loop,
            second_interval: vec![second, second_after],
        }],
    })
    .expect("report builds")
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
