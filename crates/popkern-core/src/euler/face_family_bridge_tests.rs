use crate::euler::FaceFamilyBridgeCycleRequest;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

#[test]
fn face_family_bridge_commits_disjoint_boundary_cycles_transactionally() {
    let mut brep = boundary_bridge_fixture();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let outside_loop = brep.faces[faces[2]].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let second_start = brep.loops[second_loop].half_edge;
    let outside_start = brep.loops[outside_loop].half_edge;
    let outside_next = brep.half_edges[outside_start].next;
    let outside_second =
        brep.half_edges[brep.half_edges[brep.half_edges[outside_next].next].next].next;
    let before = topology_counts(&brep);

    let plan = brep
        .preflight_face_family_bridge(vec![
            FaceFamilyBridgeCycleRequest {
                first_loop,
                first_interval: vec![first_start, brep.half_edges[first_start].next],
                second_loop: outside_loop,
                second_interval: vec![outside_next, outside_start],
            },
            FaceFamilyBridgeCycleRequest {
                first_loop: second_loop,
                first_interval: vec![second_start, brep.half_edges[second_start].next],
                second_loop: outside_loop,
                second_interval: vec![outside_second, brep.half_edges[outside_second].next],
            },
        ])
        .expect("disjoint aggregate cycles preflight");

    let result = brep
        .face_family_bridge(&plan)
        .expect("aggregate bridge commits");

    assert_eq!(result.cycles.len(), 2);
    assert_eq!(
        (
            brep.vertices.len() as i64 - before.0 as i64,
            brep.edges.len() as i64 - before.1 as i64,
            brep.faces.len() as i64 - before.2 as i64,
            brep.loops.len() as i64 - before.3 as i64,
        ),
        (0, 4, 2, 2)
    );
    brep.validate_half_edge_pairing().unwrap();
    brep.validate_half_edge_cycles().unwrap();
    brep.validate_manifold_edges().unwrap();
}

#[test]
fn face_family_bridge_refuses_overlapping_boundary_consumption() {
    let brep = boundary_bridge_fixture();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let outside_loop = brep.faces[faces[2]].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let second_start = brep.loops[second_loop].half_edge;
    let outside_start = brep.loops[outside_loop].half_edge;
    let outside_next = brep.half_edges[outside_start].next;

    let err = brep
        .preflight_face_family_bridge(vec![
            FaceFamilyBridgeCycleRequest {
                first_loop,
                first_interval: vec![first_start, brep.half_edges[first_start].next],
                second_loop: outside_loop,
                second_interval: vec![outside_next, outside_start],
            },
            FaceFamilyBridgeCycleRequest {
                first_loop: second_loop,
                first_interval: vec![second_start, brep.half_edges[second_start].next],
                second_loop: outside_loop,
                second_interval: vec![outside_next, outside_start],
            },
        ])
        .expect_err("overlapping outside boundary edge is non-manifold");

    assert!(err.to_string().contains("cycles overlap on boundary edge"));
}

fn boundary_bridge_fixture() -> BRep {
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
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(1.0, 0.0, 2.0),
            Point3::new(2.0, 0.0, 2.0),
            Point3::new(3.0, 0.0, 2.0),
            Point3::new(4.0, 0.0, 2.0),
            Point3::new(5.0, 0.0, 2.0),
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
