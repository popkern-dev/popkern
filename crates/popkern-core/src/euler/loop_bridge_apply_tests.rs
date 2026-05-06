use crate::euler::TwoLoopBridgeRequest;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

#[test]
fn two_loop_bridge_commits_between_two_boundary_face_loops() {
    let mut brep = two_boundary_quads();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let first_loop = brep.faces[faces[0]].outer_loop;
    let second_loop = brep.faces[faces[1]].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let first_end = brep.half_edges[brep.half_edges[first_start].next].next;
    let second_start = brep.loops[second_loop].half_edge;
    let second_end = brep.half_edges[brep.half_edges[second_start].next].next;
    let before_counts = topology_counts(&brep);

    let plan = brep
        .preflight_two_loop_bridge(TwoLoopBridgeRequest {
            first_loop,
            first_start,
            first_end,
            second_loop,
            second_start,
            second_end,
        })
        .expect("boundary loops can be planned");
    let result = brep.two_loop_bridge(&plan).expect("bridge commits");

    assert_eq!(topology_counts(&brep), (8, 10, 4, 4, 20));
    assert_eq!(
        (
            brep.vertices.len() as i64 - before_counts.0 as i64,
            brep.edges.len() as i64 - before_counts.1 as i64,
            brep.faces.len() as i64 - before_counts.2 as i64,
            brep.loops.len() as i64 - before_counts.3 as i64,
        ),
        (0, 2, 2, 2)
    );
    assert_eq!(
        brep.loops[result.first_bridge_loop].face,
        result.first_bridge_face
    );
    assert_eq!(
        brep.loops[result.second_bridge_loop].face,
        result.second_bridge_face
    );
    assert!(brep.boundary_edges.is_empty());
    brep.validate_half_edge_pairing().unwrap();
    brep.validate_half_edge_cycles().unwrap();
    brep.validate_manifold_edges().unwrap();
}

#[test]
fn two_loop_bridge_rolls_back_for_manifold_interior_operands() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let (first_face, second_face) = disjoint_cube_face_pair(&brep);
    let first_loop = brep.faces[first_face].outer_loop;
    let second_loop = brep.faces[second_face].outer_loop;
    let first_start = brep.loops[first_loop].half_edge;
    let first_end = brep.half_edges[first_start].next;
    let second_start = brep.loops[second_loop].half_edge;
    let second_end = brep.half_edges[second_start].next;
    let before_counts = topology_counts(&brep);

    let plan = brep
        .preflight_two_loop_bridge(TwoLoopBridgeRequest {
            first_loop,
            first_start,
            first_end,
            second_loop,
            second_start,
            second_end,
        })
        .expect("preflight accepts distinct cube face loops");
    let err = brep
        .two_loop_bridge(&plan)
        .expect_err("manifold interior operands must not commit");

    assert!(err.to_string().contains("not an unowned boundary side"));
    assert_eq!(topology_counts(&brep), before_counts);
}

fn disjoint_cube_face_pair(
    brep: &BRep,
) -> (crate::topology::FaceHandle, crate::topology::FaceHandle) {
    let faces = brep.faces.keys().collect::<Vec<_>>();
    for (i, first) in faces.iter().enumerate() {
        for second in faces.iter().skip(i + 1) {
            if !loop_edges_share_any(brep, *first, *second) {
                return (*first, *second);
            }
        }
    }
    panic!("cube fixture does not contain disjoint face pair")
}

fn loop_edges_share_any(
    brep: &BRep,
    first: crate::topology::FaceHandle,
    second: crate::topology::FaceHandle,
) -> bool {
    let first_loop = brep.faces[first].outer_loop;
    let second_loop = brep.faces[second].outer_loop;
    let mut first_edges = std::collections::BTreeSet::new();
    let mut cursor = brep.loops[first_loop].half_edge;
    loop {
        first_edges.insert(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[first_loop].half_edge {
            break;
        }
    }
    let mut cursor = brep.loops[second_loop].half_edge;
    loop {
        if first_edges.contains(&brep.half_edges[cursor].edge) {
            return true;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == brep.loops[second_loop].half_edge {
            break;
        }
    }
    false
}

fn topology_counts(brep: &BRep) -> (usize, usize, usize, usize, usize) {
    (
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
        brep.loops.len(),
        brep.half_edges.len(),
    )
}

fn two_boundary_quads() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    add_boundary_quad(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    let second_face = add_boundary_quad(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
    );
    brep.shells[shell].face = second_face;
    brep
}

fn add_boundary_quad(
    brep: &mut BRep,
    shell: crate::topology::ShellHandle,
    positions: [Point3; 4],
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

    let vertices = positions.map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let half_edges = vertices.map(|vertex| {
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
        half_edge
    });
    for i in 0..4 {
        brep.half_edges[half_edges[i]].next = half_edges[(i + 1) % 4];
        brep.half_edges[half_edges[i]].prev = half_edges[(i + 3) % 4];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}
