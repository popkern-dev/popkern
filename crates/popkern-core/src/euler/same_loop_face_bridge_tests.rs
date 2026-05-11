use crate::euler::SameLoopFaceBridgeRequest;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex};

#[test]
fn same_loop_face_bridge_commits_disjoint_boundary_intervals() {
    let mut brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let start = brep.loops[loop_].half_edge;
    let first_after = brep.half_edges[start].next;
    let second_start = nth_half_edge(&brep, start, 3);
    let second_after = brep.half_edges[second_start].next;
    let before = topology_counts(&brep);

    let plan = brep
        .preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_,
            first_interval: vec![start, first_after],
            second_interval: vec![second_start, second_after],
        })
        .expect("same-loop bridge preflight accepts disjoint boundary arcs");
    let result = brep
        .same_loop_face_bridge(&plan)
        .expect("same-loop bridge commits");

    assert_eq!(
        (
            brep.vertices.len() as i64 - before.0 as i64,
            brep.edges.len() as i64 - before.1 as i64,
            brep.faces.len() as i64 - before.2 as i64,
            brep.loops.len() as i64 - before.3 as i64,
        ),
        (0, 2, 1, 1)
    );
    assert!(brep.faces.contains_key(result.bridge_face));
    assert!(brep.loops.contains_key(result.bridge_loop));
    assert!(brep.boundary_edges.contains(&result.first_connector_edge));
    assert!(brep.boundary_edges.contains(&result.second_connector_edge));
    brep.validate_half_edge_pairing()
        .expect("pairing validates");
    brep.validate_half_edge_cycles().expect("cycles validate");
    brep.validate_manifold_edges().expect("manifold validates");
}

#[test]
fn same_loop_face_bridge_refuses_overlapping_interval_edges() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let start = brep.loops[loop_].half_edge;
    let next = brep.half_edges[start].next;
    let second_after = brep.half_edges[next].next;

    let err = brep
        .preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_,
            first_interval: vec![start, next],
            second_interval: vec![start, next, second_after],
        })
        .expect_err("overlapping intervals are non-manifold");

    assert!(err.to_string().contains("overlap"));
}

#[test]
fn same_loop_face_bridge_refuses_degenerate_connector() {
    let brep = same_loop_boundary_fixture();
    let face = brep.faces.keys().next().expect("fixture has a face");
    let loop_ = brep.faces[face].outer_loop;
    let start = brep.loops[loop_].half_edge;
    let next = brep.half_edges[start].next;
    let after_next = brep.half_edges[next].next;

    let err = brep
        .preflight_same_loop_face_bridge(SameLoopFaceBridgeRequest {
            loop_,
            first_interval: vec![start, next],
            second_interval: vec![next, after_next],
        })
        .expect_err("adjacent intervals create a zero-length connector");

    assert!(err.to_string().contains("degenerate"));
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

fn topology_counts(brep: &BRep) -> (usize, usize, usize, usize) {
    (
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
        brep.loops.len(),
    )
}

fn nth_half_edge(brep: &BRep, start: HalfEdgeHandle, offset: usize) -> HalfEdgeHandle {
    let mut cursor = start;
    for _ in 0..offset {
        cursor = brep.half_edges[cursor].next;
    }
    cursor
}
