use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Vertex};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn empty_brep() -> BRep {
    BRep::new(tol())
}

// ---------- weld_coincident_vertices ----------

#[test]
fn weld_coincident_vertices_no_op_on_empty_brep() {
    let mut brep = empty_brep();
    let n = brep.weld_coincident_vertices(&tol());
    assert_eq!(n, 0);
}

#[test]
fn weld_coincident_vertices_no_op_on_single_vertex() {
    let mut brep = empty_brep();
    brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let n = brep.weld_coincident_vertices(&tol());
    assert_eq!(n, 0);
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn weld_coincident_vertices_collapses_two_at_same_position() {
    let mut brep = empty_brep();
    let _v1 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let _v2 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let n = brep.weld_coincident_vertices(&tol());
    assert_eq!(n, 1);
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn weld_coincident_vertices_leaves_distinct_vertices_alone() {
    let mut brep = empty_brep();
    brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let n = brep.weld_coincident_vertices(&tol());
    assert_eq!(n, 0);
    assert_eq!(brep.vertices.len(), 2);
}

#[test]
fn weld_coincident_vertices_unit_cube_unaffected() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let pre = cube.vertices.len();
    let n = cube.weld_coincident_vertices(&tol());
    assert_eq!(n, 0);
    assert_eq!(cube.vertices.len(), pre);
}

#[test]
fn weld_coincident_vertices_three_at_same_position() {
    let mut brep = empty_brep();
    for _ in 0..3 {
        brep.vertices.insert(Vertex {
            position: Point3::new(0.5, 0.5, 0.5),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
    }
    let n = brep.weld_coincident_vertices(&tol());
    assert_eq!(n, 2);
    assert_eq!(brep.vertices.len(), 1);
}

// ---------- weld_vertex_pair ----------

#[test]
fn weld_vertex_pair_collapses_two_into_one() {
    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.weld_vertex_pair(v1, v2).unwrap();
    assert_eq!(brep.vertices.len(), 1);
    assert!(brep.vertices.contains_key(v1));
    assert!(!brep.vertices.contains_key(v2));
}

#[test]
fn weld_vertex_pair_preserves_primary_position() {
    // The primary's position is left unchanged even when the
    // secondary had a different position. This is by design —
    // the caller is responsible for moving both to the same
    // location before welding if exact merge is desired.
    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(7.0, 8.0, 9.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.weld_vertex_pair(v1, v2).unwrap();
    let pos = brep.vertices[v1].position;
    assert_eq!(pos.x, 1.0);
    assert_eq!(pos.y, 2.0);
    assert_eq!(pos.z, 3.0);
}

#[test]
fn weld_vertex_pair_rejects_self_weld() {
    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let err = brep.weld_vertex_pair(v1, v1).unwrap_err();
    assert!(matches!(err, CoreError::PreconditionFailed(_)));
}

#[test]
fn weld_vertex_pair_rejects_unknown_primary() {
    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let bogus = VertexHandle::default();
    let err = brep.weld_vertex_pair(bogus, v1).unwrap_err();
    assert!(matches!(err, CoreError::InvalidHandle(_)));
}

#[test]
fn weld_vertex_pair_rejects_unknown_secondary() {
    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let bogus = VertexHandle::default();
    let err = brep.weld_vertex_pair(v1, bogus).unwrap_err();
    assert!(matches!(err, CoreError::InvalidHandle(_)));
}

#[test]
fn weld_vertex_pair_updates_half_edge_origins() {
    use crate::topology::HalfEdge;

    let mut brep = empty_brep();
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    // Create a half-edge originating at v2 (no full topology, just
    // enough to verify the rewrite).
    let he = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: Default::default(),
        loop_: Default::default(),
        pcurve: None,
    });

    brep.weld_vertex_pair(v1, v2).unwrap();
    assert_eq!(brep.half_edges[he].origin, v1);
}
