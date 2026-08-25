use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, HalfEdgeHandle, VertexHandle};

fn make_cube() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    brep
}

fn valence(brep: &BRep, v: VertexHandle) -> usize {
    brep.half_edges.values().filter(|he| he.origin == v).count()
}

fn any_outgoing(brep: &BRep, v: VertexHandle) -> HalfEdgeHandle {
    brep.half_edges
        .iter()
        .find(|(_, he)| he.origin == v)
        .map(|(h, _)| h)
        .unwrap()
}

#[test]
fn split_single_edge_off_cube_corner() {
    let mut brep = make_cube();
    let v = brep.vertices.keys().next().unwrap();
    assert_eq!(valence(&brep, v), 3);
    let o = any_outgoing(&brep, v);
    let pos = Point3::new(2.0, 2.0, 2.0);

    let result = brep.split_vertex(o, o, pos).unwrap();

    assert_eq!(result.moved, vec![o]);
    assert_eq!(brep.half_edges[o].origin, result.vertex);
    assert_eq!(brep.vertices[result.vertex].position, pos);
    // Old vertex: 3 - 1 moved + connector = 3; new vertex: 1 moved + connector = 2.
    assert_eq!(valence(&brep, v), 3);
    assert_eq!(valence(&brep, result.vertex), 2);
    // Connector direction contract: half_edge runs old -> new.
    assert_eq!(brep.half_edges[result.half_edge].origin, v);
    assert_eq!(brep.half_edges[result.half_edge_twin].origin, result.vertex);
    brep.validate().unwrap();
}

#[test]
fn split_face_wedge_sector_rehomes_face_corner() {
    let mut brep = make_cube();
    let v = brep.vertices.keys().next().unwrap();
    let o = any_outgoing(&brep, v);
    // Wedge (x, o) in o's loop: x ends at v. Sector [twin(x), o] moves
    // both wedge edges, so the wedge face's corner re-homes entirely.
    let x = brep.half_edges[o].prev;
    let wedge_loop = brep.half_edges[o].loop_;
    let he_first = brep.half_edges[x].twin;
    let pos = Point3::new(-1.0, -1.0, -1.0);

    let result = brep.split_vertex(he_first, o, pos).unwrap();

    assert_eq!(result.moved, vec![he_first, o]);
    // The wedge face's loop passes through the new vertex, not v.
    let mut cursor = brep.loops[wedge_loop].half_edge;
    let start = cursor;
    loop {
        assert_ne!(brep.half_edges[cursor].origin, v);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    // Old vertex keeps 1 original edge + connector; new vertex has the
    // 2 moved edges + connector.
    assert_eq!(valence(&brep, v), 2);
    assert_eq!(valence(&brep, result.vertex), 3);
    // Connector halves live in the two cut-wedge loops, which flank the
    // moved sector (neither is the wedge loop itself).
    assert_ne!(brep.half_edges[result.half_edge].loop_, wedge_loop);
    assert_ne!(brep.half_edges[result.half_edge_twin].loop_, wedge_loop);
    brep.validate().unwrap();
}

#[test]
fn split_vertex_rejects_full_fan_sector() {
    let mut brep = make_cube();
    let v = brep.vertices.keys().next().unwrap();
    let o = any_outgoing(&brep, v);
    // Walk the fan to the last outgoing before wrapping.
    let mut last = o;
    loop {
        let succ = brep.half_edges[brep.half_edges[last].twin].next;
        if succ == o {
            break;
        }
        last = succ;
    }
    let err = brep.split_vertex(o, last, Point3::new(0.0, 0.0, 5.0));
    assert!(err.is_err());
}

#[test]
fn split_vertex_rejects_mismatched_origins() {
    let mut brep = make_cube();
    let v = brep.vertices.keys().next().unwrap();
    let o = any_outgoing(&brep, v);
    let other = brep
        .half_edges
        .iter()
        .find(|(_, he)| he.origin != v)
        .map(|(h, _)| h)
        .unwrap();
    assert!(
        brep.split_vertex(o, other, Point3::new(0.0, 0.0, 5.0))
            .is_err()
    );
}

#[test]
fn split_vertex_reanchors_old_vertex_when_anchor_moves() {
    let mut brep = make_cube();
    let v = brep.vertices.keys().next().unwrap();
    // Force the anchor onto the half-edge we are about to move.
    let o = any_outgoing(&brep, v);
    brep.vertices[v].outgoing = o;
    let result = brep.split_vertex(o, o, Point3::new(3.0, 0.0, 0.0)).unwrap();
    let anchor = brep.vertices[v].outgoing;
    assert_eq!(brep.half_edges[anchor].origin, v);
    assert_eq!(anchor, result.half_edge);
    brep.validate().unwrap();
}
