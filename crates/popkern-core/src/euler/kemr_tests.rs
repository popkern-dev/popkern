use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, HalfEdgeHandle};

fn bridge_loop_fixture() -> (BRep, crate::euler::MevResult) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    brep.mev(mev1.half_edge_twin, Point3::new(2.0, 0.0, 0.0))
        .unwrap();
    brep.mev(mev1.half_edge, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    (brep, mev1)
}

fn loop_len(brep: &BRep, start: HalfEdgeHandle) -> usize {
    let mut cursor = start;
    let mut count = 0;
    loop {
        count += 1;
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return count;
        }
        assert!(count < 100, "infinite loop detected");
    }
}

#[test]
fn kemr_creates_inner_loop() {
    let (mut brep, mev1) = bridge_loop_fixture();

    let _kemr_r = brep.kemr(mev1.half_edge).unwrap();

    assert_eq!(brep.edges.len(), 2);
    assert_eq!(brep.half_edges.len(), 4);
    assert_eq!(brep.loops.len(), 2);
    assert!(brep.faces.values().any(|f| !f.inner_loops.is_empty()));
    for loop_data in brep.loops.values() {
        assert_eq!(loop_len(&brep, loop_data.half_edge), 2);
    }
}

#[test]
fn mekr_closes_inner_loop() {
    let (mut brep, mev1) = bridge_loop_fixture();
    let kemr_r = brep.kemr(mev1.half_edge).unwrap();

    let inner_he = brep.loops[kemr_r.inner_loop].half_edge;
    let face_with_hole = brep
        .faces
        .values()
        .find(|face| !face.inner_loops.is_empty())
        .unwrap();
    let outer_he = brep.loops[face_with_hole.outer_loop].half_edge;

    let _mekr_r = brep.mekr(outer_he, inner_he).unwrap();

    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.loops.len(), 1);
    assert!(brep.faces.values().all(|face| face.inner_loops.is_empty()));
    assert_eq!(
        loop_len(&brep, brep.loops.values().next().unwrap().half_edge),
        6
    );
}

#[test]
fn kemr_precondition_different_loops() {
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

    assert!(brep.kemr(mef.half_edge).is_err());
}

#[test]
fn kemr_invalid_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    assert!(brep.kemr(Default::default()).is_err());
}

#[test]
fn mekr_invalid_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    assert!(brep.mekr(Default::default(), Default::default()).is_err());
}

#[test]
fn mekr_same_loop_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.mekr(mev1.half_edge, mev1.half_edge_twin).is_err());
}

#[test]
fn mekr_different_face_fails() {
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

    assert!(brep.mekr(mef.half_edge, mef.half_edge_twin).is_err());
}

#[test]
fn kemr_mekr_roundtrip() {
    let (mut brep, mev1) = bridge_loop_fixture();
    let edges_before = brep.edges.len();
    let loops_before = brep.loops.len();
    let kemr_r = brep.kemr(mev1.half_edge).unwrap();
    let inner_he = brep.loops[kemr_r.inner_loop].half_edge;
    let face_with_hole = brep
        .faces
        .values()
        .find(|face| !face.inner_loops.is_empty())
        .unwrap();
    let outer_he = brep.loops[face_with_hole.outer_loop].half_edge;

    brep.mekr(outer_he, inner_he).unwrap();

    assert_eq!(brep.edges.len(), edges_before);
    assert_eq!(brep.loops.len(), loops_before);
}
