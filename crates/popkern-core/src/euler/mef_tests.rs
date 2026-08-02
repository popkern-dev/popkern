use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn square_wire() -> (
    BRep,
    crate::euler::MvfsResult,
    crate::euler::MevResult,
    crate::euler::MevResult,
    crate::euler::MevResult,
) {
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
    (brep, mvfs, mev1, mev2, mev3)
}

#[test]
fn mef_splits_face_into_two() {
    let (mut brep, _mvfs, mev1, _mev2, mev3) = square_wire();

    let _mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    assert_eq!(brep.faces.len(), 2);
    assert_eq!(brep.loops.len(), 2);
    assert_eq!(brep.edges.len(), 4);
    assert_eq!(brep.half_edges.len(), 8);
}

#[test]
fn kef_merges_two_faces() {
    let (mut brep, _mvfs, mev1, _mev2, mev3) = square_wire();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    brep.kef(mef_r.half_edge).unwrap();

    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.loops.len(), 1);
    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.half_edges.len(), 6);
}

#[test]
fn mef_precondition_same_loop() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mvfs2 = brep.mvfs(Point3::new(5.0, 0.0, 0.0)).unwrap();
    let mev2 = brep
        .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(6.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.mef(mev1.half_edge, mev2.half_edge).is_err());
}

#[test]
fn mef_same_half_edge_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.mef(mev.half_edge, mev.half_edge).is_err());
}

#[test]
fn mef_invalid_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.mef(Default::default(), mev.half_edge).is_err());
    assert!(brep.mef(mev.half_edge, Default::default()).is_err());
}

#[test]
fn kef_same_loop_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.kef(mev.half_edge).is_err());
}

#[test]
fn kef_invalid_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    assert!(brep.kef(Default::default()).is_err());
}

/// PK-232: kef keeps both endpoint vertices alive, so it must
/// re-point any `outgoing` that referenced one of the removed
/// half-edges. Left stale, every downstream vertex-star walk that
/// starts at `outgoing` silently reports nothing (the oblique-apex
/// anchor collection returned false before inspecting a single
/// face and the whole smooth blend fell to Mitred).
#[test]
fn kef_repairs_stale_outgoing_on_both_endpoints() {
    let (mut brep, _mvfs, mev1, _mev2, mev3) = square_wire();
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    // Force both endpoint vertices to use the diagonal's half-edges
    // as their outgoing handles.
    let he = mef_r.half_edge;
    let twin = mef_r.half_edge_twin;
    let he_origin = brep.half_edges[he].origin;
    let twin_origin = brep.half_edges[twin].origin;
    brep.vertices[he_origin].outgoing = he;
    brep.vertices[twin_origin].outgoing = twin;

    brep.kef(he).unwrap();

    for vertex in [he_origin, twin_origin] {
        let outgoing = brep.vertices[vertex].outgoing;
        let live = brep
            .half_edges
            .get(outgoing)
            .expect("outgoing must reference a live half-edge after kef");
        assert_eq!(
            live.origin, vertex,
            "repaired outgoing must originate at its vertex"
        );
    }
}
