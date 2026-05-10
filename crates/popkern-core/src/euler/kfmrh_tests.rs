use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn kfmrh_merges_shells() {
    let mut brep = BRep::new(Tolerance::default());

    let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
    let mev4 = brep
        .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
        .unwrap();
    let mev5 = brep
        .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
        .unwrap();
    let mev6 = brep
        .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
        .unwrap();
    let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

    assert_eq!(brep.shells.len(), 2);
    assert_eq!(brep.solids.len(), 2);

    brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

    let target = &brep.faces[mvfs1.face];
    assert!(
        !target.inner_loops.is_empty(),
        "Target face should have inner loop"
    );
}

#[test]
fn kfmrh_invalid_face_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.kfmrh(Default::default(), mvfs.face);
    assert!(result.is_err());
    let result = brep.kfmrh(mvfs.face, Default::default());
    assert!(result.is_err());
}

#[test]
fn kfmrh_same_shell_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.kfmrh(mvfs1.face, mvfs1.face);
    assert!(result.is_err());
}

#[test]
fn kfmrh_face_with_inner_loops_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
    let mev4 = brep
        .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
        .unwrap();
    let mev5 = brep
        .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
        .unwrap();
    let mev6 = brep
        .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
        .unwrap();
    let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

    brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

    let mvfs3 = brep.mvfs(Point3::new(5.0, 5.0, 0.0)).unwrap();
    let result = brep.kfmrh(mvfs1.face, mvfs3.face);
    assert!(result.is_err());
}

#[test]
fn mfkrh_invalid_loop_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.mfkrh(Default::default());
    assert!(result.is_err());
}

#[test]
fn mfkrh_outer_loop_fails() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.mfkrh(mvfs.loop_);
    assert!(result.is_err());
}

#[test]
fn mfkrh_separates_shells() {
    let mut brep = BRep::new(Tolerance::default());

    let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

    let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
    let mev4 = brep
        .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
        .unwrap();
    let mev5 = brep
        .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
        .unwrap();
    let mev6 = brep
        .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
        .unwrap();
    let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

    brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

    let inner_loop = brep.faces[mvfs1.face].inner_loops[0];
    let result = brep.mfkrh(inner_loop).unwrap();

    assert!(brep.faces.contains_key(result.face));
    assert!(brep.shells.contains_key(result.shell));
    assert!(brep.faces[mvfs1.face].inner_loops.is_empty());
}
