use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn mvfs_creates_vertex_face_shell_solid() {
    let mut brep = BRep::new(Tolerance::default());
    let result = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    assert_eq!(brep.vertices.len(), 1);
    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.loops.len(), 1);
    assert_eq!(brep.shells.len(), 1);
    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.edges.len(), 0);
    assert_eq!(brep.half_edges.len(), 0);

    let face = &brep.faces[result.face];
    assert_eq!(face.outer_loop, result.loop_);
    assert!(face.inner_loops.is_empty());
    assert_eq!(brep.loops[result.loop_].face, result.face);
    assert_eq!(brep.shells[result.shell].solid, result.solid);
    assert_eq!(brep.solids[result.solid].shell, result.shell);
    assert_eq!(brep.vertices[result.vertex].position.x, 0.0);
}

#[test]
fn kvfs_removes_everything() {
    let mut brep = BRep::new(Tolerance::default());
    let result = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.kvfs(result.vertex).unwrap();

    assert_eq!(brep.vertices.len(), 0);
    assert_eq!(brep.faces.len(), 0);
    assert_eq!(brep.loops.len(), 0);
    assert_eq!(brep.shells.len(), 0);
    assert_eq!(brep.solids.len(), 0);
}

#[test]
fn kvfs_invalid_handle_fails() {
    let mut brep = BRep::new(Tolerance::default());
    assert!(brep.kvfs(Default::default()).is_err());
}

#[test]
fn two_mvfs_creates_two_solids() {
    let mut brep = BRep::new(Tolerance::default());
    let r1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let r2 = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();

    assert_eq!(brep.vertices.len(), 2);
    assert_eq!(brep.solids.len(), 2);
    assert_ne!(r1.solid, r2.solid);
}

#[test]
fn kvfs_fails_with_outgoing_edges() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.kvfs(mvfs.vertex).is_err());
}

#[test]
fn kvfs_seed_removes_one_of_two_seed_solids() {
    let mut brep = BRep::new(Tolerance::default());
    let r1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _r2 = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();

    brep.kvfs_seed(r1.into()).unwrap();

    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn kvfs_vertex_wrapper_rejects_ambiguous_seed_solids() {
    let mut brep = BRep::new(Tolerance::default());
    let r1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _r2 = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();

    let err = brep.kvfs(r1.vertex).expect_err("ambiguous seed rejected");

    assert!(err.to_string().contains("exactly one mvfs seed"));
    assert_eq!(brep.solids.len(), 2);
    assert_eq!(brep.vertices.len(), 2);
}

#[test]
fn kvfs_roundtrip() {
    let mut brep = BRep::new(Tolerance::default());
    let r = brep.mvfs(Point3::new(5.0, 5.0, 5.0)).unwrap();
    brep.kvfs(r.vertex).unwrap();

    assert_eq!(brep.vertices.len(), 0);
    assert_eq!(brep.faces.len(), 0);
    assert_eq!(brep.loops.len(), 0);
    assert_eq!(brep.shells.len(), 0);
    assert_eq!(brep.solids.len(), 0);
}

#[test]
fn kvfs_allows_seed_after_kev_returns_loop_to_empty() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    brep.kev(mev.half_edge).unwrap();

    brep.kvfs(mvfs.vertex).unwrap();

    assert_eq!(brep.vertices.len(), 0);
}

#[test]
fn kvfs_seed_with_face_with_inner_loops_rejects() {
    let mut brep = BRep::new(Tolerance::default());
    let r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.faces[r.face].inner_loops.push(Default::default());

    let err = brep
        .kvfs_seed(r.into())
        .expect_err("inner loop corruption rejected");

    assert!(err.to_string().contains("isolated mvfs face"));
    assert_eq!(brep.solids.len(), 1);
}
