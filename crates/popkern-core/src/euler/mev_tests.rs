use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn dissolve_same_loop_bridge_removes_a_spur() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap();
    let n_he_before = brep.half_edges.len();
    let n_edges_before = brep.edges.len();
    let n_vertices_before = brep.vertices.len();

    let v = brep.vertices.keys().next().unwrap();
    let v_out = brep.vertices[v].outgoing;
    let mev_result = brep
        .mev(v_out, Point3::new(0.5, 0.5, 0.5))
        .expect("mev creates spur");
    let new_vertex = mev_result.vertex;
    assert_eq!(brep.half_edges.len(), n_he_before + 2);
    assert_eq!(brep.edges.len(), n_edges_before + 1);
    assert_eq!(brep.vertices.len(), n_vertices_before + 1);

    let spur_out = mev_result.half_edge;
    let spur_in = mev_result.half_edge_twin;
    let l_out = brep.half_edges[spur_out].loop_;
    let l_in = brep.half_edges[spur_in].loop_;
    assert_eq!(l_out, l_in, "spur should be same-loop");

    brep.dissolve_same_loop_bridge(spur_out)
        .expect("dissolve same-loop spur");

    assert_eq!(brep.half_edges.len(), n_he_before);
    assert_eq!(brep.edges.len(), n_edges_before);
    assert_eq!(brep.vertices.len(), n_vertices_before + 1);
    let new_v_outgoing_count = brep
        .half_edges
        .values()
        .filter(|h| h.origin == new_vertex)
        .count();
    assert_eq!(
        new_v_outgoing_count, 0,
        "dissolved spur leaves new_vertex with valence 0"
    );
}

#[test]
fn dissolve_same_loop_bridge_rejects_cross_face() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let primary = brep.edges[edge].half_edge;
    let result = brep.dissolve_same_loop_bridge(primary);
    assert!(
        result.is_err(),
        "dissolve must reject cross-face edges; got {result:?}"
    );
}

#[test]
fn mev_at_vertex_creates_first_edge() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    let mev_r = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert_eq!(brep.vertices.len(), 2);
    assert_eq!(brep.edges.len(), 1);
    assert_eq!(brep.half_edges.len(), 2);

    let he = &brep.half_edges[mev_r.half_edge];
    let twin = &brep.half_edges[mev_r.half_edge_twin];
    assert_eq!(he.twin, mev_r.half_edge_twin);
    assert_eq!(twin.twin, mev_r.half_edge);
    assert_eq!(he.origin, mvfs_r.vertex);
    assert_eq!(twin.origin, mev_r.vertex);
    assert_eq!(he.next, mev_r.half_edge_twin);
    assert_eq!(twin.next, mev_r.half_edge);
    assert_eq!(he.prev, mev_r.half_edge_twin);
    assert_eq!(twin.prev, mev_r.half_edge);
    assert_eq!(brep.loops[mvfs_r.loop_].half_edge, mev_r.half_edge);
    assert_eq!(brep.vertices[mvfs_r.vertex].outgoing, mev_r.half_edge);
    assert_eq!(brep.vertices[mev_r.vertex].outgoing, mev_r.half_edge_twin);
}

#[test]
fn mev_extends_chain() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();

    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 2);
    assert_eq!(brep.half_edges.len(), 4);

    let he1 = &brep.half_edges[mev1.half_edge];
    assert_eq!(he1.next, mev2.half_edge);
}

#[test]
fn mev_three_vertices_in_chain() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let _mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();

    assert_eq!(brep.vertices.len(), 4);
    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.half_edges.len(), 6);
}

#[test]
fn kev_removes_last_vertex() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_r = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    brep.kev(mev_r.half_edge).unwrap();

    assert_eq!(brep.vertices.len(), 1);
    assert_eq!(brep.edges.len(), 0);
    assert_eq!(brep.half_edges.len(), 0);
}

#[test]
fn kev_from_chain_end() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();

    brep.kev(mev2.half_edge).unwrap();

    assert_eq!(brep.vertices.len(), 2);
    assert_eq!(brep.edges.len(), 1);
    assert_eq!(brep.half_edges.len(), 2);
}

#[test]
fn mev_at_vertex_invalid_vertex() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.mev_at_vertex(Default::default(), mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0));
    assert!(result.is_err());
}

#[test]
fn mev_at_vertex_invalid_loop() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.mev_at_vertex(
        mvfs_r.vertex,
        Default::default(),
        Point3::new(1.0, 0.0, 0.0),
    );
    assert!(result.is_err());
}

#[test]
fn mev_at_vertex_loop_already_has_half_edges() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _mev = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let result = brep.mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(2.0, 0.0, 0.0));
    assert!(result.is_err());
}

#[test]
fn mev_invalid_half_edge() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.mev(Default::default(), Point3::new(1.0, 0.0, 0.0));
    assert!(result.is_err());
}

#[test]
fn kev_invalid_half_edge() {
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let result = brep.kev(Default::default());
    assert!(result.is_err());
}

#[test]
fn kev_non_valence_1_vertex() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs_r = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs_r.vertex, mvfs_r.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let _mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let result = brep.kev(mev1.half_edge);
    assert!(result.is_err());
}
