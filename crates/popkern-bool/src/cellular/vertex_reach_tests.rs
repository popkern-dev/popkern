use super::*;
use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;

/// A unit cube with one corner moved `d` off its edges' ends, and the
/// corner's handle and new place.
fn cube_with_corner_off(d: f64) -> (BRep, VertexHandle, Point3) {
    let tol = Tolerance::default();
    let mut cube = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        tol,
    )
    .expect("cube");
    let v = cube.vertices.keys().next().expect("vertex");
    let p = cube.vertices[v].position + Vec3::new(d, 0.0, 0.0);
    cube.vertices[v].position = p;
    (cube, v, p)
}

/// A corner 1.5 tolerances off its edges' ends, at an anchor whose
/// group grew to three, takes that distance as its tolerance and its
/// edges at least it, and the body validates (it reads open before).
/// Past the group's tolerance, or at no anchor, nothing changes.
#[test]
fn a_vertex_at_an_anchor_records_its_edges_ends() {
    let (mut cube, v, p) = cube_with_corner_off(1.5e-6);
    assert!(cube.validate().is_err());
    record_vertex_reaches(&mut cube, &[(p, 3e-6)]);
    let t = cube.vertex_tolerance(v);
    assert!((1.5e-6..1.6e-6).contains(&t), "vertex tolerance {t}");
    for e in cube.edges.keys() {
        let (a, b) = cube.edge_vertices(e).expect("ends");
        if a == v || b == v {
            assert!(cube.edge_tolerance(e) >= t);
        }
    }
    assert!(cube.validate().is_ok(), "{:?}", cube.validate());
    for reaches in [vec![(p, 1.2e-6)], vec![(Point3::new(9.0, 9.0, 9.0), 3e-6)]] {
        let (mut cube, v, _) = cube_with_corner_off(1.5e-6);
        record_vertex_reaches(&mut cube, &reaches);
        assert!(cube.vertex_tolerance(v) <= 1e-6);
    }
}
