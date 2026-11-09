use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

/// Every vertex of a cube is a disc of three faces, every edge two fins.
#[test]
fn a_cube_is_manifold_everywhere() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    for v in cube.vertices.keys() {
        assert_eq!(
            cube.vertex_neighbourhood(v),
            VertexNeighbourhood::Disc { pinch_uses: vec![] }
        );
    }
    for e in cube.edges.keys() {
        assert_eq!(cube.edge_neighbourhood(e), EdgeNeighbourhood::Manifold);
    }
}

/// A sphere's poles close one fan round the seam: the one face meets
/// the pole once, across the seam to itself.
#[test]
fn a_spheres_poles_are_discs() {
    let ball = BRep::sphere(2.0, Tolerance::default()).expect("sphere");
    for v in ball.vertices.keys() {
        let n = ball.vertex_neighbourhood(v);
        assert!(n.is_manifold(), "{n:?}");
    }
}

/// An open box (a cube without one face) is a half disc at each of the
/// open face's four corners, its four rim edges free, and a disc at the
/// other four corners. (A constructor's two-face sheet is a disc
/// everywhere: its front and back faces meet at its rim.)
#[test]
fn an_open_box_is_a_half_disc_at_its_rim() {
    let open = BRep::open_box(Tolerance::default()).expect("open box");
    let (mut half, mut whole) = (0, 0);
    for v in open.vertices.keys() {
        match open.vertex_neighbourhood(v) {
            VertexNeighbourhood::HalfDisc { pinch_uses } if pinch_uses.is_empty() => half += 1,
            VertexNeighbourhood::Disc { pinch_uses } if pinch_uses.is_empty() => whole += 1,
            n => panic!("{n:?}"),
        }
    }
    assert_eq!((half, whole), (4, 4));
    let free = open
        .edges
        .keys()
        .filter(|&e| open.edge_neighbourhood(e) == EdgeNeighbourhood::Free)
        .count();
    assert_eq!(free, 4);
}

/// Two cubes sharing one corner, glued into one general body, meet in
/// two fans there.
#[test]
fn cubes_glued_at_a_corner_have_two_fans_there() {
    use crate::geometry::Transform3;
    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).expect("cube");
    let b = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::translation(1.0, 1.0, 1.0))
        .expect("moves");
    let mut body = a;
    body.append_body(&b);
    let glued = body.glue_coincident();
    assert!(glued.is_ok(), "{glued:?}");
    let corner = Point3::new(1.0, 1.0, 1.0);
    let v = body
        .vertices
        .iter()
        .find(|(_, x)| x.position.distance_to(&corner) < 1e-9)
        .map(|(h, _)| h)
        .expect("the shared corner");
    assert_eq!(
        body.vertex_neighbourhood(v),
        VertexNeighbourhood::Fans { count: 2 }
    );
}

/// The one-pass query answers as asking vertex by vertex, on a cube, a
/// sphere, an open box, and two cubes glued at a corner (whose glue
/// removed a vertex, leaving a gap in the vertex handles).
#[test]
fn every_vertex_at_once_answers_as_one_at_a_time() {
    use crate::geometry::Transform3;
    let tol = Tolerance::default();
    let mut glued = BRep::unit_cube(tol).expect("cube");
    glued.append_body(
        &BRep::unit_cube(tol)
            .expect("cube")
            .transform(&Transform3::translation(1.0, 1.0, 1.0))
            .expect("moves"),
    );
    assert!(glued.glue_coincident().is_ok());
    for body in [
        BRep::unit_cube(tol).expect("cube"),
        BRep::sphere(2.0, tol).expect("sphere"),
        BRep::open_box(tol).expect("open box"),
        glued,
    ] {
        let all = body.vertex_neighbourhoods();
        assert_eq!(all.len(), body.vertices.len());
        for (v, n) in all {
            assert_eq!(n, body.vertex_neighbourhood(v));
        }
    }
}

/// A cylinder's lateral face meets its rim vertices on both sides of its
/// seam: one run round the fan, crossing the seam into itself, not a
/// pinch use (read as one, every fillet and shell of a cylinder
/// declined).
#[test]
fn a_seam_crossing_is_no_pinch_use() {
    let rod = BRep::cylinder(2.0, 5.0, Tolerance::default()).expect("cylinder");
    for (v, n) in rod.vertex_neighbourhoods() {
        assert_eq!(n, VertexNeighbourhood::Disc { pinch_uses: vec![] }, "{v:?}");
    }
    assert!(rod.pinch_uses().is_empty());
}
