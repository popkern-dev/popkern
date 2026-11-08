use super::*;
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, VertexHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

fn at(body: &BRep, p: Point3) -> Vec<VertexHandle> {
    body.vertices
        .iter()
        .filter(|(_, v)| v.position.distance_to(&p) < 1e-12)
        .map(|(h, _)| h)
        .collect()
}

/// Two cubes meeting at a corner glue into one vertex with two fans (a
/// general body); ungluing gives the second fan its own vertex again.
#[test]
fn two_cubes_glue_at_a_corner_and_come_apart() {
    let mut pair = BRep::unit_cube(tol()).expect("cube");
    let other = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(1.0, 1.0, 1.0)))
        .expect("cube");
    pair.append_body(&other);
    let before = NonManifoldCounts::of(&pair);
    let corner = at(&pair, Point3::new(1.0, 1.0, 1.0));
    let glued = pair.glue_vertices(corner[0], corner[1]).expect("glue");
    assert_eq!(glued.distance, 0.0);
    assert!(glued.removed_shells.is_empty());
    assert_eq!(
        before.delta(NonManifoldCounts::of(&pair)),
        [-1, 0, 0, 0, 0, 0, 0, 0, 0]
    );
    assert!(pair.is_nonmanifold_vertex(corner[0]));
    assert_eq!(pair.body_kind(), BodyKind::General);
    assert_valid(&pair, "glued at a corner");
    let made = pair.unglue_vertex(corner[0]).expect("unglue");
    assert_eq!(made.len(), 1);
    assert_eq!(NonManifoldCounts::of(&pair), before);
    assert_eq!(pair.body_kind(), BodyKind::Solid);
    assert_valid(&pair, "apart again");
}

/// Two wires meeting end to end glue into one path in one shell, and
/// come apart into two shells again.
#[test]
fn wire_ends_glue_into_one_shell() {
    let mut body = BRep::new(tol());
    let out = body.outside;
    let p = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 0.0)];
    let v: Vec<_> = p
        .iter()
        .map(|&(x, y)| {
            body.make_acorn(out, Point3::new(x, y, 0.0))
                .expect("acorn")
                .vertex
        })
        .collect();
    body.make_wire_edge(out, v[0], v[1]).expect("first");
    body.make_wire_edge(out, v[2], v[3]).expect("second");
    let before = NonManifoldCounts::of(&body);
    assert_eq!(before.shells, 2);
    let glued = body.glue_vertices(v[1], v[2]).expect("glue");
    assert_eq!(glued.removed_shells.len(), 1);
    assert_eq!(
        before.delta(NonManifoldCounts::of(&body)),
        [-1, 0, 0, 0, 0, -1, 0, 0, 0]
    );
    assert_eq!(body.body_kind(), BodyKind::Wire);
    assert_valid(&body, "one path");
    assert_eq!(body.unglue_vertex(v[1]).expect("unglue").len(), 1);
    assert_eq!(NonManifoldCounts::of(&body), before);
    assert_valid(&body, "two wires");
}

/// An acorn glued onto a cube's corner is absorbed; an acorn of another
/// region, a far vertex, the same vertex and an unknown one decline,
/// leaving the body unchanged.
#[test]
fn acorns_absorb_and_bad_glues_decline() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let corner = at(&cube, Point3::new(1.0, 1.0, 1.0))[0];
    let solid = cube.solid_regions()[0];
    let inside = cube
        .make_acorn(solid, Point3::new(1.0, 1.0, 1.0))
        .expect("acorn in the solid")
        .vertex;
    let far = cube
        .make_acorn(cube.outside, Point3::new(5.0, 5.0, 5.0))
        .expect("far acorn")
        .vertex;
    let before = NonManifoldCounts::of(&cube);
    assert!(matches!(
        cube.glue_vertices(corner, far),
        Err(EulerError::NotCoincident { .. })
    ));
    assert_eq!(
        cube.glue_vertices(corner, corner),
        Err(EulerError::SameVertex { vertex: corner })
    );
    assert_eq!(
        cube.glue_vertices(corner, VertexHandle::default()),
        Err(EulerError::UnknownVertex(VertexHandle::default()))
    );
    assert_eq!(NonManifoldCounts::of(&cube), before);
    // The corner touches the solid, so an acorn of the solid there is
    // absorbed too.
    let glued = cube.glue_vertices(corner, inside).expect("absorbed");
    assert_eq!(glued.removed_shells.len(), 1);
    assert_eq!(
        before.delta(NonManifoldCounts::of(&cube)),
        [-1, 0, 0, 0, 0, -1, 0, 0, -1]
    );
    // Two acorns of different regions do not meet: the outside's acorn
    // finds the solid's acorn outside its region.
    let lone = cube
        .make_acorn(solid, Point3::new(5.0, 5.0, 5.0))
        .expect("lone");
    assert_eq!(
        cube.glue_vertices(far, lone.vertex),
        Err(EulerError::VertexNotInRegion {
            vertex: lone.vertex,
            region: cube.outside
        })
    );
}
