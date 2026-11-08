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

/// Counts deltas in field order: V, E, fins, L, F, shells, regions,
/// wire edges, acorns.
type Delta = [isize; 9];

#[test]
fn an_acorn_is_made_and_killed() {
    let mut body = BRep::new(tol());
    let before = NonManifoldCounts::of(&body);
    let acorn = body
        .make_acorn(body.outside, Point3::new(1.0, 2.0, 3.0))
        .expect("acorn");
    let delta: Delta = [1, 0, 0, 0, 0, 1, 0, 0, 1];
    assert_eq!(before.delta(NonManifoldCounts::of(&body)), delta);
    assert_eq!(body.body_kind(), BodyKind::Acorn);
    assert_valid(&body, "acorn");
    body.kill_acorn(acorn.vertex).expect("kill");
    assert_eq!(NonManifoldCounts::of(&body), before);
    assert_eq!(body.body_kind(), BodyKind::Empty);
}

#[test]
fn killing_a_vertex_that_is_no_acorn_declines() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let before = NonManifoldCounts::of(&cube);
    let v = cube.vertices.keys().next().expect("vertex");
    assert_eq!(
        cube.kill_acorn(v),
        Err(EulerError::NotAnAcorn { vertex: v })
    );
    assert_eq!(NonManifoldCounts::of(&cube), before);
    let gone = VertexHandle::default();
    assert_eq!(cube.kill_acorn(gone), Err(EulerError::UnknownVertex(gone)));
}

/// Acorns at the given points of the outside.
fn acorns(body: &mut BRep, points: &[(f64, f64, f64)]) -> Vec<VertexHandle> {
    points
        .iter()
        .map(|&(x, y, z)| {
            body.make_acorn(body.outside, Point3::new(x, y, z))
                .expect("acorn")
                .vertex
        })
        .collect()
}

/// A wire between two acorns absorbs both: one shell holds the edge.
/// Killing it gives both acorns back.
#[test]
fn a_wire_edge_joins_two_acorns_and_its_kill_restores_them() {
    let mut body = BRep::new(tol());
    let v = acorns(&mut body, &[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]);
    let before = NonManifoldCounts::of(&body);
    let wire = body.make_wire_edge(body.outside, v[0], v[1]).expect("wire");
    let made: Delta = [0, 1, 2, 0, 0, -1, 0, 1, -2];
    assert_eq!(before.delta(NonManifoldCounts::of(&body)), made);
    assert_eq!(body.body_kind(), BodyKind::Wire);
    assert_eq!(body.wire_edge_shell(wire.edge), Some(wire.shell));
    assert_valid(&body, "wire");
    let kill = body.kill_wire_edge(wire.edge).expect("kill");
    assert_eq!(NonManifoldCounts::of(&body), before);
    assert_eq!(kill.acorns.len(), 2);
    assert!(kill.split.is_empty());
    assert_eq!(body.body_kind(), BodyKind::Acorn);
    assert_valid(&body, "acorns again");
}

/// A path a-b-c-d cut in the middle splits its shell in two; cutting an
/// end edge leaves its far end an acorn. Three wires at one vertex are
/// a general body.
#[test]
fn killing_a_wire_edge_splits_its_shell_or_leaves_an_acorn() {
    let mut body = BRep::new(tol());
    let v = acorns(
        &mut body,
        &[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (3.0, 0.0, 0.0),
        ],
    );
    let out = body.outside;
    let ab = body.make_wire_edge(out, v[0], v[1]).expect("ab").edge;
    let bc = body.make_wire_edge(out, v[1], v[2]).expect("bc").edge;
    let cd = body.make_wire_edge(out, v[2], v[3]).expect("cd").edge;
    assert_eq!(NonManifoldCounts::of(&body).shells, 1);
    assert_valid(&body, "path");

    let mut cut = body.clone();
    let kill = cut.kill_wire_edge(bc).expect("cut");
    assert_eq!((kill.split.len(), kill.acorns.len()), (1, 0));
    assert_eq!(NonManifoldCounts::of(&cut).shells, 2);
    assert_ne!(cut.wire_edge_shell(ab), cut.wire_edge_shell(cd));
    assert_valid(&cut, "two pieces");

    let mut trimmed = body.clone();
    let kill = trimmed.kill_wire_edge(ab).expect("trim");
    assert_eq!(kill.acorns.len(), 1);
    assert_eq!(kill.acorns[0].0, v[0]);
    assert_eq!(NonManifoldCounts::of(&trimmed).shells, 2);
    assert_valid(&trimmed, "a wire and an acorn");

    let e = acorns(&mut body, &[(1.0, 1.0, 0.0)]);
    body.make_wire_edge(out, v[1], e[0]).expect("branch");
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "branching wire");
}

/// A wire standing on a cube's corner, reaching out into the outside,
/// joins the outside's face shell; killing it returns its far end as an
/// acorn.
#[test]
fn a_wire_stands_on_a_solid_and_comes_off_again() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-12)
        .map(|(h, _)| h)
        .expect("corner");
    let tip = acorns(&mut cube, &[(2.0, 2.0, 2.0)])[0];
    let before = NonManifoldCounts::of(&cube);
    let wire = cube
        .make_wire_edge(cube.outside, corner, tip)
        .expect("wire");
    assert!(cube.shell_has_faces(wire.shell));
    assert_eq!(
        before.delta(NonManifoldCounts::of(&cube)),
        [0, 1, 2, 0, 0, -1, 0, 1, -1]
    );
    assert_eq!(cube.body_kind(), BodyKind::General);
    assert_valid(&cube, "wire on a cube");
    let kill = cube.kill_wire_edge(wire.edge).expect("kill");
    assert_eq!(kill.acorns, vec![(tip, kill.acorns[0].1)]);
    assert_eq!(NonManifoldCounts::of(&cube), before);
    assert_eq!(cube.acorn_shell(corner), None);
    assert_valid(&cube, "cube and acorn");
}

/// Declines leave the body as it was: an end in another region, two
/// face shells joined by a wire, a wire from a vertex to itself, an
/// edge that is not a wire.
#[test]
fn wire_edges_decline_without_changing_the_body() {
    let mut pair = BRep::unit_cube(tol()).expect("cube");
    let far = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("far cube");
    pair.append_body(&far);
    let solid = pair.solid_regions()[0];
    let near_corner = pair
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-12)
        .map(|(h, _)| h)
        .expect("near");
    let far_corner = pair
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(3.0, 0.0, 0.0)) < 1e-12)
        .map(|(h, _)| h)
        .expect("far");
    let tip = acorns(&mut pair, &[(1.5, 5.0, 0.0)])[0];
    let before = NonManifoldCounts::of(&pair);
    let out = pair.outside;
    assert_eq!(
        pair.make_wire_edge(solid, near_corner, tip),
        Err(EulerError::VertexNotInRegion {
            vertex: tip,
            region: solid
        })
    );
    assert!(matches!(
        pair.make_wire_edge(out, near_corner, far_corner),
        Err(EulerError::WireJoinsFaceShells(..))
    ));
    assert_eq!(
        pair.make_wire_edge(out, tip, tip),
        Err(EulerError::DegenerateWireEdge { vertex: tip })
    );
    let face_edge = pair.edges.keys().next().expect("edge");
    assert_eq!(
        pair.kill_wire_edge(face_edge),
        Err(EulerError::NotAWireEdge { edge: face_edge })
    );
    assert_eq!(NonManifoldCounts::of(&pair), before);
}
