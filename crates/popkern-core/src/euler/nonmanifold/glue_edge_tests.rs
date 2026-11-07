use super::*;
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, EdgeHandle, VertexHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

fn vertices_at(body: &BRep, p: Point3) -> Vec<VertexHandle> {
    body.vertices
        .iter()
        .filter(|(_, v)| v.position.distance_to(&p) < 1e-9)
        .map(|(h, _)| h)
        .collect()
}

/// The edges from `a` to `b` (either way).
fn edges_between(body: &BRep, a: VertexHandle, b: VertexHandle) -> Vec<EdgeHandle> {
    body.edges
        .keys()
        .filter(|&e| {
            body.edge_vertices(e)
                .is_some_and(|(x, y)| (x == a && y == b) || (x == b && y == a))
        })
        .collect()
}

/// A unit cube and `other` placed by `place`, in one body, the vertices
/// at (1, 1, 0) and (1, 1, 1) glued. Returns the body, those two
/// vertices and the two edges between them.
fn sharing_the_edge(place: Transform3) -> (BRep, [VertexHandle; 2], Vec<EdgeHandle>) {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    let other = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&place))
        .expect("other");
    body.append_body(&other);
    let ends = [Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 1.0, 1.0)].map(|p| {
        let v = vertices_at(&body, p);
        assert_eq!(v.len(), 2, "both cubes have a corner at {p:?}");
        body.glue_vertices(v[0], v[1]).expect("glue corner");
        v[0]
    });
    let edges = edges_between(&body, ends[0], ends[1]);
    assert_eq!(edges.len(), 2);
    (body, ends, edges)
}

/// Two cubes meeting along an edge (one at the origin, one moved by
/// (1, 1, 0)) glue it into a ring of four, the faces alternating
/// between the cubes; the two outside shells join. Ungluing and then
/// ungluing the corners gives back the two cubes.
#[test]
fn two_cubes_glue_along_an_edge_and_come_apart() {
    let apart = {
        let mut b = BRep::unit_cube(tol()).expect("cube");
        b.append_body(
            &BRep::unit_cube(tol())
                .and_then(|c| c.transform(&Transform3::translation(1.0, 1.0, 0.0)))
                .expect("cube"),
        );
        NonManifoldCounts::of(&b)
    };
    let (mut body, ends, edges) = sharing_the_edge(Transform3::translation(1.0, 1.0, 0.0));
    let before = NonManifoldCounts::of(&body);
    let glued = body.glue_edges(edges[0], edges[1]).expect("glue");
    assert_eq!(glued.ring.len(), 4);
    assert_eq!(glued.merged_shells.len(), 1);
    assert_eq!(
        before.delta(NonManifoldCounts::of(&body)),
        [0, -1, 0, 0, 0, -1, 0, 0, 0]
    );
    // About the edge the cubes' faces come in pairs: two neighbours in
    // the ring share a cube twice round, and differ twice.
    let solid_of = |f| {
        body.fin_face(f)
            .and_then(|face| body.face_back_shell(face))
            .and_then(|s| body.shell_region(s))
    };
    let regions: Vec<_> = glued.ring.iter().map(|&f| solid_of(f)).collect();
    let same = (0..4)
        .filter(|&i| regions[i] == regions[(i + 1) % 4])
        .count();
    assert_eq!(same, 2, "{regions:?}");
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_eq!(body.solid_regions().len(), 2);
    assert_valid(&body, "glued along an edge");

    let other_cube: Vec<_> = glued
        .ring
        .iter()
        .copied()
        .filter(|&f| solid_of(f) != regions[0])
        .collect();
    let apart_edge = body.unglue_edge(edges[0], &other_cube).expect("unglue");
    assert_eq!(body.radial_fins(edges[0]).len(), 2);
    assert_eq!(body.radial_fins(apart_edge).len(), 2);
    assert_eq!(NonManifoldCounts::of(&body), before);
    for v in ends {
        body.unglue_vertex(v).expect("unglue corner");
    }
    assert_eq!(NonManifoldCounts::of(&body), apart);
    // Apart again, the two cubes touch along an edge without sharing it:
    // a body crossing itself there, as before gluing (validation reports
    // it; recorded contacts or gluing are the valid forms).
    assert_eq!(body.body_kind(), BodyKind::Solid);
    assert!(matches!(
        body.validate_all().as_slice(),
        [crate::error::ValidationError::SelfIntersection { .. }]
    ));
}

/// A cube turned so its solid wedge about the edge overlaps the
/// first's: a wedge would lie between the two solids, and the glue
/// declines naming it, the body unchanged.
#[test]
fn a_glue_whose_faces_cross_declines() {
    let turned = Transform3::translation(1.0, 1.0, 0.0)
        .compose(&Transform3::rotation_z(135.0_f64.to_radians()));
    let (mut body, _, edges) = sharing_the_edge(turned);
    let before = NonManifoldCounts::of(&body);
    match body.glue_edges(edges[0], edges[1]) {
        Err(EulerError::RegionConflict { edge, .. }) => assert_eq!(edge, edges[0]),
        other => panic!("{other:?}"),
    }
    assert_eq!(NonManifoldCounts::of(&body), before);
}

/// Edges that do not share their vertices, an edge glued to itself and
/// a run that is the whole ring decline.
#[test]
fn bad_edge_glues_decline() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let mut edges = cube.edges.keys();
    let (a, b) = (edges.next().expect("edge"), edges.next().expect("edge"));
    assert_eq!(
        cube.glue_edges(a, b),
        Err(EulerError::EdgeEndsDiffer { keep: a, gone: b })
    );
    assert_eq!(cube.glue_edges(a, a), Err(EulerError::SameEdge { edge: a }));
    let ring = cube.radial_fins(a);
    assert_eq!(
        cube.unglue_edge(a, &ring),
        Err(EulerError::NotARingRun { edge: a })
    );
}
