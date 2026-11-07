use super::super::NonManifoldCounts;
use super::*;
use crate::geometry::Vec3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, VertexUse};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

/// The unit cube's top face (z = 1).
fn top_face(body: &BRep) -> FaceHandle {
    body.faces
        .keys()
        .find(|&f| {
            body.face_point_trim_containment(f, &Point3::new(0.5, 0.5, 1.0))
                .is_some_and(|p| p.distance < 1e-9)
        })
        .expect("top face")
}

/// A vertex loop on a cube's top face is one more vertex and loop, the
/// body still a solid that validates; killing it gives the cube back.
#[test]
fn a_vertex_loop_stands_on_a_face_and_goes_again() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let before = (cube.vertices.len(), cube.loops.len());
    let counts = NonManifoldCounts::of(&cube);
    let top = top_face(&cube);
    let made = cube
        .make_vertex_loop(top, Point3::new(0.5, 0.5, 1.0))
        .expect("vertex loop");
    assert_eq!(
        (cube.vertices.len(), cube.loops.len()),
        (before.0 + 1, before.1 + 1)
    );
    assert_eq!(
        cube.vertex_uses(made.vertex),
        vec![VertexUse::Loop(made.loop_)]
    );
    assert_eq!(cube.body_kind(), BodyKind::Solid);
    assert_valid(&cube, "vertex loop on a cube");
    cube.kill_vertex_loop(made.loop_).expect("kill");
    assert_eq!((cube.vertices.len(), cube.loops.len()), before);
    assert_eq!(NonManifoldCounts::of(&cube), counts);
    assert_valid(&cube, "killed again");
}

/// Points off the surface, on an edge, outside the face and on another
/// vertex loop decline, the body unchanged; a vertex loop moved off its
/// face is reported (negative control).
#[test]
fn vertex_loops_decline_off_the_face_interior() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let top = top_face(&cube);
    cube.make_vertex_loop(top, Point3::new(0.25, 0.25, 1.0))
        .expect("first");
    let before = (cube.vertices.len(), cube.loops.len());
    for p in [
        Point3::new(0.5, 0.5, 1.1),
        Point3::new(0.5, 0.0, 1.0),
        Point3::new(2.0, 2.0, 1.0),
        Point3::new(0.25, 0.25, 1.0),
    ] {
        assert!(
            matches!(
                cube.make_vertex_loop(top, p),
                Err(EulerError::NotOnFaceInterior { .. })
            ),
            "{p:?}"
        );
        assert_eq!((cube.vertices.len(), cube.loops.len()), before);
    }
    let lone = cube.vertex_loops();
    assert_eq!(lone.len(), 1);
    cube.vertices[lone[0].1].position = Point3::new(0.25, 0.25, 1.2);
    assert!(matches!(
        cube.validate_vertex_loops_all().as_slice(),
        [crate::error::ValidationError::VertexLoop { .. }]
    ));
}

/// A cube with a second standing on its corner, the corner on the
/// first's top face centre, the face's vertex loop glued to the corner.
/// Returns the body and the shared vertex.
fn corner_on_a_face() -> (BRep, VertexHandle) {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    // The diagonal (1, 1, 1) turned to +z, the corner at the origin
    // then moved to the top face's centre: every other point lies above.
    let diagonal = Vec3::new(1.0, 1.0, 1.0);
    let axis = diagonal.cross(&Vec3::new(0.0, 0.0, 1.0));
    let angle = (1.0 / 3.0_f64.sqrt()).acos();
    let turn =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, angle).expect("turn");
    let standing = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(0.5, 0.5, 1.0).compose(&turn)))
        .expect("standing cube");
    let map = body.append_body(&standing);
    let corner = map
        .vertices
        .values()
        .copied()
        .find(|&v| {
            body.vertices[v]
                .position
                .distance_to(&Point3::new(0.5, 0.5, 1.0))
                < 1e-9
        })
        .expect("corner");
    let top = top_face(&body);
    let made = body
        .make_vertex_loop(top, Point3::new(0.5, 0.5, 1.0))
        .expect("vertex loop");
    body.glue_vertices(made.vertex, corner).expect("glue");
    (body, made.vertex)
}

/// A second cube standing on its corner, the corner on the first's top
/// face: the face's vertex loop glued to the corner makes one general
/// body (the vertex a fan and a loop), which validates and comes apart
/// again.
#[test]
fn a_corner_glued_to_a_vertex_loop_makes_a_general_body() {
    let (mut body, shared) = corner_on_a_face();
    let uses = body.vertex_uses(shared);
    assert!(
        matches!(uses.as_slice(), [VertexUse::Fan(_), VertexUse::Loop(_)]),
        "{uses:?}"
    );
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "corner on a face");
    let apart = body.unglue_vertex(shared).expect("unglue");
    assert_eq!(apart.len(), 1);
    assert_eq!(body.vertex_uses(shared).len(), 1);
}

/// Split into manifold pieces, the vertex loop goes and the contact is
/// recorded by the two faces, as the split policy records a vertex on a
/// face; gluing the contacts makes the vertex loop again, nothing kept.
#[test]
fn a_vertex_loop_splits_and_glues_again() {
    let (mut body, _) = corner_on_a_face();
    let loops = body.loops.len();
    let contacts = body.split_to_manifold().expect("split");
    assert_eq!(contacts.len(), 1);
    assert!(body.vertex_loops().is_empty());
    assert_eq!(body.loops.len(), loops - 1);
    assert_eq!(body.body_kind(), BodyKind::Solid);
    assert_valid(&body, "pieces");
    let glued = body.glue_contacts().expect("glue");
    assert_eq!((glued.glued, glued.kept.len()), (1, 0));
    assert_eq!(body.vertex_loops().len(), 1);
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "glued again");
}

/// A stored body keeps its vertex loops: the loop names its vertex
/// after a save and load, and a loop of fins writes no vertex field.
#[test]
fn vertex_loops_survive_a_save_and_load() {
    let (body, shared) = corner_on_a_face();
    let text = serde_json::to_string(&body).expect("save");
    let again: BRep = serde_json::from_str(&text).expect("load");
    assert_eq!(again.vertex_loops(), body.vertex_loops());
    assert_eq!(again.vertex_uses(shared), body.vertex_uses(shared));
    assert_valid(&again, "loaded");
    let plain = serde_json::to_string(&BRep::unit_cube(tol()).expect("cube")).expect("save");
    assert!(
        !plain.contains("\"vertex\":"),
        "loops of fins write no vertex"
    );
}
