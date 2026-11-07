use super::*;
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, FaceHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

/// A unit cube with another standing on it (moved up by one), the four
/// corners where they meet glued. Returns the body and the two faces at
/// z = 1.
fn box_on_box() -> (BRep, [FaceHandle; 2]) {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    body.append_body(
        &BRep::unit_cube(tol())
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 1.0)))
            .expect("cube"),
    );
    for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let p = Point3::new(x, y, 1.0);
        let v: Vec<_> = body
            .vertices
            .iter()
            .filter(|(_, v)| v.position.distance_to(&p) < 1e-9)
            .map(|(h, _)| h)
            .collect();
        body.glue_vertices(v[0], v[1]).expect("corner");
    }
    let at_one: Vec<FaceHandle> = body
        .faces
        .keys()
        .filter(|&f| {
            body.loop_cycle(body.faces[f].outer_loop)
                .iter()
                .all(|&h| (body.vertices[body.half_edges[h].origin].position.z - 1.0).abs() < 1e-9)
        })
        .collect();
    assert_eq!(at_one.len(), 2);
    (body, [at_one[0], at_one[1]])
}

/// The touching faces of a box standing on a box glue into one face
/// between the two solids: its edges merge, the outside's two shells
/// join, and each solid's shell holds one side of the face.
#[test]
fn a_box_on_a_box_glues_into_one_internal_face() {
    let (mut body, [keep, gone]) = box_on_box();
    let before = NonManifoldCounts::of(&body);
    let glued = body.glue_faces(keep, gone).expect("glue");
    assert_eq!(glued.merged_edges.len(), 4);
    assert_eq!(
        before.delta(NonManifoldCounts::of(&body)),
        [0, -4, -4, -1, -1, -1, 0, 0, 0]
    );
    let sides = body.face_sides(keep).expect("sides");
    let solids = body.solid_regions();
    assert_eq!(solids.len(), 2);
    let regions = [sides.back, sides.front].map(|s| body.shell_region(s).expect("region"));
    assert!(regions.iter().all(|r| solids.contains(r)));
    assert_ne!(regions[0], regions[1]);
    assert!(body.edges.keys().all(|e| body.radial_fins(e).len() <= 3));
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "internal face");
}

/// Ungluing the internal face gives it a reversed twin with a void of
/// no thickness between; gluing the two again restores the face.
#[test]
fn an_internal_face_unglues_and_glues_again() {
    let (mut body, [keep, gone]) = box_on_box();
    body.glue_faces(keep, gone).expect("glue");
    let internal = NonManifoldCounts::of(&body);
    let twin = body.unglue_face(keep).expect("unglue");
    assert_eq!(
        internal.delta(NonManifoldCounts::of(&body)),
        [0, 0, 4, 1, 1, 1, 1, 0, 0]
    );
    let gap = body.face_front_shell(keep).expect("front");
    assert_eq!(body.face_front_shell(twin), Some(gap));
    body.glue_faces(keep, twin).expect("glue again");
    assert_eq!(NonManifoldCounts::of(&body), internal);
    assert_valid(&body, "glued again");
}

/// Faces that do not share their boundary decline, the body unchanged;
/// a solid's boundary face is no internal face to unglue.
#[test]
fn faces_that_do_not_match_decline() {
    let (mut body, [keep, gone]) = box_on_box();
    let before = NonManifoldCounts::of(&body);
    let side = body
        .faces
        .keys()
        .find(|&f| f != keep && f != gone)
        .expect("another face");
    assert_eq!(
        body.glue_faces(keep, side),
        Err(EulerError::FacesDoNotMatch { keep, gone: side })
    );
    assert_eq!(
        body.glue_faces(keep, keep),
        Err(EulerError::FacesDoNotMatch { keep, gone: keep })
    );
    assert_eq!(
        body.unglue_face(keep),
        Err(EulerError::NotAnInternalFace { face: keep })
    );
    assert_eq!(NonManifoldCounts::of(&body), before);
}
