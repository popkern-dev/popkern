use super::super::NonManifoldCounts;
use super::*;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

fn placed(body: BRep, t: Transform3) -> BRep {
    body.transform(&t).expect("placed")
}

/// `body` moved by each of `steps` in turn.
fn moved(body: BRep, steps: &[Transform3]) -> BRep {
    steps.iter().fold(body, |b, t| placed(b, *t))
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    placed(
        BRep::unit_cube(tol()).expect("cube"),
        Transform3::translation(x, y, z),
    )
}

/// The unit cube with `other` appended beside it.
fn with(other: &BRep) -> BRep {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(other);
    body
}

fn region_volumes_sum(body: &BRep) -> usize {
    body.regions
        .keys()
        .filter(|&r| body.region_kind(r) == Some(crate::topology::RegionKind::Solid))
        .count()
}

/// Two cubes sharing an edge: their coincident vertices and edge glue
/// into one general body, nothing imprinted.
#[test]
fn cubes_on_an_edge_glue_from_geometry() {
    let mut body = with(&cube_at(1.0, 1.0, 0.0));
    let before = NonManifoldCounts::of(&body);
    let done = body.glue_touching().expect("glue");
    assert_eq!(done.imprinted, 0);
    assert_eq!((done.glued.vertices, done.glued.edges), (2, 1));
    let after = NonManifoldCounts::of(&body);
    assert_eq!(after.vertices, before.vertices - 2);
    assert_eq!(after.edges, before.edges - 1);
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_eq!(region_volumes_sum(&body), 2);
    assert_valid(&body, "cubes on an edge");
}

/// A cube standing on its edge across the unit cube's top face: the
/// face's edges split under the edge's ends, the face imprinted along
/// it, and the edge glued into a ring of four.
#[test]
fn an_edge_across_a_face_is_imprinted_and_glued() {
    let standing = moved(
        BRep::unit_cube(tol()).expect("cube"),
        &[
            Transform3::translation(0.0, -0.5, -0.5),
            Transform3::rotation_x(std::f64::consts::FRAC_PI_4),
            Transform3::translation(0.0, 0.5, 1.0 + std::f64::consts::FRAC_1_SQRT_2),
        ],
    );
    let mut body = with(&standing);
    let faces = body.faces.len();
    let done = body.glue_touching().expect("glue");
    assert_eq!(done.curves, 1, "{done:?}");
    assert_eq!(done.edges_split, 2, "{done:?}");
    assert_eq!(done.imprinted, 1, "{done:?}");
    assert_eq!(done.glued.edges, 1, "{done:?}");
    assert_eq!(body.faces.len(), faces + 1);
    assert!(body.edges.keys().any(|e| body.radial_fins(e).len() == 4));
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "edge across a face");
}

/// A rod lying on the top face away from its edges: the tangent ruling
/// ends on the face's interior (vertex loops joined into a slit) and on
/// the rod's cap circles (the side face split along the ruling).
#[test]
fn a_rod_lying_on_a_face_touches_along_a_slit() {
    let rod = moved(
        BRep::cylinder(0.25, 0.6, tol()).expect("rod"),
        &[
            Transform3::rotation_y(std::f64::consts::FRAC_PI_2),
            Transform3::translation(0.2, 0.5, 1.25),
        ],
    );
    let mut body = with(&rod);
    let done = body.glue_touching().expect("glue");
    assert_eq!(done.curves, 1, "{done:?}");
    assert_eq!(done.vertex_loops, 2, "{done:?}");
    assert_eq!(done.glued.edges, 1, "{done:?}");
    assert!(body.edges.keys().any(|e| body.radial_fins(e).len() == 4));
    assert!(body.faces.values().any(|f| f.inner_loops.len() == 1));
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "rod on a face");
}

/// A ball resting on the top face: a vertex loop on each face, glued.
#[test]
fn a_ball_resting_on_a_face_touches_at_a_vertex() {
    // Tilted, so the sphere's pole is off the contact and its face
    // needs a vertex loop too.
    let ball = moved(
        BRep::sphere(0.3, tol()).expect("ball"),
        &[
            Transform3::rotation_x(0.7),
            Transform3::translation(0.5, 0.4, 1.3),
        ],
    );
    let mut body = with(&ball);
    let done = body.glue_touching().expect("glue");
    assert_eq!(done.points, 1, "{done:?}");
    assert_eq!(done.glued.vertices, 1, "{done:?}");
    assert_eq!(body.vertex_loops().len(), 2);
    assert_eq!(body.body_kind(), BodyKind::General);
    assert_valid(&body, "ball on a face");
}

/// Pieces apart touch nowhere: nothing found, the body unchanged
/// (negative control).
#[test]
fn pieces_apart_are_left_alone() {
    let mut body = with(&cube_at(1.5, 0.0, 0.0));
    let before = NonManifoldCounts::of(&body);
    let done = body.glue_touching().expect("glue");
    assert_eq!(done, TouchGlue::default());
    assert_eq!(NonManifoldCounts::of(&body), before);
}

/// A rod as wide as the cube lying on its top face, its seam on the
/// contact line, near the origin, far from it and at a hundredth scale
/// (where the trim classifier's band is ten tolerances wide): one contact curve
/// between two points (the seam lies along the ruling, so no run end
/// may slide along it; the cap circles graze the side faces' top edges,
/// which is no contact curve), imprinted on the top face and glued.
#[test]
fn a_rod_on_its_seam_touches_along_one_curve_near_and_far() {
    for shift in [
        Transform3::identity(),
        Transform3::translation(1.0e4, -2.0e4, 3.0e4),
        Transform3::scale_uniform(1.0e-2),
    ] {
        let rod = moved(
            BRep::cylinder(0.5, 1.0, tol()).expect("rod"),
            &[
                Transform3::rotation_y(std::f64::consts::FRAC_PI_2),
                Transform3::translation(0.0, 0.5, 1.5),
                shift,
            ],
        );
        let mut body = placed(cube_at(0.0, 0.0, 0.0), shift);
        body.append_body(&rod);
        let done = body.glue_touching().expect("glue");
        assert_eq!((done.points, done.curves), (2, 1), "{done:?}");
        assert_eq!((done.imprinted, done.glued.edges), (1, 1), "{done:?}");
        assert_eq!(body.body_kind(), BodyKind::General);
        assert_valid(&body, "rod on its seam");
    }
}
