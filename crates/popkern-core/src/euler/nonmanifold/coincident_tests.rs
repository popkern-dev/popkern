use super::super::NonManifoldCounts;
use super::*;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// Two cubes touching along an edge, nothing shared: the glue makes one
/// vertex at each end and one ring of four, a general body that
/// validates.
#[test]
fn cubes_touching_along_an_edge_glue_into_one_body() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(1.0, 1.0, 0.0));
    let glued = body.glue_coincident().expect("glue");
    assert_eq!(
        glued,
        CoincidentGlue {
            vertices: 2,
            edges: 1
        }
    );
    assert_eq!(body.body_kind(), BodyKind::General);
    assert!(body.edges.keys().any(|e| body.radial_fins(e).len() == 4));
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// Negative control: a valid body has nothing coincident to glue.
#[test]
fn a_valid_body_has_nothing_to_glue() {
    let mut cube = cube_at(0.0, 0.0, 0.0);
    let before = NonManifoldCounts::of(&cube);
    assert_eq!(cube.glue_coincident(), Ok(CoincidentGlue::default()));
    assert_eq!(NonManifoldCounts::of(&cube), before);
}
