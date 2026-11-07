use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::options::{BoolOptions, NonManifoldPolicy};
use crate::{BoolOperation, boolean_with_options};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

fn under(policy: NonManifoldPolicy, op: BoolOperation, a: BRep, b: BRep) -> BRep {
    boolean_with_options(
        op,
        a,
        b,
        tol(),
        None,
        &BoolOptions {
            non_manifold: policy,
        },
    )
    .expect("boolean serves")
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

/// The volume the body's faces enclose (outward faces).
fn volume(body: &BRep) -> f64 {
    let faces: Vec<_> = body.faces.keys().collect();
    body.enclosed_folded_volume(&faces).expect("volume")
}

/// Two boxes meeting along an edge union into one general body: the
/// edge's ring holds both boxes' faces, the two solids stay apart, and
/// nothing is left recorded. Split again it is the split policy's
/// pieces.
#[test]
fn boxes_meeting_along_an_edge_union_into_one_general_body() {
    let general = under(
        NonManifoldPolicy::General,
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
    );
    assert_eq!(general.body_kind(), BodyKind::General);
    assert!(general.topology_contacts.is_empty());
    assert_eq!(general.solid_regions().len(), 2);
    assert!(
        general
            .edges
            .keys()
            .any(|e| general.radial_fins(e).len() == 4)
    );
    assert_valid(&general, "general union");

    let pieces = under(
        NonManifoldPolicy::SplitToManifold,
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
    );
    let mut split = general.clone();
    split.split_to_manifold().expect("split");
    let counts = |b: &BRep| (b.vertices.len(), b.edges.len(), b.faces.len());
    // The same pieces: the split policy names its contacts by the faces
    // that touch, the split by the coincident edges and vertex copies;
    // both validate.
    assert_eq!(counts(&split), counts(&pieces));
    assert_valid(&split, "split back into pieces");
    let v = volume(&split);
    assert!((v - 2.0).abs() < 1e-6, "{v}");
}

/// Two boxes sharing only a corner union into one general body with one
/// vertex of two fans.
#[test]
fn boxes_meeting_at_a_corner_union_into_one_general_body() {
    let general = under(
        NonManifoldPolicy::General,
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 1.0),
    );
    assert_eq!(general.body_kind(), BodyKind::General);
    assert!(general.topology_contacts.is_empty());
    assert_eq!(
        general
            .vertices
            .keys()
            .filter(|&v| general.is_nonmanifold_vertex(v))
            .count(),
        1
    );
    assert_valid(&general, "corner union");
}

/// A general body is an operand: uniting it with a third box beside
/// one of its boxes along another edge serves a general body of three
/// solids; cutting a notch out of one box keeps the contact edge.
#[test]
fn a_general_body_is_a_boolean_operand() {
    let general = under(
        NonManifoldPolicy::General,
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
    );
    let three = under(
        NonManifoldPolicy::General,
        BoolOperation::Union,
        general.clone(),
        cube_at(-1.0, 1.0, 0.0),
    );
    assert_eq!(three.body_kind(), BodyKind::General);
    assert_eq!(three.solid_regions().len(), 3);
    assert_valid(&three, "three boxes");

    let notched = under(
        NonManifoldPolicy::General,
        BoolOperation::Subtraction,
        general,
        cube_at(-0.5, -0.5, 0.5),
    );
    assert_eq!(notched.body_kind(), BodyKind::General);
    assert!(
        notched
            .edges
            .keys()
            .any(|e| notched.radial_fins(e).len() == 4)
    );
    assert_valid(&notched, "notched");
    let mut split = notched;
    split.split_to_manifold().expect("split");
    let v = volume(&split);
    assert!((v - (2.0 - 0.125)).abs() < 1e-6, "{v}");
}

/// The default policy still declines a general operand's boolean, and
/// the N-ary boolean serves general bodies under the general policy.
#[test]
fn the_general_policy_is_opt_in_and_reaches_nary_booleans() {
    let general = under(
        NonManifoldPolicy::General,
        BoolOperation::Union,
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
    );
    assert!(matches!(
        crate::boolean(
            BoolOperation::Union,
            general,
            cube_at(5.0, 0.0, 0.0),
            tol(),
            None
        ),
        Err(crate::BoolError::UnsupportedBodyKind { .. })
    ));
    let (a, b, c) = (
        cube_at(0.0, 0.0, 0.0),
        cube_at(1.0, 1.0, 0.0),
        cube_at(2.0, 0.0, 0.0),
    );
    let nary = crate::nary_boolean(
        &[&a, &b, &c],
        &crate::SetExpr::Any,
        tol(),
        &BoolOptions {
            non_manifold: NonManifoldPolicy::General,
        },
    )
    .expect("n-ary");
    assert_eq!(nary.body_kind(), BodyKind::General);
    assert!(nary.topology_contacts.is_empty());
    assert_valid(&nary, "n-ary general");
}
