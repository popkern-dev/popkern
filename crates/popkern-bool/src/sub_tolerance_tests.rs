use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;

/// A profile 12 long and `thickness` wide extruded 6: how a script
/// makes a saw (the frustum constructor refuses a wall under the
/// tolerance).
fn slab(thickness: f64, place: &Transform3) -> BRep {
    let pts = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(12.0, 0.0, 0.0),
        Point3::new(12.0, thickness, 0.0),
        Point3::new(0.0, thickness, 0.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 6.0, Tolerance::default())
        .and_then(|b| b.transform(place))
        .expect("slab")
}

fn places() -> [Transform3; 2] {
    [
        Transform3::identity(),
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis"),
    ]
}

#[test]
fn a_slab_a_tolerance_thick_has_no_inside() {
    let tol = Tolerance::default();
    for place in places() {
        let thin = thin_operand(&slab(1e-6, &place), &tol).expect("thin");
        assert!(thin.width <= tol.linear(), "{thin:?}");
        assert!(thin_operand(&slab(0.5e-6, &place), &tol).is_some());
        // Ten tolerances thick, and a unit cube, have an inside.
        assert!(thin_operand(&slab(1e-5, &place), &tol).is_none());
        assert!(thin_operand(&slab(1.0, &place), &tol).is_none());
    }
}

#[test]
fn a_thin_operand_resolves_to_the_other_or_empty() {
    let tol = Tolerance::default();
    let (body, thin) = (
        slab(1.0, &Transform3::identity()),
        slab(1e-6, &Transform3::identity()),
    );
    let faces = |b: &BRep| b.faces.len();
    assert_eq!(
        faces(&resolve(BoolOperation::Subtraction, &body, &thin, 1)),
        faces(&body)
    );
    assert_eq!(
        faces(&resolve(BoolOperation::Union, &body, &thin, 1)),
        faces(&body)
    );
    assert_eq!(
        faces(&resolve(BoolOperation::Union, &thin, &body, 0)),
        faces(&body)
    );
    assert_eq!(
        faces(&resolve(BoolOperation::Intersection, &body, &thin, 1)),
        0
    );
    assert_eq!(
        faces(&resolve(BoolOperation::Subtraction, &thin, &body, 0)),
        0
    );
    let _ = tol;
}
