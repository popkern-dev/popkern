//! A tool inside A that shares coplanar faces with A (PK-965).
//!
//! The coincident-containment shortcut used to serve A − B as A with
//! a cavity whenever every face of B probed Inside A. A corner notch
//! flush on three of A's faces then came back as an outer shell plus
//! a cavity shell lying on it — validator-clean, topologically wrong
//! (the notch is open to the outside). Containment decides a
//! subtraction only when A ⊂ B; these pins hold the notch, its
//! rotated and scaled images, and the cases the shortcut still owns.

use popkern_bool::{
    BoolOperation, PointClassification, boolean, face_bounding_boxes, point_in_solid,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn block(lo: (f64, f64, f64), size: f64, place: &Transform3) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(size))
        .expect("scales")
        .transform(&Transform3::translation(lo.0, lo.1, lo.2))
        .expect("translates")
        .transform(place)
        .expect("places")
}

fn classify(body: &BRep, p: Point3, tol: &Tolerance) -> PointClassification {
    point_in_solid(&p, body, &face_bounding_boxes(body), tol).expect("classifies")
}

/// The notch must be one shell with the topology of a cube missing a
/// corner octant, the notch point Outside and the remaining material
/// Inside — whatever rigid motion and scale place the pair.
fn assert_corner_notch(place: Transform3, scale: f64) {
    let tol = Tolerance::default();
    let a = block((0.0, 0.0, 0.0), 2.0 * scale, &place);
    let b = block((scale, 0.0, 0.0), scale, &place);
    let r = boolean(BoolOperation::Subtraction, a, b, tol, None).expect("notch serves");
    r.validate().expect("valid");
    assert_eq!(r.shells.len(), 1, "a flush notch is open, not a cavity");
    assert_eq!(
        (r.vertices.len(), r.edges.len(), r.faces.len()),
        (14, 21, 9),
        "cube-minus-octant topology"
    );
    let notch = place.apply_point(&Point3::new(1.5 * scale, 0.5 * scale, 0.5 * scale));
    let kept = place.apply_point(&Point3::new(0.5 * scale, 1.5 * scale, 1.5 * scale));
    assert_eq!(classify(&r, notch, &tol), PointClassification::Outside);
    assert_eq!(classify(&r, kept, &tol), PointClassification::Inside);
}

#[test]
fn three_flush_corner_notch_is_open() {
    assert_corner_notch(Transform3::identity(), 1.0);
}

#[test]
fn three_flush_corner_notch_is_open_rotated() {
    let axis = Vec3::new(1.0, 2.0, 3.0);
    for angle in [0.37, 1.9, 4.2] {
        let place = Transform3::rotation_around_axis(Point3::new(0.3, -0.2, 0.1), axis, angle)
            .expect("rotation");
        assert_corner_notch(place, 1.0);
    }
}

#[test]
fn three_flush_corner_notch_is_open_at_scale() {
    for scale in [1.0e-2, 1.0e2] {
        assert_corner_notch(Transform3::identity(), scale);
    }
}

/// The shortcut still owns the cases containment decides: B ⊂ A for
/// union (A) and intersection (B), and A ⊂ B for subtraction (empty).
#[test]
fn containment_still_decides_union_intersection_and_swallowed_subtraction() {
    let tol = Tolerance::default();
    let id = Transform3::identity();
    let a = || block((0.0, 0.0, 0.0), 2.0, &id);
    let b = || block((1.0, 0.0, 0.0), 1.0, &id);
    let u = boolean(BoolOperation::Union, a(), b(), tol, None).expect("union");
    assert_eq!((u.vertices.len(), u.faces.len(), u.shells.len()), (8, 6, 1));
    let i = boolean(BoolOperation::Intersection, a(), b(), tol, None).expect("intersection");
    assert_eq!((i.vertices.len(), i.faces.len(), i.shells.len()), (8, 6, 1));
    let s = boolean(BoolOperation::Subtraction, b(), a(), tol, None).expect("swallowed");
    assert_eq!(s.faces.len(), 0, "B − A is empty when B ⊂ A");
}
