//! Membership of a shell class no face of an operand bounds is the
//! classification of a point on one of its faces: an offset probe can
//! cross a gap thinner than the offset into the other operand.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// A 100 × 100 × 10 plate bored through by an r = 10 hole, and an
/// r = 9.998 pin through the bore from z = −10 to 20: a 0.002 radial
/// gap, so the bodies share no material. The pin's wall crosses the
/// plate's top and bottom planes inside the bore (records that split
/// only the pin).
fn plate_and_pin(tol: Tolerance) -> (BRep, BRep) {
    let plate = BRep::prism(
        &[
            Point3::new(-50.0, -50.0, 0.0),
            Point3::new(50.0, -50.0, 0.0),
            Point3::new(50.0, 50.0, 0.0),
            Point3::new(-50.0, 50.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .expect("plate");
    let bore = BRep::cylinder(10.0, 12.0, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("placed bore");
    let plate =
        crate::boolean(BoolOperation::Subtraction, plate, bore, tol, None).expect("bored plate");
    let pin = BRep::cylinder(9.998, 30.0, tol)
        .expect("pin")
        .transform(&Transform3::translation(0.0, 0.0, -10.0))
        .expect("placed pin");
    (plate, pin)
}

#[test]
fn a_pin_clear_of_its_bore_intersects_to_nothing() {
    let tol = Tolerance::default();
    let (plate, pin) = plate_and_pin(tol);
    let ix = crate::intersect::find_intersections(&plate, &pin, &tol).expect("intersections");
    let served =
        crate::complex::serve::rescue(&plate, &pin, &ix, BoolOperation::Intersection, &tol)
            .expect("manifold");
    let (body, _) = served.expect("the complex serves the intersection");
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for p in [
        Point3::new(0.0, 0.0, 5.0),
        Point3::new(9.99, 0.0, 5.0),
        Point3::new(30.0, 30.0, 5.0),
    ] {
        let got = crate::classify::point_in_solid_final(&p, &body, &boxes, &tol);
        assert!(
            !matches!(got, Ok(crate::PointClassification::Inside)),
            "no material at {p:?}: {} faces",
            body.faces.len()
        );
    }
}
