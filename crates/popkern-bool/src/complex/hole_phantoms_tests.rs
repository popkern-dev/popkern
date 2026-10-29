//! A hub standing in a bearing race's bore, their bottoms on one plane
//! (the cycloid M module's flange and inner race, PK-756): the race's
//! bottom walks a disc over its own bore, and the complex leaves it out
//! instead of declining an ambiguous radial order on the bore rim.

use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// A ring r ∈ [inner, outer], z ∈ [0, height].
fn ring(inner: f64, outer: f64, height: f64, tol: Tolerance) -> BRep {
    let body = BRep::cylinder(outer, height, tol).expect("body");
    let bore = BRep::cylinder(inner, height + 2.0, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("placed bore");
    crate::boolean(BoolOperation::Subtraction, body, bore, tol, None).expect("ring")
}

#[test]
fn a_disc_over_a_faces_own_bore_bounds_nothing() {
    let tol = Tolerance::default();
    let hub = ring(27.0, 29.0, 20.0, tol);
    let race = ring(29.0, 34.0, 9.0, tol);
    for (a, b) in [(&hub, &race), (&race, &hub)] {
        let ix = crate::intersect::find_intersections(a, b, &tol).expect("intersections");
        let served =
            crate::complex::serve::rescue(a, b, &ix, BoolOperation::Union, &tol).expect("manifold");
        let (body, _) = served.expect("the complex serves the union");
        body.validate().expect("valid");
        let boxes = crate::intersect::compute_face_aabbs(&body);
        for (p, inside) in [
            (Point3::new(28.0, 0.0, 5.0), true),
            (Point3::new(31.0, 0.0, 5.0), true),
            (Point3::new(0.0, 28.0, 15.0), true),
            (Point3::new(31.0, 0.0, 15.0), false),
            (Point3::new(0.0, 0.0, 5.0), false),
            (Point3::new(0.0, -20.0, 5.0), false),
        ] {
            let got =
                crate::classify::point_in_solid_final(&p, &body, &boxes, &tol).expect("classifies");
            assert_eq!(
                got == crate::PointClassification::Inside,
                inside,
                "at {p:?}"
            );
        }
    }
}
