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

/// The box [x0,x1]×[y0,y1]×[z0,z1].
fn block(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let profile = [
        Point3::new(lo.0, lo.1, lo.2),
        Point3::new(hi.0, lo.1, lo.2),
        Point3::new(hi.0, hi.1, lo.2),
        Point3::new(lo.0, hi.1, lo.2),
    ];
    BRep::prism(
        &profile,
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        hi.2 - lo.2,
        tol,
    )
    .expect("block")
}

/// A box with a tab on its back face (the face keeps the tab's
/// footprint as a hole), crossed by a slab whose ends lie on the box's
/// sides (a coincident interface): the slab splits the back face into
/// three, and the footprint lies outside the largest, which the walk
/// elects outer. The walk files the footprint ring apart; the complex
/// puts it back as the hole of the piece around it (measured on the
/// print-in-place hinge's swung lid: three knuckle footprints left
/// twelve free edges).
#[test]
fn a_hole_ring_filed_apart_rejoins_the_piece_around_it() {
    let tol = Tolerance::default();
    let body = block((0.0, 0.0, 0.0), (80.0, 50.0, 30.0), tol);
    let tab = block((60.0, 50.0, 27.0), (72.0, 52.0, 29.0), tol);
    let tabbed = crate::boolean(BoolOperation::Union, body, tab, tol, None).expect("tabbed box");
    let slab = block((0.0, 45.0, 19.0), (80.0, 55.0, 21.0), tol);
    for op in [BoolOperation::Intersection, BoolOperation::Union] {
        let ix = crate::intersect::find_intersections(&tabbed, &slab, &tol).expect("intersections");
        let served =
            crate::complex::serve::rescue(&tabbed, &slab, &ix, op, &tol).expect("manifold");
        let (body, _) = served.unwrap_or_else(|| panic!("the complex serves the {op:?}"));
        body.validate().expect("valid");
    }
}
