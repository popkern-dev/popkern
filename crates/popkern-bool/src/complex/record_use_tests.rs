//! A record piece outside the unused face's trim is a phantom split,
//! not a leak: the keyway slot cut into a bored block serves through
//! the complex.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// The axis-aligned box over `[x0, x1] × [y0, y1] × [z0, z1]`.
fn cuboid(x: (f64, f64), y: (f64, f64), z: (f64, f64), tol: Tolerance) -> BRep {
    let profile = [
        Point3::new(x.0, y.0, z.0),
        Point3::new(x.1, y.0, z.0),
        Point3::new(x.1, y.1, z.0),
        Point3::new(x.0, y.1, z.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), z.1 - z.0, tol).expect("cuboid")
}

/// A 1 × 1 × 0.3 block centred on the z axis, bored through by a
/// 0.15-radius hole.
fn bored_block(tol: Tolerance) -> BRep {
    let block = cuboid((-0.5, 0.5), (-0.5, 0.5), (0.0, 0.3), tol);
    let bore = BRep::cylinder(0.15, 0.5, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -0.1))
        .expect("placed bore");
    crate::boolean(BoolOperation::Subtraction, block, bore, tol, None).expect("bored block")
}

/// A 0.08 × 0.06 slot over y ∈ [0.14, 0.20], through the block: its
/// inner wall y = 0.14 crosses the bore, so the records where that wall
/// meets the block's top and bottom faces run, between x = ±0.04,
/// across those faces' holes.
fn key(tol: Tolerance) -> BRep {
    cuboid((-0.04, 0.04), (0.14, 0.2), (-0.14, 0.44), tol)
}

#[test]
fn a_record_across_a_hole_does_not_decline_the_complex() {
    let tol = Tolerance::default();
    let block = bored_block(tol);
    let key = key(tol);
    let ix = crate::intersect::find_intersections(&block, &key, &tol).expect("intersections");
    let served = crate::complex::serve::rescue(&block, &key, &ix, BoolOperation::Subtraction, &tol)
        .expect("manifold");
    let (body, _) = served.expect("the complex serves the keyway cut");
    body.validate().expect("valid");
    // Material where the block stays, none in the bore or the slot.
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for (p, inside) in [
        ((0.3, 0.0, 0.15), true),
        ((0.0, -0.3, 0.15), true),
        ((0.1, 0.17, 0.15), true),
        ((0.0, 0.0, 0.15), false),
        ((0.0, 0.17, 0.15), false),
        ((0.0, 0.17, 0.02), false),
        ((0.0, 0.17, -0.05), false),
    ] {
        let point = Point3::new(p.0, p.1, p.2);
        let got =
            crate::classify::point_in_solid_final(&point, &body, &boxes, &tol).expect("classifies");
        let want = if inside {
            crate::PointClassification::Inside
        } else {
            crate::PointClassification::Outside
        };
        assert_eq!(got, want, "at {p:?}");
    }
}
