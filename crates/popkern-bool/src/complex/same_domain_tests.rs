//! A tube flush in a D-bore (PK-926's fixture): the tube's wall and
//! the D's wall coincide below the D's top rim, and the tube's seam
//! splits its copy of that patch where the D's copy is whole. The
//! pieces join into one same-domain face with the D's fragment, so the
//! complex serves every operation instead of declining an ambiguous
//! radial order on the rim.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

fn tube_and_d(tol: Tolerance) -> (BRep, BRep) {
    let outer = BRep::cylinder(1.0, 1.0, tol).expect("outer");
    let bore = BRep::cylinder(0.9, 1.6, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -0.3))
        .expect("placed");
    let tube = crate::boolean(BoolOperation::Subtraction, outer, bore, tol, None).expect("tube");
    let disc = BRep::cylinder(1.0, 1.0, tol)
        .expect("disc")
        .transform(&Transform3::translation(0.0, 0.0, -0.5))
        .expect("placed");
    let square = [
        Point3::new(-2.0, -0.3, -1.0),
        Point3::new(2.0, -0.3, -1.0),
        Point3::new(2.0, 2.0, -1.0),
        Point3::new(-2.0, 2.0, -1.0),
    ];
    let mut slab = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("slab");
    slab.reverse_body_orientation();
    let d = crate::boolean(BoolOperation::Intersection, disc, slab, tol, None).expect("D");
    (tube, d)
}

#[test]
fn a_seam_split_copy_joins_its_coincident_partner() {
    let tol = Tolerance::default();
    let (tube, d) = tube_and_d(tol);
    // (point, in tube, in D)
    let probes = [
        (Point3::new(0.0, 0.95, 0.25), true, true),
        (Point3::new(0.0, -0.95, 0.25), true, false),
        (Point3::new(0.0, 0.95, 0.75), true, false),
        (Point3::new(0.0, 0.0, 0.0), false, true),
        (Point3::new(0.0, 0.0, 0.75), false, false),
    ];
    for op in [
        BoolOperation::Union,
        BoolOperation::Intersection,
        BoolOperation::Subtraction,
    ] {
        for swapped in [false, true] {
            let (a, b) = if swapped { (&d, &tube) } else { (&tube, &d) };
            let ix = crate::intersect::find_intersections(a, b, &tol).expect("intersections");
            let served = crate::complex::serve::rescue(a, b, &ix, op, &tol).expect("manifold");
            let (body, _) =
                served.unwrap_or_else(|| panic!("{op:?} swapped={swapped}: the complex serves"));
            body.validate()
                .unwrap_or_else(|e| panic!("{op:?} swapped={swapped}: {e}"));
            let boxes = crate::intersect::compute_face_aabbs(&body);
            for (p, in_tube, in_d) in probes {
                let (in_a, in_b) = if swapped {
                    (in_d, in_tube)
                } else {
                    (in_tube, in_d)
                };
                let want = match op {
                    BoolOperation::Union => in_a || in_b,
                    BoolOperation::Intersection => in_a && in_b,
                    BoolOperation::Subtraction => in_a && !in_b,
                };
                let got = crate::classify::point_in_solid_final(&p, &body, &boxes, &tol)
                    .expect("classifies");
                assert_eq!(
                    got == crate::PointClassification::Inside,
                    want,
                    "{op:?} swapped={swapped} at {p:?}"
                );
            }
        }
    }
}

/// A plate whose hole is most of its area, and a block flush with
/// both its faces straddling the hole's rim. On the plate's faces the
/// hole loop is rerouted round the block's footprint; as the tool that
/// big rerouted loop winding against the outer was filed as a separate
/// region, the footprint hung on the outer rim as a hole, and the
/// footprint copies never formed (the ring gear's clash lenses,
/// measured as `SameDomainMismatch` on every planet).
#[test]
fn a_tool_hole_rerouted_round_a_flush_block_stays_a_hole() {
    let tol = Tolerance::default();
    let square = [
        Point3::new(-5.0, -5.0, 0.0),
        Point3::new(5.0, -5.0, 0.0),
        Point3::new(5.0, 5.0, 0.0),
        Point3::new(-5.0, 5.0, 0.0),
    ];
    let plate = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 1.0, tol).expect("plate");
    let bore = BRep::cylinder(4.5, 3.0, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("through");
    let plate = crate::boolean(BoolOperation::Subtraction, plate, bore, tol, None).expect("ring");
    let foot = [
        Point3::new(4.0, -0.5, 0.0),
        Point3::new(4.8, -0.5, 0.0),
        Point3::new(4.8, 0.5, 0.0),
        Point3::new(4.0, 0.5, 0.0),
    ];
    let block = BRep::prism(&foot, Vec3::new(0.0, 0.0, 1.0), 1.0, tol).expect("block");
    // (point, in block, in plate): the block inside the hole, the
    // block over the plate, the plate beside the block.
    let probes = [
        (Point3::new(4.2, 0.0, 0.5), true, false),
        (Point3::new(4.7, 0.0, 0.5), true, true),
        (Point3::new(4.7, 2.0, 0.5), false, true),
    ];
    for op in [
        BoolOperation::Union,
        BoolOperation::Intersection,
        BoolOperation::Subtraction,
    ] {
        for swapped in [false, true] {
            let (a, b) = if swapped {
                (&plate, &block)
            } else {
                (&block, &plate)
            };
            let ix = crate::intersect::find_intersections(a, b, &tol).expect("intersections");
            let served = crate::complex::serve::rescue(a, b, &ix, op, &tol).expect("manifold");
            let (body, _) =
                served.unwrap_or_else(|| panic!("{op:?} swapped={swapped}: the complex serves"));
            body.validate()
                .unwrap_or_else(|e| panic!("{op:?} swapped={swapped}: {e}"));
            let boxes = crate::intersect::compute_face_aabbs(&body);
            for (p, in_block, in_plate) in probes {
                let (in_a, in_b) = if swapped {
                    (in_plate, in_block)
                } else {
                    (in_block, in_plate)
                };
                let want = match op {
                    BoolOperation::Union => in_a || in_b,
                    BoolOperation::Intersection => in_a && in_b,
                    BoolOperation::Subtraction => in_a && !in_b,
                };
                let got = crate::classify::point_in_solid_final(&p, &body, &boxes, &tol)
                    .expect("classifies");
                assert_eq!(
                    got == crate::PointClassification::Inside,
                    want,
                    "{op:?} swapped={swapped} at {p:?}"
                );
            }
        }
    }
}
