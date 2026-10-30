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

/// A unit box with an r = 0.5 rod lying along its top face's centre
/// line, the rod's end caps flush with the box's x = 0 and x = 1 faces:
/// each cap's rim circle touches the box's top edge at one point.
fn rod_on_box(place: &Transform3, tol: Tolerance) -> (BRep, BRep) {
    let block = BRep::unit_cube(tol)
        .expect("cube")
        .transform(place)
        .expect("placed cube");
    let rod = BRep::cylinder(0.5, 1.0, tol)
        .expect("rod")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("along x")
        .transform(&Transform3::translation(0.0, 0.5, 1.5))
        .expect("resting")
        .transform(place)
        .expect("placed rod");
    (block, rod)
}

/// The rim circles touch the end faces' top edges tangentially, so no
/// arc of them is a boundary of an end face: the complex intersects
/// the touching bodies to nothing and reports the union's contacts,
/// instead of declining on a face whose two sides joined (an end face
/// that took the whole rim circle as a hole, and a phantom disc of it).
#[test]
fn a_rod_resting_on_a_box_with_flush_ends_takes_the_complex() {
    let tol = Tolerance::default();
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.1),
        Vec3::new(1.0, 2.0, 3.0),
        0.37,
    )
    .expect("rotation");
    for place in [Transform3::identity(), turned] {
        let (block, rod) = rod_on_box(&place, tol);
        let ix = crate::intersect::find_intersections(&block, &rod, &tol).expect("intersections");
        let merge = |op| {
            crate::cellular::boolean::cellular_graph_merge(
                &block,
                &rod,
                &ix,
                op,
                &tol,
                None,
                crate::complex::serve::SelectionSource::Complex,
            )
        };
        let (body, _) = merge(BoolOperation::Intersection)
            .unwrap_or_else(|d| panic!("the complex intersects: {}", d.into_error()));
        assert!(body.faces.is_empty(), "{} faces", body.faces.len());
        let union = merge(BoolOperation::Union).map(|_| ());
        let why = union
            .err()
            .map(|d| d.into_error().to_string())
            .unwrap_or_default();
        assert!(why.contains("non-manifold"), "union: {why}");
    }
}

/// A 20 × 20 × 4 plate with an r = 5 round notch in its x = 10 edge,
/// and an r = 1.6 pin from z = `pin_z.0` to `pin_z.1` standing in the
/// notch `depth` proud of its wall at x = 5: their intersection is a
/// lens `depth` deep. The notch is part of the plate faces' outer
/// loops, as a cycloid disc's window flank is.
fn pin_proud_of_notch(
    depth: f64,
    pin_z: (f64, f64),
    place: &Transform3,
    tol: Tolerance,
) -> (BRep, BRep) {
    let plate = BRep::prism(
        &[
            Point3::new(-10.0, -10.0, 0.0),
            Point3::new(10.0, -10.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(-10.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        4.0,
        tol,
    )
    .expect("plate");
    let notch = BRep::cylinder(5.0, 6.0, tol)
        .expect("notch")
        .transform(&Transform3::rotation_z(std::f64::consts::PI))
        .expect("seam at the lens")
        .transform(&Transform3::translation(10.0, 0.0, -1.0))
        .expect("through");
    let plate = crate::boolean(BoolOperation::Subtraction, plate, notch, tol, None)
        .expect("notched plate")
        .transform(place)
        .expect("placed plate");
    let pin = BRep::cylinder(1.6, pin_z.1 - pin_z.0, tol)
        .expect("pin")
        .transform(&Transform3::translation(
            10.0 - (5.0 - 1.6 + depth),
            0.0,
            pin_z.0,
        ))
        .expect("proud")
        .transform(place)
        .expect("placed pin");
    (plate, pin)
}

/// A shallow lens with the notch wall's seam vertex at its middle: the
/// pin's cap circles cross the notch arc twice each and the complex
/// intersects the bodies to the lens, under a general placement too.
/// (The cycloid disc's flank, many short edges, gives the cap circles
/// crossing noise at every junction; that case is pinned in
/// `popkern-rhai`'s `a_rotor_geared_at_minus_n_interferes`, PK-1022.)
#[test]
fn a_pin_proud_of_its_notch_wall_intersects_to_the_lens() {
    let tol = Tolerance::default();
    let depth = 2.7e-3;
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.1),
        Vec3::new(1.0, 2.0, 3.0),
        0.37,
    )
    .expect("rotation");
    for place in [Transform3::identity(), turned] {
        let (plate, pin) = pin_proud_of_notch(depth, (-1.0, 5.0), &place, tol);
        let ix = crate::intersect::find_intersections(&plate, &pin, &tol).expect("intersections");
        let (mut body, _) = crate::cellular::boolean::cellular_graph_merge(
            &plate,
            &pin,
            &ix,
            BoolOperation::Intersection,
            &tol,
            None,
            crate::complex::serve::SelectionSource::Complex,
        )
        .unwrap_or_else(|d| panic!("the complex intersects: {}", d.into_error()));
        body.validate().expect("valid lens");
        let _ = body.recompute_pcurves_where_missing();
        let boxes = crate::intersect::compute_face_aabbs(&body);
        for (x, inside) in [
            (5.0 - 0.5 * depth, true),
            (5.0 + 0.5 * depth, false),
            (5.0 - 1.5 * depth, false),
        ] {
            let p = place.apply_point(&Point3::new(x, 0.0, 2.0));
            let got =
                crate::classify::point_in_solid_final(&p, &body, &boxes, &tol).expect("classifies");
            let want = if inside {
                crate::PointClassification::Inside
            } else {
                crate::PointClassification::Outside
            };
            assert_eq!(got, want, "x = {x}: {} faces", body.faces.len());
        }
    }
}
