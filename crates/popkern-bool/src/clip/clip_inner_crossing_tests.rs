//! A pin standing in a plate's enclosed window, 2.7e-3 into the window
//! wall: its cap circles cross only the window's rim on the plate's
//! faces, mostly inside the hole.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{BoolOperation, PointClassification};

const DEPTH: f64 = 2.7e-3;

fn plate_and_pin(tol: Tolerance) -> (BRep, BRep) {
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
    let window = BRep::cylinder(5.0, 6.0, tol)
        .expect("window")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("through");
    let plate =
        crate::boolean(BoolOperation::Subtraction, plate, window, tol, None).expect("window");
    let pin = BRep::cylinder(1.6, 6.0, tol)
        .expect("pin")
        .transform(&Transform3::translation(5.0 - 1.6 + DEPTH, 0.0, -1.0))
        .expect("proud of the wall");
    (plate, pin)
}

/// Every operation serves valid, and points in the lens, the pin's
/// part in the window, the plate beside the lens and away from it fall
/// where each operation puts them. All three declined: the cap circles
/// dropped and the lens's records were stranded half claimed.
#[test]
fn a_pin_proud_of_an_enclosed_window_wall_serves_every_operation() {
    let tol = Tolerance::default();
    let (plate, pin) = plate_and_pin(tol);
    // The lens, the pin inside the window and the plate beside the
    // lens; then the plate far off and the pin above and below it.
    let points = [
        Point3::new(5.0 + 0.5 * DEPTH, 0.0, 2.0),
        Point3::new(5.0 - 0.5 * DEPTH, 0.0, 2.0),
        Point3::new(5.0 + 1.5 * DEPTH, 0.0, 2.0),
        Point3::new(7.0, 0.0, 2.0),
        Point3::new(3.4, 0.0, 4.5),
        Point3::new(3.4, 0.0, -0.5),
    ];
    // Rays from the window next to the lens cross the window wall's
    // notch: those points also check the classifier skips hits inside
    // it (PK-1038).
    for (op, want) in [
        (
            BoolOperation::Intersection,
            [true, false, false, false, false, false],
        ),
        (BoolOperation::Union, [true; 6]),
        (
            BoolOperation::Subtraction,
            [false, false, true, true, false, false],
        ),
    ] {
        let body = crate::boolean(op, plate.clone(), pin.clone(), tol, None)
            .unwrap_or_else(|e| panic!("{op:?}: {e}"));
        assert!(
            body.validate_all().is_empty(),
            "{op:?}: {:?}",
            body.validate_all()
        );
        let boxes = crate::intersect::compute_face_aabbs(&body);
        for (p, inside) in points.into_iter().zip(want) {
            let got = crate::classify::point_in_solid_final(&p, &body, &boxes, &tol);
            let want = if inside {
                PointClassification::Inside
            } else {
                PointClassification::Outside
            };
            assert_eq!(got.ok(), Some(want), "{op:?} at {p:?}");
        }
    }
}
