//! A conical dimple touching a plate's edge, the scene turned and the
//! cone turned 0.7 about its axis (PK-1143): its rim grazed the edge
//! mid-record in nine reports, kept whole as a graze that does not flip
//! in/out, and the odd count flipped every clip window after it; the
//! clip kept only the 3.9e-3 between two graze groups and the plate
//! served UNCUT (volume 2000 for 1932.98), silently. A group that does
//! not flip keeps an even number of reports. The corpus row this
//! belongs to is a target (other placements decline honestly), so this
//! run is pinned on its own: it must serve the dimpled plate.

use std::f64::consts::PI;

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn the_turned_cone_dimple_cuts_the_plate() {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(20.0, 0.0, z),
            Point3::new(20.0, 20.0, z),
            Point3::new(0.0, 20.0, z),
        ]
    };
    let plate = BRep::polygon_frustum(&rect(0.0), &rect(5.0), tol).expect("plate");
    let pivot = Point3::new(16.0, 10.0, 5.0);
    let cone = BRep::cone(4.0, 4.0, tol)
        .and_then(|c| c.transform(&Transform3::rotation_x(PI)))
        .and_then(|c| c.transform(&Transform3::translation(16.0, 10.0, 5.0)))
        .and_then(|c| {
            c.transform(
                &Transform3::rotation_around_axis(pivot, Vec3::new(0.0, 0.0, 1.0), 0.7)
                    .expect("turn"),
            )
        })
        .expect("cone");
    let place = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("place");
    let body = boolean(
        BoolOperation::Subtraction,
        plate.transform(&place).expect("place"),
        cone.transform(&place).expect("place"),
        tol,
        None,
    )
    .expect("serves");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let v = popkern_measure::body_volume(&body, &tol).expect("volume");
    let want = 2000.0 - 64.0 * PI / 3.0;
    assert!((v - want).abs() < 1e-6 * want, "volume {v} for {want}");
}
