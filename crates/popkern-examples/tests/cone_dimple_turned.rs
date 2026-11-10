//! A conical dimple touching a plate's edge, the scene turned (PK-1143).
//!
//! With the cone turned 0.7 about its axis the rim grazed the edge
//! mid-record in nine reports, kept whole as a graze that does not flip
//! in/out, and the odd count flipped every clip window after it; the
//! clip kept only the 3.9e-3 between two graze groups and the plate
//! served UNCUT (volume 2000 for 1932.98), silently. A group that does
//! not flip keeps an even number of reports.
//!
//! Unturned, the rim circle's closure sat 0.24 short of the touch, and
//! the clip's far retrim moved that end onto the touch junction, which
//! the loop passes mid-curve: the rim came back as the 0.24 between
//! them and the cut declined. A piece that is its closed source whole
//! has no window end to refine.
//!
//! At kilometre scale with the cone turned pi, the line x cone solver
//! read the plate edge touching the rim as a miss (its discriminant
//! band was in squared units), so the rim's two window ends either
//! side of the touch were never refined onto it and the cut declined.
//! Tangency is decided by the line's distance to the cone.

use std::f64::consts::PI;

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn turned() -> Transform3 {
    Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("place")
}

/// The plate less the dimple, the cone turned `self_turn` about its
/// axis, both placed by `place` at `scale`: the dimpled plate's volume.
fn assert_dimple_cuts(self_turn: f64, place: Transform3, scale: f64) {
    let base = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(20.0, 0.0, z),
            Point3::new(20.0, 20.0, z),
            Point3::new(0.0, 20.0, z),
        ]
    };
    let plate = BRep::polygon_frustum(&rect(0.0), &rect(5.0), base).expect("plate");
    let pivot = Point3::new(16.0, 10.0, 5.0);
    let cone = BRep::cone(4.0, 4.0, base)
        .and_then(|c| c.transform(&Transform3::rotation_x(PI)))
        .and_then(|c| c.transform(&Transform3::translation(16.0, 10.0, 5.0)))
        .and_then(|c| {
            c.transform(
                &Transform3::rotation_around_axis(pivot, Vec3::new(0.0, 0.0, 1.0), self_turn)
                    .expect("turn"),
            )
        })
        .expect("cone");
    let tol = base.with_linear(base.linear() * scale.max(1.0));
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
    let want = (2000.0 - 64.0 * PI / 3.0) * scale.powi(3);
    assert!((v - want).abs() < 1e-6 * want, "volume {v} for {want}");
}

#[test]
fn the_turned_cone_dimple_cuts_the_plate() {
    assert_dimple_cuts(0.7, turned(), 1.0);
}

#[test]
fn a_rim_closing_short_of_the_touch_keeps_its_whole_loop() {
    assert_dimple_cuts(0.0, turned(), 1.0);
    assert_dimple_cuts(PI, turned(), 1.0);
}

#[test]
fn a_kilometre_rim_touching_the_edge_is_refined_onto_the_touch() {
    assert_dimple_cuts(PI, Transform3::scale_uniform(1e3), 1e3);
}
