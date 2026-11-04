//! A sphere centred on a cube corner, its seam meridian lying in the
//! cube's y = 0 face (PK-1036): the face's record from the equator to
//! the pole is the seam's upper piece itself. Union, subtraction and
//! intersection serve the exact octant volumes, `validate_all` clean,
//! for the scene as built and rigidly rotated and scaled (the seam
//! stays in the face).
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn scene(motion: &Transform3, scale: f64) -> (BRep, BRep) {
    let tol = Tolerance::default();
    let place = |b: BRep| {
        b.transform(&Transform3::scale(scale, scale, scale))
            .and_then(|b| b.transform(motion))
            .expect("placed")
    };
    (
        place(BRep::unit_cube(tol).expect("cube")),
        place(BRep::sphere(0.7, tol).expect("sphere")),
    )
}

#[test]
fn a_sphere_whose_seam_lies_in_a_cube_face_serves_every_boolean() {
    let tol = Tolerance::default();
    let r: f64 = 0.7;
    let sphere = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    let octant = sphere / 8.0;
    let motions = [
        ("as built", Transform3::identity()),
        ("about z", Transform3::rotation_z(0.9)),
        (
            "about x then far",
            Transform3::translation(40.0, -7.0, 3.0).compose(&Transform3::rotation_x(2.1)),
        ),
    ];
    for (what, motion) in &motions {
        for scale in [1e-2, 1.0, 1e2] {
            let cases = [
                (BoolOperation::Union, 1.0 + sphere - octant),
                (BoolOperation::Subtraction, 1.0 - octant),
                (BoolOperation::Intersection, octant),
            ];
            for (op, unit_volume) in cases {
                let (cube, ball) = scene(motion, scale);
                let at = format!("{op:?} {what} x{scale}");
                let body =
                    boolean(op, cube, ball, tol, None).unwrap_or_else(|e| panic!("{at}: {e}"));
                let issues = body.validate_all();
                assert!(issues.is_empty(), "{at}: {issues:?}");
                let v = popkern_measure::body_volume(&body, &tol).expect("volume");
                let want = unit_volume * scale.powi(3);
                assert!(
                    (v - want).abs() <= 1e-6 * want,
                    "{at}: volume {v} vs {want}"
                );
            }
        }
    }
}
