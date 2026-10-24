//! PK-933: a box intersected with a sphere that overhangs all six faces
//! measures its closed-form volume under any rigid rotation.
//!
//! Every intersection is an exact plane x sphere circle, yet one rotation
//! read 1.6e-6 high: the volume quadrature's sliver correction locates a
//! boundary chord's ends on their trim, and on a CLOSED trim it clamped
//! the Newton step at the pcurve's parameter seam, so the chords ending
//! there matched no trim and dropped their slivers. The parameter of a
//! closed trim is periodic now.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::body_volume;

#[test]
fn a_box_cut_by_a_sphere_measures_its_volume_under_rotation() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    // Sphere radius 0.7 at the box centre, less six caps of height 0.2.
    let (r, h) = (0.7_f64, 0.2_f64);
    let want = pi * (4.0 / 3.0 * r * r * r - 6.0 * h * h * (3.0 * r - h) / 3.0);
    let rotations = [
        ("unrotated", Transform3::identity()),
        ("x 1.0", Transform3::rotation_x(1.0)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("x 1.2", Transform3::rotation_x(1.2)),
        ("y 1.1", Transform3::rotation_y(1.1)),
        ("z 0.7", Transform3::rotation_z(0.7)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ];
    for (name, rot) in rotations {
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&rot))
            .expect("box");
        let sphere = BRep::sphere(r, tol)
            .and_then(|s| s.transform(&Transform3::translation(0.5, 0.5, 0.5)))
            .and_then(|s| s.transform(&rot))
            .expect("sphere");
        let body = boolean(BoolOperation::Intersection, cube, sphere, tol, None)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(body.faces.len(), 7, "{name}: six discs and the sphere");
        let v = body_volume(&body, &tol).expect("volume");
        assert!(
            (v - want).abs() < 1e-10 * want,
            "{name}: volume {v} vs {want}"
        );
    }
}
