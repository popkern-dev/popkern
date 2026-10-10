//! The ray classifier's search window on a plane face whose curved rim
//! bulges past its few vertices (PK-752).

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::classify::point_in_solid_final;
use crate::error::PointClassification;
use crate::intersect::compute_face_aabbs;

/// A disc (r 60, z 0..12) whose BOTTOM face gains notch corners from a
/// block poking in at x ≈ −55 (through the bottom only), so its loop
/// carries a few vertices while the rest of its rim bulges 110 away;
/// the top face keeps its single-vertex rim.
fn notched_disc() -> BRep {
    let tol = Tolerance::default();
    let disc = BRep::cylinder(60.0, 12.0, tol).expect("disc");
    let block = BRep::prism(
        &[
            Point3::new(-70.0, -10.0, -20.0),
            Point3::new(-50.0, -10.0, -20.0),
            Point3::new(-50.0, 10.0, -20.0),
            Point3::new(-70.0, 10.0, -20.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        26.0,
        tol,
    )
    .expect("block");
    crate::union(disc, block, tol, None).expect("union")
}

#[test]
fn a_point_below_a_notched_disc_far_from_its_vertices_is_outside() {
    let body = notched_disc();
    let aabbs = compute_face_aabbs(&body);
    let tol = Tolerance::default();
    for p in [
        Point3::new(30.0, 0.0, -1.0),
        Point3::new(-20.0, 30.0, -1.0),
        Point3::new(10.0, -40.0, -0.5),
    ] {
        let c = point_in_solid_final(&p, &body, &aabbs, &tol).expect("classified");
        assert_eq!(c, PointClassification::Outside, "{p:?}");
    }
    let inside = Point3::new(30.0, 0.0, 6.0);
    let c = point_in_solid_final(&inside, &body, &aabbs, &tol).expect("classified");
    assert_eq!(c, PointClassification::Inside);
    // The same scene turned: the window is the face's own, in any frame.
    let turned = body
        .transform(&Transform3::rotation_y(0.4).compose(&Transform3::rotation_z(0.7)))
        .expect("turned");
    let t = Transform3::rotation_y(0.4).compose(&Transform3::rotation_z(0.7));
    let aabbs = compute_face_aabbs(&turned);
    let c = point_in_solid_final(
        &t.apply_point(&Point3::new(30.0, 0.0, -1.0)),
        &turned,
        &aabbs,
        &tol,
    )
    .expect("classified");
    assert_eq!(c, PointClassification::Outside, "turned");
}
