//! A part cut by a part thinner than the linear tolerance serves whole
//! and its feature log says why (PK-1131).

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::Part;

#[test]
fn a_saw_under_the_tolerance_leaves_a_note_in_the_log() {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(10.0, 0.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(0.0, 10.0, z),
        ]
    };
    let cube = Part::from_brep(
        BRep::polygon_frustum(&rect(0.0), &rect(10.0), tol).expect("cube"),
        tol,
        "cube",
    )
    .expect("cube part");
    let saw_pts = [
        Point3::new(5.0, -1.0, 5.0),
        Point3::new(5.000001, -1.0, 5.0),
        Point3::new(5.000001, 11.0, 5.0),
        Point3::new(5.0, 11.0, 5.0),
    ];
    let saw = Part::from_brep(
        BRep::prism(&saw_pts, Vec3::new(0.0, 0.0, 1.0), 6.0, tol).expect("saw"),
        tol,
        "saw",
    )
    .expect("saw part");
    let cut = cube.subtract(&saw).expect("serves");
    assert_eq!(cut.brep().faces.len(), 6);
    let last = cut.features().entries().last().expect("entry");
    assert!(
        last.summary.contains("under the linear tolerance"),
        "{}",
        last.summary
    );
}
