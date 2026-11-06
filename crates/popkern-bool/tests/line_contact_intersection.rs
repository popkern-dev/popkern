//! A cylinder lying on a box touches it along a line (PK-966).
//!
//! The regularized intersection of bodies that only touch is empty.
//! It used to come back as a solid wrapping one shell with no faces —
//! `validate()` clean, χ-clean — so a consumer counting solids saw a
//! body that was not there. Empty results are bare bodies, and
//! `validate_all` now rejects a shell that contains no face.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn cylinder_on_box() -> (BRep, BRep) {
    let tol = Tolerance::default();
    let block = BRep::unit_cube(tol).expect("cube builds");
    // Axis along x after the rotation, resting on the top face z = 1.
    let rod = BRep::cylinder(0.5, 1.0, tol)
        .expect("cylinder builds")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("rotates")
        .transform(&Transform3::translation(0.0, 0.5, 1.5))
        .expect("translates");
    (block, rod)
}

#[test]
fn line_contact_intersection_is_an_empty_body() {
    let (block, rod) = cylinder_on_box();
    let r = boolean(
        BoolOperation::Intersection,
        block,
        rod,
        Tolerance::default(),
        None,
    )
    .expect("intersection serves");
    assert_eq!(r.faces.len(), 0);
    assert_eq!(r.back_shell_count(), 0, "no face-less shell");
    assert_eq!(r.solid_regions().len(), 0, "no solid around nothing");
    assert!(r.validate_all().is_empty());
}

/// Negative control: the validator catches the face-less shell the
/// boolean used to serve.
#[test]
fn validate_all_rejects_a_shell_without_faces() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let solid = cube.solid_regions()[0];
    cube.add_region_shell(solid, Default::default());
    let errors = cube.validate_all();
    assert!(
        errors
            .iter()
            .any(|e| e.to_string().contains("contains no face")),
        "{errors:?}"
    );
}
