//! Imprinted faces carry where they lie relative to the tool.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::FaceClass;

/// A cube imprinted by a cube overlapping its corner: the three faces
/// the tool crosses split into a piece inside the tool and one outside;
/// every face is classified, none on the tool.
#[test]
fn an_overlapping_cube_splits_faces_inside_and_outside() {
    let tol = Tolerance::default();
    let target = BRep::unit_cube(tol).expect("cube");
    let tool = BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::translation(0.5, 0.5, 0.5)))
        .expect("tool");
    let r = crate::ops::imprint::imprint(&target, &tool, &tol).expect("imprint");
    assert!(r.imprinted_target.validate().is_ok());
    assert_eq!(r.classes.len(), r.imprinted_target.faces.len());
    assert_eq!(r.faces_classified(FaceClass::Inside).count(), 3);
    assert_eq!(r.faces_classified(FaceClass::Outside).count(), 6);
    assert_eq!(r.faces_classified(FaceClass::OnAligned).count(), 0);
    assert_eq!(r.faces_classified(FaceClass::OnOpposed).count(), 0);
}

fn block(lo: [f64; 3], hi: [f64; 3]) -> BRep {
    use popkern_core::geometry::{Point3, Vec3};
    BRep::prism(
        &[
            Point3::new(lo[0], lo[1], lo[2]),
            Point3::new(hi[0], lo[1], lo[2]),
            Point3::new(hi[0], hi[1], lo[2]),
            Point3::new(lo[0], hi[1], lo[2]),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        hi[2] - lo[2],
        Tolerance::default(),
    )
    .expect("block")
}

/// A block resting on a plate: the plate's top splits into the contact
/// patch, on the block's bottom face and facing it, and the rest.
#[test]
fn a_block_resting_on_a_plate_imprints_its_contact_patch() {
    let tol = Tolerance::default();
    let plate = block([0.0, 0.0, 0.0], [4.0, 4.0, 1.0]);
    let resting = block([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
    let r = crate::ops::imprint::imprint(&plate, &resting, &tol).expect("imprint");
    assert!(r.imprinted_target.validate().is_ok());
    let contact: Vec<_> = r.faces_classified(FaceClass::OnOpposed).collect();
    assert_eq!(contact.len(), 1, "{:?}", r.classes);
    assert_eq!(r.faces_classified(FaceClass::Inside).count(), 0);
    assert_eq!(r.imprinted_target.faces.len(), 7);
}
