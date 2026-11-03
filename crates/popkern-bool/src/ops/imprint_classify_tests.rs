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
