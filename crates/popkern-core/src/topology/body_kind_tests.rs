use super::BodyKind;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn an_empty_body_is_empty() {
    assert_eq!(BRep::new(Tolerance::default()).body_kind(), BodyKind::Empty);
}

#[test]
fn a_cube_is_a_solid() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert_eq!(cube.body_kind(), BodyKind::Solid);
    assert!(cube.edges.keys().all(|e| !cube.is_free_edge(e)));
}

/// The seam holds the side face on both sides, but the cylinder closes
/// across it: not a free edge (the rule a sheet's boundary must not
/// share).
#[test]
fn a_cylinder_with_its_seam_is_a_solid() {
    let rod = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    assert_eq!(rod.body_kind(), BodyKind::Solid);
}

#[test]
fn a_cube_with_a_face_retired_is_a_sheet() {
    let mut open = BRep::unit_cube(Tolerance::default()).expect("cube");
    let face = open.faces.keys().next().expect("a face");
    let exposed = open.retire_faces_to_boundary(&[face]).expect("retires");
    assert_eq!(open.body_kind(), BodyKind::Sheet);
    for e in exposed.exposed_edges {
        assert!(open.is_free_edge(e), "{e:?}");
    }
}

#[test]
fn labels_name_the_kind_for_messages() {
    assert_eq!(BodyKind::Sheet.label(), "sheet");
    assert_eq!(BodyKind::Solid.label(), "solid");
}
