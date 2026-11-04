use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BodyKind;

#[test]
fn a_two_face_quad_becomes_a_one_face_sheet() {
    let tol = Tolerance::default();
    let quad = BRep::sheet_quad(-1.0, 2.0, -1.0, 2.0, 0.5, tol).expect("quad");
    assert_eq!(
        quad.body_kind(),
        BodyKind::Solid,
        "the two-face form reads as closed"
    );
    let sheet = single_face_sheet(&quad).expect("converted");
    assert_eq!(sheet.body_kind(), BodyKind::Sheet);
    assert_eq!((sheet.faces.len(), sheet.edges.len()), (1, 4));
    let face = sheet.faces.values().next().expect("face");
    assert_eq!(face.inner_loops.len(), 1);
    assert!(sheet.half_edges.values().all(|h| h.pcurve.is_some()));
    assert!(single_face_sheet(&sheet).is_none(), "already one face");
    assert!(single_face_sheet(&BRep::unit_cube(tol).expect("cube")).is_none());
}
