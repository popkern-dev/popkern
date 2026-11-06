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

/// A wire edge keeps two mated fins in no loop: a wire edge, neither
/// free nor a face use, with its vertices read through the mate.
#[test]
fn a_polyline_is_a_wire_of_wire_edges() {
    use crate::constructors::BRepBuilder;
    use crate::geometry::Point3;
    let mut b = BRepBuilder::new(Tolerance::default());
    let v: Vec<_> = [0.0, 1.0, 2.0]
        .iter()
        .map(|&x| b.add_vertex(Point3::new(x, 0.0, 0.0)))
        .collect();
    let (first, _, _) = b.add_edge(v[0], v[1]);
    b.add_edge(v[1], v[2]);
    let wire = b.build_unchecked();
    assert_eq!(wire.body_kind(), BodyKind::Wire);
    assert_eq!(wire.wire_edges().len(), 2);
    assert!(wire.free_edges().is_empty());
    assert_eq!(wire.edge_fin_count(first), 0);
    assert_eq!(wire.edge_vertices(first), Some((v[0], v[1])));
}

/// A solid's and a sheet's edges all bound faces: none is a wire edge.
#[test]
fn face_edges_are_not_wire_edges() {
    let mut open = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert!(open.wire_edges().is_empty());
    let face = open.faces.keys().next().expect("a face");
    open.retire_faces_to_boundary(&[face]).expect("retires");
    assert!(open.wire_edges().is_empty());
}
