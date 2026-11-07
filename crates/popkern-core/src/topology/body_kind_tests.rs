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

/// A cube with a face retired is a solid opened mid-edit (its faces
/// still bound the solid region) until it settles into the void, when
/// it is a sheet.
#[test]
fn a_cube_with_a_face_retired_is_a_sheet_once_settled() {
    let mut open = BRep::unit_cube(Tolerance::default()).expect("cube");
    let face = open.faces.keys().next().expect("a face");
    let exposed = open.retire_faces_to_boundary(&[face]).expect("retires");
    assert_eq!(open.body_kind(), BodyKind::General);
    open.settle_open_shells();
    assert_eq!(open.body_kind(), BodyKind::Sheet);
    for e in exposed.exposed_edges {
        assert!(open.is_free_edge(e), "{e:?}");
    }
}

fn cube_at(offset: f64) -> BRep {
    use crate::geometry::transform::Transform3;
    BRep::unit_cube(Tolerance::default())
        .and_then(|c| c.transform(&Transform3::translation(offset, offset, offset)))
        .expect("cube")
}

/// Two lumps: two disjoint cubes are a solid; two sharing a corner are
/// general (two fans there); a cube beside a sheet is general.
#[test]
fn lumps_are_one_kind_or_general() {
    let mut pair = cube_at(0.0);
    pair.append_body(&cube_at(3.0));
    assert_eq!(pair.body_kind(), BodyKind::Solid);

    let mut touching = cube_at(0.0);
    touching.append_body(&cube_at(1.0));
    let corner: Vec<_> = touching
        .vertices
        .iter()
        .filter(|(_, v)| {
            v.position
                .distance_to(&crate::geometry::Point3::new(1.0, 1.0, 1.0))
                < 1e-12
        })
        .map(|(h, _)| h)
        .collect();
    let (keep, gone) = (corner[0], corner[1]);
    for fin in touching.half_edges.values_mut() {
        if fin.origin == gone {
            fin.origin = keep;
        }
    }
    touching.vertices.remove(gone);
    assert_eq!(touching.body_kind(), BodyKind::General);

    let mut beside = cube_at(0.0);
    let mut sheet = cube_at(3.0);
    let face = sheet.faces.keys().next().expect("a face");
    sheet.retire_faces_to_boundary(&[face]).expect("retires");
    sheet.settle_open_shells();
    beside.append_body(&sheet);
    assert_eq!(beside.body_kind(), BodyKind::General);
}

/// A hollow solid (a cavity's faces between the solid and a bounded
/// void) is a solid.
#[test]
fn a_hollow_cube_is_a_solid() {
    use crate::geometry::transform::Transform3;
    let mut host = cube_at(0.0);
    let mut tool = BRep::unit_cube(Tolerance::default())
        .and_then(|t| t.transform(&Transform3::scale_uniform(0.5)))
        .and_then(|t| t.transform(&Transform3::translation(0.25, 0.25, 0.25)))
        .expect("tool");
    tool.reverse_body_orientation();
    let solid = host.solid_regions()[0];
    let map = host.append_body(&tool);
    host.invert_regions_into(&[map.regions[&tool.solid_regions()[0]]], solid);
    assert_eq!(host.body_kind(), BodyKind::Solid);
}

/// Isolated vertices are an acorn; a wire branching at a vertex is
/// general.
#[test]
fn acorns_and_branching_wires() {
    use crate::constructors::BRepBuilder;
    use crate::geometry::Point3;
    let mut acorn = BRepBuilder::new(Tolerance::default());
    acorn.add_vertex(Point3::new(0.0, 0.0, 0.0));
    assert_eq!(acorn.build_unchecked().body_kind(), BodyKind::Acorn);

    let mut b = BRepBuilder::new(Tolerance::default());
    let hub = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    for x in [1.0, -1.0, 2.0] {
        let tip = b.add_vertex(Point3::new(x, x * x, 0.0));
        b.add_edge(hub, tip);
    }
    assert_eq!(b.build_unchecked().body_kind(), BodyKind::General);
}

/// The entry gate passes an accepted kind and names the operation, the
/// kind found and the kinds accepted otherwise.
#[test]
fn require_kind_gates_an_operation() {
    use crate::error::CoreError;
    let cube = cube_at(0.0);
    assert_eq!(
        cube.require_kind("fillet", &[BodyKind::Solid])
            .expect("a solid"),
        BodyKind::Solid
    );
    let err = cube
        .require_kind("thicken", &[BodyKind::Sheet])
        .expect_err("not a sheet");
    assert!(matches!(
        err,
        CoreError::UnsupportedBodyKind {
            operation: "thicken",
            found: BodyKind::Solid,
            ..
        }
    ));
    assert_eq!(err.to_string(), "thicken needs sheets; the body is a solid");
}

#[test]
fn labels_name_the_kind_for_messages() {
    assert_eq!(BodyKind::Sheet.label(), "sheet");
    assert_eq!(BodyKind::Solid.label(), "solid");
    assert_eq!(
        BodyKind::list_label(&[BodyKind::Solid, BodyKind::Sheet]),
        "solids or sheets"
    );
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
