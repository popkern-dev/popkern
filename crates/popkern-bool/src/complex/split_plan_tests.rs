//! The complex splits non-manifold unions into manifold pieces meeting
//! at recorded contacts: one result edge per selected wedge, one vertex
//! per fan.

use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, TopologyContactKind};

use crate::BoolOperation;

fn cube_at(x: f64, y: f64, z: f64, tol: Tolerance) -> BRep {
    BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::translation(x, y, z))
        .expect("placed")
}

/// The union of a unit cube and a copy at `(x, y, z)`, served by the
/// complex split and accepted with its contacts recorded.
fn split_union(x: f64, y: f64, z: f64) -> BRep {
    let tol = Tolerance::default();
    let (a, b) = (cube_at(0.0, 0.0, 0.0, tol), cube_at(x, y, z, tol));
    let ix = crate::intersect::find_intersections(&a, &b, &tol).expect("intersections");
    let (body, _) = crate::complex::serve::rescue_split(&a, &b, &ix, BoolOperation::Union, &tol)
        .expect("the split is manifold")
        .expect("the complex serves the split");
    crate::postprocess::contact_record::accept_split(body, &tol).expect("the split validates")
}

fn counts(b: &BRep) -> (usize, usize, usize) {
    (b.vertices.len(), b.edges.len(), b.faces.len())
}

/// Along an edge: four kept fins pair into two result edges, and the
/// edge's end vertices split into a fan each. Every recorded contact
/// (one per touching face or vertex pair) lies on the shared line.
#[test]
fn cubes_touching_along_an_edge_split_into_two_cubes() {
    let body = split_union(1.0, 1.0, 0.0);
    assert_eq!(counts(&body), (16, 24, 12));
    let contacts = &body.topology_contacts;
    assert!(contacts.iter().any(|c| c.kind == TopologyContactKind::Edge));
    for c in contacts {
        let p = c.location;
        assert!(
            (p.x - 1.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9,
            "{c:?}"
        );
    }
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
}

/// At a corner: no edge is shared, the corner vertex splits, and every
/// contact is a vertex contact at the corner.
#[test]
fn cubes_touching_at_a_corner_split_into_two_cubes() {
    let body = split_union(1.0, 1.0, 1.0);
    assert_eq!(counts(&body), (16, 24, 12));
    let contacts = &body.topology_contacts;
    assert!(!contacts.is_empty());
    for c in contacts {
        assert_eq!(c.kind, TopologyContactKind::Vertex);
        assert!(
            c.location.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9,
            "{c:?}"
        );
    }
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
}
