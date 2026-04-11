use crate::geometry::Point3;
use crate::precision::EntityHandle;
use crate::tolerance::Tolerance;
use crate::topology::brep::BRep;
use crate::topology::{EdgeHandle, VertexHandle};

// ---------------------------------------------------------------------------
// entities_coincident
// ---------------------------------------------------------------------------

#[test]
fn entities_coincident_same_vertex() {
    let mut brep = BRep::new(Tolerance::default());
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    // Same handle — distance is 0, always coincident.
    assert!(brep.entities_coincident(EntityHandle::Vertex(va), EntityHandle::Vertex(va)));
}

#[test]
fn entities_coincident_close_vertices() {
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1e-7, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    // Distance 1e-7 < default vertex tolerance 1e-6 => coincident.
    assert!(brep.entities_coincident(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn entities_coincident_far_vertices() {
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    // Distance 1.0 >> tolerance => not coincident.
    assert!(!brep.entities_coincident(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn entities_coincident_uses_max_tolerance() {
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-2), // Wider tolerance
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(5e-3, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None, // Uses default 1e-6
        name: None,
    });
    // Distance 5e-3 < max(1e-2, 1e-6) = 1e-2 => coincident.
    assert!(brep.entities_coincident(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn entities_coincident_edge_edge_returns_false() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edges: Vec<_> = brep.edges.keys().collect();
    // Edge-edge is not yet implemented, should return false.
    assert!(!brep.entities_coincident(EntityHandle::Edge(edges[0]), EntityHandle::Edge(edges[1])));
}

#[test]
fn entities_coincident_mixed_types_returns_false() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    let eh = brep.edges.keys().next().unwrap();
    assert!(!brep.entities_coincident(EntityHandle::Vertex(vh), EntityHandle::Edge(eh)));
}

// ---------------------------------------------------------------------------
// imprecision_overlap
// ---------------------------------------------------------------------------

#[test]
fn imprecision_overlap_close_vertices() {
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-3),
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1.5e-3, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-3),
        name: None,
    });
    // Distance 1.5e-3 < 1e-3 + 1e-3 = 2e-3 => overlap.
    assert!(brep.imprecision_overlap(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn imprecision_overlap_no_overlap() {
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-3),
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-3),
        name: None,
    });
    // Distance 1.0 >> 1e-3 + 1e-3 = 2e-3 => no overlap.
    assert!(!brep.imprecision_overlap(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn imprecision_overlap_vs_coincident_boundary() {
    // Test the region between max(a,b) and a+b where overlap is true
    // but coincident is false.
    let tol = Tolerance::default().with_linear(1e-6);
    let mut brep = BRep::new(tol);
    let va = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-2),
        name: None,
    });
    let vb = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1.5e-2, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: Some(1e-2),
        name: None,
    });
    // Distance 1.5e-2:
    //   max(1e-2, 1e-2) = 1e-2 => NOT coincident (1.5e-2 >= 1e-2)
    //   1e-2 + 1e-2 = 2e-2 => overlap (1.5e-2 < 2e-2)
    assert!(!brep.entities_coincident(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
    assert!(brep.imprecision_overlap(EntityHandle::Vertex(va), EntityHandle::Vertex(vb)));
}

#[test]
fn imprecision_overlap_edge_edge_returns_false() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edges: Vec<_> = brep.edges.keys().collect();
    assert!(!brep.imprecision_overlap(EntityHandle::Edge(edges[0]), EntityHandle::Edge(edges[1])));
}

// ---------------------------------------------------------------------------
// recover_edge_consistency
// ---------------------------------------------------------------------------

#[test]
fn recover_edge_consistency_returns_current_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    let expected = brep.edge_tolerance(eh);
    let result = brep.recover_edge_consistency(eh).unwrap();
    assert!((result - expected).abs() < 1e-15);
}

#[test]
fn recover_edge_consistency_invalid_handle() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = EdgeHandle::default();
    let result = brep.recover_edge_consistency(bogus);
    assert!(result.is_err());
}

#[test]
fn recover_edge_consistency_sets_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    // Widen tolerance first.
    brep.set_entity_tolerance(EntityHandle::Edge(eh), 1e-3);
    // Recover — stub just returns the current value.
    let result = brep.recover_edge_consistency(eh).unwrap();
    assert!((result - 1e-3).abs() < 1e-15);
    // Tolerance should be set explicitly.
    assert!(brep.edges[eh].tolerance.is_some());
}

// ---------------------------------------------------------------------------
// recover_vertex_precision
// ---------------------------------------------------------------------------

#[test]
fn recover_vertex_precision_returns_current_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    let expected = brep.vertex_tolerance(vh);
    let result = brep.recover_vertex_precision(vh).unwrap();
    assert!((result - expected).abs() < 1e-15);
}

#[test]
fn recover_vertex_precision_invalid_handle() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = VertexHandle::default();
    let result = brep.recover_vertex_precision(bogus);
    assert!(result.is_err());
}

#[test]
fn recover_vertex_precision_sets_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    // Widen tolerance first.
    brep.set_entity_tolerance(EntityHandle::Vertex(vh), 1e-3);
    // Recover — stub returns current value.
    let result = brep.recover_vertex_precision(vh).unwrap();
    assert!((result - 1e-3).abs() < 1e-15);
    assert!(brep.vertices[vh].tolerance.is_some());
}
