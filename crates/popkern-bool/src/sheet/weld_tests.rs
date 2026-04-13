use super::*;
use popkern_core::geometry::Point3;

#[test]
fn weld_merges_two_vertices_at_same_position() {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1.0, 2.0, 3.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    assert_eq!(brep.vertices.len(), 2);
    let welds = weld_coincident_vertices(&mut brep, &tol);
    assert_eq!(welds, 1);
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn weld_leaves_distinct_vertices_untouched() {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let welds = weld_coincident_vertices(&mut brep, &tol);
    assert_eq!(welds, 0);
    assert_eq!(brep.vertices.len(), 2);
}

#[test]
fn weld_transitive_merge_via_middle_vertex() {
    // Three vertices: A at (0,0,0), B at (1e-7, 0, 0), C at (2e-7, 0, 0).
    // With default tol.linear() == 1e-6, A-B is within tol, B-C is
    // within tol, A-C is within tol. All three collapse into one.
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1e-7, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(2e-7, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let welds = weld_coincident_vertices(&mut brep, &tol);
    assert_eq!(welds, 2);
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn weld_rewrites_half_edge_origins_to_canonical() {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    let v1 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(5.0, 5.0, 5.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(5.0, 5.0, 5.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let he1 = brep.half_edges.insert(popkern_core::topology::HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: Default::default(),
        loop_: Default::default(),
        pcurve: None,
    });
    let he2 = brep.half_edges.insert(popkern_core::topology::HalfEdge {
        origin: v2,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: Default::default(),
        loop_: Default::default(),
        pcurve: None,
    });

    weld_coincident_vertices(&mut brep, &tol);

    let o1 = brep.half_edges[he1].origin;
    let o2 = brep.half_edges[he2].origin;
    assert_eq!(o1, o2, "both half-edges should reference the same origin");
    assert!(
        brep.vertices.get(o1).is_some(),
        "canonical vertex must still exist"
    );
}

/// Build a minimal BRep with two edges sharing the same endpoint
/// pair — simulates what sheet-sheet union produces after imprint +
/// vertex weld along a transversal crossing.
fn build_duplicate_edges_fixture() -> (BRep, VertexHandle, VertexHandle, EdgeHandle, EdgeHandle) {
    use popkern_core::topology::{Edge, HalfEdge, Vertex};
    let mut brep = BRep::new(Tolerance::default());
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let edge_a = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });
    let he_a1 = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: edge_a,
        loop_: Default::default(),
        pcurve: None,
    });
    let he_a2 = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: he_a1,
        next: Default::default(),
        prev: Default::default(),
        edge: edge_a,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.half_edges[he_a1].twin = he_a2;
    brep.edges[edge_a].half_edge = he_a1;

    let edge_b = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });
    let he_b1 = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: edge_b,
        loop_: Default::default(),
        pcurve: None,
    });
    let he_b2 = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: he_b1,
        next: Default::default(),
        prev: Default::default(),
        edge: edge_b,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.half_edges[he_b1].twin = he_b2;
    brep.edges[edge_b].half_edge = he_b1;

    (brep, v1, v2, edge_a, edge_b)
}

#[test]
fn weld_edges_collapses_duplicate_edge() {
    let (mut brep, _, _, edge_a, edge_b) = build_duplicate_edges_fixture();
    assert_eq!(brep.edges.len(), 2);
    let welds = weld_coincident_edges(&mut brep);
    assert_eq!(welds, 1, "exactly one duplicate edge should be welded");
    assert_eq!(brep.edges.len(), 1);
    let survivor = if brep.edges.contains_key(edge_a) {
        edge_a
    } else {
        edge_b
    };
    for (_, he) in &brep.half_edges {
        assert_eq!(
            he.edge, survivor,
            "every half-edge should reference the canonical edge"
        );
    }
}

#[test]
fn weld_edges_populates_radial_cycle() {
    let (mut brep, _, _, _, _) = build_duplicate_edges_fixture();
    weld_coincident_edges(&mut brep);
    assert_eq!(
        brep.nonmanifold_radial.len(),
        4,
        "radial cycle should have 4 entries for the 4 half-edges"
    );
    let start = *brep.nonmanifold_radial.keys().next().unwrap();
    let mut cursor = start;
    let mut steps = 0;
    loop {
        cursor = brep.radial_next(cursor).expect("valid handle");
        steps += 1;
        if cursor == start {
            break;
        }
        assert!(steps < 10, "cycle did not close within 10 steps");
    }
    assert_eq!(steps, 4, "4 half-edges should form a cycle of length 4");
}

#[test]
fn weld_edges_leaves_distinct_endpoints_untouched() {
    use popkern_core::topology::{Edge, HalfEdge, Vertex};
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v3 = brep.vertices.insert(Vertex {
        position: Point3::new(2.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    for (a, b) in [(v1, v2), (v2, v3)] {
        let eh = brep.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
        });
        let h1 = brep.half_edges.insert(HalfEdge {
            origin: a,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: eh,
            loop_: Default::default(),
            pcurve: None,
        });
        let h2 = brep.half_edges.insert(HalfEdge {
            origin: b,
            twin: h1,
            next: Default::default(),
            prev: Default::default(),
            edge: eh,
            loop_: Default::default(),
            pcurve: None,
        });
        brep.half_edges[h1].twin = h2;
        brep.edges[eh].half_edge = h1;
    }
    let welds = weld_coincident_edges(&mut brep);
    assert_eq!(
        welds, 0,
        "edges with different endpoint pairs must not weld"
    );
    assert_eq!(brep.edges.len(), 2);
    assert!(
        brep.nonmanifold_radial.is_empty(),
        "no radial cycles should be populated when nothing welds"
    );
}
