//! Tests for §5.3 vertex classifier.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// A cube corner with all 3 incident edges filleted: all convex,
/// uniform convexity — `SameSideConvex`.
#[test]
fn cube_corner_three_convex_edges_is_same_side_convex() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let report = classify_vertex_blend(&cube, vertex, &incident, &tol());
    assert_eq!(
        report.classification,
        VertexBlendClassification::SameSideConvex
    );
    assert_eq!(report.n_incident_filleted_edges, 3);
    assert_eq!(report.n_convex, 3);
    assert_eq!(report.n_concave, 0);
    assert_eq!(report.n_tangent, 0);
    assert!(!report.has_boundary_edge);
}

/// Vertex with no filleted incident edges → degenerate default
/// (SameSideConvex, zero counts). Callers should check
/// `n_incident_filleted_edges == 0`.
#[test]
fn no_incident_filleted_edges_returns_degenerate_default() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let report = classify_vertex_blend(&cube, vertex, &[], &tol());
    assert_eq!(report.n_incident_filleted_edges, 0);
    assert_eq!(
        report.classification,
        VertexBlendClassification::SameSideConvex
    );
}

/// Filtering: pass edges that aren't incident to the vertex. They
/// should be silently skipped.
#[test]
fn non_incident_edges_are_filtered() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let opposite = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 1.0 && v.position.y == 1.0 && v.position.z == 1.0)
        .map(|(h, _)| h)
        .unwrap();
    // Get an edge incident to the OPPOSITE vertex (not to `vertex`).
    let non_incident_edge = cube
        .half_edges
        .values()
        .find(|he| he.origin == opposite)
        .map(|he| he.edge)
        .unwrap();
    let report = classify_vertex_blend(&cube, vertex, &[non_incident_edge], &tol());
    assert_eq!(report.n_incident_filleted_edges, 0);
}

/// Open-shell boundary vertex classifies as `FreeBoundary`.
#[test]
fn open_box_boundary_vertex_is_free_boundary() {
    let ob = BRep::open_box(tol()).unwrap();
    // Find a vertex incident to at least one boundary edge.
    let (vertex, incident_boundary) = ob
        .half_edges
        .values()
        .filter(|he| ob.boundary_edges.contains(&he.edge))
        .map(|he| (he.origin, he.edge))
        .next()
        .unwrap();
    // Include ALL edges at this vertex; the boundary one triggers
    // the FreeBoundary classification regardless of other counts.
    let incident: Vec<_> = ob
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let report = classify_vertex_blend(&ob, vertex, &incident, &tol());
    assert!(report.has_boundary_edge);
    assert_eq!(
        report.classification,
        VertexBlendClassification::FreeBoundary,
        "expected FreeBoundary; report has n_convex={}, n_concave={}, \
         boundary_edge_ref={:?}",
        report.n_convex,
        report.n_concave,
        incident_boundary,
    );
}

/// Two non-boundary convex edges at a cube corner (partial fillet
/// set) still classify as SameSideConvex with N=2.
#[test]
fn partial_cube_corner_two_edges_is_same_side_convex() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let mut incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    incident.truncate(2);
    let report = classify_vertex_blend(&cube, vertex, &incident, &tol());
    assert_eq!(
        report.classification,
        VertexBlendClassification::SameSideConvex
    );
    assert_eq!(report.n_incident_filleted_edges, 2);
    assert_eq!(report.n_convex, 2);
}
