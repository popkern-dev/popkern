//! Census and rule tests for one Catmull–Clark refinement step.

use super::net::ControlNet;
use super::types::PolyhedralSplineError;
use crate::geometry::Point3;

fn cube() -> ControlNet {
    let vertices = vec![
        Point3::new(-1.0, -1.0, -1.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(1.0, 1.0, -1.0),
        Point3::new(-1.0, 1.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(1.0, -1.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(-1.0, 1.0, 1.0),
    ];
    let faces = vec![
        vec![0, 3, 2, 1],
        vec![4, 5, 6, 7],
        vec![0, 1, 5, 4],
        vec![1, 2, 6, 5],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
    ];
    ControlNet::new(vertices, faces).expect("cube cage")
}

fn valences(net: &ControlNet) -> Vec<usize> {
    let mut counts = vec![std::collections::BTreeSet::new(); net.vertices.len()];
    for face in &net.faces {
        let n = face.vertices.len();
        for corner in 0..n {
            let a = face.vertices[corner];
            let b = face.vertices[(corner + 1) % n];
            counts[a].insert(b);
            counts[b].insert(a);
        }
    }
    counts.into_iter().map(|set| set.len()).collect()
}

#[test]
fn cube_refines_to_all_quads_with_expected_census() {
    let refined = cube().catmull_clark_refined().expect("refine");
    // V + F + E = 8 + 6 + 12 = 26 vertices; sum of arities = 24 quads.
    assert_eq!(refined.vertices.len(), 26);
    assert_eq!(refined.faces.len(), 24);
    assert!(refined.faces.iter().all(|f| f.vertices.len() == 4));
    let vals = valences(&refined);
    // Original corners keep valence 3; face points valence 4; edge
    // points valence 4.
    assert!(vals[..8].iter().all(|&v| v == 3), "corners keep valence 3");
    assert!(vals[8..].iter().all(|&v| v == 4), "new points valence 4");
    // The refined net must itself be a manifold cage.
    assert!(refined.halfedges().is_ok());
}

#[test]
fn refined_cube_is_closed_and_shrinks_toward_the_limit() {
    let refined = cube().catmull_clark_refined().expect("refine");
    // Every undirected edge appears in exactly two faces (closed).
    let mut counts = std::collections::HashMap::<(usize, usize), usize>::new();
    for face in &refined.faces {
        let n = face.vertices.len();
        for corner in 0..n {
            let a = face.vertices[corner];
            let b = face.vertices[(corner + 1) % n];
            *counts.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    assert!(counts.values().all(|&c| c == 2), "refined cage is closed");
    // CC pulls the cage toward the limit surface: refined vertex
    // points sit strictly inside the original cube.
    for v in &refined.vertices[..8] {
        assert!(v.point.x.abs() < 1.0 && v.point.y.abs() < 1.0 && v.point.z.abs() < 1.0);
    }
    // Symmetry: corner points contract uniformly.
    let first = refined.vertices[0].point;
    for v in &refined.vertices[..8] {
        assert!((v.point.x.abs() - first.x.abs()).abs() < 1e-12);
    }
}

#[test]
fn ngon_face_point_keeps_ngon_valence() {
    // Triangular prism: two triangles + three quads. Face points of
    // the triangles must have valence 3, of the quads valence 4;
    // original vertices (valence 3) keep valence 3.
    let vertices = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(0.5, 1.0, 1.0),
    ];
    let faces = vec![
        vec![2, 1, 0],
        vec![3, 4, 5],
        vec![0, 1, 4, 3],
        vec![1, 2, 5, 4],
        vec![2, 0, 3, 5],
    ];
    let net = ControlNet::new(vertices, faces).expect("prism cage");
    let refined = net.catmull_clark_refined().expect("refine");
    // V + F + E = 6 + 5 + 9 = 20; faces = arity sum = 18 quads.
    assert_eq!(refined.vertices.len(), 20);
    assert_eq!(refined.faces.len(), 18);
    let vals = valences(&refined);
    assert!(vals[..6].iter().all(|&v| v == 3));
    assert_eq!(vals[6], 3, "triangle face point valence 3");
    assert_eq!(vals[7], 3, "triangle face point valence 3");
    assert!(vals[8..11].iter().all(|&v| v == 4), "quad face points");
    assert!(vals[11..].iter().all(|&v| v == 4), "edge points");
}

#[test]
fn open_strip_uses_boundary_rules() {
    // A single quad: all vertices and edges are on the boundary.
    let net = ControlNet::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![vec![0, 1, 2, 3]],
    )
    .expect("quad");
    let refined = net.catmull_clark_refined().expect("refine");
    assert_eq!(refined.vertices.len(), 4 + 1 + 4);
    assert_eq!(refined.faces.len(), 4);
    // Boundary edge points are midpoints.
    let edge_point = refined.vertices[5].point;
    assert!((edge_point.z).abs() < 1e-15);
    // Boundary vertex rule (prev + 6s + next)/8 keeps the corner in
    // plane and pulls it toward its rim neighbours.
    let corner = refined.vertices[0].point;
    assert!((corner.x - 0.125).abs() < 1e-15 && (corner.y - 0.125).abs() < 1e-15);
}

#[test]
fn boundary_pinch_declines_naming_the_vertex() {
    // Two quads sharing only vertex 2: four boundary edges meet there.
    let net = ControlNet::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
        ],
        vec![vec![0, 1, 2, 3], vec![2, 4, 5, 6]],
    )
    .expect("pinched cage");
    match net.catmull_clark_refined() {
        Err(PolyhedralSplineError::BoundaryPinch {
            vertex,
            boundary_edges,
        }) => {
            assert_eq!(vertex, 2);
            assert_eq!(boundary_edges, 4);
        }
        other => panic!("expected BoundaryPinch, got {other:?}"),
    }
}

#[test]
fn non_manifold_cage_declines_before_refining() {
    // Three faces sharing the edge 0-1: the third reuses a directed
    // edge, which the halfedge construction reports by name.
    let net = ControlNet::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
            Point3::new(0.5, 0.0, 1.0),
        ],
        vec![vec![0, 1, 2], vec![1, 0, 3], vec![0, 1, 4]],
    )
    .expect("indices are valid even when non-manifold");
    match net.catmull_clark_refined() {
        Err(PolyhedralSplineError::DuplicateDirectedEdge { from, to }) => {
            assert_eq!((from, to), (0, 1));
        }
        other => panic!("expected DuplicateDirectedEdge, got {other:?}"),
    }
}
