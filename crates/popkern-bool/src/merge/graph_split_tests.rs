use super::*;
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
use popkern_core::geometry::Point3;

fn simple_edge_graph() -> IntersectionGraph {
    let mut g = IntersectionGraph::new();
    let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = g.add_node(Point3::new(2.0, 0.0, 0.0), NodeSource::Intersection);
    let edge_h = popkern_core::topology::EdgeHandle::default();
    g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);
    g
}

#[test]
fn classify_finds_vertex() {
    let g = simple_edge_graph();
    match classify_endpoint(&g, &Point3::new(0.0, 0.0, 0.0), 1e-6) {
        EndpointClass::OnVertex(idx) => assert_eq!(idx, 0),
        other => panic!("Expected OnVertex, got {:?}", other),
    }
}

#[test]
fn classify_finds_edge_midpoint() {
    let g = simple_edge_graph();
    match classify_endpoint(&g, &Point3::new(1.0, 0.0, 0.0), 1e-6) {
        EndpointClass::OnEdge(he_idx, t) => {
            assert_eq!(he_idx, 0);
            assert!((t - 0.5).abs() < 1e-6, "t should be 0.5, got {t}");
        }
        other => panic!("Expected OnEdge, got {:?}", other),
    }
}

#[test]
fn classify_interior_point() {
    let g = simple_edge_graph();
    match classify_endpoint(&g, &Point3::new(0.0, 5.0, 0.0), 1e-6) {
        EndpointClass::OnFaceInterior => {}
        other => panic!("Expected OnFaceInterior, got {:?}", other),
    }
}

#[test]
fn split_edge_creates_new_node_and_pairs() {
    let mut g = simple_edge_graph();
    let face_ref = FaceRef {
        operand: Operand::A,
        face: popkern_core::topology::FaceHandle::default(),
    };
    // Assign face to the forward half-edge.
    g.set_half_edge_face(0, face_ref);

    let new_node = split_graph_edge(&mut g, 0, Point3::new(1.0, 0.0, 0.0), 0.5);

    // Should have 3 nodes now.
    assert_eq!(g.nodes.len(), 3);
    assert_eq!(new_node, 2);

    // Should have 4 half-edges (original pair + new pair).
    assert_eq!(g.half_edges.len(), 4);

    // Original pair: he[0] start(0) → new_node(2), he[1] new_node(2) → start(0).
    assert_eq!(g.half_edges[0].origin, 0);
    assert_eq!(g.half_edges[1].origin, 2); // twin now starts at new_node
    assert_eq!(g.half_edges[0].twin, 1);
    assert_eq!(g.half_edges[1].twin, 0);

    // New pair: he[2] new_node(2) → end(1), he[3] end(1) → new_node(2).
    assert_eq!(g.half_edges[2].origin, 2);
    assert_eq!(g.half_edges[3].origin, 1);
    assert_eq!(g.half_edges[2].twin, 3);
    assert_eq!(g.half_edges[3].twin, 2);

    // Face preserved on forward half-edges.
    assert_eq!(g.half_edges[0].face, Some(face_ref));
    assert_eq!(g.half_edges[2].face, Some(face_ref));

    // face_edges should contain both forward half-edges for that face.
    let face_hes = &g.face_edges[&face_ref];
    assert!(face_hes.contains(&0));
    assert!(face_hes.contains(&2));
}

// ── Additional classify_endpoint edge-case tests ───────────────────

#[test]
fn classify_finds_second_vertex() {
    let g = simple_edge_graph();
    // The second node is at (2, 0, 0).
    match classify_endpoint(&g, &Point3::new(2.0, 0.0, 0.0), 1e-6) {
        EndpointClass::OnVertex(idx) => assert_eq!(idx, 1),
        other => panic!("Expected OnVertex(1), got {:?}", other),
    }
}

#[test]
fn classify_vertex_within_tolerance() {
    let g = simple_edge_graph();
    // Point very close to vertex 0 but not exact.
    let near_v0 = Point3::new(1e-8, 0.0, 0.0);
    match classify_endpoint(&g, &near_v0, 1e-6) {
        EndpointClass::OnVertex(idx) => assert_eq!(idx, 0),
        other => panic!("Expected OnVertex(0), got {:?}", other),
    }
}

#[test]
fn classify_on_edge_quarter_point() {
    let g = simple_edge_graph();
    // Point at t=0.25 along the edge from (0,0,0) to (2,0,0).
    match classify_endpoint(&g, &Point3::new(0.5, 0.0, 0.0), 1e-6) {
        EndpointClass::OnEdge(he_idx, t) => {
            assert_eq!(he_idx, 0);
            assert!((t - 0.25).abs() < 1e-6, "t should be 0.25, got {t}");
        }
        other => panic!("Expected OnEdge, got {:?}", other),
    }
}

#[test]
fn classify_point_near_edge_but_off_axis() {
    let g = simple_edge_graph();
    // Point above the edge midpoint but within tolerance.
    let p = Point3::new(1.0, 1e-8, 0.0);
    match classify_endpoint(&g, &p, 1e-6) {
        EndpointClass::OnEdge(he_idx, t) => {
            assert_eq!(he_idx, 0);
            assert!((t - 0.5).abs() < 1e-6, "t should be ~0.5, got {t}");
        }
        other => panic!("Expected OnEdge, got {:?}", other),
    }
}

#[test]
fn classify_point_too_far_from_edge_is_interior() {
    let g = simple_edge_graph();
    // Point above the edge midpoint but outside tolerance.
    let p = Point3::new(1.0, 0.1, 0.0);
    match classify_endpoint(&g, &p, 1e-6) {
        EndpointClass::OnFaceInterior => {}
        other => panic!("Expected OnFaceInterior, got {:?}", other),
    }
}

#[test]
fn classify_empty_graph_returns_interior() {
    let g = IntersectionGraph::new();
    match classify_endpoint(&g, &Point3::new(0.0, 0.0, 0.0), 1e-6) {
        EndpointClass::OnFaceInterior => {}
        other => panic!("Expected OnFaceInterior, got {:?}", other),
    }
}

#[test]
fn classify_with_multiple_edges() {
    // Build a graph with two separate edges and verify classification
    // picks the correct one.
    let mut g = IntersectionGraph::new();
    let edge_h = popkern_core::topology::EdgeHandle::default();

    // Edge 1: (0,0,0) → (2,0,0)
    let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = g.add_node(Point3::new(2.0, 0.0, 0.0), NodeSource::Intersection);
    g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

    // Edge 2: (0,5,0) → (2,5,0)
    let n2 = g.add_node(Point3::new(0.0, 5.0, 0.0), NodeSource::Intersection);
    let n3 = g.add_node(Point3::new(2.0, 5.0, 0.0), NodeSource::Intersection);
    g.add_half_edge_pair(n2, n3, EdgeSource::OriginalA(edge_h), None);

    // Point on the second edge midpoint.
    match classify_endpoint(&g, &Point3::new(1.0, 5.0, 0.0), 1e-6) {
        EndpointClass::OnEdge(he_idx, t) => {
            assert_eq!(he_idx, 2, "should be on second edge pair");
            assert!((t - 0.5).abs() < 1e-6, "t should be 0.5, got {t}");
        }
        other => panic!("Expected OnEdge, got {:?}", other),
    }
}

#[test]
fn classify_near_endpoint_excluded_from_on_edge() {
    let g = simple_edge_graph();
    // A point very close to node 0 along the edge, within tolerance of
    // the endpoint exclusion zone -> should be caught as vertex instead.
    let near_start = Point3::new(1e-8, 0.0, 0.0);
    match classify_endpoint(&g, &near_start, 1e-6) {
        EndpointClass::OnVertex(idx) => assert_eq!(idx, 0),
        other => panic!(
            "Expected OnVertex (near-endpoint should match vertex), got {:?}",
            other
        ),
    }
}

// ── point_on_segment edge cases ────────────────────────────────────

#[test]
fn point_on_segment_degenerate_edge_returns_none() {
    // Degenerate edge: both endpoints the same.
    let p0 = Point3::new(1.0, 2.0, 3.0);
    let result = point_on_segment(&Point3::new(1.0, 2.0, 3.0), &p0, &p0, 1e-6);
    assert!(result.is_none(), "degenerate edge should return None");
}

#[test]
fn point_on_segment_past_end_returns_none() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(2.0, 0.0, 0.0);
    // Beyond the end of the segment.
    let result = point_on_segment(&Point3::new(3.0, 0.0, 0.0), &p0, &p1, 1e-6);
    assert!(result.is_none(), "point past end should return None");
}

#[test]
fn split_graph_edge_source_propagation() {
    // Verify that SplitA/SplitB sources are correctly derived.
    let mut g = IntersectionGraph::new();
    let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = g.add_node(Point3::new(4.0, 0.0, 0.0), NodeSource::Intersection);
    let edge_h = popkern_core::topology::EdgeHandle::default();
    g.add_half_edge_pair(n0, n1, EdgeSource::OriginalB(edge_h), None);

    let _new_node = split_graph_edge(&mut g, 0, Point3::new(2.0, 0.0, 0.0), 0.5);

    // New half-edges (indices 2,3) should have SplitB source.
    assert!(
        matches!(g.half_edges[2].source, EdgeSource::SplitB { .. }),
        "expected SplitB, got {:?}",
        g.half_edges[2].source
    );
    assert!(
        matches!(g.half_edges[3].source, EdgeSource::SplitB { .. }),
        "expected SplitB, got {:?}",
        g.half_edges[3].source
    );
}

#[test]
fn split_workhorse_resolves_the_correct_turn_on_a_coiled_curve() {
    // PK-451: the split-window/adoption workhorse kept a fixed
    // 32-sample coarse scan — the density class whose 16-sample
    // variant lost the seam pierce to the ADJACENT COIL TURN at the
    // classification site (PK-449). MEASURED boundary: 32 samples
    // survive every phase of the 6.5-turn coil seam (4.9/turn,
    // 4096-phase sweep clean) but lose phases on a 13-turn coil
    // (2.5/turn — this pin fails under the fixed count). The coarse
    // count now adapts to the NURBS knot-span count (floor 32, so
    // few-span windows scan bit-identically to the historical body);
    // a
    // point ON a many-turn curve must invert to its own turn at
    // EVERY phase, not just the lucky ones.
    use popkern_core::geometry::CurveOnEdge;
    use popkern_core::geometry::curve::{CurveEval as _, CurveGeometry};
    let turns = 13.0;
    let n_pts = 524;
    let pts: Vec<Point3> = (0..n_pts)
        .map(|i| {
            let s = turns * i as f64 / (n_pts - 1) as f64;
            let a = std::f64::consts::TAU * s;
            Point3::new(0.26 * a.cos(), 0.26 * a.sin(), 0.15 * s)
        })
        .collect();
    let params: Vec<f64> = (0..n_pts).map(|i| i as f64 / (n_pts - 1) as f64).collect();
    let nurbs = popkern_core::geometry::fitting::interpolate_curve(&pts, &params, 3)
        .expect("helix interpolation");
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    };
    for k in 0..64 {
        let t_probe = 0.02 + 0.96 * k as f64 / 63.0;
        let probe = curve.geometry.point_at(t_probe);
        let (t, d) = nearest_param_on_curve(&curve, &probe);
        assert!(
            d < 1e-9 && (t - t_probe).abs() < 1e-3,
            "phase {t_probe:.4}: inverted to t {t:.4} d {d:.3e} (adjacent-turn basin)"
        );
    }
}
