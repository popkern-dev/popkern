use super::*;

#[test]
fn rectangle_four_lines() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.nodes.len(), 4);
    assert_eq!(graph.edges.len(), 4);
}

#[test]
fn zero_length_line_filtered() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(0.0, 0.0); // same position
    let _l = sk.add_line(p0, p1);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert_eq!(graph.edges.len(), 0);
    assert_eq!(warnings.len(), 1);
    assert!(matches!(warnings[0], SketchWarning::DegenerateEdge));
}

#[test]
fn two_separate_shapes() {
    let mut sk = Sketch::new();
    // Shape 1: triangle
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(0.5, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    // Shape 2: separate line segment
    let p3 = sk.add_point(10.0, 10.0);
    let p4 = sk.add_point(11.0, 10.0);
    let _l3 = sk.add_line(p3, p4);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.edges.len(), 4); // 3 triangle + 1 separate
    assert_eq!(graph.nodes.len(), 5); // 3 triangle + 2 separate
}

#[test]
fn construction_lines_excluded() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);

    // Mark first line as construction
    sk.set_construction(ElementHandle::Line(l0), true);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.edges.len(), 1);
}

#[test]
fn arc_endpoints_computed_correctly() {
    use std::f64::consts::{FRAC_PI_2, PI};

    let mut sk = Sketch::new();
    // Arc centered at origin, radius 1, from 0 to PI/2
    let _arc = sk.add_arc(0.0, 0.0, 1.0, 0.0, FRAC_PI_2);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.edges.len(), 1);

    let edge = &graph.edges[0];
    // Start at angle 0: (1, 0)
    assert!((edge.start_pos.0 - 1.0).abs() < 1e-10);
    assert!((edge.start_pos.1 - 0.0).abs() < 1e-10);
    // End at angle PI/2: (0, 1)
    assert!((edge.end_pos.0 - 0.0).abs() < 1e-10);
    assert!((edge.end_pos.1 - 1.0).abs() < 1e-10);

    // Full circle arc (start == end) is split into two half-arcs.
    let mut sk2 = Sketch::new();
    let _arc2 = sk2.add_arc(0.0, 0.0, 1.0, 0.0, 2.0 * PI);
    let (graph2, warnings2) = build_edge_graph(&sk2, &tol);
    // Full circle should produce 2 edges (split at the midpoint)
    // and 2 nodes (start/end and midpoint).
    assert!(
        warnings2.is_empty(),
        "full circle should not produce warnings"
    );
    assert_eq!(
        graph2.edges.len(),
        2,
        "full circle should split into 2 edges"
    );
    assert_eq!(graph2.nodes.len(), 2, "full circle should have 2 nodes");
}

#[test]
fn full_circle_produces_two_edges_and_two_nodes() {
    use std::f64::consts::PI;
    let mut sk = Sketch::new();
    let _arc = sk.add_arc(0.0, 0.0, 1.0, 0.0, 2.0 * PI);

    let tol = Tolerance::default();
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.edges.len(), 2, "full circle split into 2 edges");
    assert_eq!(graph.nodes.len(), 2, "full circle has 2 nodes");

    // The two edges should form a cycle: edge0.end == edge1.start,
    // edge1.end == edge0.start.
    let e0 = &graph.edges[0];
    let e1 = &graph.edges[1];
    assert_eq!(e0.end_node, e1.start_node);
    assert_eq!(e1.end_node, e0.start_node);
}

#[test]
fn node_merging_within_tolerance() {
    let mut sk = Sketch::new();
    // Two lines sharing an endpoint, but with slight numerical noise
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0 + 1e-9, 0.0); // very close to p1
    let p3 = sk.add_point(2.0, 0.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p2, p3);

    let tol = Tolerance::default(); // linear = 1e-6
    let (graph, warnings) = build_edge_graph(&sk, &tol);

    assert!(warnings.is_empty());
    assert_eq!(graph.edges.len(), 2);
    // p1 and p2 should merge into the same node
    assert_eq!(graph.nodes.len(), 3);
    // The two edges should share a node
    assert_eq!(graph.edges[0].end_node, graph.edges[1].start_node);
}
