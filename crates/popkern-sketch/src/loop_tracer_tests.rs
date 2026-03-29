use super::*;
use crate::edge_graph::build_edge_graph;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;

#[test]
fn rectangle_single_ccw_loop() {
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
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert_eq!(
        loops.len(),
        1,
        "expected 1 loop (rectangle), got {}",
        loops.len()
    );
    assert!(
        loops[0].signed_area > 0.0,
        "rectangle loop should be CCW (positive area)"
    );
    assert!(
        (loops[0].signed_area - 1.0).abs() < 1e-10,
        "rectangle area should be 1.0"
    );
}

#[test]
fn two_separate_triangles() {
    let mut sk = Sketch::new();
    // Triangle 1
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(0.5, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    // Triangle 2 (separated)
    let p3 = sk.add_point(5.0, 0.0);
    let p4 = sk.add_point(6.0, 0.0);
    let p5 = sk.add_point(5.5, 1.0);
    let _l3 = sk.add_line(p3, p4);
    let _l4 = sk.add_line(p4, p5);
    let _l5 = sk.add_line(p5, p3);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert_eq!(
        loops.len(),
        2,
        "expected 2 loops (two triangles), got {}",
        loops.len()
    );
    for l in &loops {
        assert!(l.signed_area > 0.0, "triangle loop should be CCW");
    }
}

#[test]
fn rectangle_with_chord_two_loops() {
    let mut sk = Sketch::new();
    // Rectangle with a vertical chord splitting it.
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(2.0, 0.0);
    let p2 = sk.add_point(2.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    let p4 = sk.add_point(1.0, 0.0);
    let p5 = sk.add_point(1.0, 1.0);

    let _l0 = sk.add_line(p0, p4);
    let _l1 = sk.add_line(p4, p1);
    let _l2 = sk.add_line(p1, p2);
    let _l3 = sk.add_line(p2, p5);
    let _l4 = sk.add_line(p5, p3);
    let _l5 = sk.add_line(p3, p0);
    let _chord = sk.add_line(p4, p5);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert_eq!(
        loops.len(),
        2,
        "expected 2 loops (split rectangle), got {}",
        loops.len()
    );
    let total_area: f64 = loops.iter().map(|l| l.signed_area).sum();
    assert!(
        (total_area - 2.0).abs() < 1e-10,
        "total area should be 2.0, got {}",
        total_area
    );
}

#[test]
fn empty_graph_no_loops() {
    let sk = Sketch::new();
    let graph = EdgeGraph {
        nodes: Vec::new(),
        edges: Vec::new(),
    };
    let loops = trace_loops(&graph, &sk);
    assert!(loops.is_empty());
}

#[test]
fn single_edge_no_loop() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let _l = sk.add_line(p0, p1);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert!(loops.is_empty(), "single edge should produce no loops");
}

#[test]
fn triangle_area_correct() {
    let mut sk = Sketch::new();
    // Right triangle: (0,0), (1,0), (0,1) -- area = 0.5
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(0.0, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert_eq!(loops.len(), 1);
    assert!(
        (loops[0].signed_area - 0.5).abs() < 1e-10,
        "right triangle area should be 0.5, got {}",
        loops[0].signed_area
    );
}
