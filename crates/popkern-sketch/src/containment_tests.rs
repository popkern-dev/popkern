use super::*;
use crate::edge_graph::build_edge_graph;
use crate::loop_tracer::trace_loops;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;

/// Helper: build a rectangle loop for testing.
fn make_rectangle_loop() -> (TracedLoop, Sketch) {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(2.0, 0.0);
    let p2 = sk.add_point(2.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);
    assert_eq!(loops.len(), 1);
    (loops.into_iter().next().unwrap(), sk)
}

#[test]
fn point_inside_rectangle() {
    let (loop_, sk) = make_rectangle_loop();
    assert!(point_in_loop((1.0, 0.5), &loop_, &sk));
}

#[test]
fn point_outside_rectangle() {
    let (loop_, sk) = make_rectangle_loop();
    assert!(!point_in_loop((3.0, 0.5), &loop_, &sk));
    assert!(!point_in_loop((-1.0, 0.5), &loop_, &sk));
    assert!(!point_in_loop((1.0, 2.0), &loop_, &sk));
}

#[test]
fn point_clearly_outside_all_directions() {
    let (loop_, sk) = make_rectangle_loop();
    // Far away in all directions.
    assert!(!point_in_loop((10.0, 10.0), &loop_, &sk));
    assert!(!point_in_loop((-10.0, -10.0), &loop_, &sk));
}

#[test]
fn ray_crosses_segment_basic() {
    // Vertical segment from (1, -1) to (1, 1), ray from origin in +X.
    assert!(ray_crosses_segment((0.0, 0.0), (1.0, -1.0), (1.0, 1.0)));
    // Segment entirely above: (1, 1) to (1, 2).
    assert!(!ray_crosses_segment((0.0, 0.0), (1.0, 1.0), (1.0, 2.0)));
    // Segment to the left: (-1, -1) to (-1, 1).
    assert!(!ray_crosses_segment((0.0, 0.0), (-1.0, -1.0), (-1.0, 1.0)));
}

#[test]
fn triangle_containment() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(4.0, 0.0);
    let p2 = sk.add_point(2.0, 3.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);
    assert_eq!(loops.len(), 1);
    let loop_ = &loops[0];

    // Centroid should be inside.
    assert!(point_in_loop((2.0, 1.0), loop_, &sk));
    // Point far outside.
    assert!(!point_in_loop((10.0, 10.0), loop_, &sk));
}
