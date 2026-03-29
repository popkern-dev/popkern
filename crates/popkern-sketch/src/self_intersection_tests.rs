use super::*;
use crate::edge_graph::build_edge_graph;
use crate::loop_tracer::trace_loops;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;

#[test]
fn non_crossing_rectangle_ok() {
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

    assert!(check_self_intersection(&loops, &sk, &tol).is_ok());
}

#[test]
fn crossing_lines_detected() {
    // Create an X shape: two crossing line segments.
    // The edge graph will have 5 nodes (4 endpoints + 1 intersection-ish).
    // But since we build from raw sketch edges, these are just two
    // crossing lines that form loops via the graph.
    //
    // We build a bowtie (hourglass) shape: (0,0)-(1,1)-(1,0)-(0,1)-(0,0)
    // which creates a self-intersecting path.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 1.0);
    let p2 = sk.add_point(1.0, 0.0);
    let p3 = sk.add_point(0.0, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    // The loop tracer will find loops. At least one loop should have
    // crossing non-adjacent edges.
    let result = check_self_intersection(&loops, &sk, &tol);

    // In a bowtie, the diagonals (l0 and l2) cross. They form part
    // of the traced loops. The tracer should find a 4-edge loop where
    // edges 0 and 2 are non-adjacent and cross at (0.5, 0.5).
    assert!(result.is_err(), "bowtie should have self-intersection");
    if let Err(SketchError::SelfIntersectingProfile(x, y)) = result {
        assert!((x - 0.5).abs() < 0.1, "crossing x ~ 0.5, got {x}");
        assert!((y - 0.5).abs() < 0.1, "crossing y ~ 0.5, got {y}");
    }
}

#[test]
fn adjacent_edges_not_flagged() {
    // Triangle: all edges are adjacent -- no non-adjacent pairs for n=3.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(0.5, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert!(check_self_intersection(&loops, &sk, &tol).is_ok());
}

#[test]
fn adjacency_check_wraps_around() {
    assert!(are_adjacent(0, 1, 5));
    assert!(are_adjacent(1, 0, 5));
    assert!(are_adjacent(0, 4, 5)); // wrap-around
    assert!(are_adjacent(4, 0, 5));
    assert!(!are_adjacent(0, 2, 5));
    assert!(!are_adjacent(1, 3, 5));
}

#[test]
fn segment_intersection_basic() {
    let tol = Tolerance::default();

    // Crossing at (0.5, 0.5).
    let result = segment_intersection((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0), &tol);
    assert!(result.is_some());
    let (x, y) = result.unwrap();
    assert!((x - 0.5).abs() < 1e-10);
    assert!((y - 0.5).abs() < 1e-10);

    // Parallel lines: no intersection.
    let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), &tol);
    assert!(result.is_none());

    // T-junction (endpoint touch): should not be flagged as interior crossing.
    let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (0.5, 0.0), (0.5, 1.0), &tol);
    // The t parameter for the horizontal segment is 0.5 (interior),
    // but u for the vertical segment is 0.0 (endpoint) -- should be None.
    assert!(result.is_none());
}
