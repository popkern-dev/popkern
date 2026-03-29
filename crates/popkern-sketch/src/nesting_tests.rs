use super::*;
use crate::edge_graph::build_edge_graph;
use crate::loop_tracer::trace_loops;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;

#[test]
fn single_rectangle_one_face_root() {
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
    let nesting = classify_nesting(&loops, &sk);

    assert_eq!(nesting.len(), 1, "single rectangle = 1 root");
    assert!(nesting[0].is_face, "root should be a face");
    assert!(nesting[0].children.is_empty(), "no children");
}

#[test]
fn rectangle_with_inner_triangle_hole() {
    let mut sk = Sketch::new();
    // Outer rectangle: 0,0 to 10,10
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(10.0, 0.0);
    let p2 = sk.add_point(10.0, 10.0);
    let p3 = sk.add_point(0.0, 10.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    // Inner triangle: inside the rectangle
    let p4 = sk.add_point(3.0, 3.0);
    let p5 = sk.add_point(7.0, 3.0);
    let p6 = sk.add_point(5.0, 7.0);
    let _l4 = sk.add_line(p4, p5);
    let _l5 = sk.add_line(p5, p6);
    let _l6 = sk.add_line(p6, p4);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);

    assert_eq!(loops.len(), 2, "should have outer rect + inner triangle");

    let nesting = classify_nesting(&loops, &sk);
    assert_eq!(nesting.len(), 1, "one root (the outer rectangle)");
    assert!(nesting[0].is_face, "root is a face");
    assert_eq!(
        nesting[0].children.len(),
        1,
        "one child (the triangle hole)"
    );
    assert!(!nesting[0].children[0].is_face, "child is a hole");
}

#[test]
fn two_disjoint_rectangles_two_roots() {
    let mut sk = Sketch::new();
    // Rectangle 1
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    // Rectangle 2 (separated)
    let p4 = sk.add_point(5.0, 0.0);
    let p5 = sk.add_point(6.0, 0.0);
    let p6 = sk.add_point(6.0, 1.0);
    let p7 = sk.add_point(5.0, 1.0);
    let _l4 = sk.add_line(p4, p5);
    let _l5 = sk.add_line(p5, p6);
    let _l6 = sk.add_line(p6, p7);
    let _l7 = sk.add_line(p7, p4);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);
    let nesting = classify_nesting(&loops, &sk);

    assert_eq!(nesting.len(), 2, "two disjoint rectangles = 2 roots");
    assert!(nesting[0].is_face);
    assert!(nesting[1].is_face);
    assert!(nesting[0].children.is_empty());
    assert!(nesting[1].children.is_empty());
}

#[test]
fn three_nested_loops() {
    // Outer, middle, inner -- rectangle > hole > island.
    let mut sk = Sketch::new();

    // Outer rectangle: 0,0 to 20,20
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(20.0, 0.0);
    let p2 = sk.add_point(20.0, 20.0);
    let p3 = sk.add_point(0.0, 20.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    // Middle rectangle (hole): 3,3 to 17,17
    let p4 = sk.add_point(3.0, 3.0);
    let p5 = sk.add_point(17.0, 3.0);
    let p6 = sk.add_point(17.0, 17.0);
    let p7 = sk.add_point(3.0, 17.0);
    let _l4 = sk.add_line(p4, p5);
    let _l5 = sk.add_line(p5, p6);
    let _l6 = sk.add_line(p6, p7);
    let _l7 = sk.add_line(p7, p4);

    // Inner rectangle (island): 7,7 to 13,13
    let p8 = sk.add_point(7.0, 7.0);
    let p9 = sk.add_point(13.0, 7.0);
    let p10 = sk.add_point(13.0, 13.0);
    let p11 = sk.add_point(7.0, 13.0);
    let _l8 = sk.add_line(p8, p9);
    let _l9 = sk.add_line(p9, p10);
    let _l10 = sk.add_line(p10, p11);
    let _l11 = sk.add_line(p11, p8);

    let tol = Tolerance::default();
    let (graph, _) = build_edge_graph(&sk, &tol);
    let loops = trace_loops(&graph, &sk);
    assert_eq!(loops.len(), 3);

    let nesting = classify_nesting(&loops, &sk);
    assert_eq!(nesting.len(), 1, "one root (outermost)");

    let root = &nesting[0];
    assert!(root.is_face, "root = face (depth 0)");
    assert_eq!(root.children.len(), 1, "one child (middle = hole)");

    let hole = &root.children[0];
    assert!(!hole.is_face, "middle = hole (depth 1)");
    assert_eq!(hole.children.len(), 1, "one grandchild (inner = island)");

    let island = &hole.children[0];
    assert!(island.is_face, "inner = face (depth 2)");
    assert!(island.children.is_empty());
}

#[test]
fn compute_test_point_rectangle() {
    use crate::loop_tracer::DirectedEdge;
    use popkern_constraint::ElementHandle;

    // Manually create a mock loop with known positions.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(2.0, 0.0);
    let lh = sk.add_line(p0, p1);

    let loop_ = TracedLoop {
        edges: vec![
            DirectedEdge {
                element: ElementHandle::Line(lh),
                start: (0.0, 0.0),
                end: (2.0, 0.0),
                reversed: false,
            },
            DirectedEdge {
                element: ElementHandle::Line(lh),
                start: (2.0, 0.0),
                end: (2.0, 2.0),
                reversed: false,
            },
            DirectedEdge {
                element: ElementHandle::Line(lh),
                start: (2.0, 2.0),
                end: (0.0, 2.0),
                reversed: false,
            },
            DirectedEdge {
                element: ElementHandle::Line(lh),
                start: (0.0, 2.0),
                end: (0.0, 0.0),
                reversed: false,
            },
        ],
        signed_area: 4.0,
    };

    let pt = compute_test_point(&loop_);
    assert!((pt.0 - 1.0).abs() < 1e-10, "centroid x = 1.0");
    assert!((pt.1 - 1.0).abs() < 1e-10, "centroid y = 1.0");
}
