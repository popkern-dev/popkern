//! Tests for the BRep builder module.

use super::*;
use crate::edge_graph;
use crate::loop_tracer;
use crate::nesting;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;

/// Build a simple rectangle sketch and run the full pipeline up to
/// nesting classification.
fn rectangle_pipeline() -> (Vec<TracedLoop>, Vec<NestingNode>, Sketch) {
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
    let (graph, _warnings) = edge_graph::build_edge_graph(&sk, &tol);
    let loops = loop_tracer::trace_loops(&graph);
    let roots = nesting::classify_nesting(&loops, &sk);

    (loops, roots, sk)
}

#[test]
fn rectangle_builds_valid_brep() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, regions) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    // Should have 1 region.
    assert_eq!(regions.len(), 1, "expected 1 face region");

    // The region should reference a valid face.
    let region = &regions[0];
    assert!(
        brep.faces.contains_key(region.face),
        "face handle should be valid"
    );

    // Should have 4 vertices, 4 edges, 2 faces (front + back), 2 loops.
    assert_eq!(brep.vertices.len(), 4);
    assert_eq!(brep.edges.len(), 4);
    assert_eq!(brep.faces.len(), 2);
    assert_eq!(brep.loops.len(), 2);
}

#[test]
fn rectangle_edges_have_curves() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, _) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    for (_, edge) in &brep.edges {
        assert!(
            edge.curve.is_some(),
            "every edge should have curve geometry"
        );
    }
}

#[test]
fn rectangle_half_edges_have_pcurves() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, _) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    for (_, he) in &brep.half_edges {
        assert!(he.pcurve.is_some(), "every half-edge should have a pcurve");
    }
}

#[test]
fn rectangle_faces_have_surfaces() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, _) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    for (_, face) in &brep.faces {
        assert!(
            face.surface.is_some(),
            "every face should have surface geometry"
        );
    }
}

#[test]
fn rectangle_euler_poincare_passes() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, _) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    // V=4, E=4, F=2, L=2, S=1
    // V-E+F-(L-F) = 4-4+2-(2-2) = 2 = 2*1
    assert!(
        brep.validate_euler_poincare().is_ok(),
        "Euler-Poincare should pass for a 2-sided rectangle sheet"
    );
}

#[test]
fn triangle_builds_valid_brep() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(2.0, 0.0);
    let p2 = sk.add_point(1.0, 2.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p0);

    let tol = Tolerance::default();
    let (graph, _) = edge_graph::build_edge_graph(&sk, &tol);
    let loops = loop_tracer::trace_loops(&graph);
    let roots = nesting::classify_nesting(&loops, &sk);

    let plane = SketchPlane::xy();
    let (brep, regions) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    assert_eq!(regions.len(), 1);
    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.faces.len(), 2);
    assert!(brep.validate_euler_poincare().is_ok());
}

#[test]
fn rectangle_half_edge_cycles_are_closed() {
    let (loops, roots, sk) = rectangle_pipeline();
    let plane = SketchPlane::xy();
    let tol = Tolerance::default();

    let (brep, _) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    // Walk every loop's half-edges and verify they form closed cycles.
    for (lh, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        if start == popkern_core::topology::HalfEdgeHandle::default() {
            continue;
        }
        let mut cursor = start;
        let mut count = 0;
        loop {
            count += 1;
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
            assert!(count < 100, "loop {lh:?} did not close");
        }
        assert!(
            count >= 3,
            "loop should have at least 3 half-edges, got {count}"
        );
    }
}

#[test]
fn rectangle_with_hole_builds() {
    // Outer rectangle + inner smaller rectangle (hole).
    let mut sk = Sketch::new();
    // Outer: (0,0) - (4,0) - (4,4) - (0,4)
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(4.0, 0.0);
    let p2 = sk.add_point(4.0, 4.0);
    let p3 = sk.add_point(0.0, 4.0);
    let _l0 = sk.add_line(p0, p1);
    let _l1 = sk.add_line(p1, p2);
    let _l2 = sk.add_line(p2, p3);
    let _l3 = sk.add_line(p3, p0);

    // Inner: (1,1) - (3,1) - (3,3) - (1,3)
    let q0 = sk.add_point(1.0, 1.0);
    let q1 = sk.add_point(3.0, 1.0);
    let q2 = sk.add_point(3.0, 3.0);
    let q3 = sk.add_point(1.0, 3.0);
    let _m0 = sk.add_line(q0, q1);
    let _m1 = sk.add_line(q1, q2);
    let _m2 = sk.add_line(q2, q3);
    let _m3 = sk.add_line(q3, q0);

    let tol = Tolerance::default();
    let (graph, _) = edge_graph::build_edge_graph(&sk, &tol);
    let loops = loop_tracer::trace_loops(&graph);
    let roots = nesting::classify_nesting(&loops, &sk);

    let plane = SketchPlane::xy();
    let (brep, regions) = build_brep(&roots, &loops, &sk, &plane, &tol).unwrap();

    // Should have 1 face region with 1 hole.
    assert_eq!(regions.len(), 1, "expected 1 region for outer+hole");
    assert_eq!(regions[0].holes.len(), 1, "expected 1 hole");

    // 8 vertices (4 outer + 4 inner), 8 edges, 2 faces, 4 loops
    assert_eq!(brep.vertices.len(), 8);
    assert_eq!(brep.edges.len(), 8);
    assert_eq!(brep.faces.len(), 2);
    // 2 outer loops + 2 inner loops = 4 loops
    assert_eq!(brep.loops.len(), 4);
}
