//! Unit tests for the walk module. Split out of `mod.rs` to keep
//! the parent under the 500-line structural cap.

use super::*;
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

/// Empty BReps for tests that only use planar geometry.
fn empty_breps() -> (BRep, BRep) {
    (
        BRep::new(Tolerance::default()),
        BRep::new(Tolerance::default()),
    )
}

/// Build a triangle face in the graph with 3 nodes and 3 edge pairs.
fn triangle_graph() -> (IntersectionGraph, FaceRef) {
    let mut graph = IntersectionGraph::new();

    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

    let eh = popkern_core::topology::EdgeHandle::default();
    let face_ref = FaceRef {
        operand: Operand::A,
        face: popkern_core::topology::FaceHandle::default(),
    };

    let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
    let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
    let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

    graph.set_half_edge_face(fwd01, face_ref);
    graph.set_half_edge_face(fwd12, face_ref);
    graph.set_half_edge_face(fwd20, face_ref);

    (graph, face_ref)
}

#[test]
fn triangle_produces_one_outer_loop() {
    let (graph, _face_ref) = triangle_graph();
    let (ba, bb) = empty_breps();

    let fragments = walk_all_fragments(
        &graph,
        &ba,
        &bb,
        &popkern_core::tolerance::Tolerance::default(),
    );

    assert_eq!(fragments.len(), 1, "Expected 1 fragment for triangle");
    assert_eq!(
        fragments[0].outer_loop.len(),
        3,
        "Expected 3 half-edges in outer loop"
    );
    assert!(fragments[0].inner_loops.is_empty());
}

#[test]
fn square_face_produces_one_outer_loop() {
    let mut graph = IntersectionGraph::new();

    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    let n2 = graph.add_node(Point3::new(1.0, 1.0, 0.0), NodeSource::Intersection);
    let n3 = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);

    let eh = popkern_core::topology::EdgeHandle::default();
    let face_ref = FaceRef {
        operand: Operand::A,
        face: popkern_core::topology::FaceHandle::default(),
    };

    let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
    let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
    let (fwd23, _) = graph.add_half_edge_pair(n2, n3, EdgeSource::OriginalA(eh), None);
    let (fwd30, _) = graph.add_half_edge_pair(n3, n0, EdgeSource::OriginalA(eh), None);

    graph.set_half_edge_face(fwd01, face_ref);
    graph.set_half_edge_face(fwd12, face_ref);
    graph.set_half_edge_face(fwd23, face_ref);
    graph.set_half_edge_face(fwd30, face_ref);
    let (ba, bb) = empty_breps();

    let fragments = walk_all_fragments(
        &graph,
        &ba,
        &bb,
        &popkern_core::tolerance::Tolerance::default(),
    );

    assert_eq!(fragments.len(), 1, "Expected 1 fragment for square");
    assert_eq!(
        fragments[0].outer_loop.len(),
        4,
        "Expected 4 half-edges in outer loop"
    );
    assert!(fragments[0].inner_loops.is_empty());
}
