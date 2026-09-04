use super::*;
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle as CoreFaceHandle;

/// Build a minimal graph with one triangular face fragment and assemble it.
#[test]
fn assemble_single_triangle() {
    let mut graph = IntersectionGraph::new();

    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

    let eh = popkern_core::topology::EdgeHandle::default();
    let face_ref = FaceRef {
        operand: Operand::A,
        face: CoreFaceHandle::default(),
    };

    let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
    let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
    let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

    graph.set_half_edge_face(fwd01, face_ref);
    graph.set_half_edge_face(fwd12, face_ref);
    graph.set_half_edge_face(fwd20, face_ref);

    let fragment = Fragment {
        face: face_ref,
        outer_loop: vec![fwd01, fwd12, fwd20],
        inner_loops: vec![],
    };

    let selection = FragmentSelection {
        keep: true,
        flip: false,
    };

    // Use an empty BRep as source (no surface geometry to copy).
    let brep_a = BRep::new(Tolerance::default());
    let brep_b = BRep::new(Tolerance::default());

    let result = assemble_result(
        &graph,
        &[fragment],
        &[selection],
        &brep_a,
        &brep_b,
        Tolerance::default(),
    );

    assert!(result.is_ok(), "Assembly should succeed");
    let brep = result.unwrap_or_else(|_| panic!("unreachable"));
    assert_eq!(brep.vertices.len(), 3, "Expected 3 vertices");
    assert_eq!(brep.edges.len(), 3, "Expected 3 edges");
    assert_eq!(brep.faces.len(), 1, "Expected 1 face");
}

/// A loop whose resolved half-edges are all the wrong-direction twins in
/// the walker's original sequence (the manifold-assembly collision on a
/// multi-edge shared loop over a union-derived cap: the primary
/// half-edge was claimed by the adjacent, oppositely-oriented fragment,
/// so `resolve_loop_hes` fell back to the twins) still forms one closed
/// cycle — `reorder_into_closed_chain` re-sequences it head-to-tail
/// without fabricating or re-choosing any edge.
#[test]
fn reorder_recovers_twin_fallback_loop() {
    use popkern_core::constructors::builder::BRepBuilder;

    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));

    // Forward half-edges run v0->v1->v2->v3->v0; their twins run the
    // reverse of each edge.
    let (_, f01, t10) = b.add_edge(v0, v1);
    let (_, f12, t21) = b.add_edge(v1, v2);
    let (_, f23, t32) = b.add_edge(v2, v3);
    let (_, f30, t03) = b.add_edge(v3, v0);

    // An already-connected chain is left untouched (returns None).
    assert!(reorder_into_closed_chain(b.brep(), &[f01, f12, f23, f30]).is_none());

    // The twin-fallback delivers reversed edges in the forward walker
    // order — not a head-to-tail chain. Reorder recovers the reversed
    // cycle.
    let scrambled = [t10, t21, t32, t03];
    let chained =
        reorder_into_closed_chain(b.brep(), &scrambled).expect("twin set forms one closed cycle");
    assert_eq!(chained.len(), 4);
    let start_of = |h| b.brep().half_edges[h].origin;
    let end_of = |h| b.brep().half_edges[b.brep().half_edges[h].twin].origin;
    for i in 0..4 {
        assert_eq!(
            end_of(chained[i]),
            start_of(chained[(i + 1) % 4]),
            "reordered handles must chain head-to-tail"
        );
    }

    // A set that cannot close (a dangling spur) is declined, so the
    // caller keeps its order and declines honestly downstream.
    assert!(reorder_into_closed_chain(b.brep(), &[f01, f12]).is_none());
}

/// Verify that unselected fragments are skipped.
#[test]
fn unselected_fragments_are_skipped() {
    let mut graph = IntersectionGraph::new();

    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

    let eh = popkern_core::topology::EdgeHandle::default();
    let face_ref = FaceRef {
        operand: Operand::A,
        face: CoreFaceHandle::default(),
    };

    let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
    let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
    let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

    graph.set_half_edge_face(fwd01, face_ref);
    graph.set_half_edge_face(fwd12, face_ref);
    graph.set_half_edge_face(fwd20, face_ref);

    let fragment = Fragment {
        face: face_ref,
        outer_loop: vec![fwd01, fwd12, fwd20],
        inner_loops: vec![],
    };

    let selection = FragmentSelection {
        keep: false,
        flip: false,
    };

    let brep_a = BRep::new(Tolerance::default());
    let brep_b = BRep::new(Tolerance::default());

    let result = assemble_result(
        &graph,
        &[fragment],
        &[selection],
        &brep_a,
        &brep_b,
        Tolerance::default(),
    );

    assert!(result.is_ok());
    let brep = result.unwrap_or_else(|_| panic!("unreachable"));
    assert_eq!(brep.faces.len(), 0, "No faces should be in result");
}
