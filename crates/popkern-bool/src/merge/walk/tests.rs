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

/// PK-438: the exterior-disjoint promotion moves a multi-edge inner
/// loop only when the fragment's OUTER wall is kept (inside the opposite
/// operand) AND the inner loop is bbox-disjoint from it AND the inner
/// centroid is OUTSIDE the opposite operand (a below-body tool
/// protrusion). Everything else stays attached: a disjoint loop that
/// AGREES with the outer (both inside the operand — a real hole), a loop
/// overlapping the outer's box, a single self-loop rim, and — the holder
/// falsification — every inner loop of a fragment whose outer is OUTSIDE
/// the operand (a workpiece face carrying a genuine through-hole).
#[test]
fn exterior_disjoint_promotion_separates_protrusion_from_kept_holes() {
    use popkern_core::geometry::Vec3;

    let mut graph = IntersectionGraph::new();
    let eh = popkern_core::topology::EdgeHandle::default();
    // Add a straight triangle loop with the given node positions; return
    // its three forward half-edge indices.
    let tri = |graph: &mut IntersectionGraph, pts: [Point3; 3]| -> Vec<usize> {
        let n: Vec<usize> = pts
            .iter()
            .map(|&p| graph.add_node(p, NodeSource::Intersection))
            .collect();
        let (a, _) = graph.add_half_edge_pair(n[0], n[1], EdgeSource::OriginalA(eh), None);
        let (b, _) = graph.add_half_edge_pair(n[1], n[2], EdgeSource::OriginalA(eh), None);
        let (c, _) = graph.add_half_edge_pair(n[2], n[0], EdgeSource::OriginalA(eh), None);
        vec![a, b, c]
    };

    // Opposite solid (operand B): a cube spanning [0,1]² × [2,3].
    let cube = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(1.0, 0.0, 2.0),
            Point3::new(1.0, 1.0, 2.0),
            Point3::new(0.0, 1.0, 2.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("cube");
    let empty_a = BRep::new(Tolerance::default());
    let opposite = hole_nesting::OppositeAabbs::new(&empty_a, &cube, None);

    // Outer wall INSIDE the cube (z = 2.4) — the kept-wall configuration.
    let outer_inside = tri(
        &mut graph,
        [
            Point3::new(0.3, 0.3, 2.4),
            Point3::new(0.7, 0.3, 2.4),
            Point3::new(0.5, 0.7, 2.4),
        ],
    );
    // Protrusion: disjoint (z = 10) and centroid OUTSIDE the cube.
    let protrusion = tri(
        &mut graph,
        [
            Point3::new(0.4, 0.2, 10.0),
            Point3::new(0.6, 0.2, 10.0),
            Point3::new(0.5, 0.4, 10.0),
        ],
    );
    // Agreeing hole: disjoint from the outer (z = 2.9) but centroid also
    // INSIDE the cube — same membership as the outer, so it stays.
    let agree_hole = tri(
        &mut graph,
        [
            Point3::new(0.4, 0.4, 2.9),
            Point3::new(0.6, 0.4, 2.9),
            Point3::new(0.5, 0.6, 2.9),
        ],
    );
    // Nested hole: box overlaps the outer (z = 2.4) — not disjoint.
    let nested = tri(
        &mut graph,
        [
            Point3::new(0.35, 0.35, 2.4),
            Point3::new(0.65, 0.35, 2.4),
            Point3::new(0.5, 0.6, 2.4),
        ],
    );
    // Single self-loop rim: disjoint (z = 10) but length 1.
    let nr = graph.add_node(Point3::new(0.5, 0.5, 10.0), NodeSource::Intersection);
    let (rim, _) = graph.add_half_edge_pair(nr, nr, EdgeSource::OriginalA(eh), None);

    let mut inners = vec![
        protrusion.clone(),
        agree_hole.clone(),
        nested.clone(),
        vec![rim],
    ];
    let promoted = hole_nesting::promote_exterior_disjoint_holes(
        &graph,
        &mut inners,
        &outer_inside,
        Operand::A,
        &opposite,
        &Tolerance::default(),
    );
    assert_eq!(promoted.len(), 1, "only the protrusion promotes");
    assert_eq!(
        promoted[0], protrusion,
        "the promoted loop is the protrusion"
    );
    assert_eq!(inners.len(), 3, "three kept inner loops remain");
    assert!(
        inners.contains(&agree_hole),
        "agreeing inside-cube hole stays"
    );
    assert!(inners.contains(&nested), "box-overlapping hole stays");
    assert!(inners.contains(&vec![rim]), "single self-loop rim stays");

    // Holder falsification: a WORKPIECE face whose outer is OUTSIDE the
    // (thin) tool carries a genuine through-hole whose centroid also
    // reads Outside. Same disjoint protrusion loop, but with an outer
    // that is Outside the cube, must NOT promote anything.
    let outer_outside = tri(
        &mut graph,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
        ],
    );
    let mut inners2 = vec![protrusion.clone()];
    let promoted2 = hole_nesting::promote_exterior_disjoint_holes(
        &graph,
        &mut inners2,
        &outer_outside,
        Operand::A,
        &opposite,
        &Tolerance::default(),
    );
    assert!(
        promoted2.is_empty(),
        "no promotion when the outer wall is not kept (outside the operand)"
    );
    assert_eq!(inners2.len(), 1, "the loop stays a hole");
}
