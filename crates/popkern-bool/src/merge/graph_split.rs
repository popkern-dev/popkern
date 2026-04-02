//! Endpoint classification and graph edge splitting.
//!
//! Provides tools to determine where a 3D point falls relative to
//! existing graph topology (on a vertex, on an edge, or in a face
//! interior), and to split a graph half-edge pair at a given point.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::geometry::Point3;

use super::graph::{EdgeSource, GraphHalfEdge, IntersectionGraph, NodeSource};

/// Classification of an intersection curve endpoint relative to the graph.
#[derive(Debug)]
#[allow(clippy::enum_variant_names)]
pub(crate) enum EndpointClass {
    /// The point coincides with an existing graph node.
    OnVertex(usize),
    /// The point lies on a graph edge. Carries the half-edge index and the
    /// parametric position along the segment (0.0 = origin, 1.0 = end).
    OnEdge(usize, f64),
    /// The point is interior to a face (not on any existing edge or vertex).
    OnFaceInterior,
}

/// Classify where a point falls in the intersection graph.
///
/// Checks nodes first (within `tolerance`), then edges (point-to-segment
/// distance), and falls back to [`EndpointClass::OnFaceInterior`].
pub(crate) fn classify_endpoint(
    graph: &IntersectionGraph,
    point: &Point3,
    tolerance: f64,
) -> EndpointClass {
    // 1. Check all nodes.
    if let Some(idx) = graph.find_node(point, tolerance) {
        return EndpointClass::OnVertex(idx);
    }

    // 2. Check all half-edges (only check forward half-edges, i.e. even indices).
    //    For each edge pair (i, i+1), check if point is on the line segment.
    let mut i = 0;
    while i < graph.half_edges.len() {
        let he = &graph.half_edges[i];
        let twin = &graph.half_edges[he.twin];

        let p0 = graph.nodes[he.origin].position;
        let p1 = graph.nodes[twin.origin].position;

        if let Some(t) = point_on_segment(point, &p0, &p1, tolerance) {
            return EndpointClass::OnEdge(i, t);
        }

        // Skip the twin — we only need to check each edge once.
        i += 2;
    }

    // 3. Not on any vertex or edge.
    EndpointClass::OnFaceInterior
}

/// Check whether `point` lies on segment `p0..p1` within `tolerance`.
///
/// Returns `Some(t)` where t in (0,1) is the parameter, or `None`.
/// Endpoints (t near 0 or 1) are excluded — they would be caught by
/// the vertex check.
fn point_on_segment(point: &Point3, p0: &Point3, p1: &Point3, tolerance: f64) -> Option<f64> {
    let seg = *p1 - *p0;
    let seg_len_sq = seg.length_squared();
    if seg_len_sq < tolerance * tolerance {
        return None; // degenerate edge
    }

    let to_point = *point - *p0;
    let t = to_point.dot(&seg) / seg_len_sq;

    // Exclude endpoints (those should match as vertices).
    if t <= tolerance || t >= 1.0 - tolerance {
        return None;
    }

    // Check perpendicular distance.
    let projected = *p0 + seg * t;
    let dist = point.distance_to(&projected);
    if dist < tolerance { Some(t) } else { None }
}

/// Split a graph edge at a point, creating a new node.
///
/// The original half-edge pair (`he_idx` and its twin) is shortened to end
/// at the new node. A second pair is created from the new node to the
/// original endpoint. Face assignments are preserved on all affected
/// half-edges.
///
/// Returns the index of the new node.
pub(crate) fn split_graph_edge(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    split_point: Point3,
    _t_param: f64,
) -> usize {
    let twin_idx = graph.half_edges[he_idx].twin;

    // Original endpoints.
    let _start_node = graph.half_edges[he_idx].origin;
    let end_node = graph.half_edges[twin_idx].origin;

    // Preserve face assignments and source before mutation.
    let fwd_face = graph.half_edges[he_idx].face;
    let twin_face = graph.half_edges[twin_idx].face;
    let source = graph.half_edges[he_idx].source.clone();

    // Derive a split source from the original.
    let split_source = match &source {
        EdgeSource::OriginalA(eh) => EdgeSource::SplitA { original: *eh },
        EdgeSource::OriginalB(eh) => EdgeSource::SplitB { original: *eh },
        other => other.clone(),
    };

    // Create the new node at the split point.
    let new_node = graph.add_node(split_point, NodeSource::Intersection);

    // Create the second segment: new_node → end_node (and twin end_node → new_node).
    let new_fwd_idx = graph.half_edges.len();
    let new_twin_idx = new_fwd_idx + 1;

    graph.half_edges.push(GraphHalfEdge {
        origin: new_node,
        twin: new_twin_idx,
        face: fwd_face,
        curve: None,  // Curve splitting handled separately.
        pcurve: None, // Pcurve splitting handled separately.
        source: split_source.clone(),
    });
    graph.half_edges.push(GraphHalfEdge {
        origin: end_node,
        twin: new_fwd_idx,
        face: twin_face,
        curve: None,
        pcurve: None,
        source: split_source,
    });

    // Update outgoing lists.
    graph.nodes[new_node].outgoing.push(new_fwd_idx);
    graph.nodes[end_node].outgoing.push(new_twin_idx);

    // Remove old outgoing entry for end_node → twin_idx, since twin now
    // originates at end_node but goes to new_node instead of start_node.
    // Actually twin_idx still originates at end_node, but it now should go
    // to new_node. We update its target by re-pointing.

    // Update the original pair to end at new_node instead of end_node.
    // The forward half-edge (he_idx) keeps its origin (start_node).
    // Its twin (twin_idx) now originates at new_node instead of end_node.
    graph.half_edges[twin_idx].origin = new_node;

    // Fix outgoing: remove twin_idx from end_node, add to new_node.
    graph.nodes[end_node]
        .outgoing
        .retain(|&idx| idx != twin_idx);
    graph.nodes[new_node].outgoing.push(twin_idx);

    // Wire twins correctly:
    // he_idx ↔ twin_idx  (start_node → new_node, new_node → start_node)
    // new_fwd_idx ↔ new_twin_idx (new_node → end_node, end_node → new_node)
    // These are already set correctly by the push and original values.

    // Update face_edges: add new half-edges to the same faces.
    if let Some(face) = fwd_face {
        graph.face_edges.entry(face).or_default().push(new_fwd_idx);
    }
    if let Some(face) = twin_face {
        graph.face_edges.entry(face).or_default().push(new_twin_idx);
    }

    new_node
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};

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
}
