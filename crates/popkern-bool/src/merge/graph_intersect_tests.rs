//! Record insertion: shared edges, node reuse and the point-contact
//! and consequence rules for records whose ends share a node.

use super::*;
use popkern_core::{
    geometry::{NurbsCurve, Point3, curve_on_surface::CurveOnSurface},
    topology::FaceHandle,
};

/// Build a minimal IntersectionRecord with a line from p0 to p1.
fn make_record(
    p0: Point3,
    p1: Point3,
    face_a: FaceHandle,
    face_b: FaceHandle,
) -> IntersectionRecord {
    // Build a degree-1 NURBS line.
    let nurbs = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    };
    let cos = CurveOnSurface {
        curve_3d: nurbs.clone(),
        curve_uv1: nurbs.clone(),
        curve_uv2: nurbs,
        is_closed: false,
    };
    IntersectionRecord {
        residue_owner_a: None,
        face_a,
        face_b,
        curve: cos,
        is_tangent: false,
        is_coincident: false,
        start_boundary: Default::default(),
        end_boundary: Default::default(),
    }
}

#[test]
fn insert_shared_edge_between_two_faces() {
    let mut graph = IntersectionGraph::new();

    // Pre-create nodes at the endpoints.
    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    assert_eq!(n0, 0);
    assert_eq!(n1, 1);

    let face_a = FaceHandle::default();
    let face_b = FaceHandle::default();

    let record = make_record(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        face_a,
        face_b,
    );

    insert_intersection_curves(&mut graph, &[record], 1e-6);

    // Should still have 2 nodes (reused existing ones).
    assert_eq!(graph.nodes.len(), 2);
    // Should have 2 half-edges (one pair).
    assert_eq!(graph.half_edges.len(), 2);

    // Forward half-edge should have face A.
    let face_a_ref = FaceRef {
        operand: Operand::A,
        face: face_a,
    };
    let face_b_ref = FaceRef {
        operand: Operand::B,
        face: face_b,
    };
    assert_eq!(graph.half_edges[0].face, Some(face_a_ref));
    assert_eq!(graph.half_edges[1].face, Some(face_b_ref));

    // face_edges should list the shared edge for both faces.
    assert!(graph.face_edges[&face_a_ref].contains(&0));
    assert!(graph.face_edges[&face_b_ref].contains(&1));

    // Forward half-edge should carry the curve.
    assert!(graph.half_edges[0].curve.is_some());
}

#[test]
fn insert_creates_nodes_when_not_found() {
    let mut graph = IntersectionGraph::new();

    let face_a = FaceHandle::default();
    let face_b = FaceHandle::default();

    let record = make_record(
        Point3::new(5.0, 5.0, 5.0),
        Point3::new(6.0, 6.0, 6.0),
        face_a,
        face_b,
    );

    insert_intersection_curves(&mut graph, &[record], 1e-6);

    // Nodes should be created.
    assert_eq!(graph.nodes.len(), 2);
    assert_eq!(graph.half_edges.len(), 2);
}

/// A record whose ends share a node: the point-contact gate keeps a
/// spur longer than 500 tolerances and drops a loop inside them;
/// with the table deciding vertices, its consequence rule drops the
/// spur (it bounds nothing) and keeps the loop.
#[test]
fn the_table_keeps_loops_and_drops_spurs_whatever_their_size() {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let cubic = |pts: [Point3; 4]| {
        let mut r = make_record(pts[0], pts[3], FaceHandle::default(), FaceHandle::default());
        r.curve.curve_3d = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: pts.to_vec(),
            weights: vec![1.0; 4],
        };
        r
    };
    let records = [
        cubic([p(50.0, 0.0), p(50.002, 0.0), p(50.002, 0.0), p(50.0, 0.0)]),
        cubic([
            p(60.0, 0.0),
            p(60.0002, 2e-4),
            p(60.0002, -2e-4),
            p(60.0, 0.0),
        ]),
    ];
    let kept_at = |ends: Option<&[crate::tolerance_zone::table::RecordEnds]>| {
        let mut graph = IntersectionGraph::new();
        graph.add_node(p(0.0, 0.0), NodeSource::Intersection);
        graph.add_node(p(100.0, 0.0), NodeSource::Intersection);
        insert_intersection_curves_at(&mut graph, &records, 1e-6, ends);
        let xs: Vec<f64> = graph
            .half_edges
            .iter()
            .map(|h| graph.nodes[h.origin].position.x)
            .collect();
        xs
    };
    assert_eq!(kept_at(None), vec![50.0, 50.0]);
    let ends = [p(50.0, 0.0), p(60.0, 0.0)].map(|at| crate::tolerance_zone::table::RecordEnds {
        start: at,
        end: at,
        tolerance: 1e-6,
        paves: Vec::new(),
    });
    assert_eq!(kept_at(Some(&ends)), vec![60.0, 60.0]);
}

#[test]
fn validate_on_valid_graph() {
    let mut graph = IntersectionGraph::new();
    let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
    let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
    let edge_h = popkern_core::topology::EdgeHandle::default();
    graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

    let errors = graph.validate();
    assert!(errors.is_empty(), "Expected no errors, got: {errors:?}");
}

#[test]
fn validate_detects_orphan_node() {
    let mut graph = IntersectionGraph::new();
    // Node with no outgoing edges.
    graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);

    let errors = graph.validate();
    assert_eq!(errors.len(), 1);
    assert!(errors[0].contains("no outgoing edges"));
}

/// A record whose interior holds another record's end (a T-junction)
/// arrives split there when the table names the pave: the junction node
/// is the other record's end node, not a twin of it, and the two
/// pieces keep the record's true geometry either side of it.
#[test]
fn a_record_arrives_split_at_the_tables_paves() {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let (fa, fb) = (FaceHandle::default(), FaceHandle::default());
    let records = [
        make_record(p(0.0, 0.0), p(4.0, 0.0), fa, fb),
        make_record(p(1.0, 0.0), p(1.0, 3.0), fa, fb),
    ];
    let ends = [
        crate::tolerance_zone::table::RecordEnds {
            start: p(0.0, 0.0),
            end: p(4.0, 0.0),
            tolerance: 1e-6,
            paves: vec![(0.25, p(1.0, 0.0))],
        },
        crate::tolerance_zone::table::RecordEnds {
            start: p(1.0, 0.0),
            end: p(1.0, 3.0),
            tolerance: 1e-6,
            paves: Vec::new(),
        },
    ];
    let mut graph = IntersectionGraph::new();
    insert_intersection_curves_at(&mut graph, &records, 1e-6, Some(&ends));
    assert_eq!(graph.nodes.len(), 4, "the junction is one node");
    assert_eq!(graph.half_edges.len(), 6, "the first record is two pieces");
    let junction = graph.find_node(&p(1.0, 0.0), 1e-6).expect("junction node");
    assert_eq!(graph.nodes[junction].outgoing.len(), 3, "a T-junction");
    let windows: Vec<(f64, f64)> = [0, 4]
        .iter()
        .filter_map(|&i| graph.half_edges[i].curve.as_ref())
        .map(|c| (c.t_start, c.t_end))
        .collect();
    assert_eq!(windows.len(), 2);
    assert!((windows[0].1 - 0.25).abs() < 1e-9 && (windows[1].0 - 0.25).abs() < 1e-9);
}
