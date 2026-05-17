use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, BrepSewingProvenance, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
    Shell, ShellHandle, Solid, Vertex,
};

fn classification_options() -> BrepSewingClassificationOptions {
    BrepSewingClassificationOptions {
        exact_tolerance: 1e-9,
        snap_tolerance: 1e-3,
        samples: 3,
        require_surface_agreement: false,
        allow_nonmanifold_edges: false,
    }
}

fn provenance() -> BrepSewingProvenance {
    BrepSewingProvenance {
        first: "left test face".to_string(),
        second: "right test face".to_string(),
    }
}

fn edge_pair(first: EdgeHandle, second: EdgeHandle) -> BrepSewingEdgePair {
    BrepSewingEdgePair {
        first,
        second,
        provenance: provenance(),
    }
}

fn two_square_sheet(offset: f64, reversed_right: bool) -> (BRep, EdgeHandle, EdgeHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    let left = add_square(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    let right_points = if reversed_right {
        [
            Point3::new(2.0 + offset, 0.0, 0.0),
            Point3::new(2.0 + offset, 1.0, 0.0),
            Point3::new(1.0 + offset, 1.0, 0.0),
            Point3::new(1.0 + offset, 0.0, 0.0),
        ]
    } else {
        [
            Point3::new(2.0 + offset, 0.0, 0.0),
            Point3::new(1.0 + offset, 0.0, 0.0),
            Point3::new(1.0 + offset, 1.0, 0.0),
            Point3::new(2.0 + offset, 1.0, 0.0),
        ]
    };
    let right = add_square(&mut brep, shell, right_points);
    if brep.shells[shell].face == FaceHandle::default() {
        brep.shells[shell].face = left.face;
    }
    let right_edge = if reversed_right {
        right.edges[2]
    } else {
        right.edges[1]
    };
    (brep, left.edges[1], right_edge)
}

struct SquareHandles {
    face: FaceHandle,
    edges: [EdgeHandle; 4],
}

fn add_square(brep: &mut BRep, shell: ShellHandle, points: [Point3; 4]) -> SquareHandles {
    let vertices = points.map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        })
    });
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;

    let half_edges = vertices.map(|origin| {
        brep.half_edges.insert(HalfEdge {
            origin,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: EdgeHandle::default(),
            loop_,
            pcurve: None,
        })
    });
    let edges = half_edges.map(|half_edge| {
        let edge = brep.edges.insert(Edge {
            half_edge,
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[half_edge].edge = edge;
        brep.boundary_edges.insert(edge);
        edge
    });

    for index in 0..4 {
        let prev = (index + 3) % 4;
        let next = (index + 1) % 4;
        brep.half_edges[half_edges[index]].prev = half_edges[prev];
        brep.half_edges[half_edges[index]].next = half_edges[next];
        brep.vertices[vertices[index]].outgoing = half_edges[index];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    SquareHandles { face, edges }
}

#[test]
fn classification_identifies_exact_clean_pair_without_mutation() {
    let (brep, left_edge, right_edge) = two_square_sheet(0.0, true);
    let report = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );

    assert_eq!(brep.boundary_edges.len(), 8);
    assert_eq!(brep.vertices.len(), 8);
    assert_eq!(report.exact_clean_pairs, 1);
    assert_eq!(report.near_clean_pairs, 0);
    assert_eq!(report.rejected_pairs, 0);
    assert_bucket(&report, BrepSewingClassificationBucket::ExactClean, 1);
    assert!(matches!(
        report.pairs[0].classification,
        BrepEdgePairClassification::ExactClean
    ));
    assert!(
        !report.shared_fact_ids.is_empty(),
        "classification report should cite shared pair facts"
    );
    assert!(
        !report.pairs[0].shared_fact_ids.is_empty(),
        "classified pair should cite shared pair facts"
    );
}

#[test]
fn classification_identifies_near_pair_inside_snap_tolerance() {
    let (brep, left_edge, right_edge) = two_square_sheet(5e-5, true);
    let report = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );

    assert_eq!(brep.boundary_edges.len(), 8);
    assert_eq!(report.exact_clean_pairs, 0);
    assert_eq!(report.near_clean_pairs, 1);
    assert_eq!(report.rejected_pairs, 0);
    assert_bucket(&report, BrepSewingClassificationBucket::NearClean, 1);
    assert!(matches!(
        report.pairs[0].classification,
        BrepEdgePairClassification::NearClean
    ));
    assert!(report.max_endpoint_gap > classification_options().exact_tolerance);
}

#[test]
fn classification_rejects_pair_beyond_snap_tolerance() {
    let (brep, left_edge, right_edge) = two_square_sheet(5e-2, true);
    let report = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );

    assert_eq!(report.exact_clean_pairs, 0);
    assert_eq!(report.near_clean_pairs, 0);
    assert_eq!(report.rejected_pairs, 1);
    assert_bucket(
        &report,
        BrepSewingClassificationBucket::EndpointGapExceeded,
        1,
    );
    assert!(matches!(
        report.pairs[0].classification,
        BrepEdgePairClassification::Rejected(BrepSewingDeclineReason::EndpointGapExceeded { .. })
    ));
    assert!(
        !report.pairs[0].shared_fact_ids.is_empty(),
        "rejected pair should retain shared evidence"
    );
}

#[test]
fn classification_rejects_same_direction_candidate() {
    let (brep, left_edge, right_edge) = two_square_sheet(0.0, false);
    let report = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );

    assert_eq!(report.exact_clean_pairs, 0);
    assert_eq!(report.near_clean_pairs, 0);
    assert_eq!(report.rejected_pairs, 1);
    assert_bucket(
        &report,
        BrepSewingClassificationBucket::ForwardOrientation,
        1,
    );
    assert!(matches!(
        report.pairs[0].classification,
        BrepEdgePairClassification::Rejected(BrepSewingDeclineReason::EndpointGapExceeded { .. })
    ));
}

fn assert_bucket(
    report: &BrepSewingClassificationReport,
    bucket: BrepSewingClassificationBucket,
    count: usize,
) {
    assert_eq!(
        report
            .buckets
            .iter()
            .find(|item| item.bucket == bucket)
            .map(|item| item.count),
        Some(count)
    );
}
