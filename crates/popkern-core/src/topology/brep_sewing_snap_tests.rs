use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, BrepSewingClassificationOptions, BrepSewingProvenance, Edge, Face, FaceHandle, HalfEdge,
    HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle, Solid, Vertex,
    classify_boundary_edge_pairs,
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

fn snap_options(max_vertex_motion: f64) -> BrepSnapSewingOptions {
    BrepSnapSewingOptions {
        max_vertex_motion,
        clean_sewing: BrepSewingOptions {
            tolerance: 1e-9,
            samples: 3,
            require_surface_agreement: false,
            validate_output: false,
            allow_nonmanifold_edges: false,
        },
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
fn snap_sewing_near_pair_moves_clone_then_welds() {
    let (brep, left_edge, right_edge) = two_square_sheet(5e-5, true);
    let before_vertices = brep.vertices.len();
    let before_boundary_edges = brep.boundary_edges.len();
    let classification = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );

    let result =
        snap_and_sew_classified_edge_pairs(&brep, &classification.pairs, snap_options(1e-3))
            .unwrap();

    assert_eq!(brep.vertices.len(), before_vertices, "input is untouched");
    assert_eq!(
        brep.boundary_edges.len(),
        before_boundary_edges,
        "input boundary set is untouched"
    );
    assert_eq!(result.report.snapped_pairs, 1);
    assert_eq!(result.report.exact_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 0);
    assert!(
        !result.report.shared_fact_ids.is_empty(),
        "snap-sewing report should cite shared pair facts"
    );
    assert!(result.report.max_vertex_motion > 0.0);
    assert!(result.report.max_vertex_motion <= 2.5e-5 + 1e-12);
    assert_eq!(result.brep.vertices.len(), 6);
    assert_eq!(result.brep.boundary_edges.len(), 6);
    assert!(
        result
            .brep
            .vertices
            .values()
            .any(|vertex| vertex.tolerance.unwrap_or(0.0) >= result.report.max_vertex_motion)
    );
}

#[test]
fn snap_sewing_declines_rejected_classification() {
    let (brep, left_edge, right_edge) = two_square_sheet(5e-2, true);
    let classification = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );
    let result =
        snap_and_sew_classified_edge_pairs(&brep, &classification.pairs, snap_options(1e-3))
            .unwrap();

    assert_eq!(result.report.snapped_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert!(
        !result.report.declined_pairs[0].shared_fact_ids.is_empty(),
        "snap-sewing decline should cite shared pair facts"
    );
    assert!(matches!(
        result.report.declined_pairs[0].reason,
        BrepSnapSewingDeclineReason::RejectedClassification(_)
    ));
    assert_eq!(result.brep.boundary_edges.len(), brep.boundary_edges.len());
}

#[test]
fn snap_sewing_declines_motion_over_policy_limit() {
    let (brep, left_edge, right_edge) = two_square_sheet(5e-5, true);
    let classification = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );
    let result =
        snap_and_sew_classified_edge_pairs(&brep, &classification.pairs, snap_options(1e-6))
            .unwrap();

    assert_eq!(result.report.snapped_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert!(matches!(
        result.report.declined_pairs[0].reason,
        BrepSnapSewingDeclineReason::VertexMotionExceeded { .. }
    ));
    assert_eq!(result.brep.vertices.len(), brep.vertices.len());
    assert_eq!(result.brep.boundary_edges.len(), brep.boundary_edges.len());
}

#[test]
fn snap_sewing_rejected_orientation_does_not_mutate_input() {
    let (brep, left_edge, right_edge) = two_square_sheet(0.0, false);
    let classification = classify_boundary_edge_pairs(
        &brep,
        &[edge_pair(left_edge, right_edge)],
        classification_options(),
    );
    let result =
        snap_and_sew_classified_edge_pairs(&brep, &classification.pairs, snap_options(1e-3))
            .unwrap();

    assert_eq!(result.report.snapped_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert_eq!(brep.vertices.len(), 8);
    assert_eq!(brep.boundary_edges.len(), 8);
}
