use super::*;
use std::collections::HashSet;

use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, BrepEdgePairPolicyContext, BrepEdgePairPolicyDeclineReason, BrepEdgePairPolicyOutcome,
    BrepTopologyRepairPolicy, Edge, Face, FaceHandle, HalfEdge, Loop, LoopHandle, Shell,
    ShellHandle, Solid, Vertex, assess_edge_pair_topology_policy,
};

fn options() -> BrepSewingOptions {
    BrepSewingOptions {
        tolerance: 1e-9,
        samples: 3,
        require_surface_agreement: false,
        validate_output: false,
        allow_nonmanifold_edges: false,
    }
}

fn provenance() -> BrepSewingProvenance {
    BrepSewingProvenance {
        first: "left test face".to_string(),
        second: "right test face".to_string(),
    }
}

fn policy() -> BrepTopologyRepairPolicy {
    BrepTopologyRepairPolicy {
        allow_nonmanifold_edges: false,
    }
}

fn two_square_sheet() -> (BRep, EdgeHandle, EdgeHandle) {
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
    let right = add_square(
        &mut brep,
        shell,
        [
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
    );
    if brep.shells[shell].face == FaceHandle::default() {
        brep.shells[shell].face = left.face;
    }
    (brep, left.edges[1], right.edges[2])
}

fn policy_context(used_edges: &HashSet<EdgeHandle>) -> BrepEdgePairPolicyContext<'_> {
    BrepEdgePairPolicyContext { used_edges }
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
fn sewing_clean_non_pns_boundary_pair_is_transactional() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let result = sew_clean_boundary_edge_pairs(
        &brep,
        &[BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        }],
        options(),
    )
    .unwrap();

    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
    assert_eq!(result.report.candidate_pairs, 1);
    assert_eq!(result.report.accepted_pairs, 1);
    assert!(result.report.declined_pairs.is_empty());
    assert_eq!(result.report.remaining_boundary_edges, 6);
    assert_eq!(result.brep.vertices.len(), 6);
    assert!(result.brep.validate_vertex_consistency().is_ok());
}

#[test]
fn migration_matrix_sewing_report_contract() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let result = sew_clean_boundary_edge_pairs(
        &brep,
        &[BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        }],
        options(),
    )
    .unwrap();

    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
    assert_eq!(result.report.candidate_pairs, 1);
    assert_eq!(result.report.accepted_pairs, 1);
    assert_eq!(result.report.declined_pairs.len(), 0);
    assert_eq!(result.report.remaining_boundary_edges, 6);
    assert_eq!(result.report.max_endpoint_gap, 0.0);
    assert_eq!(result.report.max_surface_gap, 0.0);
}

#[test]
fn policy_accepts_clean_boundary_pair_before_geometry_checks() {
    let (brep, left_edge, right_edge) = two_square_sheet();
    let used_edges = HashSet::new();
    let decision = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&used_edges),
    );

    assert!(decision.is_accepted());
    assert!(matches!(
        decision.outcome,
        BrepEdgePairPolicyOutcome::Accepted
    ));
    assert!(decision.first_edge.is_boundary);
    assert!(decision.second_edge.is_boundary);
}

#[test]
fn policy_reports_same_edge_used_boundary_missing_and_ambiguous_cases() {
    let (mut brep, left_edge, right_edge) = two_square_sheet();
    let mut used_edges = HashSet::new();
    used_edges.insert(left_edge);
    let used = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&used_edges),
    );
    assert!(matches!(
        used.outcome,
        BrepEdgePairPolicyOutcome::Declined(BrepEdgePairPolicyDeclineReason::EdgeAlreadyUsed)
    ));

    let empty = HashSet::new();
    let same = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: left_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&empty),
    );
    assert!(matches!(
        same.outcome,
        BrepEdgePairPolicyOutcome::Declined(BrepEdgePairPolicyDeclineReason::SameEdge)
    ));

    brep.boundary_edges.remove(&left_edge);
    let not_boundary = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&empty),
    );
    assert!(matches!(
        not_boundary.outcome,
        BrepEdgePairPolicyOutcome::Declined(BrepEdgePairPolicyDeclineReason::NotBoundaryEdge)
    ));

    brep.boundary_edges.insert(left_edge);
    brep.edges.remove(right_edge);
    let missing = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&empty),
    );
    assert!(matches!(
        missing.outcome,
        BrepEdgePairPolicyOutcome::Declined(BrepEdgePairPolicyDeclineReason::MissingEdge)
    ));

    let other_edge = brep
        .boundary_edges
        .iter()
        .copied()
        .find(|edge| *edge != left_edge && brep.edges.contains_key(*edge))
        .unwrap();
    brep.edges[left_edge].half_edge = HalfEdgeHandle::default();
    let ambiguous = assess_edge_pair_topology_policy(
        &brep,
        &BrepSewingEdgePair {
            first: left_edge,
            second: other_edge,
            provenance: provenance(),
        },
        policy(),
        policy_context(&empty),
    );
    assert!(matches!(
        ambiguous.outcome,
        BrepEdgePairPolicyOutcome::Declined(BrepEdgePairPolicyDeclineReason::AmbiguousTopology)
    ));
}

#[test]
fn sewing_declines_endpoint_gap_without_mutating_input() {
    let (brep, left_edge, _) = two_square_sheet();
    let other_edge = brep
        .boundary_edges
        .iter()
        .copied()
        .find(|edge| {
            if *edge == left_edge {
                return false;
            }
            let Some(left) = edge_segment(&brep, left_edge) else {
                return false;
            };
            let Some(other) = edge_segment(&brep, *edge) else {
                return false;
            };
            left.0
                .distance_to(&other.1)
                .max(left.1.distance_to(&other.0))
                > 0.5
        })
        .unwrap();
    let result = sew_clean_boundary_edge_pairs(
        &brep,
        &[BrepSewingEdgePair {
            first: left_edge,
            second: other_edge,
            provenance: provenance(),
        }],
        options(),
    )
    .unwrap();

    assert_eq!(brep.boundary_edges.len(), 8);
    assert_eq!(result.report.accepted_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert!(matches!(
        result.report.declined_pairs[0].reason,
        BrepSewingDeclineReason::EndpointGapExceeded { .. }
    ));
    assert_eq!(result.report.remaining_boundary_edges, 8);
}

#[test]
fn sewing_declines_already_joined_candidate_with_explanation() {
    let (mut brep, left_edge, right_edge) = two_square_sheet();
    let left_half_edge = brep.edges[left_edge].half_edge;
    let right_half_edge = brep.edges[right_edge].half_edge;
    brep.half_edges[left_half_edge].twin = right_half_edge;

    let result = sew_clean_boundary_edge_pairs(
        &brep,
        &[BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        }],
        options(),
    )
    .unwrap();

    assert_eq!(result.report.accepted_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert!(matches!(
        result.report.declined_pairs[0].reason,
        BrepSewingDeclineReason::AlreadyJoined
    ));
    let policy = result.report.declined_pairs[0].policy.as_ref().unwrap();
    assert!(policy.explanation.contains("already joined"));
    assert!(policy.first_edge.has_existing_twin);
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
}

#[test]
fn sewing_declines_radial_nonmanifold_candidate_with_explanation() {
    let (mut brep, left_edge, right_edge) = two_square_sheet();
    let left_half_edge = brep.edges[left_edge].half_edge;
    brep.nonmanifold_radial
        .insert(left_half_edge, left_half_edge);

    let result = sew_clean_boundary_edge_pairs(
        &brep,
        &[BrepSewingEdgePair {
            first: left_edge,
            second: right_edge,
            provenance: provenance(),
        }],
        options(),
    )
    .unwrap();

    assert_eq!(result.report.accepted_pairs, 0);
    assert_eq!(result.report.declined_pairs.len(), 1);
    assert!(matches!(
        result.report.declined_pairs[0].reason,
        BrepSewingDeclineReason::NonManifoldConflict
    ));
    let policy = result.report.declined_pairs[0].policy.as_ref().unwrap();
    assert!(policy.explanation.contains("radial non-manifold"));
    assert!(policy.first_edge.has_radial_cycle);
    assert_eq!(brep.boundary_edges.len(), 8, "input BRep is untouched");
}
