use super::{BridgeResourcePairing, CylinderBridgeResourceSegment, preflight_shared_start_bridge};
use crate::euler::{
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement,
};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, HalfEdge, HalfEdgeHandle, Loop, Shell, Solid, Vertex,
};

#[test]
fn shared_start_bridge_preflight_requires_complete_candidate_gate() {
    let fixture = shared_start_cross_loop_fixture();

    let preflight = preflight_shared_start_bridge(
        &fixture.brep,
        &fixture.pairings,
        &fixture.authoring,
        &fixture.cylinder_segments,
    );

    let label = preflight.diagnostic_label();
    assert!(label.contains("status=ready"), "{label}");
    assert!(label.contains("candidate_gate_complete=true"), "{label}");
    assert!(label.contains("construction_candidates=1"), "{label}");
    assert!(
        label.contains("candidate_requirements=live_revalidation:1 postcondition_verification:1"),
        "{label}"
    );
    assert!(
        label.contains("downstream_readiness=facts:1 ready:1"),
        "{label}"
    );
}

struct SharedStartCrossLoopFixture {
    brep: BRep,
    pairings: Vec<BridgeResourcePairing>,
    authoring: Vec<CylinderMediatedSupportContinuationSegmentAuthoringPlan>,
    cylinder_segments: Vec<CylinderBridgeResourceSegment>,
}

fn shared_start_cross_loop_fixture() -> SharedStartCrossLoopFixture {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let first = add_boundary_polygon(
        &mut brep,
        shell,
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    let second = add_boundary_polygon(
        &mut brep,
        shell,
        &[
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ],
    );
    brep.shells[shell].face = first;
    let outside = first_two_half_edges(&brep, brep.faces[first].outer_loop);
    let cylinder = first_two_half_edges(&brep, brep.faces[second].outer_loop);
    let authoring = vec![authoring_segment(
        brep.faces[first].outer_loop,
        outside.clone(),
    )];
    let cylinder_segments = vec![CylinderBridgeResourceSegment {
        edges: edges_for_half_edges(&brep, &cylinder),
        half_edges: cylinder,
        request_indices: vec![0],
    }];
    SharedStartCrossLoopFixture {
        brep,
        pairings: vec![BridgeResourcePairing {
            outside_segment_index: 0,
            cylinder_segment_index: 0,
        }],
        authoring,
        cylinder_segments,
    }
}

fn authoring_segment(
    outside_loop: crate::topology::LoopHandle,
    half_edges: Vec<HalfEdgeHandle>,
) -> CylinderMediatedSupportContinuationSegmentAuthoringPlan {
    CylinderMediatedSupportContinuationSegmentAuthoringPlan {
        outside_loop,
        segment_index: 0,
        request_indices: vec![0],
        half_edges,
        boundary_bridge_half_edges: vec![],
        cylinder_owned_half_edges: vec![],
        outside_owned_half_edges: vec![],
        support_owned_half_edges: vec![],
        unassigned_twin_half_edges: vec![],
        unsupported_owned_half_edges: vec![],
        context_owned_half_edges: vec![],
        unsupported_owner_faces: vec![],
        context_owner_faces: vec![],
        owner_repair_spans: vec![],
        requirement: CylinderMediatedSupportContinuationSegmentRequirement::BoundaryBridge,
    }
}

fn first_two_half_edges(brep: &BRep, loop_: crate::topology::LoopHandle) -> Vec<HalfEdgeHandle> {
    let first = brep.loops[loop_].half_edge;
    vec![first, brep.half_edges[first].next]
}

fn edges_for_half_edges(brep: &BRep, half_edges: &[HalfEdgeHandle]) -> Vec<EdgeHandle> {
    half_edges
        .iter()
        .map(|half_edge| brep.half_edges[*half_edge].edge)
        .collect()
}

fn add_boundary_polygon(
    brep: &mut BRep,
    shell: crate::topology::ShellHandle,
    points: &[Point3],
) -> crate::topology::FaceHandle {
    let face = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: vec![],
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;
    let mut half_edges = Vec::new();
    for point in points {
        let vertex = brep.vertices.insert(Vertex {
            position: *point,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let edge = brep.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        let half_edge = brep.half_edges.insert(HalfEdge {
            origin: vertex,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });
        brep.vertices[vertex].outgoing = half_edge;
        brep.edges[edge].half_edge = half_edge;
        brep.boundary_edges.insert(edge);
        half_edges.push(half_edge);
    }
    for index in 0..half_edges.len() {
        let current = half_edges[index];
        let next = half_edges[(index + 1) % half_edges.len()];
        let prev = half_edges[(index + half_edges.len() - 1) % half_edges.len()];
        brep.half_edges[current].next = next;
        brep.half_edges[current].prev = prev;
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}
