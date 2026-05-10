use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, Shell, ShellHandle, Solid,
    Vertex, VertexHandle,
};

use super::super::repeated_apex_region_candidates::{CandidateStatus, RegionOperationKind};
use super::super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitiveDerivation;
use super::super::repeated_apex_side_support_plan::{
    SideSupportRegionPlanSelection, SideSupportRegionRefusalReason,
};
use super::super::repeated_apex_trim_subcomplex::{
    TrimRewriteSubcomplexBuilder, TrimRewriteSubcomplexPlan,
};
use super::*;

#[test]
fn read_only_snapshot_accepts_existing_support_loop_split_candidate() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let plan = trim_plan(|plan| {
        plan.push_new_generator_trim(
            0,
            cylinder_face,
            support_face,
            support_vertices[0],
            support_vertices[1],
        );
    });

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    assert!(snapshot.validate_fresh(&brep));
    assert_eq!(snapshot.candidates.len(), 1);
    let candidate = snapshot.candidates.iter().next().expect("candidate");
    assert_eq!(candidate.kind(), RegionOperationKind::SplitFaceWithEdge);
    assert!(matches!(candidate.status(), CandidateStatus::Accepted));
    assert!(matches!(
        snapshot.selection,
        SideSupportRegionPlanSelection::Selected(_)
    ));
    assert!(snapshot.diagnostic_label().contains("selected_plan"));
    assert!(snapshot.diagnostic_label().contains("ready_for_executor"));
    assert!(
        snapshot
            .diagnostic_label()
            .contains("total_delta=dV=0 dE=1 dF=1 dL=1")
    );
}

#[test]
fn executable_continuation_refuses_stale_knowledge_graph() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let plan = trim_plan(|plan| {
        plan.push_new_generator_trim(
            0,
            cylinder_face,
            support_face,
            support_vertices[0],
            support_vertices[1],
        );
    });
    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let mut stale_brep = brep.clone();
    stale_brep.vertices.insert(Vertex {
        position: Point3::new(10.0, 10.0, 10.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    assert_eq!(
        snapshot.apply_executable_continuation(&mut stale_brep),
        Err(SideSupportRegionContinuationFailure::StaleKnowledgeGraph)
    );
}

#[test]
fn read_only_snapshot_blocks_bridge_plan_that_leaves_support_family() {
    let (brep, plan) = shared_edge_bridge_fixture(false);

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let labels = snapshot.diagnostic_label();
    assert!(labels.contains("outside-support-family"));
    assert!(labels.contains("bridge-support-region"));
    assert!(labels.contains("merge-adjacent-support-region"));
    assert!(labels.contains("total_delta=dV=0 dE=2 dF=2 dL=2"));
    assert!(labels.contains("dV=0 dE=2 dF=2 dL=2"));
    assert!(labels.contains("one_endpoint_on_support_loop=true"));
    assert!(labels.contains("MissingValidatedFaceFamilyBridgePrimitive"));
    assert!(labels.contains("face-family-bridge primitive_ready=false"));
    assert!(labels.contains("shared-bridge-operand-edges"), "{labels}");
    assert!(labels.contains("support-region-merge primitive_ready=false"));
    assert!(labels.contains("core-merge-refused="), "{labels}");
    assert!(labels.contains("MissingValidatedSupportRegionMergePrimitive"));
    assert!(labels.contains("MissingValidatedFaceFamilyBridgePrimitive"));
    assert!(labels.contains("execution refused stage=Preflight"));
    assert_eq!(snapshot.candidates.len(), 6);
    assert!(snapshot.candidates.iter().any(|candidate| {
        candidate.kind() == RegionOperationKind::BridgeSupportRegion
            && matches!(candidate.status(), CandidateStatus::Accepted)
    }));
    assert!(snapshot.candidates.iter().any(|candidate| {
        candidate.kind() == RegionOperationKind::MergeAdjacentSupportRegion
            && matches!(candidate.status(), CandidateStatus::Accepted)
    }));
    match &snapshot.selection {
        SideSupportRegionPlanSelection::Blocked(_, refusal) => assert_eq!(
            refusal.reason(),
            SideSupportRegionRefusalReason::MissingValidatedFaceFamilyBridgePrimitive
        ),
        selection => panic!("expected blocked bridge selection, got {selection:?}"),
    }
}

#[test]
fn read_only_snapshot_selects_merge_when_shared_edge_is_mergeable() {
    let (mut brep, plan) = shared_edge_bridge_fixture(true);
    let before_counts = (
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
        brep.loops.len(),
    );

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let labels = snapshot.diagnostic_label();
    assert!(labels.contains("shared-bridge-operand-edges"), "{labels}");
    assert!(
        labels.contains("support-region-merge primitive_ready=true"),
        "{labels}"
    );
    assert!(labels.contains("core_merge_probe=ready"), "{labels}");
    assert!(
        labels.contains("post_merge_side_trims_ready=true"),
        "{labels}"
    );
    assert!(
        labels.contains("combined_delta=dV=0 dE=1 dF=1 dL=1"),
        "{labels}"
    );
    assert!(labels.contains("exact_curve_authoring=refused"), "{labels}");
    assert!(
        labels.contains("MissingExactSideGeneratorCurveAuthoring"),
        "{labels}"
    );
    assert_eq!(
        apply_selected_side_support_continuation(&mut brep, &plan)
            .expect("not-ready continuation is a valid non-mutating outcome"),
        SideSupportRegionContinuationOutcome::NotReady
    );
    assert_eq!(
        (
            brep.vertices.len(),
            brep.edges.len(),
            brep.faces.len(),
            brep.loops.len()
        ),
        before_counts
    );
}

#[test]
fn read_only_snapshot_ready_for_boundary_bridge_plan() {
    let brep = two_boundary_quads();
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let outside_face = faces[1];
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let outside_vertices = outer_loop_vertices(&brep, outside_face);
    let plan = trim_plan(|plan| {
        plan.push_new_generator_trim(
            0,
            outside_face,
            support_face,
            support_vertices[0],
            outside_vertices[0],
        );
        plan.push_new_generator_trim(
            1,
            outside_face,
            support_face,
            outside_vertices[2],
            support_vertices[2],
        );
    });

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let labels = snapshot.diagnostic_label();
    let bridge = FaceFamilyBridgePrimitiveDerivation::build(&snapshot.selection, &brep);
    let topology_report = bridge.topology_operator_report(&brep);
    let topology_label = topology_report.diagnostic_label();
    assert!(bridge.is_complete(), "{labels}");
    assert!(topology_report.is_ready(), "{topology_label}");
    assert!(
        bridge
            .diagnostic_label()
            .contains("boundary_readiness=boundary-ready"),
        "{}",
        bridge.diagnostic_label()
    );
    assert!(labels.contains("execution ready_for_executor"), "{labels}");

    let mut applied = brep.clone();
    bridge
        .apply_topology_operator(&mut applied)
        .expect("boundary-ready bridge applies");
    assert_eq!(applied.edges.len(), brep.edges.len() + 2);
    assert_eq!(applied.faces.len(), brep.faces.len() + 2);
    assert_eq!(applied.loops.len(), brep.loops.len() + 2);
    assert!(applied.boundary_edges.is_empty());
}

#[test]
fn read_only_snapshot_refuses_bridge_region_with_one_trim() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let cap_face = face_with_vertex_outside_support_loop(&brep, support_face);
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let cap_vertices = outer_loop_vertices(&brep, cap_face);
    let outside_vertex = cap_vertices
        .into_iter()
        .find(|vertex| !support_vertices.contains(vertex))
        .expect("fixture has a cap vertex outside support loop");
    let plan = trim_plan(|plan| {
        plan.push_new_generator_trim(
            0,
            cylinder_face,
            support_face,
            support_vertices[0],
            outside_vertex,
        );
    });

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let labels = snapshot.diagnostic_label();
    assert!(labels.contains("UnsupportedBridgeTrimCount"));
    assert!(!labels.contains("MissingValidatedFaceFamilyBridgePrimitive"));
    assert_eq!(snapshot.candidates.len(), 3);
    match &snapshot.selection {
        SideSupportRegionPlanSelection::Refused(refusal) => assert_eq!(
            refusal.reason(),
            SideSupportRegionRefusalReason::UnsupportedBridgeTrimCount
        ),
        selection => panic!("expected refused bridge selection, got {selection:?}"),
    }
}

#[test]
fn executable_continuation_refuses_outside_outside_plan_before_mutation() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let cap_face = face_with_vertex_outside_support_loop(&brep, support_face);
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let outside_vertices = outer_loop_vertices(&brep, cap_face)
        .into_iter()
        .filter(|vertex| !support_vertices.contains(vertex))
        .collect::<Vec<_>>();
    assert!(
        outside_vertices.len() >= 2,
        "fixture has at least two vertices outside support loop"
    );
    let plan = trim_plan(|plan| {
        plan.push_new_generator_trim(
            0,
            cylinder_face,
            support_face,
            outside_vertices[0],
            outside_vertices[1],
        );
    });
    let mut candidate = brep.clone();

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);
    let labels = snapshot.diagnostic_label();

    assert!(
        labels.contains("outside-outside-support-continuation"),
        "{labels}"
    );
    assert!(
        labels.contains("MissingValidatedOutsideOutsideSupportContinuationPrimitive"),
        "{labels}"
    );
    assert!(
        labels.contains("outside-outside-support primitive_ready=false"),
        "{labels}"
    );
    assert!(labels.contains("class=cylinder-mediated"), "{labels}");
    assert!(labels.contains("cylinder_intervals=["), "{labels}");
    assert!(
        labels.contains(
            "delta=unknown reason=outside-outside-support-continuation-operator-required"
        ),
        "{labels}"
    );
    let refusal = apply_selected_side_support_continuation(&mut candidate, &plan)
        .expect_err("outside/outside continuation must refuse before mutation");
    let ApexBuildError::PlanMismatch { reason } = refusal else {
        panic!("expected PlanMismatch refusal");
    };
    assert!(
        reason.contains(
            "planner-refused stage=SelectPlan \
                 reason=MissingValidatedOutsideOutsideSupportContinuationPrimitive"
        ),
        "{reason}"
    );
    assert!(reason.contains(&labels), "{reason}");
    assert_eq!(candidate.vertices.len(), brep.vertices.len());
    assert_eq!(candidate.edges.len(), brep.edges.len());
    assert_eq!(candidate.faces.len(), brep.faces.len());
    assert_eq!(candidate.loops.len(), brep.loops.len());
}

fn trim_plan(build: impl FnOnce(&mut TrimRewriteSubcomplexBuilder)) -> TrimRewriteSubcomplexPlan {
    let mut builder = TrimRewriteSubcomplexBuilder::new();
    build(&mut builder);
    builder.finish()
}

fn shared_edge_bridge_fixture(clear_surfaces: bool) -> (BRep, TrimRewriteSubcomplexPlan) {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    if clear_surfaces {
        for face in brep.faces.values_mut() {
            face.surface = None;
        }
    }
    let faces = brep.faces.keys().collect::<Vec<_>>();
    let support_face = faces[0];
    let cylinder_face = faces[1];
    let cap_face = face_with_vertex_outside_support_loop(&brep, support_face);
    let support_vertices = outer_loop_vertices(&brep, support_face);
    let cap_vertices = outer_loop_vertices(&brep, cap_face);
    let outside_vertices = cap_vertices
        .into_iter()
        .filter(|vertex| !support_vertices.contains(vertex))
        .collect::<Vec<_>>();
    assert!(
        outside_vertices.len() >= 2,
        "fixture has at least two cap vertices outside support loop"
    );
    let outside_start = outside_vertices[0];
    let outside_end = outside_vertices[1];
    let plan = trim_plan(|plan| {
        plan.push_existing_seam(0, cylinder_face, cap_face, outside_start, outside_end);
        plan.push_new_generator_trim(
            0,
            cylinder_face,
            support_face,
            support_vertices[0],
            outside_start,
        );
        plan.push_new_generator_trim(
            1,
            cylinder_face,
            support_face,
            outside_end,
            support_vertices[1],
        );
    });
    (brep, plan)
}

fn face_with_vertex_outside_support_loop(brep: &BRep, support_face: FaceHandle) -> FaceHandle {
    let support_vertices = outer_loop_vertices(brep, support_face);
    brep.faces
        .keys()
        .find(|face| {
            *face != support_face
                && outer_loop_vertices(brep, *face)
                    .iter()
                    .any(|vertex| !support_vertices.contains(vertex))
        })
        .expect("fixture has a face with a vertex outside support loop")
}

fn outer_loop_vertices(brep: &BRep, face: FaceHandle) -> Vec<VertexHandle> {
    let loop_handle = brep.faces[face].outer_loop;
    let start = brep.loops[loop_handle].half_edge;
    let mut vertices = Vec::new();
    let mut cursor = start;
    loop {
        let half_edge = &brep.half_edges[cursor];
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    vertices
}

fn two_boundary_quads() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    let first_face = add_boundary_quad(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    add_boundary_quad(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
    );
    brep.shells[shell].face = first_face;
    brep
}

fn add_boundary_quad(brep: &mut BRep, shell: ShellHandle, positions: [Point3; 4]) -> FaceHandle {
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
    let vertices = positions.map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let half_edges = vertices.map(|vertex| {
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
            twin: HalfEdgeHandle::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });
        brep.edges[edge].half_edge = half_edge;
        brep.vertices[vertex].outgoing = half_edge;
        brep.boundary_edges.insert(edge);
        half_edge
    });
    for i in 0..4 {
        brep.half_edges[half_edges[i]].next = half_edges[(i + 1) % 4];
        brep.half_edges[half_edges[i]].prev = half_edges[(i + 3) % 4];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}
