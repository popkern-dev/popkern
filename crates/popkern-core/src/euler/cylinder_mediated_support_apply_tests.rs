use crate::constructors::builder::BRepBuilder;
use crate::euler::CylinderMediatedSupportContinuationRequest;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

#[test]
fn mutator_boundary_refuses_without_mutating_topology() {
    let fixture = long_interval_cylinder_mediated_fixture();
    let mut brep = fixture.brep.clone();
    let before = TopologySnapshot::capture(&brep);
    let batch = brep
        .preflight_cylinder_mediated_support_continuation_batch(vec![
            fixture.first_request(),
            fixture.second_request(),
        ])
        .expect("batch preflight validates independently valid requests");

    let err = brep
        .cylinder_mediated_support_continuation_batch(&batch)
        .expect_err("mutator boundary refuses until topology rewiring is implemented");

    let detail = err.to_string();
    assert!(
        detail.contains(
            "cylinder-mediated support continuation split cylinder resource merge construction not implemented"
        ),
        "{err}"
    );
    assert!(detail.contains("segment_authoring=segments=3"), "{err}");
    assert!(
        detail.contains("authored_delta=dV=0 dE=8 dF=3 dL=3"),
        "{err}"
    );
    assert!(
        detail.contains("execution_delta=dV=0 dE=10 dF=3 dL=3"),
        "{err}"
    );
    assert!(
        detail.contains("topology_delta_gap=dV=0 dE=4 dF=0 dL=0"),
        "{err}"
    );
    assert!(detail.contains("boundary_compatible=1"), "{err}");
    assert!(
        detail.contains("owned_reclassification_segments=2"),
        "{err}"
    );
    assert!(detail.contains("unassigned_repair_segments=0"), "{err}");
    assert!(detail.contains("cylinder_owned_half_edges=2"), "{err}");
    assert!(detail.contains("outside_owned_half_edges=0"), "{err}");
    assert!(detail.contains("support_owned_half_edges=0"), "{err}");
    assert!(detail.contains("context_owned_half_edges=0"), "{err}");
    assert!(detail.contains("unassigned_twin_half_edges=0"), "{err}");
    assert!(detail.contains("unsupported_owned_half_edges=0"), "{err}");
    assert!(detail.contains("owner_repair_spans=2"), "{err}");
    assert!(detail.contains("owner_repair_half_edges=2"), "{err}");
    assert!(
        detail.contains("owner_repair_candidate=exposed_spans=2 exposed_half_edges=2"),
        "{err}"
    );
    assert!(
        detail.contains(
            "segment_bridge_candidate=segments=3 boundary_half_edges=15 \
             ownership=shared-support-transition-ready \
             shared_support_segments=1 missing_support_segments=0 \
             max_support_faces=2 degenerate_endpoint_segments=0 max_cylinder_faces=1 \
             opposite_cylinder_half_edges=7 duplicate_opposite_cylinder_half_edges=0 \
             missing_opposite_path_segments=1"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "cylinder_resource_partition=segments=2 half_edges=7 \
             shared_segments=0 max_request_owners=1"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "cylinder_resource_blocker_exposure=spans=2 half_edges=2 \
             exposed_spans=2 exposed_half_edges=2"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "cylinder_resource_remap=segments=4 half_edges=7 reversed_segments=0 \
             unavailable_edges=0 unavailable_not_boundary=0 unavailable_live_twin=0 \
             unavailable_desired_live=0 unavailable_desired_retired=0 \
             unavailable_owner_roles=cylinder:0 support:0 context:0 \
             unavailable_blocking_roles=cylinder:0 support:0 context:0 \
             blocking_spans=0 blocking_half_edges=0 non_contiguous_blocking_spans=0 \
             unavailable_owner_faces=[] non_contiguous_segments=0"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "bridge_resource_reconciliation=outside_segments=3 cylinder_segments=4 \
             exact_matches=2 missing_cylinder_matches=1 extra_cylinder_resources=2 \
             duplicate_cylinder_request_sets=2 resource_status=split-cylinder-resources-need-merge \
             construction_actions=exact_pair:2 shared_outside:1 single_covered:0 \
             unsupported:0 \
             action_sequence=[0:exact:0:0,1:shared-outside:-:0+1,2:exact:1:1] \
             outside_request_sets=[0,0+1,1] cylinder_request_sets=[0,1,1,0] \
             missing_request_sets=[0+1] \
             split_resource_merge_candidates=missing_sets=1 coverable_sets=1 cylinder_runs=4"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "split_resource_merge_operands=[outside=1 requests=0+1 \
             endpoint_shape="
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "split_resource_endpoint_shapes=shared_start:0 shared_end:0 aligned:0 mixed:1"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "cylinder_resource_construction_candidate=family=split-cylinder-resource-merge \
             feasibility=requires-resource-merge reason=split-cylinder-resource-merge-required \
             action_candidates=exact_pair:2 shared_outside:1 single_covered:0 unsupported:0"
        ),
        "{err}"
    );
    assert!(
        detail.contains(
            "shared_start_bridge_preflight=status=not-applicable pairings=0 cycles=0 \
             cross_loop_cycles=0 same_loop_cycles=0 one_edge_interior_overlap_candidates=0 \
             repeated_one_edge_overlap_occurrences=0 one_edge_interior_overlap_plans=0 \
             one_edge_interior_overlap_boundary_edges=0 \
             one_edge_interior_overlap_normalizations=0 \
             bridge_diagnostics=operations=0 resources=0 conflicts=0 readiness=0 construction_candidates=0 postconditions=0 \
             conflicts_by_kind=duplicate_half_edge:0 duplicate_boundary_edge:0 shared_start:0 \
             shared_end:0 crossed_endpoint:0 aligned:0 nested:0 interior:0 repeated_edge:0 \
             already_boundary:0 paired_exposure:0 unsupported_owner:0 stale_handle:0 unknown:0 \
             readiness_by_status=ready:0 ready_after_normalization:0 missing_operator:0 \
             invalid_input:0 unknown:0 \
             candidate_families=disjoint_cross_loop:0 disjoint_same_loop:0 \
             paired_same_loop_twin_exposure:0 boundary_side_separator:0 unknown:0 \
             candidate_feasibility=executable_existing:0 executable_after_normalization:0 \
             requires_separator_subcomplex:0 requires_custom_aggregate_authoring:0 \
             requires_cycle_regrouping:0 invalid_non_manifold:0 unsupported_unknown_family:0 unknown:0 \
             candidate_requirements=live_revalidation:0 postcondition_verification:0 unknown:0 \
             postcondition_preview_candidates=0 predicted_delta=dV=0 dE=0 dF=0 dL=0 \
             expected_delta=dV=0 dE=0 dF=0 dL=0 reason=none"
        ),
        "{err}"
    );
    assert!(detail.contains("outside_endpoints="), "{err}");
    assert!(detail.contains("context_owner_faces=[]"), "{err}");
    assert!(detail.contains("unsupported_owner_faces=[]"), "{err}");
    assert_eq!(TopologySnapshot::capture(&brep), before);
}

#[test]
fn mutator_boundary_refuses_stale_plan_without_mutating_topology() {
    let fixture = long_interval_cylinder_mediated_fixture();
    let mut brep = fixture.brep.clone();
    let before = TopologySnapshot::capture(&brep);
    let mut batch = brep
        .preflight_cylinder_mediated_support_continuation_batch(vec![
            fixture.first_request(),
            fixture.second_request(),
        ])
        .expect("batch preflight validates independently valid requests");
    batch.plans[0].outside_interval.pop();

    let err = brep
        .cylinder_mediated_support_continuation_batch(&batch)
        .expect_err("mutator boundary revalidates stale plans");

    assert!(err.to_string().contains("plan is stale"), "{err}");
    assert_eq!(TopologySnapshot::capture(&brep), before);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TopologySnapshot {
    vertices: usize,
    edges: usize,
    half_edges: usize,
    loops: usize,
    faces: usize,
    shells: usize,
    solids: usize,
    boundary_edges: usize,
}

impl TopologySnapshot {
    fn capture(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            shells: brep.shells.len(),
            solids: brep.solids.len(),
            boundary_edges: brep.boundary_edges.len(),
        }
    }
}

struct LongIntervalCylinderMediatedFixture {
    brep: BRep,
    outside_face: FaceHandle,
    cylinder_face: FaceHandle,
    outside_loop: LoopHandle,
    outside_half_edges: Vec<HalfEdgeHandle>,
    support_faces: [FaceHandle; 2],
    support_loops: [LoopHandle; 2],
    support_cylinder_half_edges: [HalfEdgeHandle; 2],
    outside_cylinder_half_edges: [HalfEdgeHandle; 2],
    cylinder_intervals: [Vec<HalfEdgeHandle>; 2],
}

impl LongIntervalCylinderMediatedFixture {
    fn first_request(&self) -> CylinderMediatedSupportContinuationRequest {
        self.request(0, 0, 12)
    }

    fn second_request(&self) -> CylinderMediatedSupportContinuationRequest {
        self.request(1, 5, 15)
    }

    fn request(
        &self,
        support_index: usize,
        start_index: usize,
        end_index: usize,
    ) -> CylinderMediatedSupportContinuationRequest {
        CylinderMediatedSupportContinuationRequest {
            support_face: self.support_faces[support_index],
            outside_face: self.outside_face,
            cylinder_face: self.cylinder_face,
            start: self.brep.half_edges[self.outside_half_edges[start_index]].origin,
            end: self.brep.half_edges[self.outside_half_edges[end_index]].origin,
            outside_loop: self.outside_loop,
            outside_start_half_edge: self.outside_half_edges[start_index],
            outside_end_half_edge: self.outside_half_edges[end_index],
            support_loop: self.support_loops[support_index],
            support_cylinder_half_edge: self.support_cylinder_half_edges[support_index],
            outside_cylinder_half_edge: self.outside_cylinder_half_edges[support_index],
            cylinder_interval: self.cylinder_intervals[support_index].clone(),
            context_owner_faces: Vec::new(),
        }
    }
}

fn long_interval_cylinder_mediated_fixture() -> LongIntervalCylinderMediatedFixture {
    let mut builder = BRepBuilder::new(Tolerance::default());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    let outside_face = builder.add_face(shell);
    let cylinder_face = builder.add_face(shell);
    let first_support_face = builder.add_face(shell);
    let second_support_face = builder.add_face(shell);

    let outside_vertices = (0..20)
        .map(|index| builder.add_vertex(Point3::new(index as f64, 0.0, 0.0)))
        .collect::<Vec<_>>();
    let mut outside_half_edges = Vec::new();
    let mut outside_twins = Vec::new();
    for index in 0..outside_vertices.len() {
        let next = (index + 1) % outside_vertices.len();
        let (_, half_edge, twin) =
            builder.add_edge(outside_vertices[index], outside_vertices[next]);
        outside_half_edges.push(half_edge);
        outside_twins.push(twin);
    }
    let outside_loop = builder
        .add_loop(&outside_half_edges, outside_face)
        .expect("outside loop wires");

    let (first_support_loop, first_support_cylinder_half_edge, first_support_cylinder_twin) =
        add_support_square(&mut builder, first_support_face, 100.0);
    let (second_support_loop, second_support_cylinder_half_edge, second_support_cylinder_twin) =
        add_support_square(&mut builder, second_support_face, 200.0);

    let cylinder_connectors = (0..8)
        .map(|index| {
            let v0 = builder.add_vertex(Point3::new(300.0 + index as f64, 0.0, 0.0));
            let v1 = builder.add_vertex(Point3::new(300.5 + index as f64, 0.0, 0.0));
            let (_, half_edge, _) = builder.add_edge(v0, v1);
            half_edge
        })
        .collect::<Vec<_>>();
    builder
        .add_loop(
            &[
                outside_twins[3],
                cylinder_connectors[0],
                cylinder_connectors[1],
                cylinder_connectors[2],
                first_support_cylinder_twin,
                cylinder_connectors[3],
                outside_twins[7],
                second_support_cylinder_twin,
                cylinder_connectors[4],
                cylinder_connectors[5],
                cylinder_connectors[6],
                cylinder_connectors[7],
            ],
            cylinder_face,
        )
        .expect("cylinder loop wires");

    LongIntervalCylinderMediatedFixture {
        brep: boundary_normalized_test_brep(builder.build_unchecked()),
        outside_face,
        cylinder_face,
        outside_loop,
        outside_half_edges: outside_half_edges.clone(),
        support_faces: [first_support_face, second_support_face],
        support_loops: [first_support_loop, second_support_loop],
        support_cylinder_half_edges: [
            first_support_cylinder_half_edge,
            second_support_cylinder_half_edge,
        ],
        outside_cylinder_half_edges: [outside_half_edges[3], outside_half_edges[7]],
        cylinder_intervals: [
            vec![
                outside_twins[3],
                cylinder_connectors[0],
                cylinder_connectors[1],
                cylinder_connectors[2],
                first_support_cylinder_twin,
            ],
            vec![outside_twins[7], second_support_cylinder_twin],
        ],
    }
}

fn boundary_normalized_test_brep(mut brep: BRep) -> BRep {
    let edges = brep.edges.keys().collect::<Vec<_>>();
    let mut removed_half_edges = Vec::new();
    for edge in edges {
        let incident_half_edges = brep
            .half_edges
            .iter()
            .filter_map(|(half_edge, data)| (data.edge == edge).then_some(half_edge))
            .collect::<Vec<_>>();
        if incident_half_edges.len() != 2 {
            continue;
        }
        let live_half_edges = incident_half_edges
            .iter()
            .copied()
            .filter(|half_edge| brep.loops.contains_key(brep.half_edges[*half_edge].loop_))
            .collect::<Vec<_>>();
        if live_half_edges.len() != 1 {
            continue;
        }
        let boundary_half_edge = live_half_edges[0];
        let residual_half_edge = incident_half_edges
            .into_iter()
            .find(|half_edge| *half_edge != boundary_half_edge)
            .expect("two incident half-edges include one residual");
        brep.half_edges[boundary_half_edge].twin = HalfEdgeHandle::default();
        brep.edges[edge].half_edge = boundary_half_edge;
        brep.boundary_edges.insert(edge);
        brep.half_edges.remove(residual_half_edge);
        removed_half_edges.push(residual_half_edge);
    }
    repair_test_vertex_outgoing_handles(&mut brep, &removed_half_edges);
    brep
}

fn repair_test_vertex_outgoing_handles(brep: &mut BRep, removed_half_edges: &[HalfEdgeHandle]) {
    let vertices = brep.vertices.keys().collect::<Vec<_>>();
    for vertex in vertices {
        if !removed_half_edges.contains(&brep.vertices[vertex].outgoing) {
            continue;
        }
        brep.vertices[vertex].outgoing = brep
            .half_edges
            .iter()
            .find_map(|(half_edge, data)| (data.origin == vertex).then_some(half_edge))
            .unwrap_or_default();
    }
}

fn add_support_square(
    builder: &mut BRepBuilder,
    support_face: FaceHandle,
    offset: f64,
) -> (LoopHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let vertices = [
        builder.add_vertex(Point3::new(offset, 0.0, 0.0)),
        builder.add_vertex(Point3::new(offset + 1.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(offset + 1.0, 1.0, 0.0)),
        builder.add_vertex(Point3::new(offset, 1.0, 0.0)),
    ];
    let (_, support_cylinder_half_edge, support_cylinder_twin) =
        builder.add_edge(vertices[0], vertices[1]);
    let (_, side_1, _) = builder.add_edge(vertices[1], vertices[2]);
    let (_, side_2, _) = builder.add_edge(vertices[2], vertices[3]);
    let (_, side_3, _) = builder.add_edge(vertices[3], vertices[0]);
    let support_loop = builder
        .add_loop(
            &[support_cylinder_half_edge, side_1, side_2, side_3],
            support_face,
        )
        .expect("support loop wires");
    (
        support_loop,
        support_cylinder_half_edge,
        support_cylinder_twin,
    )
}
