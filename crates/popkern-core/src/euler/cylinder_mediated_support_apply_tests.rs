#[path = "cylinder_mediated_support_apply_test_fixtures.rs"]
mod fixtures;

use fixtures::{TopologySnapshot, long_interval_cylinder_mediated_fixture};

#[test]
fn mutator_boundary_commits_multi_run_connector_authoring() {
    let fixture = long_interval_cylinder_mediated_fixture();
    let mut brep = fixture.brep.clone();
    let before = TopologySnapshot::capture(&brep);
    let batch = brep
        .preflight_cylinder_mediated_support_continuation_batch(vec![
            fixture.first_request(),
            fixture.second_request(),
        ])
        .expect("batch preflight validates independently valid requests");

    let segment_authoring = brep.segment_authoring_plan(&batch).unwrap();
    assert!(
        segment_authoring
            .iter()
            .all(|segment| !segment.source_fact_ids.is_empty()),
        "side-support segment authoring should cite shared BRep repair facts"
    );

    let result = brep
        .cylinder_mediated_support_continuation_batch(&batch)
        .expect("multi-run connector authoring commits");

    assert_eq!(result.bridge_faces.len(), 3);
    assert_eq!(result.bridge_loops.len(), 3);
    assert_eq!(
        result.postconditions,
        super::CylinderMediatedSupportContinuationPostconditions::Satisfied
    );
    brep.validate_half_edge_pairing()
        .expect("pairing validates");
    brep.validate_half_edge_cycles().expect("cycles validate");
    brep.validate_manifold_edges().expect("manifold validates");
    let delta = before.delta_to(TopologySnapshot::capture(&brep));
    assert_eq!(delta.vertices, 0);
    assert_eq!(delta.edges, 13);
    assert_eq!(delta.faces, 3);
    assert_eq!(delta.loops, 3);
    assert_eq!(delta.shells, 0);
    assert_eq!(delta.solids, 0);
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
