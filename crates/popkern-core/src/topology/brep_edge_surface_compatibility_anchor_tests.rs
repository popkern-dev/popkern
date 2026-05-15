use super::*;

#[test]
fn repair_declines_reanchor_when_vertex_motion_exceeds_policy() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    let half_edge = brep.edges[edge].half_edge;
    let start_vertex = brep.half_edges[half_edge].origin;
    let end_vertex = brep.half_edges[brep.half_edges[half_edge].next].origin;
    let offset = Vec3::new(0.0, 0.0, 0.02);
    brep.vertices[start_vertex].position = brep.vertices[start_vertex].position + offset;
    brep.vertices[end_vertex].position = brep.vertices[end_vertex].position + offset;
    corrupt_edge_curve(&mut brep, edge, offset);
    let mut options = options();
    options.max_anchor_motion = 0.005;

    let failure =
        repair_brep_edge_surface_compatibility_transactionally(&brep, options).unwrap_err();

    assert_eq!(failure.report.repaired_edges, 0);
    assert_eq!(failure.report.declined_edges, 1);
    assert!(matches!(
        failure.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
            BrepEdgeSurfaceCompatibilityDeclineReason::AnchorMotionExceeded
        )
    ));
    assert!(failure.report.decisions[0].evidence.attempted_anchor_repair);
    assert!(failure.report.decisions[0].evidence.anchor_motion > options.max_anchor_motion);
    assert_eq!(
        failure.report.decisions[0]
            .evidence
            .max_allowed_anchor_motion,
        options.max_anchor_motion
    );
    assert_eq!(failure.report.vertex_stars.inspected_vertices, 2);
    assert_eq!(failure.report.vertex_stars.coherent_vertices, 0);
    assert_eq!(failure.report.vertex_stars.declined_vertices, 2);
    assert!(
        failure
            .report
            .vertex_stars
            .decisions
            .iter()
            .all(|decision| {
                matches!(
            decision.outcome,
            crate::topology::BrepEdgeSurfaceVertexStarOutcome::Declined(
                crate::topology::BrepEdgeSurfaceVertexStarDeclineReason::AnchorMotionExceeded
            )
        )
            })
    );
}
