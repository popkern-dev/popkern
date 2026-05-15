use super::*;
use crate::topology::{
    BrepEdgeSurfaceVertexCandidateOutcome, BrepEdgeSurfaceVertexCandidateSource,
};

#[test]
fn vertex_star_accepts_noncanonical_trim_endpoint_candidates() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    let he_b = brep.half_edges[he_a].twin;
    let va = brep.half_edges[he_a].origin;
    let vb = brep.half_edges[he_b].origin;
    let start = Point3::new(0.25, 0.0, 0.0);
    let end = Point3::new(0.25, 1.0, 0.0);
    brep.vertices[va].position = start;
    brep.vertices[vb].position = end;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));
    brep.half_edges[he_b].pcurve = Some(Pcurve::linear((0.25, 1.0), (0.25, 0.0), 0.0, 1.0));
    brep.edges[edge].curve = CurveOnEdge::line_segment(
        start + Vec3::new(0.0, 0.0, 0.1),
        end + Vec3::new(0.0, 0.0, 0.1),
        1e-12,
    );

    let report = match repair_brep_edge_surface_compatibility_transactionally(
        &brep,
        no_validation_options(),
    ) {
        Ok(result) => result.report,
        Err(failure) => failure.report,
    };

    assert_eq!(report.vertex_stars.inspected_vertices, 2);
    assert_eq!(report.vertex_stars.coherent_vertices, 2);
    assert!(
        report
            .vertex_stars
            .decisions
            .iter()
            .all(|decision| decision.max_anchor_motion <= no_validation_options().tolerance)
    );
    assert!(report.vertex_stars.decisions.iter().all(|decision| {
        decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.source == BrepEdgeSurfaceVertexCandidateSource::TrimEndpoint
                && diagnostic.outcome == BrepEdgeSurfaceVertexCandidateOutcome::Accepted
        })
    }));
}
