use super::*;
use crate::euler::DiagnosticHandle;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

#[test]
fn live_pyramid_apex_reports_high_valence_patch_candidate() {
    let pyramid = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, vertex)| vertex.position.z > 0.0)
        .map(|(handle, _)| handle)
        .expect("apex");
    let edges = incident_edges(&pyramid, apex);
    assert_eq!(edges.len(), 4);

    let report = pyramid.diagnose_corner_apex(request(apex, &edges));

    assert_eq!(
        report.classification.valence,
        BlendValence::HighValenceCorner
    );
    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::HighValencePatch
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("n-sided-corner-patch")
    );
    assert_eq!(
        report.strategy.patch_candidate.applicability,
        PatchCandidateApplicability::Applicable
    );
    assert_eq!(
        report.strategy.patch_candidate.n_sided_feasibility,
        PatchCandidateFeasibility::MissingPatchPrimitive("n-sided-corner-patch")
    );
    assert_eq!(
        report.strategy.patch_candidate.boundary_compatibility,
        PatchCandidateFeasibility::MissingBoundaryPrimitive("patch-boundary-compatibility")
    );
    assert_eq!(
        report.strategy.patch_candidate.quality,
        PatchCandidateQuality::ManagedRisk("high-valence-patch-quality")
    );
    assert_eq!(
        report.strategy.patch_candidate.downgrade_policy,
        PatchCandidateDowngradePolicy::DowngradeAllowed("G1-continuity")
    );
}

fn request(vertex: VertexHandle, edges: &[EdgeHandle]) -> CornerApexDiagnosticRequest {
    CornerApexDiagnosticRequest {
        scope: LocalCornerScope {
            vertex: DiagnosticHandle::Vertex(vertex),
        },
        selected_edges: edges
            .iter()
            .map(|edge| LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(*edge),
                blend_kind: BlendKind::FilletConstantRadius,
                radius_nanometers: Some(10_000_000),
            })
            .collect(),
        selection_intent: "unit-test",
        operation_context: "live-corner-apex-diagnostics",
    }
}

fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges = std::collections::BTreeSet::new();
    for (_, half_edge) in &brep.half_edges {
        if half_edge.origin == vertex
            || brep
                .half_edges
                .get(half_edge.twin)
                .is_some_and(|twin| twin.origin == vertex)
        {
            edges.insert(half_edge.edge);
        }
    }
    edges.into_iter().collect()
}
