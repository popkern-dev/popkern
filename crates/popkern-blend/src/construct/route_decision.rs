use popkern_core::topology::{FaceHandle, VertexHandle};

use super::EdgeGeometry;

pub(super) struct ApexRouteDecisionInput<'a> {
    pub vertex: VertexHandle,
    pub edge_indices: &'a [usize],
    pub edge_geometries: &'a [EdgeGeometry],
    pub route: crate::types::ApexConstructionRoute,
    pub built_faces: Vec<FaceHandle>,
    pub corner_patch_count: usize,
    pub exact_boundary_sector_count: usize,
    pub exact_boundary_prepared_curve_count: usize,
    pub exact_boundary_inner_spoke_count: usize,
    pub exact_boundary_ribbon_count: usize,
    pub exact_boundary_solved_surface_count: usize,
    pub exact_boundary_pcurve_pair_count: usize,
    pub exact_boundary_performance: crate::types::ExactBoundaryCandidatePerformanceReport,
    pub exact_boundary_report: crate::types::ExactBoundaryCapReport,
    pub fallback_reason: Option<String>,
}

pub(super) fn apex_route_decision(
    input: ApexRouteDecisionInput<'_>,
) -> crate::types::ApexConstructionRouteDecision {
    crate::types::ApexConstructionRouteDecision {
        vertex: input.vertex,
        incident_edges: input
            .edge_indices
            .iter()
            .filter_map(|&index| {
                input
                    .edge_geometries
                    .get(index)
                    .map(|geometry| geometry.analysis.edge)
            })
            .collect(),
        route: input.route,
        built_faces: input.built_faces,
        corner_patch_count: input.corner_patch_count,
        exact_boundary_sector_count: input.exact_boundary_sector_count,
        exact_boundary_prepared_curve_count: input.exact_boundary_prepared_curve_count,
        exact_boundary_inner_spoke_count: input.exact_boundary_inner_spoke_count,
        exact_boundary_ribbon_count: input.exact_boundary_ribbon_count,
        exact_boundary_solved_surface_count: input.exact_boundary_solved_surface_count,
        exact_boundary_pcurve_pair_count: input.exact_boundary_pcurve_pair_count,
        exact_boundary_performance: input.exact_boundary_performance,
        exact_boundary_report: input.exact_boundary_report,
        fallback_reason: input.fallback_reason,
    }
}

/// Route record for a corner dispatched straight to the vertex-blend path
/// (no direct-apex attempt was made).
pub(super) fn vertex_blend_dispatch_route(
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
) -> crate::types::ApexConstructionRouteDecision {
    apex_route_decision(ApexRouteDecisionInput {
        vertex,
        edge_indices,
        edge_geometries,
        route: crate::types::ApexConstructionRoute::VertexBlend,
        built_faces: Vec::new(),
        corner_patch_count: 0,
        exact_boundary_sector_count: 0,
        exact_boundary_prepared_curve_count: 0,
        exact_boundary_inner_spoke_count: 0,
        exact_boundary_ribbon_count: 0,
        exact_boundary_solved_surface_count: 0,
        exact_boundary_pcurve_pair_count: 0,
        exact_boundary_performance: Default::default(),
        exact_boundary_report: crate::types::ExactBoundaryCapReport::not_applicable("vertex_blend"),
        fallback_reason: None,
    })
}

/// Route record for a direct-apex attempt that fell back to the
/// vertex-blend path, carrying the decline reason.
pub(super) fn direct_apex_fallback_route(
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    reason: String,
) -> crate::types::ApexConstructionRouteDecision {
    apex_route_decision(ApexRouteDecisionInput {
        vertex,
        edge_indices,
        edge_geometries,
        route: crate::types::ApexConstructionRoute::DirectApexFallbackToVertexBlend,
        built_faces: Vec::new(),
        corner_patch_count: 0,
        exact_boundary_sector_count: 0,
        exact_boundary_prepared_curve_count: 0,
        exact_boundary_inner_spoke_count: 0,
        exact_boundary_ribbon_count: 0,
        exact_boundary_solved_surface_count: 0,
        exact_boundary_pcurve_pair_count: 0,
        exact_boundary_performance: Default::default(),
        exact_boundary_report: crate::types::ExactBoundaryCapReport::not_applicable(
            "direct_apex_fallback_to_vertex_blend",
        ),
        fallback_reason: Some(reason),
    })
}

/// PK-258: a completed vertex-blend corner serve records its cap faces on
/// the vertex's route (so the serving-retry keep-better contract can
/// measure BOTH candidates uniformly) and, when the corner FELL BACK
/// from direct-apex, promotes a cap/support continuity breach of the
/// commercial bars into the shared diagnostics surface.
pub(super) fn record_vertex_blend_serve_on_route(
    result: &popkern_core::topology::BRep,
    repair_report: &mut crate::types::BlendRepairReport,
    vertex: VertexHandle,
    cap_faces: Vec<FaceHandle>,
    fell_back_from_direct_apex: bool,
) {
    if fell_back_from_direct_apex
        && let Some(diagnostic) = super::cap_outer_acceptance::cap_support_continuity_diagnostic(
            result, &cap_faces, 5.0e-3, 5.0e-2,
        )
    {
        repair_report.apex_diagnostics.push(diagnostic);
    }
    if let Some(route) = repair_report
        .apex_routes
        .iter_mut()
        .rev()
        .find(|route| route.vertex == vertex)
    {
        route.built_faces = cap_faces;
    }
}
